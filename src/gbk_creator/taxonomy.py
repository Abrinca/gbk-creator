"""NCBI taxonomy helpers (E-utilities, stdlib HTTP, cached on disk).

1. ``name_to_taxid``: latin name -> taxid (esearch; exactly one hit required)
2. ``fetch_taxon``:   taxid -> scientific name + rank + lineage (efetch)

Credentials come from arguments or the env vars ``NCBI_EMAIL`` / ``NCBI_API_KEY``
(``.env`` files are never read here). The API key is sent in the POST body only and is
never cached, logged or put into exception messages.

Cache directory: ``$GBK_CREATOR_CACHE``, else ``$XDG_CACHE_HOME/gbk-creator``, else
``~/.cache/gbk-creator``. Entries never expire; delete the directory to refresh.
"""
from __future__ import annotations

import hashlib
import http.client
import json
import os
import socket
import ssl
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass
from pathlib import Path

from .errors import TaxonomyError

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOOL = "gbk-creator"
_RETRY_STATUS = {429, 500, 502, 503, 504}

_rate_lock = threading.Lock()
_last_request = 0.0


@dataclass(frozen=True)
class Taxon:
    taxid: int
    scientific_name: str
    rank: str
    lineage: tuple[str, ...]  # GenBank style: no "cellular organisms", excludes the taxon itself


def default_cache_dir() -> Path:
    env = os.environ.get("GBK_CREATOR_CACHE")
    if env:
        return Path(env).expanduser()
    base = os.environ.get("XDG_CACHE_HOME") or str(Path.home() / ".cache")
    return Path(base) / "gbk-creator"


def _classify_transport_error(e: BaseException) -> str:
    """Safe, fixed-vocabulary description of a transport failure (never the exception text, which may
    contain the request URL or body and thus the credentials)."""
    if isinstance(e, urllib.error.URLError):
        reason = e.reason
        if isinstance(reason, TimeoutError):
            return "timeout"
        if isinstance(reason, socket.gaierror):
            return "DNS lookup failed"
        if isinstance(reason, ssl.SSLError):
            return "TLS error"
        if isinstance(reason, ConnectionError):
            return "connection failed"
        return "network error"
    if isinstance(e, TimeoutError):
        return "timeout"
    if isinstance(e, ConnectionError):
        return "connection failed"
    if isinstance(e, http.client.HTTPException):
        return "malformed HTTP response"
    return f"network error ({type(e).__name__})"


class NcbiClient:
    def __init__(self, email: str | None = None, api_key: str | None = None,
                 cache_dir: str | os.PathLike | None = None):
        self.email = email or os.environ.get("NCBI_EMAIL") or None
        self.api_key = api_key or os.environ.get("NCBI_API_KEY") or None
        self.cache_dir = Path(cache_dir).expanduser() if cache_dir else default_cache_dir()
        if not self.email or "@" not in self.email:
            raise TaxonomyError(
                "An e-mail address is required for NCBI E-utilities: pass ncbi_email=... "
                "or set the environment variable NCBI_EMAIL (an API key via NCBI_API_KEY "
                "is optional but raises the rate limit)."
            )

    # -- cache -------------------------------------------------------------------
    def _cache_path(self, kind: str, key: str) -> Path:
        digest = hashlib.sha256(key.encode()).hexdigest()[:24]
        return self.cache_dir / "taxonomy" / f"{kind}-{digest}.json"

    def _cache_get(self, kind: str, key: str):
        path = self._cache_path(kind, key)
        try:
            return json.loads(path.read_text())
        except (OSError, ValueError):
            return None

    def _cache_put(self, kind: str, key: str, value) -> None:
        path = self._cache_path(kind, key)
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            tmp = path.with_name(f"{path.name}.{os.getpid()}.tmp")
            tmp.write_text(json.dumps(value))
            os.replace(tmp, path)
        except OSError:
            pass  # cache is best effort

    # -- HTTP --------------------------------------------------------------------
    def _http_post(self, endpoint: str, params: dict) -> bytes:
        """POST to an E-utility with rate limiting and retries. Overridable in tests."""
        global _last_request
        params = dict(params, tool=TOOL, email=self.email)
        if self.api_key:
            params["api_key"] = self.api_key
        data = urllib.parse.urlencode(params).encode()
        interval = 0.11 if self.api_key else 0.34
        last_error = "unknown error"
        attempts = 5
        for attempt in range(attempts):
            with _rate_lock:
                wait = _last_request + interval - time.monotonic()
                if wait > 0:
                    time.sleep(wait)
                _last_request = time.monotonic()
            retry = True
            try:
                with urllib.request.urlopen(EUTILS + endpoint, data=data, timeout=30) as resp:
                    return resp.read()
            except urllib.error.HTTPError as e:
                last_error = f"HTTP {int(e.code)}"
                retry = e.code in _RETRY_STATUS
            except Exception as e:  # noqa: BLE001 - only a safe classification of it is kept
                last_error = _classify_transport_error(e)
            if not retry:
                break
            if attempt < attempts - 1:
                time.sleep(min(2 ** attempt, 8))
        # raised outside the except blocks: no __cause__/__context__ that could carry request details
        raise TaxonomyError(self._redact(f"NCBI E-utilities request to {endpoint} failed ({last_error})."))

    def _redact(self, text: str) -> str:
        """Defense in depth: remove credentials (and their URL-encoded forms) from any message."""
        for secret in (self.api_key, self.email):
            if secret:
                for variant in {secret, urllib.parse.quote(secret), urllib.parse.quote_plus(secret),
                                urllib.parse.quote(secret, safe="")}:
                    text = text.replace(variant, "<redacted>")
        return text

    def __repr__(self) -> str:
        return f"NcbiClient(cache_dir={str(self.cache_dir)!r}, credentials=<redacted>)"

    # -- API ---------------------------------------------------------------------
    def name_to_taxid(self, name: str) -> int:
        name = name.strip()
        if not name:
            raise TaxonomyError("Empty organism name.")
        key = name.lower()
        cached = self._cache_get("name", key)
        if cached is not None:
            return int(cached["taxid"])
        raw = self._http_post("esearch.fcgi", {"db": "taxonomy", "term": name, "retmode": "json"})
        try:
            result = json.loads(raw)["esearchresult"]
        except (ValueError, KeyError):
            raise TaxonomyError(f"Unexpected esearch response for {name!r}.") from None
        ids = result.get("idlist", [])
        if "ERROR" in result or "errorlist" in result or len(ids) != 1 or not ids[0].isdigit():
            raise TaxonomyError(
                f"Organism name {name!r} must match exactly one NCBI taxonomy entry, "
                f"found {len(ids)}: {ids[:10]}. Use a more specific name or pass the taxid directly."
            )
        taxid = int(ids[0])
        self._cache_put("name", key, {"name": name, "taxid": taxid})
        return taxid

    def fetch_taxon(self, taxid: int) -> Taxon:
        if isinstance(taxid, bool) or not isinstance(taxid, int) or taxid <= 0:
            raise TaxonomyError(f"taxid must be a positive integer, got {taxid!r}.")
        cached = self._cache_get("taxid", str(taxid))
        if cached is not None:
            return Taxon(cached["taxid"], cached["scientific_name"], cached["rank"],
                         tuple(cached["lineage"]))
        raw = self._http_post("efetch.fcgi", {"db": "taxonomy", "id": str(taxid), "retmode": "xml"})
        taxon = parse_taxon_xml(raw, taxid)
        self._cache_put("taxid", str(taxid), asdict(taxon) | {"lineage": list(taxon.lineage)})
        return taxon


def parse_taxon_xml(raw: bytes | str, requested_taxid: int) -> Taxon:
    try:
        root = ET.fromstring(raw)
    except ET.ParseError:
        raise TaxonomyError(f"Unparsable efetch response for taxid {requested_taxid}.") from None
    nodes = root.findall("Taxon")
    if len(nodes) != 1:
        raise TaxonomyError(f"{requested_taxid} is not a valid NCBI taxid (efetch returned no taxon).")
    node = nodes[0]
    taxid = int((node.findtext("TaxId") or "0").strip())
    name = (node.findtext("ScientificName") or "").strip()
    if not name or not taxid:
        raise TaxonomyError(f"Incomplete efetch response for taxid {requested_taxid}.")
    if taxid != requested_taxid:
        raise TaxonomyError(
            f"taxid {requested_taxid} has been merged into {taxid} ({name!r}); please use {taxid}.")
    rank = (node.findtext("Rank") or "no rank").strip()
    lineage = tuple(p.strip() for p in (node.findtext("Lineage") or "").split(";") if p.strip())
    lineage = tuple(p for p in lineage if p != "cellular organisms")
    return Taxon(taxid=taxid, scientific_name=name, rank=rank, lineage=lineage)


def name_to_taxid(name: str, *, email: str | None = None, api_key: str | None = None,
                  cache_dir: str | os.PathLike | None = None) -> int:
    """Latin name -> NCBI taxid. Raises unless NCBI finds exactly one taxonomy entry."""
    return NcbiClient(email, api_key, cache_dir).name_to_taxid(name)


def fetch_taxon(taxid: int, *, email: str | None = None, api_key: str | None = None,
                cache_dir: str | os.PathLike | None = None) -> Taxon:
    """taxid -> scientific name, rank and lineage."""
    return NcbiClient(email, api_key, cache_dir).fetch_taxon(taxid)
