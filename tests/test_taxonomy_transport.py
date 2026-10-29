"""NcbiClient._http_post for real (urlopen mocked): credential hygiene, retries, rate limiting. Offline.

The credentials below are synthetic canaries; they must never show up in exception messages, tracebacks,
logs, reprs or cache files, in any failure mode.
"""
import http.client
import io
import logging
import socket
import ssl
import traceback
import urllib.error
import urllib.parse

import pytest

from gbk_creator import TaxonomyError
from gbk_creator import taxonomy
from gbk_creator.taxonomy import NcbiClient

EMAIL = "leak.canary+x@example.org"
KEY = "KEYCANARY0123456789abcdef&=%"
SECRETS = [EMAIL, KEY, urllib.parse.quote(EMAIL), urllib.parse.quote_plus(EMAIL),
           urllib.parse.quote(KEY), urllib.parse.quote_plus(KEY), "leak.canary", "KEYCANARY"]
URL_WITH_SECRETS = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?email={EMAIL}&api_key={KEY}"

XML = b"""<?xml version="1.0"?><TaxaSet><Taxon><TaxId>1308</TaxId><ScientificName>Streptococcus thermophilus</ScientificName>
<Rank>species</Rank><Lineage>cellular organisms; Bacteria</Lineage></Taxon></TaxaSet>"""

pytestmark = pytest.mark.transport


class Resp:
    def __init__(self, body=XML):
        self.body = body

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def read(self):
        return self.body


def http_error(code):
    return urllib.error.HTTPError(URL_WITH_SECRETS, code, f"msg {KEY} {EMAIL}", {}, io.BytesIO(b""))


class Net:
    """Fake urlopen + fake clock. `script` items: an Exception (raised) or bytes (returned)."""

    def __init__(self, monkeypatch, script):
        self.script, self.calls, self.sleeps, self.now = list(script), [], [], 1000.0
        monkeypatch.setattr(taxonomy.urllib.request, "urlopen", self.urlopen)
        monkeypatch.setattr(taxonomy.time, "sleep", self.sleep)
        monkeypatch.setattr(taxonomy.time, "monotonic", lambda: self.now)
        monkeypatch.setattr(taxonomy, "_last_request", 0.0)

    def sleep(self, seconds):
        self.sleeps.append(seconds)
        self.now += seconds

    def urlopen(self, url, data=None, timeout=None):
        self.calls.append((url, data, timeout))
        item = self.script.pop(0) if len(self.script) > 1 else self.script[0]
        if isinstance(item, BaseException):
            raise item
        return Resp(item)


@pytest.fixture
def client(tmp_path):
    return NcbiClient(email=EMAIL, api_key=KEY, cache_dir=tmp_path / "cache")


def everything_about(exc, caplog, cache_dir) -> str:
    parts = [str(exc), repr(exc), repr(exc.args), "".join(traceback.format_exception(exc)), caplog.text]
    parts += [p.read_text() for p in cache_dir.rglob("*") if p.is_file()]
    parts += [str(exc.__cause__), str(exc.__context__)]
    return "\n".join(parts)


FAILURES = {
    "urlerror text": urllib.error.URLError(f"<urlopen error {URL_WITH_SECRETS}>"),
    "dns": urllib.error.URLError(socket.gaierror(-2, f"Name or service not known {KEY} {EMAIL}")),
    "tls": urllib.error.URLError(ssl.SSLError(1, f"bad certificate {URL_WITH_SECRETS}")),
    "timeout": TimeoutError(f"timed out {KEY}"),
    "urlerror timeout": urllib.error.URLError(TimeoutError(EMAIL)),
    "connection reset": ConnectionResetError(f"reset {URL_WITH_SECRETS}"),
    "remote disconnected": http.client.RemoteDisconnected(f"closed {KEY}"),
    "incomplete read": http.client.IncompleteRead(KEY.encode(), 5),
    "oserror": OSError(f"disk {EMAIL}"),
    "valueerror": ValueError(f"unknown url type: {URL_WITH_SECRETS}"),
    "http 400": http_error(400),
    "http 404": http_error(404),
    "http 429": http_error(429),
    "http 500": http_error(500),
}


@pytest.mark.parametrize("name", FAILURES)
def test_credentials_never_leak_on_any_failure(name, client, monkeypatch, caplog, tmp_path):
    Net(monkeypatch, [FAILURES[name]])
    with caplog.at_level(logging.DEBUG):
        with pytest.raises(TaxonomyError) as info:
            client.fetch_taxon(1308)
    text = everything_about(info.value, caplog, tmp_path / "cache")
    for secret in SECRETS:
        assert secret not in text, (name, secret)
    assert info.value.__cause__ is None and info.value.__context__ is None
    assert "efetch.fcgi" in str(info.value)  # still useful: says which request failed


@pytest.mark.parametrize("exc,expected", [
    (urllib.error.URLError("x"), "network error"),
    (urllib.error.URLError(socket.gaierror(-2, "x")), "DNS lookup failed"),
    (urllib.error.URLError(ssl.SSLError("x")), "TLS error"),
    (TimeoutError("x"), "timeout"),
    (ConnectionResetError("x"), "connection failed"),
    (http.client.BadStatusLine("x"), "malformed HTTP response"),
    (ValueError("x"), "network error (ValueError)"),
    (http_error(404), "HTTP 404"),
])
def test_failures_are_classified_not_copied(exc, expected, client, monkeypatch):
    Net(monkeypatch, [exc])
    with pytest.raises(TaxonomyError) as info:
        client.fetch_taxon(1308)
    assert str(info.value) == f"NCBI E-utilities request to efetch.fcgi failed ({expected})."


def test_redaction_is_a_second_line_of_defense(client):
    msg = client._redact(f"a {EMAIL} b {urllib.parse.quote(EMAIL)} c {KEY} d {urllib.parse.quote_plus(KEY)}")
    assert msg == "a <redacted> b <redacted> c <redacted> d <redacted>"
    assert EMAIL not in repr(client) and KEY not in repr(client)


def test_credentials_travel_in_the_post_body_only(client, monkeypatch, tmp_path):
    net = Net(monkeypatch, [XML])
    assert client.fetch_taxon(1308).scientific_name == "Streptococcus thermophilus"
    url, data, timeout = net.calls[0]
    assert url == taxonomy.EUTILS + "efetch.fcgi" and timeout == 30
    assert "api_key" not in url and "@" not in url
    body = urllib.parse.parse_qs(data.decode())
    assert body["email"] == [EMAIL] and body["api_key"] == [KEY] and body["tool"] == ["gbk-creator"]
    assert body["db"] == ["taxonomy"] and body["id"] == ["1308"]
    # nothing secret is cached
    cached = [p.read_text() for p in (tmp_path / "cache").rglob("*") if p.is_file()]
    assert cached and not any(s in c for c in cached for s in SECRETS)


def test_no_api_key_means_no_api_key_param(monkeypatch, tmp_path):
    net = Net(monkeypatch, [XML])
    NcbiClient(email=EMAIL, api_key=None, cache_dir=tmp_path).fetch_taxon(1308)
    assert "api_key" not in urllib.parse.parse_qs(net.calls[0][1].decode())


# ----------------------------------------------------------------------------- retries
@pytest.mark.parametrize("status", [429, 500, 502, 503, 504])
def test_retryable_http_errors_are_retried_with_backoff(status, client, monkeypatch):
    net = Net(monkeypatch, [http_error(status), http_error(status), XML])
    assert client.fetch_taxon(1308).taxid == 1308
    assert len(net.calls) == 3
    assert net.sleeps == [1, 2]


@pytest.mark.parametrize("status", [400, 401, 403, 404])
def test_other_http_errors_are_not_retried(status, client, monkeypatch):
    net = Net(monkeypatch, [http_error(status), XML])
    with pytest.raises(TaxonomyError, match=f"HTTP {status}"):
        client.fetch_taxon(1308)
    assert len(net.calls) == 1 and net.sleeps == []


def test_gives_up_after_five_attempts(client, monkeypatch):
    net = Net(monkeypatch, [http_error(503)])
    with pytest.raises(TaxonomyError, match="HTTP 503"):
        client.fetch_taxon(1308)
    assert len(net.calls) == 5
    assert net.sleeps == [1, 2, 4, 8]  # no pointless sleep after the last attempt


def test_backoff_is_capped(client, monkeypatch):
    net = Net(monkeypatch, [http_error(503)] * 6 + [XML])
    monkeypatch.setattr(taxonomy.time, "monotonic", lambda: net.now)
    with pytest.raises(TaxonomyError):
        client.fetch_taxon(1308)
    assert max(net.sleeps) <= 8


@pytest.mark.parametrize("exc", [urllib.error.URLError("x"), TimeoutError("x"), ConnectionResetError("x"),
                                 http.client.RemoteDisconnected("x")])
def test_transport_errors_are_retried(exc, client, monkeypatch):
    net = Net(monkeypatch, [exc, XML])
    assert client.fetch_taxon(1308).taxid == 1308
    assert len(net.calls) == 2 and net.sleeps == [1]


def test_failed_lookups_are_not_cached(client, monkeypatch, tmp_path):
    Net(monkeypatch, [http_error(404)])
    with pytest.raises(TaxonomyError):
        client.fetch_taxon(1308)
    assert not list((tmp_path / "cache").rglob("*.json"))


# ----------------------------------------------------------------------------- rate limiting
@pytest.mark.parametrize("api_key,interval", [(KEY, 0.11), (None, 0.34)])
def test_requests_are_spaced(api_key, interval, monkeypatch, tmp_path):
    net = Net(monkeypatch, [XML])
    c = NcbiClient(email=EMAIL, api_key=api_key, cache_dir=tmp_path)
    for taxid in range(1, 11):
        c._http_post("efetch.fcgi", {"id": str(taxid)})
    assert len(net.calls) == 10
    assert net.sleeps == [pytest.approx(interval)] * 9  # first request is free, then one interval each


def test_no_extra_wait_when_enough_time_has_passed(monkeypatch, tmp_path):
    net = Net(monkeypatch, [XML])
    c = NcbiClient(email=EMAIL, api_key=KEY, cache_dir=tmp_path)
    c._http_post("efetch.fcgi", {})
    net.now += 5
    c._http_post("efetch.fcgi", {})
    assert net.sleeps == []


def test_rate_limit_is_shared_between_clients(monkeypatch, tmp_path):
    net = Net(monkeypatch, [XML])
    NcbiClient(email=EMAIL, api_key=KEY, cache_dir=tmp_path)._http_post("efetch.fcgi", {})
    NcbiClient(email=EMAIL, api_key=KEY, cache_dir=tmp_path)._http_post("efetch.fcgi", {})
    assert net.sleeps == [pytest.approx(0.11)]


def test_esearch_goes_through_the_same_path(client, monkeypatch):
    net = Net(monkeypatch, [b'{"esearchresult": {"count": "1", "idlist": ["1308"]}}'])
    assert client.name_to_taxid("Streptococcus thermophilus") == 1308
    assert net.calls[0][0].endswith("esearch.fcgi")
    assert urllib.parse.parse_qs(net.calls[0][1].decode())["term"] == ["Streptococcus thermophilus"]
