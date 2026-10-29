"""NCBI taxonomy helpers: offline tests with a faked HTTP layer + marked network tests."""
import json

import pytest

from conftest import ncbi_available
from gbk_creator import GenBankBuilder, Taxon, TaxonomyError, fetch_taxon, name_to_taxid
from gbk_creator.taxonomy import NcbiClient, default_cache_dir, parse_taxon_xml

EFETCH_XML = """<?xml version="1.0" ?>
<!DOCTYPE TaxaSet PUBLIC "-//NLM//DTD Taxon, 14th January 2002//EN" "https://www.ncbi.nlm.nih.gov/entrez/query/DTD/taxon.dtd">
<TaxaSet><Taxon>
    <TaxId>1308</TaxId>
    <ScientificName>Streptococcus thermophilus</ScientificName>
    <OtherNames><Synonym>Streptococcus salivarius subsp. thermophilus</Synonym></OtherNames>
    <ParentTaxId>1307</ParentTaxId>
    <Rank>species</Rank>
    <Division>Bacteria</Division>
    <Lineage>cellular organisms; Bacteria; Bacillati; Bacillota; Bacilli; Lactobacillales; Streptococcaceae; Streptococcus</Lineage>
    <LineageEx>
      <Taxon><TaxId>131567</TaxId><ScientificName>cellular organisms</ScientificName><Rank>cellular root</Rank></Taxon>
      <Taxon><TaxId>2</TaxId><ScientificName>Bacteria</ScientificName><Rank>domain</Rank></Taxon>
    </LineageEx>
</Taxon></TaxaSet>"""


class FakeNcbi(NcbiClient):
    def __init__(self, responses, **kw):
        super().__init__(**kw)
        self.responses, self.calls = responses, []

    def _http_post(self, endpoint, params):
        self.calls.append((endpoint, params))
        return self.responses[endpoint]


def esearch(ids):
    return json.dumps({"esearchresult": {"count": str(len(ids)), "idlist": ids}}).encode()


@pytest.fixture
def client_factory(tmp_path):
    def make(responses):
        return FakeNcbi(responses, email="me@example.org", api_key="SECRET", cache_dir=tmp_path / "cache")
    return make


def test_parse_efetch(taxon):
    t = parse_taxon_xml(EFETCH_XML, 1308)
    assert t == Taxon(1308, "Streptococcus thermophilus", "species", taxon.lineage)


def test_fetch_taxon_is_cached(client_factory):
    c = client_factory({"efetch.fcgi": EFETCH_XML.encode()})
    assert c.fetch_taxon(1308).scientific_name == "Streptococcus thermophilus"
    assert c.fetch_taxon(1308).rank == "species"
    assert len(c.calls) == 1
    # a fresh client with the same cache dir does not touch the network either
    c2 = client_factory({})
    assert c2.fetch_taxon(1308).lineage[0] == "Bacteria"
    assert c2.calls == []


def test_api_key_is_never_cached(client_factory, tmp_path):
    c = client_factory({"efetch.fcgi": EFETCH_XML.encode(), "esearch.fcgi": esearch(["1308"])})
    c.fetch_taxon(1308)
    c.name_to_taxid("Streptococcus thermophilus")
    files = list((tmp_path / "cache").rglob("*.json"))
    assert len(files) == 2
    assert all("SECRET" not in f.read_text() and "example.org" not in f.read_text() for f in files)


def test_name_to_taxid_exactly_one_hit(client_factory):
    c = client_factory({"esearch.fcgi": esearch(["1308"])})
    assert c.name_to_taxid(" Streptococcus thermophilus ") == 1308
    assert c.name_to_taxid("streptococcus thermophilus") == 1308  # cached, case-insensitive
    assert len(c.calls) == 1
    assert c.calls[0][1]["db"] == "taxonomy"


@pytest.mark.parametrize("ids", [[], ["1", "2"]])
def test_name_to_taxid_rejects_zero_or_many_hits(client_factory, ids):
    c = client_factory({"esearch.fcgi": esearch(ids)})
    with pytest.raises(TaxonomyError, match="exactly one"):
        c.name_to_taxid("Bacillus")


def test_failures_are_not_cached(client_factory):
    c = client_factory({"esearch.fcgi": esearch([])})
    for _ in range(2):
        with pytest.raises(TaxonomyError):
            c.name_to_taxid("nonsense")
    assert len(c.calls) == 2


def test_invalid_taxid(client_factory):
    c = client_factory({"efetch.fcgi": b"<?xml version='1.0'?><TaxaSet></TaxaSet>"})
    with pytest.raises(TaxonomyError, match="not a valid NCBI taxid"):
        c.fetch_taxon(999999999)
    for bad in (0, -1, "1308", None):
        with pytest.raises(TaxonomyError):
            c.fetch_taxon(bad)


def test_merged_taxid():
    xml = EFETCH_XML.replace("<TaxId>1308</TaxId>", "<TaxId>1309</TaxId>", 1)
    with pytest.raises(TaxonomyError, match="merged into 1309"):
        parse_taxon_xml(xml, 1308)


def test_garbage_response():
    with pytest.raises(TaxonomyError, match="Unparsable"):
        parse_taxon_xml(b"<html>", 1)


def test_email_is_required(monkeypatch):
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    with pytest.raises(TaxonomyError, match="NCBI_EMAIL"):
        NcbiClient()
    with pytest.raises(TaxonomyError, match="NCBI_EMAIL"):
        NcbiClient(email="not-an-address")


def test_credentials_from_env(monkeypatch, tmp_path):
    monkeypatch.setenv("NCBI_EMAIL", "env@example.org")
    monkeypatch.setenv("NCBI_API_KEY", "envkey")
    c = NcbiClient(cache_dir=tmp_path)
    assert (c.email, c.api_key) == ("env@example.org", "envkey")
    c = NcbiClient(email="arg@example.org", api_key="argkey", cache_dir=tmp_path)
    assert (c.email, c.api_key) == ("arg@example.org", "argkey")


def test_cache_dir_resolution(monkeypatch, tmp_path):
    monkeypatch.setenv("GBK_CREATOR_CACHE", str(tmp_path / "c"))
    assert default_cache_dir() == tmp_path / "c"
    monkeypatch.delenv("GBK_CREATOR_CACHE")
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "xdg"))
    assert default_cache_dir() == tmp_path / "xdg" / "gbk-creator"


def test_builder_uses_lookup_and_header_uses_scientific_name(monkeypatch, tmp_path, taxon):
    seen = {}

    def fake_fetch(self, taxid):
        seen["taxid"], seen["email"], seen["key"] = taxid, self.email, self.api_key
        return taxon

    monkeypatch.setattr(NcbiClient, "fetch_taxon", fake_fetch)
    b = GenBankBuilder("X-1", 1308, strain="X", ncbi_email="a@b.org", ncbi_api_key="k")
    assert b.taxon is taxon
    assert seen == {"taxid": 1308, "email": "a@b.org", "key": "k"}


@pytest.mark.network
@pytest.mark.skipif(not ncbi_available(), reason="NCBI_EMAIL not set")
def test_network_roundtrip(tmp_path):
    assert name_to_taxid("Streptococcus thermophilus", cache_dir=tmp_path) == 1308
    t = fetch_taxon(1308, cache_dir=tmp_path)
    assert t.scientific_name == "Streptococcus thermophilus"
    assert t.rank == "species"
    assert t.lineage[0] == "Bacteria" and t.lineage[-1] == "Streptococcus"
    assert "cellular organisms" not in t.lineage
    assert fetch_taxon(1308, cache_dir=tmp_path) == t  # served from the cache
    with pytest.raises(TaxonomyError):
        name_to_taxid("Definitely notanorganism xyzzy", cache_dir=tmp_path)
    with pytest.raises(TaxonomyError):
        fetch_taxon(2_000_000_000, cache_dir=tmp_path)
