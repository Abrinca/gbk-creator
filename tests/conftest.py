import datetime
import os
import random
from pathlib import Path

import pytest

from gbk_creator import GenBankBuilder, Taxon
from gbk_creator.taxonomy import NcbiClient

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
ID = "FAM13496-i1-1.1"
BUILD_DATE = datetime.datetime(2026, 1, 2, 3, 4, 5)

TAXON = Taxon(
    taxid=1308, scientific_name="Streptococcus thermophilus", rank="species",
    lineage=("Bacteria", "Bacillati", "Bacillota", "Bacilli", "Lactobacillales", "Streptococcaceae",
             "Streptococcus"))


def random_dna(n: int, seed: int = 1) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


@pytest.fixture(autouse=True)
def _no_network(request, monkeypatch):
    """Unit tests must run offline: any real E-utilities request fails the test."""
    if request.node.get_closest_marker("network"):
        return
    # offline tests never see real credentials (they could otherwise show up in assertion output)
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    if request.node.get_closest_marker("transport"):
        return  # these mock urllib.request.urlopen themselves, one level below _http_post

    def refuse(self, endpoint, params):
        raise AssertionError(f"unexpected network access ({endpoint}) in an offline test")

    monkeypatch.setattr(NcbiClient, "_http_post", refuse)


@pytest.fixture
def taxon():
    return TAXON


@pytest.fixture
def make_builder():
    def make(identifier=ID, taxid=1308, strain="FAM13496", **kw):
        kw.setdefault("taxon", TAXON)
        kw.setdefault("build_date", BUILD_DATE)
        return GenBankBuilder(identifier, taxid, strain=strain, **kw)
    return make


@pytest.fixture
def builder(make_builder):
    """Builder with two linear contigs of 2000 nt each, no genes yet."""
    b = make_builder()
    b.add_contig(f"{ID}_scf1", random_dna(2000, 1), "linear", "scaffold")
    b.add_contig(f"{ID}_scf2", random_dna(2000, 2), "linear", "scaffold")
    return b


def tag(n: int) -> str:
    return f"{ID}_{n:06d}"


def sample_dirs() -> list[str]:
    if not DATA.is_dir():
        return []
    return sorted(p.name for p in DATA.iterdir() if (p / "2_cds").is_dir())


def require_sample(name: str) -> Path:
    path = DATA / name
    if not (path / "2_cds").is_dir():
        pytest.skip(f"test data {path} not available")
    return path


def ncbi_available() -> bool:
    return bool(os.environ.get("NCBI_EMAIL"))


# ----------------------------------------------------------------------------- real-data helpers
CORE_SAMPLES = ["FAM13496-i1-1.1", "FAM14217-p1-1.1"]  # linear scaffolds / circular chromosome (-p1-)


# Per-sample items the frontend must drop in lenient mode; samples not listed are built strictly (any
# drop would raise). Entries: (kind, locus tag number, detail). All of them are PGAP quirks:
# * invalid_annotation: PGAP's preliminary EC numbers (".nN") do not match Arx's EC regex
# * orphan_gene: PGAP emitted a `gene` without any CDS/RNA feature
EXPECTED_LENIENT_DROPS = {
    "FAM1079-i1-1.1": [('invalid_annotation', '000226', 'EC:3.6.5.n1')],
    "FAM13496-i1-1.1": [('invalid_annotation', '001078', 'EC:3.6.5.n1')],
    "FAM14217-p1-1.1": [('invalid_annotation', '001792', 'EC:3.6.5.n1')],
    "FAM17891-i1-1.1": [('invalid_annotation', '001437', 'EC:3.6.5.n1')],
    "FAM17927-i1-1.1": [('invalid_annotation', '001437', 'EC:3.6.5.n1')],
    "FAM18356-i1-1.1": [('invalid_annotation', '001384', 'EC:3.6.5.n1')],
    "FAM18815-i1-1.1": [('orphan_gene', '001616', ''), ('invalid_annotation', '001275', 'EC:3.6.5.n1')],
    "FAM18969-i1-1.1": [('invalid_annotation', '000249', 'EC:3.6.5.n1')],
    "FAM19038-p1-1.1": [('invalid_annotation', '002598', 'EC:3.6.5.n1')],
    "FAM19144-i1-1.1": [('invalid_annotation', '001400', 'EC:3.6.5.n1')],
    "FAM19471-i1-1.1": [('invalid_annotation', '001414', 'EC:3.6.5.n1')],
    "FAM20446-i1-1.1": [('invalid_annotation', '000144', 'EC:3.6.5.n1')],
    "FAM20558-i1-1.1": [('invalid_annotation', '001347', 'EC:3.6.5.n1')],
    "FAM20833-i1-1.1": [('invalid_annotation', '001050', 'EC:3.6.5.n1'), ('invalid_annotation', '002467', 'EC:2.5.1.n9')],
    "FAM21277-i1-1.1": [('invalid_annotation', '000720', 'EC:3.6.5.n1')],
    "FAM21784-i1-1.1": [('invalid_annotation', '001223', 'EC:3.6.5.n1')],
    "FAM21789-i1-1.1": [('invalid_annotation', '000250', 'EC:3.6.5.n1')],
    "FAM22472-i1-1.1": [('invalid_annotation', '000715', 'EC:3.6.5.n1')],
    "FAM23169-i1-1.1": [('invalid_annotation', '000105', 'EC:3.6.5.n1')],
    "FAM23217-i1-1.1": [('invalid_annotation', '000062', 'EC:3.6.5.n1')],
    "FAM24227-i1-1.1": [('invalid_annotation', '000183', 'EC:3.6.5.n1')],
    "FAM24235-i1-1.1": [('invalid_annotation', '000854', 'EC:3.6.5.n1')],
    "FAM3257-i1-1.1": [('invalid_annotation', '001493', 'EC:3.6.5.n1')],
}


def expected_drops(sample: str) -> list[tuple[str, str, str]]:
    return [(kind, f"{sample}_{num}", detail) for kind, num, detail in EXPECTED_LENIENT_DROPS.get(sample, [])]


def oracle_taxon(sample: str) -> Taxon:
    """Taxon rebuilt from the oracle gbk itself, so real-data tests run offline."""
    from Bio import SeqIO
    rec = next(SeqIO.parse(DATA / sample / "2_cds" / f"{sample}.gbk", "genbank"))
    taxid = int(rec.features[0].qualifiers["db_xref"][0].split(":")[1])
    return Taxon(taxid, rec.annotations["organism"], "species", tuple(rec.annotations["taxonomy"]))


@pytest.fixture(scope="session")
def build_sample(tmp_path_factory):
    """``build_sample(name) -> (builder, out_dir)``: run the GFF frontend once per sample and session."""
    from gbk_creator.frontends.arx_gff import build_from_arx_gff, find_sample_files
    cache = {}

    def build(name):
        if name not in cache:
            src = require_sample(name)
            files = find_sample_files(src)
            dropped = []
            b = build_from_arx_gff(files["gff"], files["fna"], faa=files["faa"], ffn=files["ffn"],
                                   eggnog=files["eggnog"], taxon=oracle_taxon(name), build_date=BUILD_DATE,
                                   lenient=name in EXPECTED_LENIENT_DROPS, dropped=dropped)
            assert [(d.kind, d.locus_tag, d.detail) for d in dropped] == expected_drops(name)
            out = tmp_path_factory.mktemp(name) / name
            b.write(out)
            cache[name] = (b, out)
        return cache[name]

    return build
