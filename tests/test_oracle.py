"""Frontend output vs the PGAP oracle (``data/{sample}/2_cds/{sample}.gbk``); skipped without data."""
import re

import pytest
from Bio import SeqIO

from conftest import CORE_SAMPLES, DATA, expected_drops, ncbi_available, oracle_taxon, require_sample, sample_dirs
from gbk_creator import load_categories
from oracle_diff import compare_gbk, oracle_ec

ALL = sample_dirs() or CORE_SAMPLES  # without data the tests below skip themselves


@pytest.mark.parametrize("sample", ALL)
def test_semantic_diff_against_oracle(sample, build_sample):
    _, out = build_sample(sample)
    stats = {}
    diffs = compare_gbk(DATA / sample / "2_cds" / f"{sample}.gbk", out / f"{sample}.gbk", stats)
    assert not diffs, f"{len(diffs)} differences, first: {diffs[:5]}"


def test_diff_tool_detects_differences(build_sample, tmp_path):
    """The comparison must not be vacuous."""
    sample = CORE_SAMPLES[0]
    _, out = build_sample(sample)
    text = (out / f"{sample}.gbk").read_text()
    oracle = DATA / sample / "2_cds" / f"{sample}.gbk"
    for name, mutated in {
        "product": text.replace('/product="permease"', '/product="permeasE"', 1),
        "location": text.replace("complement(146..586)", "complement(146..587)", 1),
        "qualifier removed": text.replace('/transl_table=11\n', '', 1),
        "EC_number in output": text.replace('/transl_table=11\n', '/EC_number="1.1.1.1"\n                     /transl_table=11\n', 1),
        "partial marker dropped": text.replace("<1..100", "1..100"),
        "strain": text.replace('/strain="FAM13496"', '/strain="FAM13497"', 1),
    }.items():
        assert mutated != text, name
        p = tmp_path / "mutated.gbk"
        p.write_text(mutated)
        assert compare_gbk(oracle, p), name


@pytest.mark.parametrize("sample", CORE_SAMPLES)
def test_ec_sidecar_contains_oracle_ec_numbers(sample, build_sample):
    _, out = build_sample(sample)
    regex = load_categories()["EC"]
    got = {}
    for line in (out / f"{sample}.EC").read_text().splitlines():
        lt, annos = line.split("\t")
        got[lt] = set(annos.split(", "))
    expected = {lt: {e for e in ecs if regex.match(e)} for lt, ecs in oracle_ec(DATA / sample / "2_cds" / f"{sample}.gbk").items()}
    assert expected
    missing = {lt: e - got.get(lt, set()) for lt, e in expected.items() if e - got.get(lt, set())}
    assert not missing
    skipped = [e for ecs in oracle_ec(DATA / sample / "2_cds" / f"{sample}.gbk").values() for e in ecs if not regex.match(e)]
    assert all(re.search(r"\.n\d+$", e) for e in skipped), skipped  # only PGAP's preliminary EC numbers


@pytest.mark.parametrize("sample", CORE_SAMPLES)
def test_eggnog_mapping_in_sidecars(sample, build_sample):
    """Independent re-implementation of Arx's load_eggnog_file column mapping for a handful of rows."""
    src = require_sample(sample)
    path = next((src / "3_annotation" / "eggnog").glob("*.emapper.annotations"))
    sidecars = {}
    for p in (build_sample(sample)[1]).glob(f"{sample}.*"):
        if p.suffix != ".gbk":
            sidecars[p.suffix[1:]] = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1].split(", ") for l in p.read_text().splitlines()}
    n = 0
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        c = line.split("\t")
        lt = c[0].split("|")[-1]
        want = {
            "KG": [x[3:] for x in c[8].split(",")] if c[8] else [],
            "KR": c[11].split(",") if c[11] else [],
            "GO": c[6].split(",") if c[6] else [],
            "EP": [f"EP:{c[5]}"] if c[5] else [],
            "EO": [f"EO:{c[18].split('@')[0]}"] if c[18] else [],
        }
        for t, ids in want.items():
            for i in ids:
                assert i in sidecars[t][lt], (t, lt, i)
        n += 1
    assert n > 100
    # prefixes the Arx loader adds
    for lt, ids in sidecars["EP"].items():
        assert all(i.startswith("EP:") for i in ids)
    for ids in sidecars["ED"].values():
        assert all(i.startswith("ED:") for i in ids)


@pytest.mark.parametrize("sample", CORE_SAMPLES)
def test_sidecar_rules_on_real_data(sample, build_sample):
    b, out = build_sample(sample)
    cats = load_categories()
    cds = {g.locus_tag for g in b.genes if g.feature_type == "CDS"}
    assert not (out / f"{sample}.GC").exists() and not (out / f"{sample}.GP").exists()
    for p in out.glob(f"{sample}.*"):
        if p.suffix == ".gbk":
            continue
        lines = p.read_text().splitlines(keepends=True)
        assert lines and all(l.endswith("\n") and l.strip() for l in lines)
        tags = [l.split("\t")[0] for l in lines]
        assert tags == sorted(tags) and len(set(tags)) == len(tags) and set(tags) <= cds
        for l in lines:
            ids = l.rstrip("\n").split("\t")[1].split(", ")
            assert ids == sorted(set(ids))
            assert all(cats[p.suffix[1:]].match(i) for i in ids)


@pytest.mark.parametrize("sample", CORE_SAMPLES)
def test_header_counts_match_oracle(sample, build_sample):
    _, out = build_sample(sample)
    oracle = next(SeqIO.parse(DATA / sample / "2_cds" / f"{sample}.gbk", "genbank"))
    ours = next(SeqIO.parse(out / f"{sample}.gbk", "genbank"))
    o = oracle.annotations["structured_comment"]["Genome-Annotation-Data"]
    g = ours.annotations["structured_comment"]["Genome-Annotation-Data"]
    for key in ("Genes (total)", "CDSs (total)", "Genes (coding)", "CDSs (with protein)", "Genes (RNA)",
                "tRNAs", "Pseudo Genes (total)", "CDSs (without protein)"):
        assert int(g[key]) == int(o[key].replace(",", "")), key


def test_circular_sample_topology_and_wraparound(build_sample):
    sample = "FAM14217-p1-1.1"
    b, out = build_sample(sample)
    rec = next(SeqIO.parse(out / f"{sample}.gbk", "genbank"))
    assert rec.annotations["topology"] == "circular"
    assert "chromosome" in rec.description
    first_gene = rec.features[1]
    assert first_gene.type == "gene" and len(first_gene.location.parts) == 2  # join across the origin
    assert b.contigs[0].kind == "chromosome"


@pytest.mark.network
@pytest.mark.skipif(not ncbi_available(), reason="NCBI_EMAIL not set")
@pytest.mark.parametrize("sample", ["FAM13496-i1-1.1", "FAM19471-i1-1.1"])  # 2nd: PGAP name differs from NCBI today
def test_with_real_ncbi_taxonomy(sample, tmp_path):
    """Full path including the E-utilities lookup (cached in the default cache dir)."""
    from gbk_creator.frontends.arx_gff import build_from_arx_gff, find_sample_files
    src = require_sample(sample)
    f = find_sample_files(src)
    dropped = []
    b = build_from_arx_gff(f["gff"], f["fna"], faa=f["faa"], ffn=f["ffn"], eggnog=f["eggnog"],
                           lenient=True, dropped=dropped)
    assert [(d.kind, d.locus_tag, d.detail) for d in dropped] == expected_drops(sample)
    b.write(tmp_path)
    oracle = DATA / sample / "2_cds" / f"{sample}.gbk"
    assert b.taxon.taxid == oracle_taxon(sample).taxid
    assert compare_gbk(oracle, tmp_path / f"{sample}.gbk", ignore_source_organism=True) == []
    rec = next(SeqIO.parse(tmp_path / f"{sample}.gbk", "genbank"))
    assert rec.annotations["organism"] == b.taxon.scientific_name
    assert rec.annotations["taxonomy"][0] == "Bacteria"


@pytest.mark.parametrize("sample", ALL)
def test_lenient_drops_are_exactly_the_known_pgap_quirks(sample):
    """Cross-check of EXPECTED_LENIENT_DROPS against the oracle itself (independent of the frontend)."""
    src = require_sample(sample)
    oracle = src / "2_cds" / f"{sample}.gbk"
    regex = load_categories()["EC"]
    want = [("invalid_annotation", lt, ec) for lt, ecs in oracle_ec(oracle).items() for ec in sorted(ecs) if not regex.match(ec)]
    for rec in SeqIO.parse(oracle, "genbank"):
        tags = {f.qualifiers["locus_tag"][0] for f in rec.features if f.type not in ("gene", "source") and "locus_tag" in f.qualifiers}
        want += [("orphan_gene", f.qualifiers["locus_tag"][0], "") for f in rec.features
                 if f.type == "gene" and f.qualifiers["locus_tag"][0] not in tags]
    assert sorted(want) == sorted(expected_drops(sample))
