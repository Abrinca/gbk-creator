"""Format of the generated gbk + sidecars (offline, synthetic genome)."""
import re

import pytest
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation, SimpleLocation

from conftest import BUILD_DATE, random_dna
from gbk_creator import GbkCreatorError, GenBankBuilder, __version__
from gbk_creator.builder import classify_contig

LONG_ID = "S" * 20 + "-" + "9" * 29  # 50 characters


def loc(a, b, strand=1):
    return SimpleLocation(a, b, strand=strand)


@pytest.fixture
def genome(taxon, tmp_path):
    """chromosome (circular) + plasmid + 12 scaffolds, long identifier, all gene kinds, annotations."""
    ident = LONG_ID
    b = GenBankBuilder(ident, 1308, strain="S" * 20, taxon=taxon, build_date=BUILD_DATE)
    kinds = [("circular", "chromosome"), ("circular", "plasmid")] + [("linear", "scaffold")] * 10
    for i, (topo, kind) in enumerate(kinds, 1):
        b.add_contig(f"{ident}_scf{i}", random_dna(3000, i), topo, kind)
    c1, c2, c12 = (f"{ident}_scf{i}" for i in (1, 2, 12))
    lt = lambda n: f"{ident}_{n:06d}"
    wrap = CompoundLocation([loc(2900, 3000, -1), loc(0, 50, -1)][::-1])
    b.add_gene(lt(1), c1, wrap, "CDS", product="wrap protein", translation="M" * 50)
    b.add_gene(lt(2), c1, loc(100, 400, -1), "CDS", product="p2", gene="cas9", translation="MKV",
               inference=["a", "b"], protein_id=f"extdb:{lt(2)}", note="n")
    slip = CompoundLocation([loc(500, 800), loc(799, 1100)])
    b.add_gene(lt(3), c1, slip, "CDS", product="slip", translation="MKV", ribosomal_slippage=True)
    b.add_gene(lt(4), c1, SimpleLocation(1200, 1300, strand=1), "CDS", product="pseudo", pseudo=True)
    b.add_feature(c1, "repeat_region", loc(1400, 1500), rpt_family="CRISPR", rpt_type="direct",
                  rpt_unit_range="1405..1440", rpt_unit_seq="acgt")
    b.add_gene(lt(5), c1, loc(1600, 1676), "tRNA", product="tRNA-Gly", anticodon="(pos:1620..1622,aa:Gly,seq:tcc)")
    b.add_gene(lt(6), c1, loc(1700, 1800), "rRNA", product="16S ribosomal RNA")
    b.add_gene(lt(7), c1, loc(1900, 2000), "ncRNA", product="RNase P", gene="rnpB", ncRNA_class="RNase_P_RNA",
               db_xref="RFAM:RF00010")
    b.add_gene(lt(8), c1, loc(2100, 2200), "tmRNA", product="tmRNA", gene="ssrA", db_xref="RFAM:RF00023")
    b.add_feature(c1, "regulatory", loc(2300, 2400), regulatory_class="riboswitch", bound_moiety="x",
                  db_xref="RFAM:RF00174")
    b.add_gene(lt(9), c2, loc(10, 100), "CDS", product="plasmid protein", translation="MM")
    b.add_gene(lt(10), c12, loc(10, 100), "CDS", product="last", translation="MM")
    b.add_annotation(lt(10), "KG", "K00002")
    b.add_annotation(lt(2), "KG", "K00002")
    b.add_annotation(lt(2), "KG", "K00001")
    b.add_annotation(lt(4), "EC", "EC:1.1.1.-")
    b.add_annotation(lt(1), "EC", "EC:2.7.7.7")
    out = tmp_path / "out"
    b.write(out)
    return b, out, ident


def read(out, ident):
    return list(SeqIO.parse(out / f"{ident}.gbk", "genbank"))


def test_long_names_round_trip(genome):
    b, out, ident = genome
    recs = read(out, ident)
    assert [r.id for r in recs] == [f"{ident}_scf{i}" for i in range(1, 13)]
    assert all(r.name == r.id for r in recs)
    text = (out / f"{ident}.gbk").read_text()
    assert f"LOCUS       {ident}_scf12 3000 bp" in text  # widened LOCUS line (name + space + length)
    one = SeqIO.read(__import__("io").StringIO(text.split("//\n")[11] + "//\n"), "genbank")
    assert one.id == f"{ident}_scf12"
    assert [str(r.seq) for r in recs] == [c.sequence for c in b.contigs]


def test_header(genome):
    b, out, ident = genome
    text = (out / f"{ident}.gbk").read_text()
    for forbidden in ("KEYWORDS", "REFERENCE", "DBLINK", "PGAP"):
        assert forbidden not in text
    first = text.split("//\n")[0]
    assert re.search(r"^LOCUS .* 3000 bp    DNA     circular BCT 02-JAN-2026$", first, re.M)
    assert re.search(rf"^ACCESSION   {ident}_scf1$", first, re.M)
    assert re.search(rf"^VERSION     {ident}_scf1$", first, re.M)
    assert re.search(r"^SOURCE      Streptococcus thermophilus$", first, re.M)
    assert "Bacteria; Bacillati; Bacillota; Bacilli; Lactobacillales;" in first.replace("\n            ", " ")
    recs = read(out, ident)
    assert [r.annotations["topology"] for r in recs[:3]] == ["circular", "circular", "linear"]
    assert recs[0].annotations["taxonomy"][:3] == ["Bacteria", "Bacillati", "Bacillota"]
    assert recs[0].annotations["organism"] == "Streptococcus thermophilus"


def test_definition_drives_arx_classification(genome):
    b, out, ident = genome
    recs = read(out, ident)
    assert "chromosome" in recs[0].description and "plasmid" not in recs[0].description
    assert "plasmid" in recs[1].description
    assert recs[2].description == "Streptococcus thermophilus strain SSSSSSSSSSSSSSSSSSSS scaffold 3, whole genome shotgun sequence"
    kinds = [classify_contig(r.description, r.annotations["topology"]) for r in recs]
    assert kinds == ["chromosome", "plasmid"] + ["scaffold"] * 10


def test_comment_and_structured_comment(genome):
    b, out, ident = genome
    rec = read(out, ident)[0]
    assert rec.annotations["comment"] == f"Made by Abrinca GenBankBuilder {__version__}"
    sc = rec.annotations["structured_comment"]["Genome-Annotation-Data"]
    assert sc["Annotation Pipeline"] == "Abrinca GenBankBuilder"
    assert sc["Annotation Software revision"] == __version__
    assert sc["Annotation Date"] == "01/02/2026 03:04:05"
    assert {k: v for k, v in sc.items() if k not in ("Annotation Pipeline", "Annotation Software revision", "Annotation Date")} == {
        "Genes (total)": "10", "CDSs (total)": "6", "Genes (coding)": "5", "CDSs (with protein)": "5",
        "Genes (RNA)": "4", "rRNAs": "1", "tRNAs": "1", "ncRNAs": "2",
        "Pseudo Genes (total)": "1", "CDSs (without protein)": "1"}
    # same block in every record, as PGAP does
    assert all(r.annotations["structured_comment"] == rec.annotations["structured_comment"] for r in read(out, ident))
    assert (out / f"{ident}.gbk").read_text().count("##Genome-Annotation-Data-START##") == 12


def test_source_feature(genome):
    b, out, ident = genome
    for r in read(out, ident):
        src = r.features[0]
        assert src.type == "source" and (int(src.location.start), int(src.location.end)) == (0, len(r))
        assert src.qualifiers == {"organism": ["Streptococcus thermophilus"], "mol_type": ["genomic DNA"],
                                  "strain": ["S" * 20], "db_xref": ["taxon:1308"]}


def test_feature_order_and_parents(genome):
    b, out, ident = genome
    rec = read(out, ident)[0]
    seq = [(f.type, f.qualifiers.get("locus_tag", [None])[0]) for f in rec.features[1:]]
    lt = lambda n: f"{ident}_{n:06d}"
    assert seq == [
        ("gene", lt(1)), ("CDS", lt(1)),  # wrap-around starts at 0
        ("gene", lt(2)), ("CDS", lt(2)), ("gene", lt(3)), ("CDS", lt(3)), ("gene", lt(4)), ("CDS", lt(4)),
        ("repeat_region", None), ("gene", lt(5)), ("tRNA", lt(5)), ("gene", lt(6)), ("rRNA", lt(6)),
        ("gene", lt(7)), ("ncRNA", lt(7)), ("gene", lt(8)), ("tmRNA", lt(8)), ("regulatory", None)]


def test_gene_parent_locations_and_qualifiers(genome):
    b, out, ident = genome
    feats = {(f.type, f.qualifiers.get("locus_tag", [None])[0]): f for f in read(out, ident)[0].features}
    lt = lambda n: f"{ident}_{n:06d}"
    # wrap-around keeps the join on the gene, slippage joins give the gene the span (PGAP behaviour)
    assert str(feats["gene", lt(1)].location) == str(feats["CDS", lt(1)].location)
    assert len(feats["gene", lt(1)].location.parts) == 2
    assert str(feats["gene", lt(3)].location) == "[500:1100](+)"
    assert str(feats["CDS", lt(3)].location) == "join{[500:800](+), [799:1100](+)}"
    assert feats["CDS", lt(3)].qualifiers["ribosomal_slippage"] == [""]
    # pseudo: /pseudo on gene + CDS, no translation; qualifiers
    assert feats["gene", lt(4)].qualifiers == {"locus_tag": [lt(4)], "pseudo": [""]}
    assert "translation" not in feats["CDS", lt(4)].qualifiers and feats["CDS", lt(4)].qualifiers["pseudo"] == [""]
    cds2 = feats["CDS", lt(2)].qualifiers
    assert cds2["gene"] == ["cas9"] and cds2["inference"] == ["a", "b"] and cds2["transl_table"] == ["11"]
    assert cds2["codon_start"] == ["1"] and cds2["translation"] == ["MKV"]
    assert feats["gene", lt(2)].qualifiers == {"gene": ["cas9"], "locus_tag": [lt(2)]}
    assert "EC_number" not in cds2
    # raw text: transl_table/codon_start unquoted like PGAP, flags without value
    text = (out / f"{ident}.gbk").read_text()
    assert "/transl_table=11\n" in text and "/codon_start=1\n" in text
    assert re.search(r"^ +/pseudo$", text, re.M) and re.search(r"^ +/ribosomal_slippage$", text, re.M)


def test_sidecars(genome):
    b, out, ident = genome
    lt = lambda n: f"{ident}_{n:06d}"
    assert sorted(p.name for p in out.iterdir()) == sorted([f"{ident}.{e}" for e in ("gbk", "KG", "EC")] + ["organism.json", "genome.json"])
    assert (out / f"{ident}.KG").read_text() == f"{lt(2)}\tK00001, K00002\n{lt(10)}\tK00002\n"
    assert (out / f"{ident}.EC").read_text() == f"{lt(1)}\tEC:2.7.7.7\n{lt(4)}\tEC:1.1.1.-\n"
    for p in out.glob(f"{ident}.[A-Z]*"):
        text = p.read_text()
        assert text.endswith("\n") and "\n\n" not in text and not text.startswith("\n")
        assert all(len(line.split("\t")) == 2 for line in text.splitlines())


def test_output_is_deterministic(genome, tmp_path):
    b, out, ident = genome
    b.write(tmp_path / "again")
    for p in out.iterdir():
        assert p.read_text() == (tmp_path / "again" / p.name).read_text()


@pytest.mark.parametrize("note", [
    "word " * 60,                      # wrapped at spaces
    "x" * 200,                         # one long word: never broken
    'with "quotes" and ; semicolons',
    "a  double  space",
])
def test_qualifier_values_survive_wrapping(taxon, tmp_path, note):
    b = GenBankBuilder("G-1", 1308, strain="G", taxon=taxon)
    b.add_contig("G-1_scf1", random_dna(500))
    b.add_gene("G-1_000001", "G-1_scf1", loc(0, 90), "CDS", product="p " * 20, translation="M" * 30, note=note)
    b.write(tmp_path)
    cds = [f for f in SeqIO.read(tmp_path / "G-1.gbk", "genbank").features if f.type == "CDS"][0]
    assert cds.qualifiers["note"] == [note.strip()] or cds.qualifiers["note"] == [note]
    assert cds.qualifiers["translation"] == ["M" * 30]


def test_round_trip_check_catches_unreadable_output(taxon, tmp_path, monkeypatch):
    """With stock Biopython qualifier wrapping a 50-char identifier corrupts /locus_tag on re-read."""
    from Bio.SeqIO.InsdcIO import GenBankWriter
    from gbk_creator import writer
    monkeypatch.setattr(writer._GenBankWriter, "_write_feature_qualifier", GenBankWriter._write_feature_qualifier)
    b = GenBankBuilder(LONG_ID, 1308, strain="S", taxon=taxon)
    b.add_contig(f"{LONG_ID}_scf1", random_dna(500))
    b.add_gene(f"{LONG_ID}_000001", f"{LONG_ID}_scf1", loc(0, 90), "CDS", product="p", translation="MKV")
    with pytest.raises(GbkCreatorError, match="round-trip"):
        b.write(tmp_path / "out")
    assert not (tmp_path / "out").exists()
