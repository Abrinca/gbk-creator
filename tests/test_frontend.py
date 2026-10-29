"""The .arx.gff test frontend (offline, synthetic input)."""
import logging

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation, SimpleLocation

from conftest import ID, TAXON, random_dna
from gbk_creator import GbkCreatorError, main
from gbk_creator.frontends.arx_gff import (Dropped, FrontendError, build_from_arx_gff, parse_attributes, read_gff, reconstruct_join)
from gbk_creator.taxonomy import NcbiClient

SCF1, SCF2 = f"{ID}_scf1", f"{ID}_scf2"


def gff_line(seqid, type_, start, end, strand, attrs):
    return "\t".join([seqid, "feature", type_, str(start), str(end), ".", strand, ".", attrs])


def emapper_line(name, **cols):
    c = [""] * 22
    c[0] = name
    for k, v in cols.items():
        c[int(k[1:])] = v
    return "\t".join(c)


@pytest.fixture
def files(tmp_path):
    seq1, seq2 = random_dna(1000, 1), random_dna(500, 2)
    (tmp_path / "a.fna").write_text(f">{SCF1}\n{seq1}\n>{SCF2}\n{seq2}\n")
    lt = lambda n: f"{ID}_{n:06d}"
    rows = [
        f"{SCF1}\tannotation\tremark\t1\t1000\t.\t.\t.\ttopology=linear;accessions={SCF1}",
        f"{SCF2}\tannotation\tremark\t1\t500\t.\t.\t.\ttopology=circular",
        gff_line(SCF1, "source", 1, 1000, "+", "db_xref=taxon:1308;mol_type=genomic DNA;organism=x;strain=FAM13496"),
        gff_line(SCF1, "gene", "<1", 90, "+", f"ID={lt(1)}_gene;locus_tag={lt(1)}"),
        gff_line(SCF1, "CDS", "<1", 90, "+",
                 f"ID={lt(1)};Name=p;Parent={lt(1)}_gene;codon_start=2;EC_number=1.1.1.1,2.2.2.-;gene=abcd;"
                 f"inference=a%2C b,c;locus_tag={lt(1)};note=x%3B y;product=prod;protein_id=extdb:{lt(1)};"
                 f"transl_table=11;translation=MKV"),
        gff_line(SCF1, "gene", 100, ">200", "-", f"ID={lt(2)}_gene;locus_tag={lt(2)};pseudo="),
        gff_line(SCF1, "CDS", 100, ">200", "-", f"ID={lt(2)};locus_tag={lt(2)};product=pp;pseudo=;transl_table=11"),
        gff_line(SCF1, "tRNA", 300, 376, "+", f"ID={lt(3)};locus_tag={lt(3)};product=tRNA-Gly;anticodon=(pos:320..322%2Caa:Gly%2Cseq:tcc)"),
        gff_line(SCF1, "repeat_region", 400, 500, "+", "ID=x;rpt_family=CRISPR;rpt_type=direct;rpt_unit_range=1..2;rpt_unit_seq=ac"),
        gff_line(SCF2, "gene", 10, 40, "+", f"ID={lt(4)}_gene;locus_tag={lt(4)}"),
        gff_line(SCF2, "CDS", 10, 40, "+", f"ID={lt(4)};locus_tag={lt(4)};product=q;transl_table=11;translation=MM"),
    ]
    (tmp_path / "a.arx.gff").write_text("##gff-version 3\n" + "\n".join(rows) + "\n")
    emapper = [
        "#query\tseed...",
        emapper_line(f"gnl|extdb|{lt(1)}", c5="dnaA", c6="GO:0000001,GO:0000002", c7="1.1.1.1", c8="ko:K00001,ko:K00002",
                     c11="R00001", c18="COG0592@1|root,COG0592@2|Bacteria", c21="Some description, with comma; more"),
        emapper_line(f"gnl|extdb|{lt(4)}", c8="ko:K00003"),
        "## 2 queries",
    ]
    (tmp_path / "a.emapper.annotations").write_text("\n".join(emapper) + "\n")
    (tmp_path / "a.faa").write_text(f">gnl|extdb|{lt(1)} prod [x]\nMKV\n>gnl|extdb|{lt(4)} q [x]\nMM\n")
    return tmp_path


def build(files, **kw):
    kw.setdefault("taxon", TAXON)
    return build_from_arx_gff(files / "a.arx.gff", files / "a.fna", **kw)


def test_parse_attributes():
    a = parse_attributes("ID=x;note=a%3B b%2C c,d;pseudo=;anticodon=(pos:1..3%2Caa:G)")
    assert a == {"ID": ["x"], "note": ["a; b, c", "d"], "pseudo": [""], "anticodon": ["(pos:1..3,aa:G)"]}


def test_coordinates_and_partials(files):
    rows = read_gff(files / "a.arx.gff")
    cds = [r for r in rows if r.type == "CDS"]
    assert (cds[0].start, cds[0].end, cds[0].start_partial, cds[0].end_partial) == (0, 90, True, False)
    assert (cds[1].start, cds[1].end, cds[1].start_partial, cds[1].end_partial) == (99, 200, False, True)


def test_build_from_synthetic_gff(files, tmp_path):
    b = build(files, faa=files / "a.faa", eggnog=files / "a.emapper.annotations", build_date=None)
    assert (b.identifier, b.strain, b.taxid) == (ID, "FAM13496", 1308)  # all derived from the GFF
    assert [(c.name, c.topology, c.kind) for c in b.contigs] == [(SCF1, "linear", "scaffold"), (SCF2, "circular", "chromosome")]
    out = tmp_path / "out"
    b.write(out)
    rec = next(SeqIO.parse(out / f"{ID}.gbk", "genbank"))
    cds = [f for f in rec.features if f.type == "CDS"]
    assert str(cds[0].location) == "[<0:90](+)"
    assert str(cds[1].location) == "[99:>200](-)" and "pseudo" in cds[1].qualifiers and "translation" not in cds[1].qualifiers
    q = cds[0].qualifiers
    assert q["codon_start"] == ["2"] and q["inference"] == ["a, b", "c"] and q["gene"] == ["abcd"]
    assert q["protein_id"] == [f"{ID}_000001".join(["extdb:", ""])] and "EC_number" not in q
    assert [f.type for f in rec.features].count("gene") == 3
    lt = lambda n: f"{ID}_{n:06d}"
    assert (out / f"{ID}.EC").read_text() == f"{lt(1)}\tEC:1.1.1.1, EC:2.2.2.-\n"
    assert (out / f"{ID}.KG").read_text() == f"{lt(1)}\tK00001, K00002\n{lt(4)}\tK00003\n"
    assert (out / f"{ID}.GO").read_text() == f"{lt(1)}\tGO:0000001, GO:0000002\n"
    assert (out / f"{ID}.KR").read_text() == f"{lt(1)}\tR00001\n"
    assert (out / f"{ID}.EP").read_text() == f"{lt(1)}\tEP:dnaA\n"
    assert (out / f"{ID}.EO").read_text() == f"{lt(1)}\tEO:COG0592\n"
    assert (out / f"{ID}.ED").read_text() == f"{lt(1)}\tED:Some description\n"


def test_gff_translation_used_without_faa_and_must_agree_with_faa(files):
    build(files)  # translation attribute of the GFF
    (files / "bad.faa").write_text(f">{ID}_000001 x\nMKW\n>{ID}_000004 q\nMM\n")
    with pytest.raises(FrontendError, match="differ"):
        build(files, faa=files / "bad.faa")


def test_ec_from_gff_and_eggnog_are_merged_without_duplicates(files, tmp_path):
    (files / "a.emapper.annotations").write_text(emapper_line(f"gnl|extdb|{ID}_000001", c7="1.1.1.1,3.3.3.3") + "\n")
    b = build(files, eggnog=files / "a.emapper.annotations")
    b.write(tmp_path)
    assert (tmp_path / f"{ID}.EC").read_text() == f"{ID}_000001\tEC:1.1.1.1, EC:2.2.2.-, EC:3.3.3.3\n"


def test_invalid_annotation_ids_raise_by_default_and_are_dropped_when_lenient(files, caplog):
    gff = files / "a.arx.gff"
    gff.write_text(gff.read_text().replace("EC_number=1.1.1.1,2.2.2.-", "EC_number=3.6.5.n1"))
    with pytest.raises(FrontendError, match=r"EC annotation 'EC:3.6.5.n1' of .*_000001 does not match.*lenient"):
        build(files)
    dropped = []
    with caplog.at_level(logging.WARNING):
        b = build(files, lenient=True, dropped=dropped)
    assert dropped == [Dropped("invalid_annotation", f"{ID}_000001", "EC:3.6.5.n1")]
    assert "dropping EC annotation 'EC:3.6.5.n1'" in caplog.text
    assert "EC" not in b._annotations


def test_unknown_things_fail_loudly_even_when_lenient(files):
    gff = files / "a.arx.gff"
    original = gff.read_text()
    gff.write_text(original + gff_line(SCF1, "mobile_element", 600, 700, "+", "ID=m") + "\n")
    for lenient in (False, True):
        with pytest.raises(FrontendError, match="unsupported feature type"):
            build(files, lenient=lenient)
    gff.write_text(original.replace("ID=x;rpt_family", "ID=x;weird_qualifier=1;rpt_family"))
    with pytest.raises(GbkCreatorError, match="weird_qualifier"):
        build(files, lenient=True)


def test_leftover_faa_and_ffn_records_raise_by_default(files):
    (files / "left.faa").write_text((files / "a.faa").read_text() + f">gnl|extdb|{ID}_000077 stray\nMKV\n")
    with pytest.raises(FrontendError, match=r"_000077: faa record is not used.*lenient"):
        build(files, faa=files / "left.faa")
    dropped = []
    build(files, faa=files / "left.faa", lenient=True, dropped=dropped)
    assert dropped == [Dropped("unused_faa", f"{ID}_000077", "")]

    (files / "left.ffn").write_text(f">{ID}_000078 stray\nACGT\n")
    with pytest.raises(FrontendError, match=r"_000078: ffn record is not used"):
        build(files, ffn=files / "left.ffn")
    dropped = []
    build(files, ffn=files / "left.ffn", lenient=True, dropped=dropped)
    assert dropped == [Dropped("unused_ffn", f"{ID}_000078", "")]


def test_translation_on_non_cds_is_always_rejected(files):
    lt3 = f"{ID}_000003"  # the tRNA
    (files / "rna.faa").write_text((files / "a.faa").read_text() + f">gnl|extdb|{lt3} x\nMKV\n")
    for lenient in (False, True):
        with pytest.raises(FrontendError, match="translation is given for a tRNA"):
            build(files, faa=files / "rna.faa", lenient=lenient)
    gff = files / "a.arx.gff"
    gff.write_text(gff.read_text().replace(f"locus_tag={lt3};product=tRNA-Gly", f"locus_tag={lt3};product=tRNA-Gly;translation=MKV"))
    for lenient in (False, True):
        with pytest.raises(FrontendError, match="translation is given for a tRNA"):
            build(files, lenient=lenient)


def test_other_inconsistent_inputs_fail_loudly(files):
    gff = files / "a.arx.gff"
    original = gff.read_text()
    gff.write_text(original + gff_line("nope", "tRNA", 1, 40, "+", f"locus_tag={ID}_000009;product=t") + "\n")
    with pytest.raises(FrontendError, match="not in"):
        build(files)
    gff.write_text(original)
    (files / "a.emapper.annotations").write_text(emapper_line(f"gnl|extdb|{ID}_000099", c8="ko:K00001") + "\n")
    with pytest.raises(GbkCreatorError, match="unknown locus tag"):
        build(files, eggnog=files / "a.emapper.annotations", lenient=True)
    (files / "a.emapper.annotations").write_text("a\tb\n")
    with pytest.raises(FrontendError, match="22 columns"):
        build(files, eggnog=files / "a.emapper.annotations")


def test_eggnog_on_rna_is_rejected(files):
    (files / "a.emapper.annotations").write_text(emapper_line(f"gnl|extdb|{ID}_000003", c8="ko:K00001") + "\n")
    with pytest.raises(GbkCreatorError, match="only allowed on CDS"):
        build(files, eggnog=files / "a.emapper.annotations")


def test_orphan_gene_rows_raise_by_default_and_are_dropped_when_lenient(files, caplog):
    gff = files / "a.arx.gff"
    gff.write_text(gff.read_text() + gff_line(SCF1, "gene", 800, 850, "+", f"ID={ID}_000008_gene;locus_tag={ID}_000008") + "\n")
    with pytest.raises(FrontendError, match=r"_000008: gene feature without CDS/RNA feature.*lenient"):
        build(files)
    dropped = []
    with caplog.at_level(logging.WARNING):
        b = build(files, lenient=True, dropped=dropped)
    assert "gene feature without CDS/RNA" in caplog.text
    assert dropped == [Dropped("orphan_gene", f"{ID}_000008", "")]
    assert f"{ID}_000008" not in b._genes


def test_identifier_strain_taxid_overrides_and_plasmids(files):
    b = build(files, identifier=ID, strain=ID, plasmids=["2"])
    assert b.strain == ID and b.contigs[1].kind == "plasmid"
    with pytest.raises(GbkCreatorError, match="must start with the strain"):
        build(files, strain="OTHER")
    (files / "a.arx.gff").write_text((files / "a.arx.gff").read_text().replace("db_xref=taxon:1308;", ""))
    with pytest.raises(FrontendError, match="taxid"):
        build(files)
    assert build(files, taxid=1308).taxid == 1308


# ----------------------------------------------------------------------------- join reconstruction
CONTIG = random_dna(400, 7)


def rc(s):
    return str(Seq(s).reverse_complement())


def test_reconstruct_plain_span_unchanged():
    loc = SimpleLocation(10, 100, 1)
    assert reconstruct_join(loc, CONTIG[10:100], CONTIG, False) is loc
    loc = SimpleLocation(10, 100, -1)
    assert reconstruct_join(loc, rc(CONTIG[10:100]), CONTIG, False) is loc


@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("gap", [-1, 1, 2, -2])
def test_reconstruct_slippage(strand, gap):
    a, m, b = 50, 150, 240
    plus = CONTIG[a:m] + CONTIG[m + gap:b]
    nt = plus if strand == 1 else rc(plus)
    loc = reconstruct_join(SimpleLocation(a, b, strand), nt, CONTIG, False)
    assert isinstance(loc, CompoundLocation)
    assert str(loc.extract(Seq(CONTIG))) == nt  # sequence reproduced (junction may shift inside homopolymers)
    assert (int(loc.start), int(loc.end), loc.strand) == (a, b, strand)
    starts = [int(p.start) for p in loc.parts]
    assert starts == sorted(starts, reverse=strand == -1)


@pytest.mark.parametrize("strand", [1, -1])
def test_reconstruct_origin_spanning(strand):
    plus = CONTIG[-60:] + CONTIG[:90]
    nt = plus if strand == 1 else rc(plus)
    loc = reconstruct_join(SimpleLocation(0, len(CONTIG), strand), nt, CONTIG, True)
    assert str(loc.extract(Seq(CONTIG))) == nt
    expected = [(len(CONTIG) - 60, len(CONTIG)), (0, 90)]
    got = [(int(p.start), int(p.end)) for p in loc.parts]
    assert got == (expected if strand == 1 else expected[::-1])
    assert reconstruct_join(SimpleLocation(0, len(CONTIG), strand), nt, CONTIG, False) is None  # linear: no wrap


def test_reconstruct_keeps_partial_markers_and_gives_up_on_garbage():
    plus = CONTIG[50:150] + CONTIG[149:240]
    loc = reconstruct_join(SimpleLocation(BeforePosition(50), AfterPosition(240), 1), plus, CONTIG, False)
    assert isinstance(loc.parts[0].start, BeforePosition) and isinstance(loc.parts[-1].end, AfterPosition)
    assert reconstruct_join(SimpleLocation(50, 240, 1), "ACGT" * 40, CONTIG, False) is None


def test_mismatching_ffn_fails_loudly(files):
    (files / "a.ffn").write_text(f">{ID}_000001 prod\n{'ACGT' * 5}\n")
    with pytest.raises(FrontendError, match="does not reproduce the .ffn"):
        build(files, ffn=files / "a.ffn")


# ----------------------------------------------------------------------------- CLI
def test_cli_lenient_flag(files, tmp_path, monkeypatch, capsys, caplog):
    monkeypatch.setattr(NcbiClient, "fetch_taxon", lambda self, taxid: TAXON)
    monkeypatch.setenv("NCBI_EMAIL", "me@example.org")
    gff = files / "a.arx.gff"
    gff.write_text(gff.read_text().replace("EC_number=1.1.1.1,2.2.2.-", "EC_number=3.6.5.n1"))
    args = ["--gff", str(gff), "--fna", str(files / "a.fna"), "--out", str(tmp_path / "o")]
    assert main(args) == 1
    assert "--lenient" in capsys.readouterr().err
    assert main(args + ["--lenient"]) == 0
    assert (tmp_path / "o" / f"{ID}.gbk").is_file()


def test_cli(files, tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(NcbiClient, "fetch_taxon", lambda self, taxid: TAXON)
    monkeypatch.setenv("NCBI_EMAIL", "me@example.org")
    out = tmp_path / "cli-out"
    rc_ = main(["--gff", str(files / "a.arx.gff"), "--fna", str(files / "a.fna"), "--faa", str(files / "a.faa"),
                "--eggnog", str(files / "a.emapper.annotations"), "--out", str(out)])
    assert rc_ == 0
    printed = capsys.readouterr().out.split()
    assert str(out / f"{ID}.gbk") in printed and str(out / f"{ID}.KG") in printed
    assert (out / f"{ID}.gbk").is_file()


def test_cli_reports_errors_without_traceback(files, tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(NcbiClient, "fetch_taxon", lambda self, taxid: TAXON)
    monkeypatch.setenv("NCBI_EMAIL", "me@example.org")
    rc_ = main(["--gff", str(files / "a.arx.gff"), "--fna", str(files / "a.fna"), "--strain", "NOPE", "--out", str(tmp_path / "o")])
    assert rc_ == 1
    assert "must start with the strain" in capsys.readouterr().err
    assert not (tmp_path / "o").exists()
