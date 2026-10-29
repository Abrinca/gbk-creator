"""Acceptance: arx-tools (the Arx importer) must accept the output without renaming."""
import json
import shutil

import pytest
from Bio import SeqIO

from conftest import CORE_SAMPLES, DATA
from arx_tools.genbank_to_fasta import GenBankToFasta
from arx_tools.import_genome import import_genome
from arx_tools.rename_custom_annotations import CustomAnnotationFile
from arx_tools.rename_genbank import GenBankFile

SAMPLES = CORE_SAMPLES + ["FAM18815-i1-1.1"]  # + a genome with an orphan PGAP gene


def fasta(path):
    return {r.id.split("|")[-1]: str(r.seq) for r in SeqIO.parse(path, "fasta")}


@pytest.mark.parametrize("sample", SAMPLES)
def test_arx_tools_file_checks(sample, build_sample):
    b, out = build_sample(sample)
    gbk = GenBankFile(str(out / f"{sample}.gbk"))
    gbk.validate_contig_ids(genome_id=sample, contig_format="_scf{n}")
    gbk.validate_locus_tags(locus_tag_prefix=f"{sample}_")
    assert gbk.detect_strain_locus_tag_prefix() == (b.strain, f"{sample}_")
    assert gbk.taxid() == b.taxid
    for p in out.glob(f"{sample}.*"):
        if p.suffix != ".gbk":
            CustomAnnotationFile(str(p), custom_annotation_type=p.suffix[1:]).validate_locus_tags(locus_tag_prefix=f"{sample}_")


@pytest.mark.parametrize("sample", SAMPLES)
def test_genbank_to_fasta_strict_matches_oracle(sample, build_sample, tmp_path):
    """strict=True raises if `Genes (total)` / `CDSs (with protein)` disagree with the features."""
    _, out = build_sample(sample)
    gbk = str(out / f"{sample}.gbk")
    assert GenBankToFasta._get_total_genes(gbk) is not None
    assert GenBankToFasta._get_total_proteins(gbk) is not None
    GenBankToFasta.convert(gbk, str(tmp_path / "x.faa"), "faa", strict=True)
    GenBankToFasta.convert(gbk, str(tmp_path / "x.ffn"), "ffn", strict=True)
    src = DATA / sample / "2_cds"
    faa, ffn = fasta(tmp_path / "x.faa"), fasta(tmp_path / "x.ffn")
    assert len(faa) == GenBankToFasta._get_total_proteins(gbk)
    assert len(ffn) == GenBankToFasta._get_total_genes(gbk)
    oracle_faa, oracle_ffn = fasta(next(src.glob("*.faa"))), fasta(next(src.glob("*.ffn")))
    assert faa == oracle_faa
    assert ffn == oracle_ffn


def test_strict_check_is_a_real_check(build_sample, tmp_path):
    sample = CORE_SAMPLES[0]
    _, out = build_sample(sample)
    text = (out / f"{sample}.gbk").read_text()
    n = int(GenBankToFasta._get_total_genes(str(out / f"{sample}.gbk")))
    bad = tmp_path / "bad.gbk"
    bad.write_text(text.replace(f"Genes (total)                :: {n}", f"Genes (total)                :: {n + 1}"))
    with pytest.raises(AssertionError, match="Claims to have"):
        GenBankToFasta.convert(str(bad), str(tmp_path / "bad.ffn"), "ffn", strict=True)


@pytest.mark.parametrize("sample", SAMPLES)
def test_import_genome(sample, build_sample, tmp_path):
    b, out = build_sample(sample)
    folder = tmp_path / "folder_structure"
    (folder / "organisms").mkdir(parents=True)
    (folder / "version.json").write_text(json.dumps({"folder_structure_version": 3}))
    import_dir = tmp_path / "import"
    shutil.copytree(out, import_dir)

    import_genome(import_dir=str(import_dir), folder_structure_dir=str(folder), organism=b.strain, genome=sample)

    genome_dir = folder / "organisms" / b.strain / "genomes" / sample
    genome = json.loads((genome_dir / "genome.json").read_text())
    assert genome["identifier"] == sample
    assert genome["cds_tool_gbk_file"] == f"{sample}.gbk"
    for key in ("assembly_fasta_file", "cds_tool_faa_file", "cds_tool_ffn_file", "cds_tool_gff_file"):
        assert (genome_dir / genome[key]).is_file(), key
    assert {ca["type"] for ca in genome["custom_annotations"]} == {p.suffix[1:] for p in out.glob(f"{sample}.*") if p.suffix != ".gbk"}
    assert json.loads((folder / "organisms" / b.strain / "organism.json").read_text())["taxid"] == b.taxid
    # nothing was renamed: the gbk is byte-identical
    assert (genome_dir / f"{sample}.gbk").read_text() == (out / f"{sample}.gbk").read_text()
    # the fna arx-tools derives from our gbk is the input assembly
    fna = fasta(genome_dir / genome["assembly_fasta_file"])
    assert fna == {c.name: c.sequence for c in b.contigs}


# ----------------------------------------------------------------------------- new annotation category
def _new_category_genome(tmp_path):
    from Bio.SeqFeature import SimpleLocation
    from conftest import BUILD_DATE, TAXON, random_dna
    from gbk_creator import GenBankBuilder

    cats = tmp_path / "categories.json"
    cats.write_text(json.dumps({"KG": "^K[0-9]{5}$", "ZZ": "^ZZ[0-9]{3}$"}))
    b = GenBankBuilder("NEW-1.1", 1308, strain="NEW", taxon=TAXON, build_date=BUILD_DATE, annotation_regex=cats)
    b.add_contig("NEW-1.1_scf1", random_dna(600))
    for n, start in enumerate((0, 200), 1):
        b.add_gene(f"NEW-1.1_{n:06d}", "NEW-1.1_scf1", SimpleLocation(start, start + 90, 1), "CDS",
                   product="p", translation="M" * 30)
    b.add_annotation("NEW-1.1_000001", "KG", "K00001")
    b.add_annotation("NEW-1.1_000001", "ZZ", "ZZ001")
    b.add_annotation("NEW-1.1_000002", "ZZ", "ZZ002")
    out = tmp_path / "gen"
    b.write(out)
    assert (out / "NEW-1.1.ZZ").read_text() == "NEW-1.1_000001\tZZ001\nNEW-1.1_000002\tZZ002\n"
    return out


def _import(out, tmp_path, name, import_settings=None):
    folder = tmp_path / name
    (folder / "organisms").mkdir(parents=True)
    (folder / "version.json").write_text(json.dumps({"folder_structure_version": 3}))
    import_dir = tmp_path / f"{name}-in"
    shutil.copytree(out, import_dir)
    import_genome(import_dir=str(import_dir), folder_structure_dir=str(folder), organism="NEW", genome="NEW-1.1",
                  import_settings=import_settings)
    genome_dir = folder / "organisms" / "NEW" / "genomes" / "NEW-1.1"
    return genome_dir, json.loads((genome_dir / "genome.json").read_text())


def test_new_two_letter_category_needs_custom_import_settings(tmp_path):
    """arx-tools only globs a fixed list of types; a new category must be added to ARX_IMPORT_SETTINGS."""
    from arx_tools.import_genome import ImportSettings

    out = _new_category_genome(tmp_path)

    # 1. default settings: the file is copied along but NOT registered as a custom annotation
    genome_dir, genome = _import(out, tmp_path, "default")
    assert (genome_dir / "NEW-1.1.ZZ").is_file()
    assert {ca["type"] for ca in genome["custom_annotations"]} == {"KG"}

    # 2. documented extra file_finder entry: discovered, with the right type, and locus tags validated
    custom = [dict(e) for e in ImportSettings.default_settings["file_finder"]["custom_annotations"]]
    custom.append({"glob": "*.ZZ", "anno_type": "ZZ", "expected": False})
    settings = {"file_finder": {"custom_annotations": custom}}
    genome_dir, genome = _import(out, tmp_path, "custom", import_settings=settings)
    assert {(ca["type"], ca["file"]) for ca in genome["custom_annotations"]} == {
        ("KG", "NEW-1.1.KG"), ("ZZ", "NEW-1.1.ZZ")}


def test_category_names_must_be_two_letters(tmp_path):
    from gbk_creator import AnnotationError, load_categories
    for bad in ("Z", "ZZZ", "Z1", "zz", "K G", "1K"):
        p = tmp_path / "c.json"
        p.write_text(json.dumps({bad: "^x$"}))
        with pytest.raises(AnnotationError, match="exactly two upper-case letters"):
            load_categories(p)


# ----------------------------------------------------------------------------- organism.json
def _organism_genome(tmp_path, name, **organism_kwargs):
    from Bio.SeqFeature import SimpleLocation
    from conftest import BUILD_DATE, TAXON, random_dna
    from gbk_creator import GenBankBuilder

    b = GenBankBuilder("NEW-1.1", 1308, strain="NEW", taxon=TAXON, build_date=BUILD_DATE, **organism_kwargs)
    b.add_contig("NEW-1.1_scf1", random_dna(600))
    b.add_gene("NEW-1.1_000001", "NEW-1.1_scf1", SimpleLocation(0, 90, 1), "CDS", product="p", translation="M" * 30)
    out = tmp_path / name
    b.write(out)
    return out


def _import_organism(out, tmp_path, name, existing_organism_json=None):
    folder = tmp_path / f"fs-{name}"
    (folder / "organisms" / "NEW").mkdir(parents=True)
    (folder / "version.json").write_text(json.dumps({"folder_structure_version": 3}))
    if existing_organism_json is not None:
        (folder / "organisms" / "NEW" / "organism.json").write_text(json.dumps(existing_organism_json))
    import_dir = tmp_path / f"in-{name}"
    shutil.copytree(out, import_dir)
    import_genome(import_dir=str(import_dir), folder_structure_dir=str(folder), organism="NEW", genome="NEW-1.1")
    return json.loads((folder / "organisms" / "NEW" / "organism.json").read_text())


@pytest.mark.parametrize("kwargs,expected", [
    (dict(restricted=True, tags=["probiotic", "dairy"]), dict(restricted=True, tags=["probiotic", "dairy"])),
    (dict(restricted=True), dict(restricted=True, tags=[])),
    ({}, dict(restricted=False, tags=[])),
])
def test_import_genome_picks_up_our_organism_json(tmp_path, kwargs, expected):
    from arx_tools.metadata_schemas import organism_json_schema

    organism = _import_organism(_organism_genome(tmp_path, "gen", **kwargs), tmp_path, "x")
    assert organism == {"name": "NEW", "alternative_name": None, "taxid": 1308, "representative": "NEW-1.1", **expected}  # arx-tools fills alternative_name itself
    organism_json_schema.validate(organism)  # arx-tools' own schema


def test_our_organism_json_overrides_an_existing_one_where_it_has_values(tmp_path):
    from arx_tools.metadata_schemas import organism_json_schema

    existing = {"name": "NEW", "alternative_name": "kept", "taxid": 999, "restricted": True, "tags": ["old"],
                "representative": "NEW-0.9"}
    # values we set win (taxid, restricted - also False -, tags); keys we never write are untouched
    out = _organism_genome(tmp_path, "g1", restricted=False, tags=["new"])
    organism = _import_organism(out, tmp_path, "a", existing)
    assert organism == {"name": "NEW", "alternative_name": "kept", "taxid": 1308, "restricted": False,
                        "tags": ["new"], "representative": "NEW-1.1"}
    organism_json_schema.validate(organism)
    # omitted tags / empty list never erase existing values; restricted=False does overwrite
    out = _organism_genome(tmp_path, "g2")
    organism = _import_organism(out, tmp_path, "b", existing)
    assert organism["alternative_name"] == "kept" and organism["tags"] == ["old"]
    assert organism["restricted"] is False and organism["taxid"] == 1308
    out = _organism_genome(tmp_path, "g3", tags=[])
    assert _import_organism(out, tmp_path, "c", existing)["tags"] == ["old"]


# ----------------------------------------------------------------------------- genome.json
def _import_genome_json(out, tmp_path, name):
    folder = tmp_path / f"fs-{name}"
    (folder / "organisms").mkdir(parents=True)
    (folder / "version.json").write_text(json.dumps({"folder_structure_version": 3}))
    import_dir = tmp_path / f"in-{name}"
    shutil.copytree(out, import_dir)
    import_genome(import_dir=str(import_dir), folder_structure_dir=str(folder), organism="NEW", genome="NEW-1.1")
    return json.loads((folder / "organisms" / "NEW" / "genomes" / "NEW-1.1" / "genome.json").read_text())


@pytest.mark.parametrize("kwargs,old", [({}, None), (dict(old_identifier="NEW-legacy"), "NEW-legacy")])
def test_import_genome_picks_up_our_genome_json(tmp_path, kwargs, old):
    from arx_tools.metadata_schemas import genome_json_schema
    from gbk_creator import __version__

    out = _organism_genome(tmp_path, "gen", **kwargs)
    assert sorted(json.loads((out / "genome.json").read_text())) == (["cds_tool"] if old is None else ["cds_tool", "old_identifier"])
    genome = _import_genome_json(out, tmp_path, "x")
    genome_json_schema.validate(genome)  # arx-tools' own schema
    assert genome["cds_tool"] == f"Abrinca GenBankBuilder v{__version__}"
    assert genome["old_identifier"] == old
    assert genome["identifier"] == "NEW-1.1" and "contig_format" not in genome
    assert genome["cds_tool_gbk_file"] == "NEW-1.1.gbk"  # arx-tools' own file bookkeeping is intact


def test_genome_json_does_not_interfere_with_contig_format(tmp_path):
    """arx-tools reads `contig_format` from import_dir/genome.json: we must not write it (default `_scf{n}`)."""
    out = _organism_genome(tmp_path, "gen", old_identifier="x")
    assert "contig_format" not in json.loads((out / "genome.json").read_text())
    _import_genome_json(out, tmp_path, "x")  # check_files validates the `_scf{n}` contig names


def test_our_cds_tool_wins_the_merge_with_arx_tools_detection(tmp_path, monkeypatch):
    """gather_metadata: gbk-derived cds_tool* first, import_dir/genome.json merged afterwards (ours wins)."""
    from arx_tools.rename_genbank import GenBankFile

    out = _organism_genome(tmp_path, "gen")
    # today arx-tools does not recognise our COMMENT: it sets no cds_tool and would only warn
    _, genome_gbk = GenBankFile(str(out / "NEW-1.1.gbk")).metadata()
    assert "cds_tool" not in genome_gbk
    # ...and if a future arx-tools did detect something, our genome.json is merged later and still wins
    detected = {"cds_tool": "SomethingElse", "cds_tool_version": "9", "cds_tool_date": "2020-01-01"}
    monkeypatch.setattr(GenBankFile, "metadata", lambda self: ({"taxid": 1308}, dict(detected)))
    genome = _import_genome_json(out, tmp_path, "x")
    assert genome["cds_tool"].startswith("Abrinca GenBankBuilder v")
    assert genome["cds_tool_version"] == "9" and genome["cds_tool_date"] == "2020-01-01"  # not ours: untouched


# ----------------------------------------------------------------------------- old_locus_tag
def test_old_locus_tag_survives_import_and_generated_files(tmp_path):
    from Bio.SeqFeature import SimpleLocation
    from conftest import BUILD_DATE, TAXON, random_dna
    from gbk_creator import GenBankBuilder

    b = GenBankBuilder("NEW-1.1", 1308, strain="NEW", taxon=TAXON, build_date=BUILD_DATE)
    b.add_contig("NEW-1.1_scf1", random_dna(600))
    old = {1: "old_gene_A", 2: "LEGACY 0002", 3: "x" * 100}
    b.add_gene("NEW-1.1_000001", "NEW-1.1_scf1", SimpleLocation(0, 90, 1), "CDS", product="p1", translation="M" * 30,
               old_locus_tag=old[1])
    b.add_gene("NEW-1.1_000002", "NEW-1.1_scf1", SimpleLocation(100, 190, -1), "CDS", product="p2", pseudo=True,
               old_locus_tag=old[2])
    b.add_gene("NEW-1.1_000003", "NEW-1.1_scf1", SimpleLocation(200, 276, 1), "tRNA", product="tRNA-Gly",
               old_locus_tag=old[3])
    b.add_gene("NEW-1.1_000004", "NEW-1.1_scf1", SimpleLocation(300, 390, 1), "CDS", product="p4", translation="M" * 30)
    out = tmp_path / "gen"
    b.write(out)

    folder = tmp_path / "fs"
    (folder / "organisms").mkdir(parents=True)
    (folder / "version.json").write_text(json.dumps({"folder_structure_version": 3}))
    import_genome(import_dir=str(out), folder_structure_dir=str(folder), organism="NEW", genome="NEW-1.1")
    genome_dir = folder / "organisms" / "NEW" / "genomes" / "NEW-1.1"

    # the gbk is untouched and parses: both features of every locus carry the first /old_locus_tag Arx reads
    assert (genome_dir / "NEW-1.1.gbk").read_text() == (out / "NEW-1.1.gbk").read_text()
    rec = next(SeqIO.parse(genome_dir / "NEW-1.1.gbk", "genbank"))
    seen = {}
    for f in rec.features:
        if "locus_tag" in f.qualifiers:
            seen.setdefault(f.qualifiers["locus_tag"][0], []).append(f.qualifiers.get("old_locus_tag"))
    assert seen == {"NEW-1.1_000001": [[old[1]]] * 2, "NEW-1.1_000002": [[old[2]]] * 2,
                    "NEW-1.1_000003": [[old[3]]] * 2, "NEW-1.1_000004": [None] * 2}
    # arx-tools' generated faa/ffn (strict) are unaffected and keyed by the new locus tags
    assert set(fasta(genome_dir / "NEW-1.1.faa")) == {"NEW-1.1_000001", "NEW-1.1_000004"}
    assert set(fasta(genome_dir / "NEW-1.1.ffn")) == {f"NEW-1.1_00000{n}" for n in (1, 2, 3, 4)}
    assert not any(o in (genome_dir / "NEW-1.1.faa").read_text() for o in old.values())
    # the gff derived from the gbk carries it too
    gff = (genome_dir / "NEW-1.1.arx.gff").read_text()
    assert "old_locus_tag=old_gene_A" in gff and "old_locus_tag=LEGACY 0002" in gff
    # Arx reads the first /old_locus_tag of a locus: the `gene` feature comes first and has it
    gbk = GenBankFile(str(genome_dir / "NEW-1.1.gbk"))
    gbk.validate_locus_tags(locus_tag_prefix="NEW-1.1_")
