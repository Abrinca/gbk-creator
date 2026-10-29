"""organism.json / genome.json written next to the gbk (offline)."""
import json

import pytest
from Bio.SeqFeature import SimpleLocation

from conftest import ID, TAXON, random_dna
from gbk_creator import GbkCreatorError, GenBankBuilder, __version__, main
from gbk_creator.taxonomy import NcbiClient


def make(**kw):
    b = GenBankBuilder(ID, 1308, strain="FAM13496", taxon=TAXON, **kw)
    b.add_contig(f"{ID}_scf1", random_dna(300))
    b.add_gene(f"{ID}_000001", f"{ID}_scf1", SimpleLocation(0, 90, 1), "CDS", product="p", translation="M" * 30)
    return b


def written(tmp_path, **kw):
    make(**kw).write(tmp_path)
    return (tmp_path / "organism.json").read_text()


def test_defaults(tmp_path):
    text = written(tmp_path)
    assert json.loads(text) == {"name": "FAM13496", "taxid": 1308, "restricted": False}
    assert text == '{\n    "name": "FAM13496",\n    "taxid": 1308,\n    "restricted": false\n}\n'  # key order, indent=4, trailing newline; no "tags" key


def test_all_fields(tmp_path):
    text = written(tmp_path, restricted=True, tags=["probiotic", "dairy ä"])
    data = json.loads(text)
    assert list(data) == ["name", "taxid", "restricted", "tags"]
    assert data == {"name": "FAM13496", "taxid": 1308, "restricted": True, "tags": ["probiotic", "dairy ä"]}
    assert text.endswith("\n") and not text.endswith("\n\n") and '\n    "tags": [\n        "probiotic",' in text


def test_empty_tags_list_is_written_but_none_is_omitted(tmp_path):
    assert json.loads(written(tmp_path / "a", tags=[]))["tags"] == []
    assert "tags" not in json.loads(written(tmp_path / "b", tags=None))


def test_name_is_the_strain_and_representative_is_never_written(tmp_path):
    for strain in ("FAM13496", ID):
        b = GenBankBuilder(ID, 1308, strain=strain, taxon=TAXON, tags=["x"])
        b.add_contig(f"{ID}_scf1", random_dna(300))
        b.add_gene(f"{ID}_000001", f"{ID}_scf1", SimpleLocation(0, 90, 1), "CDS", product="p", translation="MM")
        b.write(tmp_path / strain)
        data = json.loads((tmp_path / strain / "organism.json").read_text())
        assert data["name"] == strain and "representative" not in data


def test_write_returns_the_json_files_and_writes_nothing_when_invalid(tmp_path):
    b = make()
    assert [p.name for p in b.write(tmp_path)] == [f"{ID}.gbk", "organism.json", "genome.json"]
    broken = make()
    broken.add_gene(f"{ID}_000000", f"{ID}_scf1", SimpleLocation(100, 190, 1), "CDS", product="p", translation="MM")
    with pytest.raises(GbkCreatorError):
        broken.write(tmp_path / "x")
    assert not (tmp_path / "x").exists()


@pytest.mark.parametrize("bad", ["False", "", 0, 1, None, "true"])
def test_restricted_must_be_a_real_bool(bad):
    with pytest.raises(GbkCreatorError, match="restricted must be True or False"):
        make(restricted=bad)


@pytest.mark.parametrize("bad,match", [
    ("probiotic", "list of strings"), (5, "list of strings"), ({"a"}, "list of strings"),
    ([""], "non-empty"), (["  "], "non-empty"), ([1], "non-empty"), (["a", None], "non-empty"), (["a\nb"], "single-line"),
    (["a", "b", "a"], "duplicates"),
])
def test_tags_validation(bad, match):
    with pytest.raises(GbkCreatorError, match=match):
        make(tags=bad)


def test_tuple_of_tags_is_accepted_and_copied(tmp_path):
    tags = ("a", "b")
    assert json.loads(written(tmp_path, tags=tags))["tags"] == ["a", "b"]


# ----------------------------------------------------------------------------- CLI
def cli_files(tmp_path):
    seq = random_dna(300)
    scf = f"{ID}_scf1"
    (tmp_path / "a.fna").write_text(f">{scf}\n{seq}\n")
    (tmp_path / "a.arx.gff").write_text(
        f"{scf}\tannotation\tremark\t1\t300\t.\t.\t.\ttopology=linear\n"
        f"{scf}\tfeature\tsource\t1\t300\t.\t+\t.\tdb_xref=taxon:1308;strain=FAM13496\n"
        f"{scf}\tfeature\tgene\t1\t90\t.\t+\t.\tlocus_tag={ID}_000001\n"
        f"{scf}\tfeature\tCDS\t1\t90\t.\t+\t0\tlocus_tag={ID}_000001;product=p;transl_table=11;translation=MM\n")
    return ["--gff", str(tmp_path / "a.arx.gff"), "--fna", str(tmp_path / "a.fna")]


def test_cli_options(tmp_path, monkeypatch):
    monkeypatch.setattr(NcbiClient, "fetch_taxon", lambda self, taxid: TAXON)
    monkeypatch.setenv("NCBI_EMAIL", "me@example.org")
    base = cli_files(tmp_path)
    assert main(base + ["--out", str(tmp_path / "plain")]) == 0
    assert json.loads((tmp_path / "plain" / "organism.json").read_text()) == {
        "name": "FAM13496", "taxid": 1308, "restricted": False}
    assert json.loads((tmp_path / "plain" / "genome.json").read_text()).keys() == {"cds_tool"}
    assert main(base + ["--out", str(tmp_path / "full"), "--old-identifier", "Old-1", "--restricted",
                        "--tag", "a", "--tag", "b"]) == 0
    assert json.loads((tmp_path / "full" / "organism.json").read_text()) == {
        "name": "FAM13496", "taxid": 1308, "restricted": True, "tags": ["a", "b"]}
    assert json.loads((tmp_path / "full" / "genome.json").read_text())["old_identifier"] == "Old-1"
    assert main(base + ["--out", str(tmp_path / "same"), "--old-identifier", ID]) == 1
    assert main(base + ["--out", str(tmp_path / "dup"), "--tag", "a", "--tag", "a"]) == 1
    assert not (tmp_path / "dup").exists()


# ----------------------------------------------------------------------------- genome.json
def test_genome_json_default(tmp_path):
    make().write(tmp_path)
    text = (tmp_path / "genome.json").read_text()
    assert json.loads(text) == {"cds_tool": f"Abrinca GenBankBuilder v{__version__}"}
    assert text == '{\n    "cds_tool": "Abrinca GenBankBuilder v%s"\n}\n' % __version__


def test_genome_json_with_old_identifier(tmp_path):
    make(old_identifier="FAM13496-old ä").write(tmp_path)
    text = (tmp_path / "genome.json").read_text()
    data = json.loads(text)
    assert list(data) == ["cds_tool", "old_identifier"]
    assert data == {"cds_tool": f"Abrinca GenBankBuilder v{__version__}", "old_identifier": "FAM13496-old ä"}
    assert text.endswith("}\n") and '\n    "old_identifier"' in text


def test_genome_json_never_carries_arx_tools_input_only_or_derived_keys(tmp_path):
    make(old_identifier="x").write(tmp_path)
    data = json.loads((tmp_path / "genome.json").read_text())
    assert set(data) == {"cds_tool", "old_identifier"}  # no identifier, no contig_format, no file names


def test_cds_tool_fits_arx_field(tmp_path):
    make().write(tmp_path)
    assert len(json.loads((tmp_path / "genome.json").read_text())["cds_tool"]) <= 50


def test_cds_tool_too_long_is_rejected(tmp_path, monkeypatch):
    monkeypatch.setattr("gbk_creator.builder.package_version", lambda: "1" * 40)
    with pytest.raises(GbkCreatorError, match="50 character"):
        make().write(tmp_path / "out")
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("bad", ["", "  ", "a\nb", 5, ["x"], b"x"])
def test_old_identifier_must_be_a_non_empty_single_line_string(bad):
    with pytest.raises(GbkCreatorError, match="old_identifier must be None or a non-empty"):
        make(old_identifier=bad)


def test_old_identifier_max_50_chars():
    make(old_identifier="o" * 50)
    with pytest.raises(GbkCreatorError, match="max 50"):
        make(old_identifier="o" * 51)


def test_old_identifier_must_differ_from_identifier():
    with pytest.raises(GbkCreatorError, match="must differ"):
        make(old_identifier=ID)
