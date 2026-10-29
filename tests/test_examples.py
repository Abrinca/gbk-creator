"""The example script and the README must not rot (offline)."""
import json
import os
import re
import shutil
import subprocess
import sys

import pytest
from Bio import SeqIO
from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation, SimpleLocation

from conftest import BUILD_DATE, ROOT, TAXON, ncbi_available, random_dna
from gbk_creator import GenBankBuilder, load_categories
from gbk_creator.taxonomy import NcbiClient
from gbk_creator.writer import render_genbank
from gbk_creator.builder import ALLOWED_GENE_QUALIFIERS, ALLOWED_UNTAGGED_QUALIFIERS

README = (ROOT / "README.md").read_text()
EXAMPLE = ROOT / "examples" / "from_database.py"


def readme_section(title):
    start = README.index(title)
    nxt = re.search(r"\n#{1,3} ", README[start + len(title):])
    return README[start:start + len(title) + (nxt.start() if nxt else len(README))]


def table_rows(section):
    return [[c.strip() for c in re.split(r"(?<!\\)\|", line.strip().strip("|"))]
            for line in section.splitlines() if line.startswith("|") and not re.match(r"\|[-| ]+\|$", line)][1:]


# ----------------------------------------------------------------------------- example script
BLOCK_NETWORK = (
    "import runpy, sys, urllib.request\n"
    "def refuse(*args, **kwargs):\n"
    "    raise RuntimeError('network access attempted')\n"
    "urllib.request.urlopen = refuse\n"
    "sys.argv = sys.argv[1:]\n"
    "runpy.run_path(sys.argv[0], run_name='__main__')\n"
)


def run_example(out, cache, email="dummy@example.org", block_network=True, extra_env=None):
    """Run the example script as a user would (taxid only), optionally with every HTTP request refused."""
    env = {k: v for k, v in os.environ.items() if not k.startswith("NCBI_")}
    env["GBK_CREATOR_CACHE"] = str(cache)
    if email:
        env["NCBI_EMAIL"] = email
    env.update(extra_env or {})
    cmd = [sys.executable, "-c", BLOCK_NETWORK, str(EXAMPLE), str(out)] if block_network \
        else [sys.executable, str(EXAMPLE), str(out)]
    return subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=120)


def seed_cache(cache):
    """Pre-seed the taxon cache with the library's own cache writer (same format as after a real lookup)."""
    NcbiClient(email="dummy@example.org", cache_dir=cache).cache_taxon(TAXON)


@pytest.fixture(scope="module")
def example_out(tmp_path_factory):
    base = tmp_path_factory.mktemp("example")
    cache, out = base / "cache", base / "out" / "ACME0042-1.1"
    seed_cache(cache)
    done = run_example(out, cache)  # any network access fails the run: the taxon must come from the cache
    assert done.returncode == 0, done.stderr
    return out, done.stdout


def test_example_needs_the_cache_or_the_network(tmp_path):
    """Control for the fixture above: without the seeded cache the blocked network is hit and the run fails."""
    done = run_example(tmp_path / "out", tmp_path / "empty-cache")
    assert done.returncode != 0
    assert "NCBI E-utilities request" in done.stderr and not (tmp_path / "out").exists()


def test_example_needs_an_email(tmp_path):
    cache = tmp_path / "cache"
    seed_cache(cache)
    done = run_example(tmp_path / "out", cache, email=None)
    assert done.returncode != 0 and "NCBI_EMAIL" in done.stderr


@pytest.mark.network
@pytest.mark.skipif(not ncbi_available(), reason="NCBI_EMAIL not set")
def test_example_for_real(tmp_path):
    """The same script against NCBI (one efetch for taxid 1308, cached in a temp dir)."""
    done = run_example(tmp_path / "out", tmp_path / "cache", email=None, block_network=False,
                       extra_env={k: os.environ[k] for k in ("NCBI_EMAIL", "NCBI_API_KEY") if k in os.environ})
    assert done.returncode == 0, done.stderr
    rec = next(SeqIO.parse(tmp_path / "out" / "ACME0042-1.1.gbk", "genbank"))
    assert rec.annotations["organism"] == "Streptococcus thermophilus" and rec.annotations["taxonomy"][0] == "Bacteria"
    assert any((tmp_path / "cache").rglob("*.json"))


def test_example_writes_files_and_prints_import_command(example_out):
    out, stdout = example_out
    names = sorted(p.name for p in out.iterdir())
    assert names == sorted(["ACME0042-1.1.gbk", "ACME0042-1.1.KG", "ACME0042-1.1.EC", "ACME0042-1.1.GO",
                            "organism.json", "genome.json"])
    assert f"import_genome {out} --folder_structure_dir /path/to/folder_structure --organism ACME0042 --genome ACME0042-1.1" in stdout


def test_example_content(example_out):
    out, _ = example_out
    recs = list(SeqIO.parse(out / "ACME0042-1.1.gbk", "genbank"))
    assert [r.id for r in recs] == ["ACME0042-1.1_scf1", "ACME0042-1.1_scf2"]
    assert [r.annotations["topology"] for r in recs] == ["linear", "circular"]
    assert "plasmid" in recs[1].description
    feats = [(f.type, f.qualifiers.get("locus_tag", ["-"])[0].rsplit("_", 1)[-1], str(f.location))
             for r in recs for f in r.features if f.type != "source"]
    cds = {tag: loc for t, tag, loc in feats if t == "CDS"}
    assert cds["000001"] == "[<0:90](+)"                                         # incomplete start
    assert cds["000003"] == "[500:800](-)"                                       # minus strand
    assert cds["000007"] == "join{[1800:1900](+), [1902:2105](+)}"               # join
    assert cds["000008"] == "[2700:>3000](+)"                                    # incomplete end
    assert cds["000009"] == "join{[1120:1200](+), [0:100](+)}"                   # across the plasmid origin
    # renumbering: locus tags follow contig order and position, old IDs are kept
    old = {f.qualifiers["locus_tag"][0].rsplit("_", 1)[-1]: f.qualifiers["old_locus_tag"][0]
           for r in recs for f in r.features if f.type == "CDS"}
    assert old["000001"] == "gene_00007" and old["000009"] == "gene_00031" and old["000010"] == "gene_00032"
    pseudo = [f for r in recs for f in r.features if f.type == "CDS" and "pseudo" in f.qualifiers]
    assert len(pseudo) == 1 and "translation" not in pseudo[0].qualifiers
    trna = [f for r in recs for f in r.features if f.type == "tRNA"]
    assert trna[0].qualifiers["anticodon"] == ["(pos:1235..1237,aa:Gly,seq:tcc)"]
    assert [t for t, _, _ in feats].count("repeat_region") == 1 and [t for t, _, _ in feats].count("rRNA") == 1
    # annotations use the new tags; the duplicate KG line is ignored
    assert (out / "ACME0042-1.1.KG").read_text() == (
        "ACME0042-1.1_000002\tK02338\nACME0042-1.1_000003\tK02056\nACME0042-1.1_000010\tK03196\n")
    assert (out / "ACME0042-1.1.EC").read_text() == "ACME0042-1.1_000002\tEC:2.7.7.7\nACME0042-1.1_000004\tEC:2.7.7.-\n"
    assert json.loads((out / "genome.json").read_text())["old_identifier"] == "legacy_genome_42"
    assert json.loads((out / "organism.json").read_text())["tags"] == ["metagenomic-sample-42"]


def test_example_output_imports_with_arx_tools(example_out, tmp_path):
    from arx_tools.genbank_to_fasta import GenBankToFasta
    from arx_tools.import_genome import import_genome

    out, _ = example_out
    folder = tmp_path / "folder_structure"
    (folder / "organisms").mkdir(parents=True)
    (folder / "version.json").write_text(json.dumps({"folder_structure_version": 3}))
    import_dir = tmp_path / "in"
    shutil.copytree(out, import_dir)
    import_genome(import_dir=str(import_dir), folder_structure_dir=str(folder), organism="ACME0042", genome="ACME0042-1.1")
    genome_dir = folder / "organisms" / "ACME0042" / "genomes" / "ACME0042-1.1"
    genome = json.loads((genome_dir / "genome.json").read_text())
    assert genome["old_identifier"] == "legacy_genome_42" and genome["identifier"] == "ACME0042-1.1"
    assert {ca["type"] for ca in genome["custom_annotations"]} == {"KG", "EC", "GO"}
    GenBankToFasta.convert(str(out / "ACME0042-1.1.gbk"), str(tmp_path / "x.faa"), "faa", strict=True)


# ----------------------------------------------------------------------------- README stays true
def test_readme_embeds_the_current_conversion_code():
    source = EXAMPLE.read_text()
    excerpt = re.search(r"# --- BEGIN excerpt[^\n]*\n(.*?)# --- END excerpt", source, re.S).group(1).rstrip("\n")
    assert f"```python\n{excerpt}\n```" in README


def test_coordinate_table_matches_what_is_written():
    rows = table_rows(readme_section("## Coordinates"))
    assert len(rows) == 9
    namespace = dict(SimpleLocation=SimpleLocation, CompoundLocation=CompoundLocation,
                     BeforePosition=BeforePosition, AfterPosition=AfterPosition)
    for db, python, text in rows:
        location = eval(python.strip("`"), namespace)  # noqa: S307 - our own README
        b = GenBankBuilder("T-1", 1308, strain="T", taxon=TAXON, build_date=BUILD_DATE)
        b.add_contig("T-1_scf1", random_dna(3000))
        b.add_gene("T-1_000001", "T-1_scf1", location, "CDS", product="p", pseudo=True)
        gbk = render_genbank(b)
        lines = [l.split(None, 1)[1].strip() for l in gbk.splitlines() if re.match(r"     (gene|CDS) +\S", l)]
        assert text.strip("`") in lines, (db, text, lines)


def test_feature_tables_match_the_code():
    def qualifiers(section):
        return {row[0].strip("`"): set(re.findall(r"`([^`]+)`", row[1])) for row in table_rows(section)}

    ref = readme_section("### Feature types")
    gene_part, feature_part = ref.split("`add_feature(")
    assert qualifiers(gene_part) == ALLOWED_GENE_QUALIFIERS
    assert qualifiers(feature_part) == ALLOWED_UNTAGGED_QUALIFIERS


def test_category_table_matches_the_default_regexes():
    rows = table_rows(readme_section("Default categories"))
    documented = {row[0].strip("`"): row[2].strip("`").replace("\\|", "|") for row in rows}
    assert documented == {t: r.pattern for t, r in load_categories().items()}


def test_readme_is_for_customers():
    for word in ("Plan.md", "PGAP", "sidecar", "oracle", "frontend"):
        assert word not in README, word
    assert "DEVELOPMENT.md" in README
    assert "Plan.md" in (ROOT / "DEVELOPMENT.md").read_text()
