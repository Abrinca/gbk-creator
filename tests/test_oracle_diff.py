"""The oracle-diff tool itself: it must accept only the documented differences (synthetic, offline)."""
import pytest
from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation, SimpleLocation

from conftest import BUILD_DATE, TAXON, random_dna
from gbk_creator import GenBankBuilder
from oracle_diff import compare_gbk

SEQ = list(random_dna(400, 3))
SEQ[99] = SEQ[100] = "A"  # homopolymer at the slippage junction: shifting it by one base keeps the sequence
SEQ = "".join(SEQ)


def part(a, b, strand=1, start=None, end=None):
    return SimpleLocation(start if start is not None else a, end if end is not None else b, strand)


def write(tmp_path, name, location, **kw):
    b = GenBankBuilder("G-1", 1308, strain="G", taxon=TAXON, build_date=BUILD_DATE)
    b.add_contig("G-1_scf1", SEQ)
    b.add_gene("G-1_000001", "G-1_scf1", location, "CDS", product="p", translation="MKV", **kw)
    b.write(tmp_path / name)
    return tmp_path / name / "G-1.gbk"


def mutated(path, old, new, tmp_path, name="mutated.gbk"):
    text = path.read_text()
    assert old in text, old
    out = tmp_path / name
    out.write_text(text.replace(old, new))
    return out


JOIN = CompoundLocation([part(50, 100), part(99, 200)])


@pytest.fixture
def oracle(tmp_path):
    return write(tmp_path, "oracle", JOIN)


def test_identical_files_have_no_differences(oracle, tmp_path):
    assert compare_gbk(oracle, write(tmp_path, "same", JOIN)) == []


def test_junction_shift_inside_homopolymer_is_the_only_join_exception(oracle, tmp_path):
    stats = {}
    shifted = write(tmp_path, "shifted", CompoundLocation([part(50, 101), part(100, 200)]))
    assert compare_gbk(oracle, shifted, stats) == []
    assert stats == {"joins_rebuilt_differently": 1}


@pytest.mark.parametrize("name,location,reason", [
    ("junction changes the sequence", CompoundLocation([part(50, 100), part(100, 200)]), "extracted sequence differs"),
    ("different span", CompoundLocation([part(50, 101), part(100, 201)]), "overall span differs"),
    ("more parts", CompoundLocation([part(50, 75), part(75, 100), part(99, 200)]), "2 parts != 3 parts"),
    ("not compound", SimpleLocation(50, 200, 1), "not both compound"),
    ("partial start added", CompoundLocation([part(50, 100, start=BeforePosition(50)), part(99, 200)]), "outer ends differ"),
    ("partial end added", CompoundLocation([part(50, 100), part(99, 200, end=AfterPosition(200))]), "outer ends differ"),
    ("fuzzy internal junction", CompoundLocation([part(50, 101), part(100, 200, start=BeforePosition(100))]), "internal junction"),
])
def test_join_exception_is_tight(oracle, tmp_path, name, location, reason):
    diffs = compare_gbk(oracle, write(tmp_path, "built", location))
    assert diffs and any(reason in d for d in diffs), diffs


def test_dropping_partial_markers_is_reported(tmp_path):
    partial = CompoundLocation([part(50, 100, start=BeforePosition(50)), part(99, 200, end=AfterPosition(200))])
    oracle = write(tmp_path, "oracle", partial)
    assert compare_gbk(oracle, write(tmp_path, "same", partial)) == []
    gbk = oracle.read_text()
    assert "join(<51..100,100..>200)" in gbk
    for old, new in (("<51..100", "51..100"), ("..>200", "..200")):
        diffs = compare_gbk(oracle, mutated(oracle, old, new, tmp_path))
        assert any("fuzzy (</>) outer ends differ" in d for d in diffs), diffs
    # the same junction shift is fine when the markers are kept
    shifted = write(tmp_path, "shifted", CompoundLocation([part(50, 101, start=BeforePosition(50)), part(100, 200, end=AfterPosition(200))]))
    assert compare_gbk(oracle, shifted) == []


def test_operator_and_strand_changes_are_reported(oracle, tmp_path):
    diffs = compare_gbk(oracle, mutated(oracle, "join(51..100,100..200)", "order(51..100,100..200)", tmp_path))
    assert any("operator" in d for d in diffs), diffs
    diffs = compare_gbk(oracle, mutated(oracle, "join(51..100,100..200)", "complement(join(51..100,100..200))", tmp_path))
    assert diffs


def test_ec_number_is_stripped_from_the_oracle_side_only(oracle, tmp_path):
    ec = '/EC_number="1.1.1.1"\n                     /transl_table=11'
    with_ec = mutated(oracle, "/transl_table=11", ec, tmp_path, "with_ec.gbk")
    plain = write(tmp_path, "plain", JOIN)
    assert compare_gbk(with_ec, plain) == []  # oracle-only /EC_number: allowed
    diffs = compare_gbk(plain, with_ec)  # generated output carrying /EC_number: reported
    assert any("qualifier /EC_number" in d for d in diffs), diffs
    assert any("qualifier /EC_number" in d for d in compare_gbk(with_ec, with_ec))  # not symmetric on purpose


def test_other_oracle_only_qualifiers_are_still_reported(oracle, tmp_path):
    extra = mutated(oracle, "/transl_table=11", '/note="x"\n                     /transl_table=11', tmp_path)
    assert compare_gbk(extra, oracle) and compare_gbk(oracle, extra)
