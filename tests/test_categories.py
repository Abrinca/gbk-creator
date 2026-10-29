
import pytest

from conftest import ROOT
from gbk_creator import AnnotationError, load_categories
from gbk_creator.categories import default_categories_path


def test_default_is_packaged_copy_of_repo_file():
    assert default_categories_path().read_text() == (ROOT / "annotation_regex.json").read_text()
    cats = load_categories()
    assert {"EC", "KG", "KR", "GO", "OL", "GC", "GP", "EP", "EO", "ED", "SL", "AR"} <= set(cats)


def test_both_formats(tmp_path):
    plain = tmp_path / "a.json"
    plain.write_text('{"KG": "^K[0-9]{5}$"}')
    arx = tmp_path / "b.json"
    arx.write_text('{"KG": {"name": "KEGG", "color": "x", "regex": "^K[0-9]{5}$", "hyperlinks": []}}')
    assert load_categories(plain)["KG"].pattern == load_categories(arx)["KG"].pattern


def test_arx_example_file_if_present():
    path = ROOT / "data" / "annotations.json"
    if not path.is_file():
        pytest.skip("data/annotations.json not available")
    cats = load_categories(path)
    assert cats["KG"].match("K00001") and not cats["KG"].match("K0001")


@pytest.mark.parametrize("content", ["not json", "[]", "{}", '{"KG": 5}', '{"KG": {"name": "x"}}',
                                     '{"KG": "("}', '{"kg": "^k$"}', '{"TOOLONGTYPE1": "x"}',
                                     '{"K": "x"}', '{"KGX": "x"}', '{"K1": "x"}'])
def test_invalid_files(tmp_path, content):
    p = tmp_path / "x.json"
    p.write_text(content)
    with pytest.raises(AnnotationError):
        load_categories(p)


def test_missing_file(tmp_path):
    with pytest.raises(AnnotationError, match="Cannot read"):
        load_categories(tmp_path / "nope.json")
