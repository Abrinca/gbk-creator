"""Custom-annotation categories (``{TYPE: regex}``) loaded from JSON.

Two JSON layouts are accepted:

* plain:     ``{"KG": "^K[0-9]{5}$", ...}``
* Arx style: ``{"KG": {"name": ..., "color": ..., "regex": "^K[0-9]{5}$", "hyperlinks": [...]}, ...}``
"""
from __future__ import annotations

import json
import os
import re
from importlib import resources

from .errors import AnnotationError

# Arx derives these from /gene and /product; sidecars for them are forbidden.
FORBIDDEN_SIDECAR_TYPES = frozenset({"GC", "GP"})
# Arx stores Annotation.anno_type in CharField(max_length=2); arx-tools globs `*.{TYPE}`
_TYPE_RE = re.compile(r"[A-Z]{2}")


def default_categories_path():
    return resources.files("gbk_creator") / "data" / "annotation_regex.json"


def parse_categories(data: object, source: str = "<dict>") -> dict[str, re.Pattern]:
    if not isinstance(data, dict) or not data:
        raise AnnotationError(f"{source}: annotation categories must be a non-empty JSON object.")
    categories = {}
    for anno_type, spec in data.items():
        if not _TYPE_RE.fullmatch(anno_type):
            raise AnnotationError(
                f"{source}: invalid annotation type {anno_type!r}: category names must be exactly two "
                "upper-case letters (e.g. 'KG'), because Arx stores them in a 2-character field.")
        regex = spec.get("regex") if isinstance(spec, dict) else spec
        if not isinstance(regex, str):
            raise AnnotationError(
                f"{source}: category {anno_type!r} needs a regex string, either directly or as "
                '{"regex": "..."}.')
        try:
            categories[anno_type] = re.compile(regex)
        except re.error as e:
            raise AnnotationError(f"{source}: bad regex for category {anno_type!r}: {e}") from None
    return categories


def load_categories(path: str | os.PathLike | None = None) -> dict[str, re.Pattern]:
    """Load categories from ``path`` (default: the packaged ``annotation_regex.json``)."""
    if path is None:
        text, source = default_categories_path().read_text(), "packaged annotation_regex.json"
    else:
        try:
            with open(path, encoding="utf-8") as f:
                text = f.read()
        except OSError as e:
            raise AnnotationError(f"Cannot read annotation categories file {path!r}: {e}") from None
        source = str(path)
    try:
        data = json.loads(text)
    except ValueError as e:
        raise AnnotationError(f"{source}: invalid JSON: {e}") from None
    return parse_categories(data, source)
