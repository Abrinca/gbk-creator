"""Exceptions raised by gbk-creator.

Everything derives from :class:`GbkCreatorError` (a ``ValueError``), so callers can
catch one class. The library follows the "fail loudly" rule: anything that cannot be
written out faithfully, or that violates the PGAP-canonical rules, raises.
"""


class GbkCreatorError(ValueError):
    """Base class: invalid input or a violated GenBank/Arx rule."""


class NamingError(GbkCreatorError):
    """Identifier, strain, contig name or locus tag does not follow the naming rules."""


class FeatureError(GbkCreatorError):
    """A gene/feature, its location or one of its qualifiers is invalid."""


class AnnotationError(GbkCreatorError):
    """A custom annotation (sidecar) is invalid or its category is unknown."""


class TaxonomyError(GbkCreatorError):
    """NCBI taxonomy lookup failed (bad name/taxid, no credentials, network error)."""
