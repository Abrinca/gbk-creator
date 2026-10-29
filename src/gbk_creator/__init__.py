"""gbk-creator: build PGAP-canonical GenBank files + custom-annotation sidecars for Arx."""
from .builder import GenBankBuilder, package_version
from .categories import load_categories
from .errors import (AnnotationError, FeatureError, GbkCreatorError, NamingError,
                     TaxonomyError)
from .taxonomy import Taxon, fetch_taxon, name_to_taxid

__version__ = package_version()

__all__ = [
    "GenBankBuilder", "Taxon", "fetch_taxon", "name_to_taxid", "load_categories", "main",
    "GbkCreatorError", "NamingError", "FeatureError", "AnnotationError", "TaxonomyError",
    "__version__",
]


def main(argv: list[str] | None = None) -> int:
    """Entry point of the ``gbk-creator`` command (test frontend: .arx.gff + .fna -> gbk)."""
    from .frontends.arx_gff import main as _main
    return _main(argv)
