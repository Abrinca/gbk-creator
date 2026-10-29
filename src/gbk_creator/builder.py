"""GenBankBuilder: collects genome data, enforces the rules, writes gbk + sidecars."""
from __future__ import annotations

import datetime as _dt
import os
import re
from collections import defaultdict
from pathlib import Path

from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation, ExactPosition, SimpleLocation

from .categories import FORBIDDEN_SIDECAR_TYPES, load_categories
from .errors import AnnotationError, FeatureError, GbkCreatorError, NamingError
from .models import Contig, Gene, Location, Qualifiers, UntaggedFeature
from .taxonomy import NcbiClient, Taxon

MAX_IDENTIFIER_LENGTH = 50
MAX_STRAIN_LENGTH = 40
MAX_OLD_LOCUS_TAG_LENGTH = 100  # Arx Gene.old_locus_tag CharField(max_length=100)
# Whitespace, ':' '/' '\' as per Plan; '|' because Arx/arx-tools split locus tags on it.
_FORBIDDEN_IDENTIFIER_CHARS = set(' \t\r\n:/\\|')

CONTIG_KINDS = ("scaffold", "chromosome", "plasmid")
TOPOLOGIES = ("linear", "circular")
GENE_FEATURE_TYPES = ("CDS", "tRNA", "rRNA", "ncRNA", "tmRNA")
UNTAGGED_FEATURE_TYPES = ("repeat_region", "regulatory", "misc_feature")
VALID_TRANSL_TABLES = frozenset(list(range(1, 7)) + list(range(9, 17)) + list(range(21, 34)))

# Qualifiers (besides locus_tag/product/gene/translation/pseudo, which are explicit
# arguments of add_gene) that occur in the PGAP oracle, per feature type.
ALLOWED_GENE_QUALIFIERS = {
    "CDS": {"gene_synonym", "inference", "note", "codon_start", "transl_table", "protein_id",
            "ribosomal_slippage", "db_xref"},
    "tRNA": {"gene_synonym", "inference", "note", "anticodon", "db_xref"},
    "rRNA": {"gene_synonym", "inference", "note", "db_xref"},
    "tmRNA": {"gene_synonym", "inference", "note", "db_xref"},
    "ncRNA": {"gene_synonym", "inference", "note", "db_xref", "ncRNA_class"},
}
ALLOWED_UNTAGGED_QUALIFIERS = {
    "repeat_region": {"inference", "note", "rpt_family", "rpt_type", "rpt_unit_range", "rpt_unit_seq"},
    "regulatory": {"inference", "note", "regulatory_class", "bound_moiety", "db_xref"},
    "misc_feature": {"note"},
}
_FLAG_QUALIFIERS = {"ribosomal_slippage"}
_INT_QUALIFIERS = {"codon_start", "transl_table"}
_FORBIDDEN_QUALIFIER_HINTS = {
    "EC_number": "EC numbers must not be written to the gbk; use add_annotation(locus_tag, 'EC', 'EC:1.2.3.4')",
    "pseudo": "use the pseudo=True argument of add_gene",
    "locus_tag": "locus tags are set via the locus_tag argument of add_gene (untagged features have none)",
    "old_locus_tag": "use the old_locus_tag argument of add_gene (untagged features have no locus tag)",
    "product": "use the product argument",
    "gene": "use the gene argument",
    "translation": "use the translation argument",
    "go_function": "GO terms must not be written to the gbk; use add_annotation(locus_tag, 'GO', 'GO:0000000')",
    "go_process": "GO terms must not be written to the gbk; use add_annotation(locus_tag, 'GO', 'GO:0000000')",
    "go_component": "GO terms must not be written to the gbk; use add_annotation(locus_tag, 'GO', 'GO:0000000')",
}
_DB_XREF_RE = re.compile(r"RFAM:RF\d{5}")
_TRANSLATION_RE = re.compile(r"[A-Za-z*]+")
_DNA_RE = re.compile(r"[ACGTRYSWKMBDHVN]+")


def classify_contig(definition: str, topology: str) -> str:
    """How Arx (GenomeContent.py) classifies a contig from its DEFINITION line and topology."""
    d = definition.lower()
    if "plasmid" in d:
        return "plasmid"
    if "chromosome" in d and topology == "circular":
        return "chromosome"
    return "scaffold"


def package_version() -> str:
    from importlib.metadata import PackageNotFoundError, version
    try:
        return version("gbk-creator")
    except PackageNotFoundError:
        return "0+unknown"


def _has_control_chars(text: str) -> bool:
    return any(ord(c) < 32 or ord(c) == 127 for c in text)


class GenBankBuilder:
    def __init__(
        self,
        identifier: str,
        taxid: int,
        strain: str | None = None,
        transl_table: int = 11,
        annotation_regex: str | os.PathLike | None = None,
        ncbi_email: str | None = None,
        ncbi_api_key: str | None = None,
        *,
        taxon: Taxon | None = None,
        cache_dir: str | os.PathLike | None = None,
        build_date: _dt.datetime | None = None,
        restricted: bool = False,
        tags: list[str] | None = None,
        old_identifier: str | None = None,
    ):
        """
        :param identifier: genome identifier (Arx "genome"); max 50 chars, no space : / \\ |
        :param taxid: NCBI taxid; scientific name and lineage are fetched from NCBI (cached)
        :param strain: Arx "organism"; identifier must start with it; default = identifier; max 40 chars
        :param transl_table: default translation table for CDS features
        :param annotation_regex: JSON file with the annotation categories (plain or Arx format)
        :param ncbi_email, ncbi_api_key: default to env vars NCBI_EMAIL / NCBI_API_KEY
        :param taxon: pre-resolved :class:`Taxon` (skips the NCBI lookup; must match ``taxid``)
        :param build_date: date written to LOCUS and the annotation comment (default: now)
        :param restricted: Arx organism "restricted" flag (a real bool), for organism.json
        :param tags: Arx organism tags (non-empty, unique strings); None = key omitted from organism.json
        :param old_identifier: previous name of this genome, written to genome.json to link the new Arx genome
            to it (non-empty single line, max 50 chars, different from ``identifier``); None = key omitted
        """
        self.identifier = self._check_identifier(identifier)
        self.strain = self._check_strain(strain, self.identifier)
        if isinstance(taxid, bool) or not isinstance(taxid, int) or taxid <= 0:
            raise GbkCreatorError(f"taxid must be a positive integer, got {taxid!r}.")
        self.taxid = taxid
        if transl_table not in VALID_TRANSL_TABLES:
            raise GbkCreatorError(f"Invalid transl_table {transl_table!r}.")
        self.transl_table = transl_table
        self.categories = load_categories(annotation_regex)
        if taxon is None:
            taxon = NcbiClient(ncbi_email, ncbi_api_key, cache_dir).fetch_taxon(taxid)
        elif taxon.taxid != taxid:
            raise GbkCreatorError(f"taxon.taxid ({taxon.taxid}) does not match taxid ({taxid}).")
        self.taxon = taxon
        self.build_date = build_date or _dt.datetime.now()
        self.restricted, self.tags = self._check_organism_fields(restricted, tags)
        self.old_identifier = self._check_old_identifier(old_identifier, self.identifier)

        self._contigs: list[Contig] = []
        self._contigs_by_name: dict[str, Contig] = {}
        self._genes: dict[str, Gene] = {}
        self._old_locus_tags: dict[str, str] = {}  # old locus tag -> new locus tag
        self._features: list[UntaggedFeature] = []
        self._annotations: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
        self._order = 0
        # the GC regex validates /gene; fall back to the packaged one if the user's JSON lacks it
        self._gc_regex = self.categories.get("GC") or load_categories()["GC"]

    # ------------------------------------------------------------------ checks
    @staticmethod
    def _check_identifier(identifier: str) -> str:
        if not isinstance(identifier, str) or not identifier:
            raise NamingError("identifier must be a non-empty string.")
        if len(identifier) > MAX_IDENTIFIER_LENGTH:
            raise NamingError(
                f"identifier {identifier!r} is {len(identifier)} characters long; max {MAX_IDENTIFIER_LENGTH}.")
        bad = sorted({c for c in identifier if c in _FORBIDDEN_IDENTIFIER_CHARS}
                     | {c for c in identifier if not c.isprintable() or ord(c) > 126})
        if bad:
            raise NamingError(
                f"identifier {identifier!r} contains forbidden characters {bad!r} "
                "(no whitespace, ':', '/', '\\', '|' or non-ASCII characters).")
        return identifier

    @staticmethod
    def _check_strain(strain: str | None, identifier: str) -> str:
        if strain is None:
            strain = identifier
        if not isinstance(strain, str) or not strain:
            raise NamingError("strain must be a non-empty string.")
        if len(strain) > MAX_STRAIN_LENGTH:
            raise NamingError(f"strain {strain!r} is {len(strain)} characters long; max {MAX_STRAIN_LENGTH}.")
        if _has_control_chars(strain):
            raise NamingError(f"strain {strain!r} contains control characters.")
        if not identifier.startswith(strain):
            raise NamingError(
                f"identifier {identifier!r} must start with the strain (Arx organism) {strain!r}.")
        return strain

    @staticmethod
    def _check_old_identifier(old_identifier, identifier):
        if old_identifier is None:
            return None
        if (not isinstance(old_identifier, str) or not old_identifier.strip()
                or _has_control_chars(old_identifier)):
            raise GbkCreatorError(
                f"old_identifier must be None or a non-empty single-line string, got {old_identifier!r}.")
        if len(old_identifier) > MAX_IDENTIFIER_LENGTH:
            raise GbkCreatorError(
                f"old_identifier {old_identifier!r} is {len(old_identifier)} characters long; "
                f"max {MAX_IDENTIFIER_LENGTH} (Arx limit).")
        if old_identifier == identifier:
            raise GbkCreatorError(f"old_identifier must differ from the identifier ({identifier!r}).")
        return old_identifier

    @staticmethod
    def _check_organism_fields(restricted, tags):
        if not isinstance(restricted, bool):
            raise GbkCreatorError(f"restricted must be True or False, got {restricted!r}.")
        if tags is not None:
            if not isinstance(tags, (list, tuple)):
                raise GbkCreatorError(f"tags must be a list of strings or None, got {tags!r}.")
            for tag in tags:
                if not isinstance(tag, str) or not tag.strip() or _has_control_chars(tag):
                    raise GbkCreatorError(f"tags must be non-empty single-line strings, got {tag!r}.")
            if len(set(tags)) != len(tags):
                dup = sorted({t for t in tags if list(tags).count(t) > 1})
                raise GbkCreatorError(f"tags must not contain duplicates: {dup}.")
            tags = list(tags)
        return restricted, tags

    def _contig(self, name: str) -> Contig:
        try:
            return self._contigs_by_name[name]
        except (KeyError, TypeError):
            raise FeatureError(
                f"Unknown contig {name!r}; add it first with add_contig(). "
                f"Known contigs: {list(self._contigs_by_name)[:5]}") from None

    # ------------------------------------------------------------------ contigs
    def add_contig(self, name: str, sequence, topology: str = "linear", kind: str = "scaffold") -> None:
        """Add a contig. Names are ``{identifier}_scf{n}`` and must be added in the order n = 1, 2, 3, ..."""
        expected = f"{self.identifier}_scf{len(self._contigs) + 1}"
        if name != expected:
            hint = ""
            if isinstance(name, str) and re.fullmatch(re.escape(self.identifier) + r"_scf0\d*", name):
                hint = " (no leading zeros)"
            elif name in self._contigs_by_name:
                hint = " (duplicate contig)"
            raise NamingError(
                f"Contig name {name!r} is invalid: contigs must be named '{{identifier}}_scf{{n}}' with "
                f"n = 1..N in the order of add_contig() calls, no leading zeros{hint}. Expected {expected!r}.")
        if topology not in TOPOLOGIES:
            raise GbkCreatorError(f"topology must be one of {TOPOLOGIES}, got {topology!r}.")
        if kind not in CONTIG_KINDS:
            raise GbkCreatorError(f"kind must be one of {CONTIG_KINDS}, got {kind!r}.")
        if kind == "chromosome" and topology != "circular":
            raise GbkCreatorError(
                f"Contig {name!r}: a chromosome must be circular (Arx classifies a contig as chromosome "
                "only if it is circular). Use kind='scaffold' for linear sequences.")
        seq = str(sequence).upper()
        if not seq:
            raise GbkCreatorError(f"Contig {name!r} has an empty sequence.")
        if not _DNA_RE.fullmatch(seq):
            bad = sorted(set(seq) - set("ACGTRYSWKMBDHVN"))
            raise GbkCreatorError(f"Contig {name!r} contains invalid nucleotide characters: {bad[:10]}")
        contig = Contig(name=name, number=len(self._contigs) + 1, sequence=seq,
                        topology=topology, kind=kind, index=len(self._contigs))
        self._contigs.append(contig)
        self._contigs_by_name[name] = contig

    # ------------------------------------------------------------------ locations
    @staticmethod
    def _check_location(location, contig: Contig, what: str) -> Location:
        if not isinstance(location, (SimpleLocation, CompoundLocation)):
            raise FeatureError(
                f"{what}: location must be a Biopython SimpleLocation/CompoundLocation, got {type(location).__name__}.")
        if isinstance(location, CompoundLocation) and location.operator != "join":
            raise FeatureError(
                f"{what}: only join() compound locations are supported, got operator {location.operator!r}.")
        parts = list(location.parts)
        strands = set()
        for part in parts:
            if part.strand not in (1, -1):
                raise FeatureError(f"{what}: strand must be +1 or -1 (got {part.strand!r}).")
            strands.add(part.strand)
            if part.ref is not None or part.ref_db is not None:
                raise FeatureError(f"{what}: remote locations (ref/ref_db) are not supported.")
            if type(part.start) not in (ExactPosition, BeforePosition):
                raise FeatureError(
                    f"{what}: unsupported start position {part.start!r}; only exact positions and "
                    "'<' (BeforePosition) partial starts are allowed.")
            if type(part.end) not in (ExactPosition, AfterPosition):
                raise FeatureError(
                    f"{what}: unsupported end position {part.end!r}; only exact positions and "
                    "'>' (AfterPosition) partial ends are allowed.")
            if not 0 <= int(part.start) < int(part.end) <= len(contig.sequence):
                raise FeatureError(
                    f"{what}: location part {int(part.start)}..{int(part.end)} (0-based, end exclusive) is "
                    f"outside contig {contig.name!r} of length {len(contig.sequence)}.")
        if len(strands) != 1:
            raise FeatureError(f"{what}: all parts of a join must be on the same strand.")
        return location

    # ------------------------------------------------------------------ qualifiers
    @staticmethod
    def _normalize_qualifiers(feature_type: str, qualifiers: dict, allowed: set[str], what: str) -> Qualifiers:
        out: Qualifiers = {}
        for key, value in qualifiers.items():
            if key not in allowed:
                hint = _FORBIDDEN_QUALIFIER_HINTS.get(key)
                if hint:
                    raise FeatureError(f"{what}: qualifier /{key} is not allowed here: {hint}.")
                raise FeatureError(
                    f"{what}: qualifier /{key} is not allowed on {feature_type} (not used by the PGAP "
                    f"oracle). Allowed: {sorted(allowed)}.")
            if key in _FLAG_QUALIFIERS:
                if value is True:
                    out[key] = [None]
                elif value in (False, None):
                    continue
                else:
                    raise FeatureError(f"{what}: /{key} is a flag; pass {key}=True.")
                continue
            values = list(value) if isinstance(value, (list, tuple)) else [value]
            if not values:
                raise FeatureError(f"{what}: /{key} has no value.")
            norm = []
            for v in values:
                if key in _INT_QUALIFIERS:
                    if isinstance(v, bool) or not isinstance(v, int):
                        raise FeatureError(f"{what}: /{key} must be an integer, got {v!r}.")
                    norm.append(v)
                    continue
                if not isinstance(v, str) or not v.strip():
                    raise FeatureError(f"{what}: /{key} must be a non-empty string, got {v!r}.")
                if _has_control_chars(v):
                    raise FeatureError(f"{what}: /{key} contains control characters / newlines: {v!r}.")
                if key == "db_xref" and not _DB_XREF_RE.fullmatch(v):
                    raise FeatureError(f"{what}: /db_xref must look like 'RFAM:RF00010', got {v!r}.")
                norm.append(v)
            out[key] = norm
        return out

    # ------------------------------------------------------------------ genes
    def add_gene(
        self,
        locus_tag: str,
        contig: str,
        location: Location,
        feature_type: str,
        product: str,
        gene: str | None = None,
        pseudo: bool = False,
        translation: str | None = None,
        codon_start: int | None = None,
        old_locus_tag: str | None = None,
        **extra_qualifiers,
    ) -> None:
        """Add a locus-tagged gene. The `gene` parent feature is created automatically (same location,
        except for joins that do not wrap around the contig origin: there it spans all parts, as in PGAP).

        ``old_locus_tag`` (your previous locus tag; non-empty single line, no ``"``, no surrounding whitespace,
        max 100 chars, unique, different from ``locus_tag``) is written to both the gene and the child feature.

        ``extra_qualifiers`` may hold only qualifiers that occur in the PGAP oracle (see
        ``ALLOWED_GENE_QUALIFIERS``); valueless flags such as ribosomal_slippage take ``True``.
        """
        what = f"Gene {locus_tag!r}"
        m = re.fullmatch(re.escape(self.identifier) + r"_(\d{6})", locus_tag) if isinstance(locus_tag, str) else None
        if not m:
            raise NamingError(
                f"{what}: locus tags must be '{self.identifier}_NNNNNN' with exactly 6 digits.")
        if locus_tag in self._genes:
            raise NamingError(f"{what}: duplicate locus tag.")
        contig_obj = self._contig(contig)
        if feature_type not in GENE_FEATURE_TYPES:
            raise FeatureError(f"{what}: feature_type must be one of {GENE_FEATURE_TYPES}, got {feature_type!r}.")
        location = self._check_location(location, contig_obj, what)
        if not isinstance(product, str) or not product.strip() or _has_control_chars(product):
            raise FeatureError(f"{what}: product is required (non-empty, single line).")
        if gene is not None:
            if not isinstance(gene, str) or _has_control_chars(gene) or not self._gc_regex.match(f"GC:{gene}"):
                raise FeatureError(
                    f"{what}: /gene {gene!r} does not match the GC regex {self._gc_regex.pattern!r} "
                    "(3-11 characters of [0-9a-zA-Z_/-' ]).")
        self._check_old_locus_tag(old_locus_tag, locus_tag, what)
        is_cds = feature_type == "CDS"
        if not isinstance(pseudo, bool):
            raise FeatureError(f"{what}: pseudo must be True or False, got {pseudo!r}.")
        if pseudo and not is_cds:
            raise FeatureError(f"{what}: only CDS genes can be pseudo (as in the PGAP oracle).")
        if pseudo and translation is not None:
            raise FeatureError(f"{what}: a pseudo gene must not have a translation.")
        if is_cds and not pseudo and translation is None:
            raise FeatureError(f"{what}: a non-pseudo CDS requires a translation (use pseudo=True otherwise).")
        if translation is not None:
            if not is_cds:
                raise FeatureError(f"{what}: only CDS features can have a translation.")
            if not isinstance(translation, str) or not _TRANSLATION_RE.fullmatch(translation):
                raise FeatureError(f"{what}: translation must be a string of amino-acid letters.")
        if not is_cds and codon_start is not None:
            raise FeatureError(f"{what}: codon_start is only valid for CDS features.")

        if is_cds:
            extra_qualifiers.setdefault("transl_table", self.transl_table)
            extra_qualifiers["codon_start"] = 1 if codon_start is None else codon_start
        q = self._normalize_qualifiers(feature_type, extra_qualifiers, ALLOWED_GENE_QUALIFIERS[feature_type], what)
        if is_cds:
            if q["codon_start"][0] not in (1, 2, 3):
                raise FeatureError(f"{what}: codon_start must be 1, 2 or 3.")
            if len(q["transl_table"]) != 1 or q["transl_table"][0] not in VALID_TRANSL_TABLES:
                raise FeatureError(f"{what}: invalid transl_table {q['transl_table']!r}.")

        gene_q: Qualifiers = {}
        if gene is not None:
            gene_q["gene"] = [gene]
        gene_q["locus_tag"] = [locus_tag]
        if "gene_synonym" in q:
            gene_q["gene_synonym"] = list(q["gene_synonym"])
        if pseudo:
            gene_q["pseudo"] = [None]

        child_q = _order_child_qualifiers(feature_type, gene, locus_tag, product, pseudo, translation, q)
        if old_locus_tag is not None:
            # like arx-tools' rename mode: appended as the last qualifier of the gene and of its child feature
            gene_q["old_locus_tag"] = [old_locus_tag]
            child_q["old_locus_tag"] = [old_locus_tag]
            self._old_locus_tags[old_locus_tag] = locus_tag
        gene_location = location
        starts = [int(p.start) for p in location.parts]
        in_order = starts == sorted(starts, reverse=location.strand == -1)
        if len(starts) > 1 and in_order:
            # PGAP: the gene of a join (slippage, intron) spans all parts (CDS: join(a..b,b..c), gene: a..c);
            # only joins across the origin of a circular contig keep the join on the gene
            low = min(location.parts, key=lambda p: int(p.start))
            high = max(location.parts, key=lambda p: int(p.end))
            gene_location = SimpleLocation(low.start, high.end, location.strand)
        self._genes[locus_tag] = Gene(
            locus_tag=locus_tag, number=int(m.group(1)), contig=contig_obj, location=location,
            gene_location=gene_location,
            feature_type=feature_type, pseudo=pseudo, gene_qualifiers=gene_q,
            child_qualifiers=child_q, order=self._next_order())

    def _check_old_locus_tag(self, old, locus_tag, what):
        if old is None:
            return
        if not isinstance(old, str) or not old or _has_control_chars(old):
            raise NamingError(f"{what}: old_locus_tag must be None or a non-empty single-line string, got {old!r}.")
        if old != old.strip():
            raise NamingError(f"{what}: old_locus_tag {old!r} must not have leading or trailing whitespace.")
        if '"' in old:
            raise NamingError(f"{what}: old_locus_tag {old!r} must not contain a double quote.")
        if len(old) > MAX_OLD_LOCUS_TAG_LENGTH:
            raise NamingError(
                f"{what}: old_locus_tag is {len(old)} characters long; max {MAX_OLD_LOCUS_TAG_LENGTH} (Arx limit).")
        if old == locus_tag:
            raise NamingError(f"{what}: old_locus_tag must differ from the new locus tag.")
        if old in self._old_locus_tags:
            raise NamingError(
                f"{what}: old_locus_tag {old!r} is already used by {self._old_locus_tags[old]}; "
                "old locus tags must be unique.")

    def add_feature(self, contig: str, feature_type: str, location: Location, **qualifiers) -> None:
        """Add a feature without locus tag (``repeat_region``, ``regulatory``, ``misc_feature``)."""
        what = f"{feature_type} feature"
        if feature_type not in UNTAGGED_FEATURE_TYPES:
            hint = " Use add_gene() for locus-tagged features." if feature_type in GENE_FEATURE_TYPES + ("gene",) else ""
            raise FeatureError(f"feature_type must be one of {UNTAGGED_FEATURE_TYPES}, got {feature_type!r}.{hint}")
        contig_obj = self._contig(contig)
        location = self._check_location(location, contig_obj, what)
        q = self._normalize_qualifiers(
            feature_type, qualifiers, ALLOWED_UNTAGGED_QUALIFIERS[feature_type], what)
        if feature_type == "regulatory" and "regulatory_class" not in q:
            raise FeatureError(f"{what}: /regulatory_class is required.")
        if feature_type == "misc_feature" and "note" not in q:
            raise FeatureError(f"{what}: /note is required.")
        order = ["regulatory_class", "inference", "rpt_family", "rpt_type", "rpt_unit_range", "rpt_unit_seq",
                 "note", "bound_moiety", "db_xref"]
        q = {k: q[k] for k in sorted(q, key=order.index)}
        self._features.append(UntaggedFeature(contig_obj, feature_type, location, q, self._next_order()))

    def _next_order(self) -> int:
        self._order += 1
        return self._order

    # ------------------------------------------------------------------ annotations
    def add_annotation(self, locus_tag: str, anno_type: str, anno_id: str) -> None:
        """Attach a custom annotation (e.g. ``'KG', 'K00001'``) to a CDS gene. Duplicates are ignored."""
        what = f"Annotation {anno_type}/{anno_id!r} for {locus_tag!r}"
        gene = self._genes.get(locus_tag)
        if gene is None:
            raise AnnotationError(f"{what}: unknown locus tag; add the gene first with add_gene().")
        if anno_type in FORBIDDEN_SIDECAR_TYPES:
            raise AnnotationError(
                f"{what}: GC and GP sidecars are forbidden; Arx derives them from /gene and /product. "
                "Set gene= / product= on add_gene() instead.")
        regex = self.categories.get(anno_type)
        if regex is None:
            raise AnnotationError(
                f"{what}: unknown annotation category {anno_type!r}. Known: {sorted(self.categories)}. "
                "To use a new category, add it to your annotation JSON file and pass it as annotation_regex=.")
        if gene.feature_type != "CDS":
            raise AnnotationError(
                f"{what}: annotations are only allowed on CDS genes, but {locus_tag} is a {gene.feature_type}.")
        if not isinstance(anno_id, str) or not anno_id:
            raise AnnotationError(f"{what}: annotation id must be a non-empty string.")
        if _has_control_chars(anno_id) or "," in anno_id or anno_id != anno_id.strip():
            raise AnnotationError(
                f"{what}: annotation ids must not contain commas, tabs, newlines or surrounding whitespace "
                "(they would corrupt the sidecar file).")
        if not regex.match(anno_id):
            raise AnnotationError(
                f"{what}: does not match the regex of category {anno_type}: {regex.pattern!r}.")
        self._annotations[anno_type][locus_tag].add(anno_id)

    # ------------------------------------------------------------------ views
    @property
    def contigs(self) -> list[Contig]:
        return list(self._contigs)

    @property
    def genes(self) -> list[Gene]:
        """Genes sorted by contig order and position (the order locus tags must follow)."""
        return sorted(self._genes.values(), key=lambda g: (g.contig.index, g.start, g.order))

    def definition(self, contig: Contig) -> str:
        t = self.taxon
        base = f"{t.scientific_name} strain {self.strain}"
        same_kind = [c for c in self._contigs if c.kind == contig.kind]
        k = same_kind.index(contig) + 1
        if contig.kind == "scaffold":
            return f"{base} scaffold {contig.number}, whole genome shotgun sequence"
        if contig.kind == "chromosome":
            if len(same_kind) == 1:
                return f"{base} chromosome, complete genome"
            return f"{base} chromosome {k}, complete sequence"
        return f"{base} plasmid unnamed{k}, complete sequence"

    # ------------------------------------------------------------------ validation
    def validate(self) -> None:
        """Check all global rules. Called by :meth:`write` before anything is written."""
        if not self._contigs:
            raise GbkCreatorError("No contigs: add at least one contig with add_contig().")
        if not self._genes:
            raise GbkCreatorError("No genes: add at least one gene with add_gene().")
        previous = None
        for g in self.genes:
            if previous is not None and g.number <= previous.number:
                raise NamingError(
                    f"Locus tags must increase strictly along contig order and position: {g.locus_tag} "
                    f"({g.contig.name}:{g.start + 1}) comes after {previous.locus_tag} "
                    f"({previous.contig.name}:{previous.start + 1}).")
            previous = g
        for c in self._contigs:
            definition = self.definition(c)
            if classify_contig(definition, c.topology) != c.kind:
                raise GbkCreatorError(
                    f"Contig {c.name!r}: the generated DEFINITION {definition!r} would make Arx classify it "
                    f"as {classify_contig(definition, c.topology)!r} instead of {c.kind!r} (do the strain or "
                    "scientific name contain 'plasmid' or 'chromosome'?).")
        for anno_type in self._annotations:
            if anno_type in FORBIDDEN_SIDECAR_TYPES or anno_type not in self.categories:
                raise AnnotationError(f"Annotation type {anno_type!r} is not allowed.")

    # ------------------------------------------------------------------ output
    def write(self, out_dir: str | os.PathLike) -> list[Path]:
        """Validate, then write ``{out_dir}/{identifier}.gbk``, ``organism.json``, ``genome.json`` and
        one ``{identifier}.{TYPE}`` per used type.

        Stale sidecars ``{identifier}.{TYPE}`` of known categories that are not written now are removed.
        Returns the written paths.
        """
        from .writer import render_genbank, render_genome_json, render_organism_json, render_sidecars

        self.validate()
        gbk_text = render_genbank(self)
        sidecars = render_sidecars(self)
        organism_text = render_organism_json(self)
        genome_text = render_genome_json(self)

        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        written = []
        path = out / f"{self.identifier}.gbk"
        path.write_text(gbk_text, encoding="utf-8")
        written.append(path)
        path = out / "organism.json"
        path.write_text(organism_text, encoding="utf-8")
        written.append(path)
        path = out / "genome.json"
        path.write_text(genome_text, encoding="utf-8")
        written.append(path)
        for anno_type, text in sidecars.items():
            path = out / f"{self.identifier}.{anno_type}"
            path.write_text(text, encoding="utf-8")
            written.append(path)
        for anno_type in set(self.categories) | FORBIDDEN_SIDECAR_TYPES:
            if anno_type not in sidecars:
                (out / f"{self.identifier}.{anno_type}").unlink(missing_ok=True)
        return written


_CDS_ORDER = ["gene", "locus_tag", "gene_synonym", "inference", "ribosomal_slippage", "note", "pseudo",
              "codon_start", "transl_table", "product", "protein_id", "db_xref", "translation"]
_RNA_ORDER = ["ncRNA_class", "gene", "locus_tag", "gene_synonym", "product", "inference", "note",
              "anticodon", "db_xref"]


def _order_child_qualifiers(feature_type, gene, locus_tag, product, pseudo, translation, q) -> Qualifiers:
    merged: Qualifiers = {"locus_tag": [locus_tag], "product": [product]}
    if gene is not None:
        merged["gene"] = [gene]
    if pseudo:
        merged["pseudo"] = [None]
    if translation is not None:
        merged["translation"] = [translation]
    merged.update(q)
    order = _CDS_ORDER if feature_type == "CDS" else _RNA_ORDER
    return {k: merged[k] for k in sorted(merged, key=lambda k: (order.index(k) if k in order else len(order), k))}
