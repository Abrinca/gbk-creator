"""Test frontend: ``.arx.gff`` + ``.fna`` (+ ``.faa``/``.ffn``) + emapper annotations -> GenBankBuilder.

The library's primary interface is the Python API (:class:`gbk_creator.GenBankBuilder`); this module
only exists to drive it with real data (the Arx-style GFF3 that arx-tools derives from a PGAP gbk).

Known GFF losses, handled here as far as possible:

* ``join(...)`` locations collapse to one span. With a ``.ffn`` they are reconstructed for
  ribosomal-slippage CDSs and tRNAs with introns (two parts overlapping by, or separated by, a few
  bases) and for genes spanning the origin of a circular contig. Inside homopolymers the junction is
  ambiguous (the extracted sequence is the same; ~27 of ~60 slippage joins differ from the oracle in
  that way). Without ``.ffn`` the collapsed span is used as is.
* Partial ends ``<`` / ``>`` survive (they are kept in the coordinate columns).
* Translations come from the GFF ``translation`` attribute and/or ``.faa`` (they must agree).
"""
from __future__ import annotations

import argparse
import logging
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

from Bio.Seq import Seq
from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation, ExactPosition, SimpleLocation

from ..builder import GENE_FEATURE_TYPES, UNTAGGED_FEATURE_TYPES, GenBankBuilder
from ..errors import GbkCreatorError

_IGNORED_ATTRIBUTES = {"ID", "Name", "Parent", "locus_tag", "product", "gene", "pseudo", "codon_start",
                       "EC_number", "Ontology_term", "translation"}
MAX_SLIPPAGE_GAP = 10
_SKIPPED_TYPES = {"remark", "source", "gene"}  # `gene` parents are created by the builder


log = logging.getLogger(__name__)


class FrontendError(GbkCreatorError):
    """Input files cannot be converted faithfully."""


# --------------------------------------------------------------------------- parsers
def read_fasta(path) -> list[tuple[str, str]]:
    """``[(first header token, sequence), ...]`` in file order."""
    records, name, chunks = [], None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name, chunks = line[1:].split(None, 1)[0], []
            elif line.strip():
                if name is None:
                    raise FrontendError(f"{path}: sequence before first FASTA header.")
                chunks.append(line.strip())
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def read_locus_fasta(path) -> dict[str, str]:
    """FASTA keyed by locus tag; handles ``>gnl|extdb|TAG product`` and ``>TAG product`` headers."""
    out = {}
    for name, seq in read_fasta(path):
        tag = name.rsplit("|", 1)[-1]
        if tag in out:
            raise FrontendError(f"{path}: duplicate record {tag!r}.")
        out[tag] = seq
    return out


def parse_attributes(text: str) -> dict[str, list[str]]:
    attrs: dict[str, list[str]] = {}
    for item in text.strip().split(";"):
        if not item:
            continue
        key, sep, value = item.partition("=")
        if not sep:
            raise FrontendError(f"Bad GFF attribute {item!r}.")
        attrs[unquote(key)] = [unquote(v) for v in value.split(",")] if value != "" else [""]
    return attrs


class GffRow:
    def __init__(self, seqid, type_, start, end, strand, attrs, line_no):
        self.seqid, self.type, self.strand, self.attrs, self.line_no = seqid, type_, strand, attrs, line_no
        self.start_partial = start.startswith("<")
        self.end_partial = end.startswith(">")
        try:
            self.start = int(start.lstrip("<>")) - 1  # 0-based
            self.end = int(end.lstrip("<>"))
        except ValueError:
            raise FrontendError(f"GFF line {line_no}: bad coordinates {start!r}..{end!r}.") from None
        if start.startswith(">") or end.startswith("<"):
            raise FrontendError(f"GFF line {line_no}: unsupported partial markers {start!r}..{end!r}.")

    def first(self, key, default=None):
        values = self.attrs.get(key)
        return values[0] if values else default


def read_gff(path) -> list[GffRow]:
    rows = []
    with open(path) as f:
        for n, line in enumerate(f, 1):
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                raise FrontendError(f"{path}:{n}: expected 9 tab-separated columns, got {len(cols)}.")
            rows.append(GffRow(cols[0], cols[2], cols[3], cols[4], cols[6], parse_attributes(cols[8]), n))
    return rows


def read_emapper(path) -> dict[str, dict[str, set[str]]]:
    """``{locus_tag: {TYPE: {ids}}}`` from an emapper 1.0.3 file (22 columns), as Arx's ``load_eggnog_file``."""
    result: dict[str, dict[str, set[str]]] = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 22:
                raise FrontendError(f"{path}: expected 22 columns (emapper 1.0.3), got {len(cols)}.")
            annos = result.setdefault(cols[0].rsplit("|", 1)[-1], {})

            def add(anno_type, values):
                annos.setdefault(anno_type, set()).update(values)

            if cols[6] != "":
                add("GO", cols[6].split(","))
            if cols[7] != "":
                add("EC", [f"EC:{x}" for x in cols[7].split(",")])
            if cols[8] != "":
                add("KG", [x[3:] for x in cols[8].split(",")])
            if cols[11] != "":
                add("KR", cols[11].split(","))
            if cols[5] != "":
                add("EP", [f"EP:{cols[5]}"])
            if cols[18] != "":
                add("EO", [f"EO:{cols[18].split('@', 1)[0]}"])
            if cols[21] != "":
                add("ED", [f"ED:{cols[21].split(',', 1)[0].split(';', 1)[0].strip()}"])
    return result


# --------------------------------------------------------------------------- locations
def _location(row: GffRow):
    if row.strand not in ("+", "-"):
        raise FrontendError(f"GFF line {row.line_no}: strand must be + or -, got {row.strand!r}.")
    strand = 1 if row.strand == "+" else -1
    start = BeforePosition(row.start) if row.start_partial else ExactPosition(row.start)
    end = AfterPosition(row.end) if row.end_partial else ExactPosition(row.end)
    return SimpleLocation(start, end, strand=strand)


def _revcomp(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def reconstruct_join(loc: SimpleLocation, nt: str, contig: str, circular: bool):
    """Rebuild a collapsed ``join(...)`` from the CDS nucleotide sequence; None if impossible.

    Handles ribosomal slippage (parts overlapping by one base) and features spanning the origin of
    a circular contig (``join(a..N, 1..b)``).
    """
    start, end, strand = int(loc.start), int(loc.end), loc.strand
    plus = nt if strand == 1 else _revcomp(nt)
    span = contig[start:end]
    if plus == span:
        return loc
    gap = len(span) - len(plus)  # >0: bases skipped between the parts; <0: parts overlap
    if 0 < abs(gap) <= MAX_SLIPPAGE_GAP and not (circular and start == 0 and end == len(contig)):
        # slippage: contig[a:a+k] + contig[a+k+gap:b]; k is ambiguous inside homopolymers (first one wins)
        for k in range(1, len(plus)):
            if plus[:k] == contig[start:start + k] and plus[k:] == contig[start + k + gap:end]:
                parts = [SimpleLocation(start, start + k, 1), SimpleLocation(start + k + gap, end, 1)]
                break
        else:
            return None
    elif circular and start == 0 and end == len(contig) and len(plus) < len(contig):
        for k in range(1, len(plus)):  # plus = contig[a:] + contig[:b]
            if contig.endswith(plus[:k]) and contig.startswith(plus[k:]):
                parts = [SimpleLocation(len(contig) - k, len(contig), 1),
                         SimpleLocation(0, len(plus) - k, 1)]
                break
        else:
            return None
    else:
        return None
    if strand == -1:
        parts = [SimpleLocation(p.start, p.end, -1) for p in reversed(parts)]
    # keep partial markers on the outer ends
    if loc.start.__class__ is BeforePosition or loc.end.__class__ is AfterPosition:
        low = min(parts, key=lambda p: int(p.start))
        high = max(parts, key=lambda p: int(p.end))
        if loc.start.__class__ is BeforePosition:
            parts[parts.index(low)] = SimpleLocation(BeforePosition(int(low.start)), low.end, low.strand)
        if loc.end.__class__ is AfterPosition:
            parts[parts.index(high)] = SimpleLocation(high.start, AfterPosition(int(high.end)), high.strand)
    return CompoundLocation(parts)


# --------------------------------------------------------------------------- conversion
@dataclass(frozen=True)
class Dropped:
    """Something ``lenient=True`` threw away instead of raising."""
    kind: str  # invalid_annotation | orphan_gene | unused_faa | unused_ffn
    locus_tag: str
    detail: str = ""


def _translation(lt: str, row: GffRow, proteins: dict[str, str]) -> str | None:
    """Protein from the .faa and/or the GFF ``translation`` attribute (they must agree if both exist)."""
    from_faa, from_gff = proteins.get(lt), row.first("translation")
    if row.type != "CDS":
        if from_faa is not None or from_gff is not None:
            raise FrontendError(f"{lt}: a translation is given for a {row.type} feature (only CDS can have one).")
        return None
    if from_faa is not None and from_gff is not None and from_faa != from_gff:
        raise FrontendError(f"{lt}: translation in .faa and in the GFF differ.")
    return from_faa if from_faa is not None else from_gff


def build_from_arx_gff(
    gff, fna, *, faa=None, ffn=None, eggnog=None, identifier=None, strain=None, taxid=None,
    plasmids=(), lenient=False, dropped: list | None = None, **builder_kwargs,
) -> GenBankBuilder:
    """Drive a :class:`GenBankBuilder` with an ``.arx.gff`` and friends (does not write anything).

    Contig order and sequences come from ``fna``; topology from the GFF ``remark`` rows; circular
    contigs become chromosomes unless listed in ``plasmids``. Missing ``identifier``/``strain``/``taxid``
    are taken from the GFF (contig names, ``source`` row ``strain`` and ``db_xref=taxon:``).

    Everything given must be consumed, otherwise this raises (:class:`FrontendError` or the builder's
    error): leftover ``.faa``/``.ffn`` records, translations on non-CDS features, annotation ids that do
    not match their category regex, ``gene`` rows without CDS/RNA child, unknown qualifiers, ...

    ``lenient=True`` is the single, explicit lossy mode: the cases above marked "lenient" (leftover
    faa/ffn records, regex-invalid annotation ids such as PGAP's preliminary ``EC:3.6.5.n1``, orphan
    ``gene`` rows) are dropped with a logged warning instead. Every drop is appended to ``dropped`` (a list
    of :class:`Dropped`) if given, so callers can assert exactly what was lost.
    """
    def drop(kind, lt, detail, message):
        if not lenient:
            raise FrontendError(f"{message} (use lenient=True / --lenient to drop it with a warning)")
        log.warning("dropping %s", message)
        if dropped is not None:
            dropped.append(Dropped(kind, lt, detail))

    rows = read_gff(gff)
    contigs = read_fasta(fna)
    if not contigs:
        raise FrontendError(f"{fna}: no sequences.")
    contig_seq = dict(contigs)
    if len(contig_seq) != len(contigs):
        raise FrontendError(f"{fna}: duplicate contig names.")

    topology = {r.seqid: r.first("topology", "linear") for r in rows if r.type == "remark"}
    source = next((r for r in rows if r.type == "source"), None)
    if identifier is None:
        m = re.fullmatch(r"(.+)_scf\d+", contigs[0][0])
        if not m:
            raise FrontendError(f"Cannot derive the identifier from contig name {contigs[0][0]!r}; pass identifier=.")
        identifier = m.group(1)
    if strain is None and source is not None:
        strain = source.first("strain")
    if taxid is None and source is not None:
        xrefs = [x for x in source.attrs.get("db_xref", []) if x.startswith("taxon:")]
        if len(xrefs) != 1:
            raise FrontendError("Cannot derive taxid: the GFF source row needs exactly one db_xref=taxon:N; pass taxid=.")
        taxid = int(xrefs[0][len("taxon:"):])
    if taxid is None:
        raise FrontendError("taxid is required.")

    b = GenBankBuilder(identifier=identifier, taxid=taxid, strain=strain, **builder_kwargs)
    plasmids = set(plasmids)
    for name, seq in contigs:
        topo = topology.get(name, "linear")
        if topo not in ("linear", "circular"):
            raise FrontendError(f"Contig {name!r}: unknown topology {topo!r}.")
        is_plasmid = name in plasmids or name.rsplit("_scf", 1)[-1] in plasmids
        kind = "plasmid" if is_plasmid else "chromosome" if topo == "circular" else "scaffold"
        b.add_contig(name, seq, topology=topo, kind=kind)

    proteins = read_locus_fasta(faa) if faa else {}
    nucleotides = read_locus_fasta(ffn) if ffn else {}
    sidecar: list[tuple[str, str, str]] = []

    tagged_children = {r.first("locus_tag") for r in rows if r.type in GENE_FEATURE_TYPES}
    used = set()
    for row in rows:
        if row.type == "gene" and row.first("locus_tag") not in tagged_children:
            # PGAP occasionally emits a gene without any CDS/RNA child; the builder only knows gene+child
            lt = row.first("locus_tag")
            drop("orphan_gene", lt, "", f"{lt}: gene feature without CDS/RNA feature")
        if row.type in _SKIPPED_TYPES:
            continue
        if row.seqid not in contig_seq:
            raise FrontendError(f"GFF line {row.line_no}: contig {row.seqid!r} is not in {fna}.")
        loc = _location(row)
        extra = {k: v for k, v in row.attrs.items() if k not in _IGNORED_ATTRIBUTES}
        if "ribosomal_slippage" in extra:
            extra["ribosomal_slippage"] = True
        if row.type in GENE_FEATURE_TYPES:
            lt = row.first("locus_tag")
            if lt is None:
                raise FrontendError(f"GFF line {row.line_no}: {row.type} without locus_tag.")
            used.add(lt)
            if lt in nucleotides:
                new = reconstruct_join(loc, nucleotides[lt], contig_seq[row.seqid],
                                       topology.get(row.seqid) == "circular")
                if new is None:
                    raise FrontendError(
                        f"{lt}: location {row.seqid}:{row.start + 1}..{row.end} does not reproduce the .ffn sequence.")
                loc = new
            if "transl_table" in extra:
                extra["transl_table"] = int(extra["transl_table"][0])
            b.add_gene(
                locus_tag=lt, contig=row.seqid, location=loc, feature_type=row.type,
                product=row.first("product"), gene=row.first("gene"), pseudo="pseudo" in row.attrs,
                translation=_translation(lt, row, proteins),
                codon_start=int(row.first("codon_start")) if "codon_start" in row.attrs else None,
                **extra)
            for ec in row.attrs.get("EC_number", []):
                sidecar.append((lt, "EC", f"EC:{ec}"))
            for go in row.attrs.get("Ontology_term", []):
                sidecar.append((lt, "GO", go))
        elif row.type in UNTAGGED_FEATURE_TYPES:
            b.add_feature(row.seqid, row.type, loc, **extra)
        else:
            raise FrontendError(f"GFF line {row.line_no}: unsupported feature type {row.type!r}.")

    for kind, records in (("unused_faa", proteins), ("unused_ffn", nucleotides)):
        for lt in sorted(set(records) - used):
            drop(kind, lt, "", f"{lt}: {kind[-3:]} record is not used by any GFF CDS/RNA feature")

    def annotate(lt, anno_type, anno_id):
        regex = b.categories.get(anno_type)
        if regex is not None and not regex.match(anno_id):
            drop("invalid_annotation", lt, anno_id,
                 f"{anno_type} annotation {anno_id!r} of {lt} does not match {regex.pattern!r}")
            return
        b.add_annotation(lt, anno_type, anno_id)

    for lt, anno_type, anno_id in sidecar:
        annotate(lt, anno_type, anno_id)
    if eggnog:
        for lt, by_type in sorted(read_emapper(eggnog).items()):
            for anno_type, ids in by_type.items():
                for anno_id in sorted(ids):
                    annotate(lt, anno_type, anno_id)
    return b


def convert(gff, fna, out_dir, **kwargs) -> list[Path]:
    """:func:`build_from_arx_gff` + ``write(out_dir)``; returns the written paths."""
    return build_from_arx_gff(gff, fna, **kwargs).write(out_dir)


def find_sample_files(sample_dir) -> dict[str, Path | None]:
    """Locate the files of an Arx-style sample directory (``1_assembly``, ``2_cds``, ``3_annotation``)."""
    root = Path(sample_dir)

    def one(pattern, required=True):
        hits = sorted(root.glob(pattern))
        if len(hits) > 1 or (required and not hits):
            raise FrontendError(f"Expected exactly one {pattern} in {root}, found {len(hits)}.")
        return hits[0] if hits else None

    return {
        "gff": one("2_cds/*.arx.gff"), "fna": one("1_assembly/*.fna"),
        "faa": one("2_cds/*.faa", False), "ffn": one("2_cds/*.ffn", False),
        "eggnog": one("3_annotation/eggnog/*.emapper.annotations", False),
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="gbk-creator",
        description="Test frontend: build a PGAP-canonical gbk + annotation sidecars from an .arx.gff, "
                    ".fna, optional .faa/.ffn and emapper annotations.")
    p.add_argument("sample_dir", nargs="?", help="Arx-style sample dir (1_assembly/, 2_cds/, 3_annotation/eggnog/)")
    p.add_argument("--gff"), p.add_argument("--fna"), p.add_argument("--faa"), p.add_argument("--ffn")
    p.add_argument("--eggnog", help="*.emapper.annotations (emapper 1.0.3)")
    p.add_argument("--out", required=True, help="output directory (receives {identifier}.gbk, ...)")
    p.add_argument("--identifier"), p.add_argument("--strain")
    p.add_argument("--taxid", type=int)
    p.add_argument("--old-identifier", help="genome.json old_identifier (previous name of this genome)")
    p.add_argument("--restricted", action="store_true", help="organism.json restricted=true")
    p.add_argument("--tag", action="append", default=None, help="organism.json tag (repeatable)")
    p.add_argument("--annotation-regex", help="annotation categories JSON (default: packaged)")
    p.add_argument("--plasmid", action="append", default=[], help="contig name or number that is a plasmid")
    p.add_argument("--ncbi-email", help="default: $NCBI_EMAIL"), p.add_argument("--ncbi-api-key", help="default: $NCBI_API_KEY")
    p.add_argument("--lenient", action="store_true",
                   help="lossy mode: drop (with a warning) leftover faa/ffn records, regex-invalid annotation ids "
                        "and orphan gene rows instead of failing")
    a = p.parse_args(argv)
    logging.basicConfig(level=logging.WARNING, format="warning: %(message)s")

    files = find_sample_files(a.sample_dir) if a.sample_dir else {}
    for key in ("gff", "fna", "faa", "ffn", "eggnog"):
        if getattr(a, key):
            files[key] = getattr(a, key)
    if not files.get("gff") or not files.get("fna"):
        p.error("need a sample_dir or both --gff and --fna")
    try:
        paths = convert(
            files["gff"], files["fna"], a.out, faa=files.get("faa"), ffn=files.get("ffn"),
            eggnog=files.get("eggnog"), identifier=a.identifier, strain=a.strain, taxid=a.taxid,
            plasmids=a.plasmid, lenient=a.lenient,
            old_identifier=a.old_identifier, restricted=a.restricted, tags=a.tag,
            annotation_regex=a.annotation_regex,
            ncbi_email=a.ncbi_email, ncbi_api_key=a.ncbi_api_key)
    except GbkCreatorError as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    for path in paths:
        print(path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
