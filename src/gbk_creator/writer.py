"""Render a validated GenBankBuilder to GenBank text and annotation sidecars."""
from __future__ import annotations

import io
import json
import re
import warnings

from Bio import SeqIO
from Bio.SeqIO.InsdcIO import GenBankWriter
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from .errors import GbkCreatorError

_MONTHS = ("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")
STRUCTURED_COMMENT = "Genome-Annotation-Data"
_RNA_TYPES = ("tRNA", "rRNA", "ncRNA", "tmRNA")


def _division(lineage: tuple[str, ...]) -> str:
    domain = lineage[0] if lineage else ""
    if domain in ("Bacteria", "Archaea"):
        return "BCT"
    if domain == "Viruses":
        return "VRL"
    if domain == "Eukaryota" and any(x in lineage for x in ("Fungi", "Viridiplantae")):
        return "PLN"
    return "UNK"


def genome_counts(builder) -> dict[str, int]:
    genes = list(builder.genes)
    cds = [g for g in genes if g.feature_type == "CDS"]
    count = lambda t: sum(g.feature_type == t for g in genes)
    return {
        "Genes (total)": len(genes),
        "CDSs (total)": len(cds),
        "Genes (coding)": sum(not g.pseudo for g in cds),
        "CDSs (with protein)": sum(g.has_protein for g in cds),
        "Genes (RNA)": sum(g.feature_type in _RNA_TYPES for g in genes),
        "rRNAs": count("rRNA"),
        "tRNAs": count("tRNA"),
        "ncRNAs": count("ncRNA") + count("tmRNA"),  # PGAP counts tmRNA as ncRNA
        "Pseudo Genes (total)": sum(g.pseudo for g in genes),
        "CDSs (without protein)": sum(not g.has_protein for g in cds),
    }


def comment_text(builder) -> str:
    from .builder import package_version
    version = package_version()
    d = builder.build_date
    rows = {
        "Annotation Pipeline": "Abrinca GenBankBuilder",
        "Annotation Software revision": version,
        "Annotation Date": d.strftime("%m/%d/%Y %H:%M:%S"),
        # plain digits: arx-tools does int(...) on these values, "2,033" would not parse
        **{k: str(v) for k, v in genome_counts(builder).items()},
    }
    width = max(len(k) for k in rows)
    lines = [f"Made by Abrinca GenBankBuilder {version}", "", f"##{STRUCTURED_COMMENT}-START##"]
    lines += [f"{k.ljust(width)} :: {v}" for k, v in rows.items()]
    lines.append(f"##{STRUCTURED_COMMENT}-END##")
    return "\n".join(lines)


def _feature(location, ftype, qualifiers) -> SeqFeature:
    return SeqFeature(location, type=ftype, qualifiers={k: list(v) for k, v in qualifiers.items()})


def build_records(builder) -> list[SeqRecord]:
    d = builder.build_date
    date = f"{d.day:02d}-{_MONTHS[d.month - 1]}-{d.year}"
    taxon = builder.taxon
    comment = comment_text(builder)

    units = {c.index: [] for c in builder.contigs}  # (start, order, [features])
    for g in builder.genes:
        units[g.contig.index].append((g.start, g.order, [
            _feature(g.gene_location, "gene", g.gene_qualifiers),
            _feature(g.location, g.feature_type, g.child_qualifiers)]))
    for f in builder._features:
        units[f.contig.index].append((f.start, f.order, [_feature(f.location, f.feature_type, f.qualifiers)]))

    records = []
    for c in builder.contigs:
        rec = SeqRecord(Seq(c.sequence), id=c.name, name=c.name, description=builder.definition(c))
        rec.annotations = {
            "molecule_type": "DNA",
            "topology": c.topology,
            "data_file_division": _division(taxon.lineage),
            "date": date,
            "source": taxon.scientific_name,
            "organism": taxon.scientific_name,
            "taxonomy": list(taxon.lineage),
            "comment": comment,
        }
        rec.features.append(SeqFeature(
            SimpleLocation(0, len(c.sequence), strand=1), type="source",
            qualifiers={"organism": [taxon.scientific_name], "mol_type": ["genomic DNA"],
                        "strain": [builder.strain], "db_xref": [f"taxon:{taxon.taxid}"]}))
        for _, _, feats in sorted(units[c.index], key=lambda u: (u[0], u[1])):
            rec.features.extend(feats)
        records.append(rec)
    return records


class _GenBankWriter(GenBankWriter):
    """Biopython's writer, except that qualifier values are never broken inside a word.

    Biopython wraps a long word (e.g. ``/locus_tag="<50 char identifier>_000001"``) at column 80 and its
    reader then re-joins the pieces *with a space*, corrupting the value. Here such lines simply stay
    longer than 80 columns. ``/translation`` keeps Biopython's mid-word wrapping (its reader removes the
    line breaks without inserting spaces).
    """

    def _write_feature_qualifier(self, key, value=None, quote=None):
        if key == "translation" or value is None:
            return super()._write_feature_qualifier(key, value, quote)
        if isinstance(value, str):
            value = value.replace('"', '""')
        if quote is None:
            quote = not (isinstance(value, int) or key in self.FTQUAL_NO_QUOTE)
        line = f'/{key}="{value}"' if quote else f"/{key}={value}"
        indent = self.QUALIFIER_INDENT_STR
        words = line.split(" ")
        lines, current = [], indent + words[0]
        for word in words[1:]:
            if len(current) + 1 + len(word) <= self.MAX_WIDTH:
                current += " " + word
            else:
                lines.append(current)
                current = indent + word
        lines.append(current)
        self.handle.write("\n".join(lines) + "\n")


def render_genbank(builder) -> str:
    records = build_records(builder)
    buf = io.StringIO()
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # any other Biopython warning means something is not writable
        # Biopython >= 1.88 widens the LOCUS line for long names (valid per GenBank release 229)
        warnings.filterwarnings("ignore", message="Increasing length of locus line")
        _GenBankWriter(buf).write_file(records)
    text = buf.getvalue()
    text = re.sub(r"^KEYWORDS    \.\n", "", text, flags=re.MULTILINE)  # Biopython always emits one
    _verify_roundtrip(builder, records, text)
    return text


def _loc_str(loc) -> str:
    return str(loc)


def _verify_roundtrip(builder, records: list[SeqRecord], text: str) -> None:
    """Re-read the generated text with Biopython and compare it with what we meant to write."""
    try:
        parsed = list(SeqIO.parse(io.StringIO(text), "genbank"))
    except Exception as e:  # pragma: no cover - indicates a bug or unwritable content
        raise GbkCreatorError(f"Generated GenBank text cannot be parsed again: {e}") from e
    if [r.id for r in parsed] != [r.id for r in records]:
        raise GbkCreatorError("Generated GenBank text does not round-trip: contig ids differ.")
    for want, got in zip(records, parsed):
        problems = []
        if str(want.seq).lower() != str(got.seq).lower():
            problems.append("sequence")
        if got.annotations.get("topology") != want.annotations["topology"]:
            problems.append("topology")
        if got.description != want.description:
            problems.append("definition")
        if len(got.features) != len(want.features):
            problems.append(f"feature count ({len(got.features)} != {len(want.features)})")
        else:
            for w, g in zip(want.features, got.features):
                if (w.type, _loc_str(w.location)) != (g.type, _loc_str(g.location)):
                    problems.append(f"feature {w.type} {w.location} read back as {g.type} {g.location}")
                    break
                want_q = {k: ["" if v is None else str(v) for v in vals] for k, vals in w.qualifiers.items()}
                if want_q != dict(g.qualifiers):
                    bad = sorted(k for k in set(want_q) | set(g.qualifiers) if want_q.get(k) != g.qualifiers.get(k))
                    problems.append(f"qualifiers {bad} of {w.type} {w.location} are altered when read back")
                    break
        if problems:
            raise GbkCreatorError(f"Generated GenBank record {want.id!r} does not round-trip: {problems[:3]}")


def render_sidecars(builder) -> dict[str, str]:
    """``{TYPE: file text}``: ``locus_tag<TAB>anno1, anno2`` per annotated gene, no header, no blank lines."""
    out = {}
    for anno_type in sorted(builder._annotations):
        by_gene = builder._annotations[anno_type]
        lines = [f"{lt}\t{', '.join(sorted(by_gene[lt]))}\n" for lt in sorted(by_gene) if by_gene[lt]]
        if lines:
            out[anno_type] = "".join(lines)
    return out


def render_organism_json(builder) -> str:
    """Arx organism.json (arx-tools merges it over the organism's existing one; it sets `representative` itself)."""
    data = {
        "name": builder.strain,
        "taxid": builder.taxid,
        "restricted": builder.restricted,
    }
    if builder.tags is not None:
        data["tags"] = list(builder.tags)
    return json.dumps(data, indent=4, ensure_ascii=False) + "\n"


MAX_CDS_TOOL_LENGTH = 50  # Arx Genome.cds_tool CharField(max_length=50)


def render_genome_json(builder) -> str:
    """Arx genome.json (merged by arx-tools over the metadata it gathers; only the keys we know)."""
    from .builder import package_version
    cds_tool = f"Abrinca GenBankBuilder v{package_version()}"
    if len(cds_tool) > MAX_CDS_TOOL_LENGTH:
        raise GbkCreatorError(f"cds_tool {cds_tool!r} exceeds Arx's {MAX_CDS_TOOL_LENGTH} character limit.")
    data = {"cds_tool": cds_tool}
    if builder.old_identifier is not None:
        data["old_identifier"] = builder.old_identifier
    return json.dumps(data, indent=4, ensure_ascii=False) + "\n"
