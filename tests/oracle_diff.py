"""Semantic diff between a PGAP oracle gbk and a gbk produced by gbk-creator.

Both files are parsed with Biopython and compared per contig: id, length, sequence, topology and every
feature (type, location, qualifiers). Differences are returned as human-readable strings; an empty list
means "semantically identical up to the ALLOWLIST below".

ALLOWLIST (every entry is a known, accepted difference; keep it tight):

* Header (not compared at all): DEFINITION, KEYWORDS, REFERENCE, COMMENT, LOCUS date, ORGANISM lineage,
  dbxrefs. Reason: gbk-creator writes its own header (current NCBI lineage, build date, own COMMENT, no
  KEYWORDS/REFERENCE) - see Plan.md "Header".
* ``EC_number`` on CDS features: oracle only (stripped from the oracle side only; the output must not have
  it). EC numbers go to the ``.EC`` sidecar instead (Plan.md
  "No /EC_number or GO in the gbk"); the sidecar is verified separately.
* ``join()`` locations: the GFF collapses joins; the frontend rebuilds them from the .ffn. The rebuilt
  parts can differ from the oracle's inside homopolymers (the junction is ambiguous), so a differing
  location is accepted iff both are compound with the same operator, part count, strands, fuzzy outer
  ends (</>) and overall span, and the extracted nucleotide sequence is identical; only internal junction
  coordinates may differ.
* ``source`` /organism (only with ``ignore_source_organism=True``, used by the NCBI-backed network test): the
  oracle carries the PGAP-era name (e.g. "Lactobacillus casei"), gbk-creator the current NCBI scientific name
  ("Lacticaseibacillus casei"); db_xref taxon, strain and mol_type are still compared. The offline tests
  rebuild the Taxon from the oracle header and therefore compare /organism too.
* Oracle-only ``gene`` features that have no CDS/RNA child (a PGAP quirk; gbk-creator only knows gene+child).
"""
from __future__ import annotations

from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation, ExactPosition

ALLOWED_ORACLE_ONLY_QUALIFIERS = {"CDS": {"EC_number"}}


def _read(path):
    with open(path) as f:
        return list(SeqIO.parse(f, "genbank"))


def _key(rec, feat, n_untagged):
    lt = feat.qualifiers.get("locus_tag", [None])[0]
    if lt is not None:
        return (feat.type, lt)
    return (feat.type, "untagged", int(feat.location.start), int(feat.location.end), feat.location.strand,
            n_untagged)


def _join_difference(orec, ofeat, brec, bfeat) -> str | None:
    """None if two differing locations are the same join up to *internal junction coordinates*, else why not.

    Required: both compound, same operator, same number of parts, same strand, same outer ends including
    their fuzziness (``<``/``>``), no fuzzy positions elsewhere, same overall span and the same extracted
    nucleotide sequence.
    """
    la, lb = ofeat.location, bfeat.location
    if not (isinstance(la, CompoundLocation) and isinstance(lb, CompoundLocation)):
        return "not both compound locations"
    if la.operator != lb.operator:
        return f"operator {la.operator} != {lb.operator}"
    if len(la.parts) != len(lb.parts):
        return f"{len(la.parts)} parts != {len(lb.parts)} parts"
    if la.strand != lb.strand or [p.strand for p in la.parts] != [p.strand for p in lb.parts]:
        return "strand differs"
    if (int(la.start), int(la.end)) != (int(lb.start), int(lb.end)):
        return "overall span differs"
    outer = lambda loc: (type(min(loc.parts, key=lambda p: int(p.start)).start),
                         type(max(loc.parts, key=lambda p: int(p.end)).end))
    if outer(la) != outer(lb):
        return "fuzzy (</>) outer ends differ"
    for loc in (la, lb):
        low = min(loc.parts, key=lambda p: int(p.start))
        high = max(loc.parts, key=lambda p: int(p.end))
        for p in loc.parts:
            if (p is not low and type(p.start) is not ExactPosition) or (p is not high and type(p.end) is not ExactPosition):
                return "fuzzy position at an internal junction"
    if _extract(orec, ofeat) != _extract(brec, bfeat):
        return "extracted sequence differs"
    return None


def _extract(rec, feat) -> str:
    return str(feat.location.extract(rec.seq)).upper()


def compare_gbk(oracle_path, built_path, stats: dict | None = None, ignore_source_organism: bool = False) -> list[str]:
    oracle, built = _read(oracle_path), _read(built_path)
    diffs: list[str] = []
    if [r.id for r in oracle] != [r.id for r in built]:
        return [f"contig ids differ: {[r.id for r in oracle][:5]} vs {[r.id for r in built][:5]}"]
    stats = {} if stats is None else stats
    for o, b in zip(oracle, built):
        where = o.id
        if len(o) != len(b) or str(o.seq).upper() != str(b.seq).upper():
            diffs.append(f"{where}: sequence differs")
        if o.annotations.get("topology") != b.annotations.get("topology"):
            diffs.append(f"{where}: topology {o.annotations.get('topology')} != {b.annotations.get('topology')}")
        if o.annotations.get("molecule_type") != b.annotations.get("molecule_type"):
            diffs.append(f"{where}: molecule_type differs")
        if o.annotations.get("data_file_division") != b.annotations.get("data_file_division"):
            diffs.append(f"{where}: division differs")

        def index(rec):
            out, n = {}, 0
            for f in rec.features:
                k = _key(rec, f, n)
                if k[1] == "untagged":
                    n += 1
                if k in out:
                    diffs.append(f"{where}: duplicate feature key {k}")
                out[k] = f
            return out

        fo, fb = index(o), index(b)
        children = {k[1] for k in fo if k[0] not in ("gene", "source") and k[1] != "untagged"}
        for k in [k for k in fo if k[0] == "gene" and k[1] not in children]:
            if k not in fb:
                del fo[k]
                stats["orphan_genes_ignored"] = stats.get("orphan_genes_ignored", 0) + 1
        if set(fo) != set(fb):
            diffs.append(f"{where}: feature sets differ: oracle-only {sorted(set(fo) - set(fb))[:3]}, "
                         f"built-only {sorted(set(fb) - set(fo))[:3]}")
        # same order along the contig (gene followed by child, sorted by position)
        if [k[:2] if k[1] != "untagged" else k[:4] for k in fo] != [k[:2] if k[1] != "untagged" else k[:4] for k in fb]:
            diffs.append(f"{where}: feature order differs")
        for k in set(fo) & set(fb):
            a, c = fo[k], fb[k]
            label = f"{where} {k[0]} {k[1]}"
            if str(a.location) != str(c.location):
                problem = _join_difference(o, a, b, c)
                if problem is None:
                    stats["joins_rebuilt_differently"] = stats.get("joins_rebuilt_differently", 0) + 1
                else:
                    diffs.append(f"{label}: location {a.location} != {c.location} ({problem})")
            # oracle-only qualifiers are stripped from the oracle side ONLY: output containing them is reported
            drop_both = {"organism"} if a.type == "source" and ignore_source_organism else set()
            drop_oracle = drop_both | ALLOWED_ORACLE_ONLY_QUALIFIERS.get(a.type, set())
            qa = {q: v for q, v in a.qualifiers.items() if q not in drop_oracle}
            qc = {q: v for q, v in c.qualifiers.items() if q not in drop_both}
            if qa != qc:
                for q in sorted(set(qa) | set(qc)):
                    if qa.get(q) != qc.get(q):
                        diffs.append(f"{label}: qualifier /{q} {qa.get(q)!r} != {qc.get(q)!r}")
    return diffs


def oracle_ec(oracle_path) -> dict[str, set[str]]:
    """``{locus_tag: {'EC:x.y.z.w', ...}}`` from the oracle's /EC_number qualifiers."""
    out: dict[str, set[str]] = {}
    for rec in _read(oracle_path):
        for f in rec.features:
            if f.type == "CDS" and "EC_number" in f.qualifiers:
                out.setdefault(f.qualifiers["locus_tag"][0], set()).update(f"EC:{e}" for e in f.qualifiers["EC_number"])
    return out
