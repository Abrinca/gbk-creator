#!/usr/bin/env python
"""Export one genome from database-shaped data into an Arx-importable folder with gbk-creator.

Part 1 (``fake_database``) is a stand-in for your SQL queries: it returns plain lists of dicts shaped like
table rows (1-based, end-inclusive coordinates, your own IDs) and invents a tiny genome so that the sequences,
coordinates and proteins agree. Replace it with real queries.

Part 2 (``export_genome``) is the actual conversion.

The organism name and lineage are looked up from the taxid at NCBI, so NCBI_EMAIL must be set (NCBI_API_KEY is
optional); answers are cached on disk.

    NCBI_EMAIL=you@example.org python examples/from_database.py [OUT_DIR]
"""
import random
import sys
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation, SimpleLocation

from gbk_creator import GenBankBuilder

IDENTIFIER = "ACME0042-1.1"  # new Arx genome name; must start with the strain
STRAIN = "ACME0042"


# =====================================================================================================
# Part 1: stand-in for the database
# =====================================================================================================
def fake_database():
    """Return (contigs, genes, other_features, annotations) as lists of dicts, like SQL rows."""
    rng = random.Random(42)
    codons = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in ("TAA", "TAG", "TGA")]

    def orf(n_codons, start=True, stop=True):
        return ("ATG" if start else "") + "".join(rng.choice(codons) for _ in range(n_codons)) + ("TAA" if stop else "")

    def noise(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    contigs = [
        {"contig_id": "ctg_A", "sequence": list(noise(3000)), "circular": False, "replicon_type": "scaffold"},
        {"contig_id": "ctg_B", "sequence": list(noise(1200)), "circular": True, "replicon_type": "plasmid"},
    ]
    by_id = {c["contig_id"]: c for c in contigs}

    # (gene_id, contig, start, end, strand, type, product, gene_name, extras); parts = join in transcription order
    rows = [
        ("gene_00031", "ctg_B", 1, 1200, "+", "CDS", "replication initiator protein", None,
         dict(parts=[(1121, 1200), (1, 100)], seq=orf(58))),                       # crosses the origin
        ("gene_00032", "ctg_B", 301, 600, "+", "CDS", "conjugal transfer protein", None, dict(seq=orf(98))),
        ("gene_00007", "ctg_A", 1, 90, "+", "CDS", "hypothetical protein", None,
         dict(start_partial=True, codon_start=2, seq="C" + orf(29, start=False, stop=False) + "GT")),  # incomplete 5' end
        ("gene_00002", "ctg_A", 101, 400, "+", "CDS", "DNA polymerase III subunit beta", "dnaN", dict(seq=orf(98))),
        ("gene_00003", "ctg_A", 501, 800, "-", "CDS", "ABC transporter ATP-binding protein", None, dict(seq=orf(98))),
        ("gene_00004", "ctg_A", 901, 1100, "+", "CDS", "IS3 family transposase", None,
         dict(is_pseudo=True, seq=noise(200))),                                    # broken gene: no protein
        ("gene_00005", "ctg_A", 1201, 1277, "+", "tRNA", "tRNA-Gly", None,
         dict(anticodon="(pos:1235..1237,aa:Gly,seq:tcc)", seq=noise(77))),
        ("gene_00006", "ctg_A", 1401, 1700, "+", "rRNA", "16S ribosomal RNA", None, dict(seq=noise(300))),
        ("gene_00008", "ctg_A", 1801, 2105, "+", "CDS", "IS3 family transposase", None,
         dict(parts=[(1801, 1900), (1903, 2105)], slippage=True, seq=orf(99))),   # two parts, 2 bases skipped
        ("gene_00009", "ctg_A", 2701, 3000, "+", "CDS", "hypothetical protein", None,
         dict(end_partial=True, seq=orf(99, stop=False))),                         # incomplete 3' end
    ]
    genes = []
    for gene_id, contig, start, end, strand, ftype, product, name, x in rows:
        parts = x.get("parts") or [(start, end)]
        seq = x["seq"]
        pos = 0
        for p_start, p_end in parts:  # implant the gene so that the coordinates really contain it
            chunk = seq[pos:pos + p_end - p_start + 1]
            pos += len(chunk)
            chunk = chunk if strand == "+" else str(Seq(chunk).reverse_complement())
            by_id[contig]["sequence"][p_start - 1:p_end] = chunk
        assert pos == len(seq), gene_id
        cs = x.get("codon_start", 1)
        usable = seq[cs - 1:len(seq) - (len(seq) - cs + 1) % 3]
        protein = None
        if ftype == "CDS" and not x.get("is_pseudo"):
            protein = str(Seq(usable).translate(table=11)).rstrip("*")
        genes.append({
            "gene_id": gene_id, "contig": contig, "start": start, "end": end, "strand": strand,
            "feature_type": ftype, "product": product, "gene_name": name,
            "parts": x.get("parts"), "start_partial": x.get("start_partial", False),
            "end_partial": x.get("end_partial", False), "codon_start": cs,
            "is_pseudo": x.get("is_pseudo", False), "protein": protein,
            "anticodon": x.get("anticodon"), "ribosomal_slippage": x.get("slippage", False),
        })
    for c in contigs:
        c["sequence"] = "".join(c["sequence"])
    other_features = [
        {"contig": "ctg_A", "start": 2401, "end": 2500, "strand": "+", "feature_type": "repeat_region",
         "qualifiers": {"rpt_family": "CRISPR", "rpt_type": "direct", "rpt_unit_range": "2401..2430",
                        "rpt_unit_seq": "gttgtacagttacttaaatcttgagagtac"}},
    ]
    annotations = [  # (your gene id, category, value); duplicates are fine, they are ignored
        {"gene_id": "gene_00002", "type": "KG", "value": "K02338"},
        {"gene_id": "gene_00002", "type": "EC", "value": "EC:2.7.7.7"},
        {"gene_id": "gene_00002", "type": "KG", "value": "K02338"},
        {"gene_id": "gene_00003", "type": "KG", "value": "K02056"},
        {"gene_id": "gene_00003", "type": "GO", "value": "GO:0005524"},
        {"gene_id": "gene_00004", "type": "EC", "value": "EC:2.7.7.-"},   # pseudo genes may carry annotations
        {"gene_id": "gene_00032", "type": "KG", "value": "K03196"},
    ]
    return contigs, genes, other_features, annotations


# =====================================================================================================
# Part 2: the conversion (this is what you adapt)
# =====================================================================================================
# --- BEGIN excerpt (embedded in README.md; tests keep both in sync)
def to_location(row):
    """DB columns -> Biopython location. DB: 1-based, end-inclusive. Biopython: 0-based, end-exclusive."""
    strand = 1 if row["strand"] == "+" else -1
    parts = row.get("parts") or [(row["start"], row["end"])]  # a join, listed in transcription order
    first_low = min(start for start, _ in parts)
    last_high = max(end for _, end in parts)
    locations = []
    for start, end in parts:
        begin = BeforePosition(start - 1) if row.get("start_partial") and start == first_low else start - 1
        finish = AfterPosition(end) if row.get("end_partial") and end == last_high else end
        locations.append(SimpleLocation(begin, finish, strand))
    return locations[0] if len(locations) == 1 else CompoundLocation(locations)


def export_genome(contigs, genes, other_features, annotations, out_dir):
    builder = GenBankBuilder(
        identifier=IDENTIFIER, taxid=1308, strain=STRAIN,
        old_identifier="legacy_genome_42",  # keeps the link to your old genome name
        tags=["metagenomic-sample-42"],     # optional organism tags in Arx
    )

    # 1. Contigs are renamed {identifier}_scf1, _scf2, ... in the order they are added.
    contig_name = {}
    for n, contig in enumerate(contigs, start=1):
        contig_name[contig["contig_id"]] = f"{IDENTIFIER}_scf{n}"
        builder.add_contig(
            contig_name[contig["contig_id"]], contig["sequence"],
            topology="circular" if contig["circular"] else "linear", kind=contig["replicon_type"])

    # 2. Locus tags must increase along contig order, then position: renumber, keep your IDs as old_locus_tag.
    contig_rank = {c["contig_id"]: n for n, c in enumerate(contigs)}
    lowest = lambda g: min(start for start, _ in (g["parts"] or [(g["start"], g["end"])]))
    new_tag = {}
    for n, g in enumerate(sorted(genes, key=lambda g: (contig_rank[g["contig"]], lowest(g))), start=1):
        new_tag[g["gene_id"]] = f"{IDENTIFIER}_{n:06d}"
        extra = {}
        if g["anticodon"]:
            extra["anticodon"] = g["anticodon"]
        if g["ribosomal_slippage"]:
            extra["ribosomal_slippage"] = True
        builder.add_gene(
            locus_tag=new_tag[g["gene_id"]], contig=contig_name[g["contig"]], location=to_location(g),
            feature_type=g["feature_type"], product=g["product"], gene=g["gene_name"],
            pseudo=g["is_pseudo"], translation=g["protein"],
            codon_start=g["codon_start"] if g["feature_type"] == "CDS" else None,
            old_locus_tag=g["gene_id"], **extra)

    # 3. Features without a locus tag (CRISPR arrays, riboswitches, ...)
    for f in other_features:
        builder.add_feature(contig_name[f["contig"]], f["feature_type"], to_location(f), **f["qualifiers"])

    # 4. Annotations (KEGG, EC, GO, ...): use the new locus tags
    for a in annotations:
        builder.add_annotation(new_tag[a["gene_id"]], a["type"], a["value"])

    return builder.write(out_dir)  # checks everything first, then writes the files
# --- END excerpt


def main(argv):
    out_dir = Path(argv[1]) if len(argv) > 1 else Path("out") / IDENTIFIER
    paths = export_genome(*fake_database(), out_dir)
    for path in paths:
        print(path)
    print("\nNext: import it into Arx with arx-tools:\n"
          f"  import_genome {out_dir} --folder_structure_dir /path/to/folder_structure "
          f"--organism {STRAIN} --genome {IDENTIFIER}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
