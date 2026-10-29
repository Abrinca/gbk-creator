# gbk-creator

gbk-creator turns genome data that lives in your own database into the files Arx imports: one GenBank file
(`.gbk`, the standard text format for sequences plus gene annotations) and a few small metadata files. You hand it
contigs, genes and annotations from Python; it checks everything against the naming rules Arx needs and writes
the folder for you. Anything it cannot write faithfully raises an error instead of producing a subtly wrong file.

```
out/ACME0042-1.1/
├── ACME0042-1.1.gbk      # contigs, genes, proteins
├── ACME0042-1.1.KG       # one small text file per annotation type (KEGG, EC, GO, ...):
├── ACME0042-1.1.EC       #   "locus_tag<TAB>annotation1, annotation2" per annotated gene
├── organism.json         # strain name, taxid, restricted flag, optional tags
└── genome.json           # annotation tool, optional old genome name
```

## Install

```sh
pip install git+https://github.com/Abrinca/gbk-creator
```

Requires Python 3.12+. The organism name and lineage in the GenBank header come from NCBI Taxonomy, looked up from
your taxid. Provide an e-mail address (NCBI asks for one) and, optionally, an API key (raises the rate limit; create one
under *Account settings* on ncbi.nlm.nih.gov):

```sh
export NCBI_EMAIL=you@example.org
export NCBI_API_KEY=...        # optional
```

or pass `ncbi_email=` / `ncbi_api_key=` to `GenBankBuilder`. Answers are cached in `~/.cache/gbk-creator`
(override with `GBK_CREATOR_CACHE`; delete the folder to refresh). The library never reads `.env` files and never stores
your key.


## A complete example

[`examples/from_database.py`](examples/from_database.py) exports a tiny invented genome. Plain lists of dicts stand
in for your SQL rows (contigs, genes with 1-based coordinates and your own IDs, other features, annotations); the
conversion renumbers genes in position order, builds, writes and prints the follow-up `import_genome` command.
It covers a normal CDS, a minus-strand CDS, a pseudo gene, a tRNA with anticodon, an rRNA, an incomplete CDS at
either end, a join, a gene across the origin of a plasmid, a CRISPR repeat and annotations.

```sh
export NCBI_EMAIL=you@example.org        # see Install
python examples/from_database.py out/ACME0042-1.1
```

The conversion part (replace `fake_database()` in the script with your own queries):

```python
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
```

## The rules

Arx expects a fixed naming scheme, so **you have to renumber contigs and genes** when you export. Keep your own IDs
with `old_locus_tag=` (per gene) and `old_identifier=` (per genome); Arx then shows them next to the new names.

| What | Rule |
|---|---|
| `identifier` (the new Arx genome name) | at most 50 characters; no spaces and none of `: / \ \|` |
| `strain` (Arx "organism") | at most 40 characters; the identifier must **start with** the strain (default: strain = identifier) |
| Contig names | `{identifier}_scf1`, `{identifier}_scf2`, ... with no leading zeros, **added in exactly that order** |
| Contig `kind` / `topology` | `scaffold`, `chromosome` or `plasmid`; `linear` or `circular`; a chromosome must be circular |
| Locus tags | `{identifier}_000001`: exactly 6 digits, unique, **increasing by contig order, then by start position** (gaps are fine) |
| Old names | `old_locus_tag=` on a gene (unique, at most 100 characters, no `"`); `old_identifier=` on the builder (at most 50 characters, different from the identifier) |
| Strand | `+1` or `-1` on every part of a location |
| Genes | every gene needs a `product` and lies inside its contig; a CDS needs a `translation` unless it is `pseudo` (then it must not have one). Sequences are trusted, not re-translated |
| Annotations | only on CDS genes (pseudo is fine, tRNA/rRNA/... are not); the value must match the category's regex; duplicates are ignored |
| Forbidden | `/EC_number` and GO terms in the gbk (use annotations instead), `GC`/`GP` annotation files (Arx derives them from `gene` and `product`) |

## Coordinates

> **Biopython coordinates are 0-based and end-exclusive. Databases and GenBank text are 1-based and end-inclusive.**
> A feature stored as `start..end` becomes `SimpleLocation(start - 1, end, strand)`.

| Database (1-based) | Python | Written to the gbk |
|---|---|---|
| `101..400`, `+` | `SimpleLocation(100, 400, +1)` | `101..400` |
| `501..800`, `-` | `SimpleLocation(500, 800, -1)` | `complement(501..800)` |
| `1..90`, start incomplete | `SimpleLocation(BeforePosition(0), 90, +1)` | `<1..90` |
| `2701..3000`, end incomplete | `SimpleLocation(2700, AfterPosition(3000), +1)` | `2701..>3000` |
| `1..90`, `-`, low end incomplete | `SimpleLocation(BeforePosition(0), 90, -1)` | `complement(<1..90)` |
| join: `1801..1900` + `1903..2105`, `+` | `CompoundLocation([SimpleLocation(1800, 1900, +1), SimpleLocation(1902, 2105, +1)])` | `join(1801..1900,1903..2105)` |
| same join, `-` (parts in reading order, i.e. high to low) | `CompoundLocation([SimpleLocation(1902, 2105, -1), SimpleLocation(1800, 1900, -1)])` | `complement(join(1801..1900,1903..2105))` |
| across the origin of a circular contig of 1200 bp, `+`: `1121..1200` then `1..100` | `CompoundLocation([SimpleLocation(1120, 1200, +1), SimpleLocation(0, 100, +1)])` | `join(1121..1200,1..100)` |
| same, `-` | `CompoundLocation([SimpleLocation(0, 100, -1), SimpleLocation(1120, 1200, -1)])` | `complement(join(1121..1200,1..100))` |

`<` and `>` always sit at the **lower** and **higher** coordinate, whatever the strand. For a join, list the parts in the
order the gene is read (5' to 3' on its own strand). The `gene` record gbk-creator adds for every CDS/RNA uses the same
location, except for ordinary joins (not across the origin), where it spans all parts.


## Reference

### Feature types

`add_gene(locus_tag, contig, location, feature_type, product, gene=None, pseudo=False, translation=None,
codon_start=None, old_locus_tag=None, **qualifiers)` writes a `gene` record plus one child record. The keyword
`qualifiers` may only be the ones in the table (they are the ones the reference annotations use); anything else raises.
`gene` must be 3-11 characters from `[0-9a-zA-Z_/' -]` (the same regex Arx applies to gene names).

| Feature type | Extra qualifiers (keyword arguments) |
|---|---|
| `CDS` | `codon_start`, `db_xref`, `gene_synonym`, `inference`, `note`, `protein_id`, `ribosomal_slippage`, `transl_table` |
| `tRNA` | `anticodon`, `db_xref`, `gene_synonym`, `inference`, `note` |
| `rRNA` | `db_xref`, `gene_synonym`, `inference`, `note` |
| `ncRNA` | `db_xref`, `gene_synonym`, `inference`, `ncRNA_class`, `note` |
| `tmRNA` | `db_xref`, `gene_synonym`, `inference`, `note` |

* `transl_table` defaults to the builder's `transl_table` (11); `codon_start` defaults to 1 for a CDS.
* `db_xref` must look like `RFAM:RF00010`. `ribosomal_slippage` is a flag: `ribosomal_slippage=True`.
* Values may be a string or a list of strings (e.g. several `inference` lines); single line only.

`add_feature(contig, feature_type, location, **qualifiers)` adds a record without locus tag:

| Feature type | Qualifiers |
|---|---|
| `repeat_region` | `inference`, `note`, `rpt_family`, `rpt_type`, `rpt_unit_range`, `rpt_unit_seq` |
| `regulatory` | `bound_moiety`, `db_xref`, `inference`, `note`, `regulatory_class` |
| `misc_feature` | `note` |

`regulatory` requires `regulatory_class`; `misc_feature` requires `note`.

**Pseudo genes** (`pseudo=True`, CDS only) get `/pseudo` on both records and no translation.

### Annotations

```python
builder.add_annotation(locus_tag, "KG", "K00001")   # file {identifier}.KG, one line per annotated gene
```

Default categories (`annotation_regex.json`; pass your own file as `annotation_regex="my.json"`):

| Type | Content | Regex |
|---|---|---|
| `EC` | enzyme number | `^EC:[0-9\.-]{1,12}$` |
| `KG` | KEGG ortholog | `^K[0-9]{5}$` |
| `KR` | KEGG reaction | `^R[0-9]{5}$` |
| `GO` | Gene Ontology term | `^GO:[0-9]{7}$` |
| `OL` | ortholog group | `^((N[0-9]+\.)?H)?OG[0-9]{7}$` |
| `GC` | gene name (derived by Arx, do not write) | `^GC:[0-9a-zA-Z\_\/\-\ \']{3,11}$` |
| `GP` | product (derived by Arx, do not write) | `^GP:.*$` |
| `EP` | eggNOG preferred name | `^EP:.*$` |
| `EO` | eggNOG ortholog group | `^EO:.*$` |
| `ED` | eggNOG description | `^ED:.*$` |
| `SL` | subcellular location | `^SL-[0-9]{4}$` |
| `AR` | free-form annotation | `^AR:.*$` |

`GC` (gene name) and `GP` (product) must not be written: Arx derives them from `gene` and `product`.
The JSON may be a plain `{"KG": "^K[0-9]{5}$"}` mapping or Arx's `annotations.json` format (`{"KG": {"regex": ...}}`).
Category names are **exactly two upper-case letters** (Arx stores them in a 2-character field).

**Adding your own category (say `ZZ`):**

1. Put it in your annotation JSON (`{"ZZ": "^ZZ[0-9]{3}$"}`) and pass it as `annotation_regex=`; `add_annotation(tag, "ZZ", "ZZ001")` then writes `{identifier}.ZZ`.
2. Add the category to **Arx's own `annotations.json`** (in the Arx folder structure), otherwise Arx rejects the file.
3. arx-tools only looks for a fixed list of types (`GC GP EP ED EO EC KG KR GO SL OL`). Point `ARX_IMPORT_SETTINGS`
   (or `import_genome --import_settings`) at a JSON file that repeats those default entries and adds yours:
   ```json
   {"file_finder": {"custom_annotations": [
       {"glob": "*.GC", "anno_type": "GC", "expected": false},
       "... one entry per default type (GP EP ED EO EC KG KR GO SL OL) ...",
       {"glob": "*.ZZ", "anno_type": "ZZ", "expected": false}
   ]}}
   ```

### Taxonomy helpers

```python
from gbk_creator import name_to_taxid, fetch_taxon
name_to_taxid("Streptococcus thermophilus")   # 1308; raises unless NCBI finds exactly one entry
fetch_taxon(1308)                              # Taxon(taxid=1308, scientific_name=..., rank='species', lineage=(...))
```

### Other builder options

`transl_table=11`, `restricted=False` and `tags=[...]` (both go to `organism.json`; `tags` is left out when `None`),
`old_identifier=...` (goes to `genome.json`), `build_date=` (default: now).

## Errors

Every problem raises `gbk_creator.GbkCreatorError` (a `ValueError`; subclasses `NamingError`, `FeatureError`,
`AnnotationError`, `TaxonomyError`) with a message that says what to change. Nothing is written when anything is wrong.

| Message (shortened) | Fix |
|---|---|
| `Contig name 'contig_7' is invalid: ... Expected 'ACME-1.1_scf2'` | name contigs `{identifier}_scf1`, `_scf2`, ... in the order you add them |
| `locus tags must be 'ACME-1.1_NNNNNN' with exactly 6 digits` | format the new tag as `f"{identifier}_{n:06d}"`; keep your ID in `old_locus_tag=` |
| `Locus tags must increase strictly along contig order and position: ..._000001 (...:201) comes after ..._000002 (...:101)` | number genes after sorting by contig order, then start position |
| `strand must be +1 or -1 (got 0)` | every part of a location needs a strand |
| `location part 0..3000 (0-based, end exclusive) is outside contig ... of length 300` | check the `start - 1` conversion and the contig sequence |
| `qualifier /EC_number is not allowed here: ... use add_annotation(...)` | EC numbers and GO terms are annotations, not qualifiers |
| `/gene 'x' does not match the GC regex` | gene names are 3-11 characters; leave `gene=None` otherwise |
| `unknown annotation category 'XX'. ... add it to your annotation JSON file` | see "Adding your own category" |
| `does not match the regex of category KG` | fix the value (`K` + 5 digits) |

`builder.validate()` runs all the whole-genome checks without writing; `builder.write(out_dir)` calls it first.

## Importing into Arx

```sh
import_genome OUT_DIR --folder_structure_dir /path/to/folder_structure --organism ACME0042 --genome ACME0042-1.1
```

`--organism` is your `strain`, `--genome` your `identifier`. arx-tools finds the `.gbk` and annotation files, creates
the sequence, protein and gene files itself and merges the two json files:

* `organism.json`: `name` (= strain), `taxid`, `restricted` and `tags` (if given). **Its values override an existing
  organism's `organism.json`** (including `restricted: false`); omitted tags leave existing ones alone.
* `genome.json`: `cds_tool` = `Abrinca GenBankBuilder v<version>` and `old_identifier` if you set it.

arx-tools will log `Failed to discover annotation information from ...gbk`; this is expected and harmless (it only looks
for the annotation tools it knows).

## Development

See [DEVELOPMENT.md](DEVELOPMENT.md).
