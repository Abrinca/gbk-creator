# gbk-creator — Plan

## Goal
A small Python **library** (`gbk_creator.GenBankBuilder`) that lets an Arx customer turn genome data
(stored in their own Postgres, etc.) into a **"PGAP-canonical" GenBank file + custom-annotation sidecars**
that import cleanly into Arx via `arx-tools import_genome`.

- Primary interface: **Python API**. Parsers (gff/faa/ffn/eggnog) exist only as a test frontend.
- Output: `out/{identifier}/{identifier}.gbk` + `out/{identifier}/{identifier}.{TYPE}` sidecars +
  `out/{identifier}/organism.json` (`name`=strain, `taxid`, `restricted`, optional `tags`; never `representative`,
  arx-tools sets that itself) + `out/{identifier}/genome.json` (`cds_tool`="Abrinca GenBankBuilder v{version}",
  optional `old_identifier`). Nothing else (arx-tools generates fna/faa/ffn/gff and the rest of genome.json).
  arx-tools merges both json files found in the import dir over what it gathers/has (organism.json over the
  organism's existing one, genome.json after its own gbk detection): our non-null values win.
- Reference/oracle: `./data/{sample}/2_cds/{sample}.gbk` (PGAP 2019/2020 output). `data/` is gitignored; tests skip when absent.
- Dependencies managed by uv. `requires-python >=3.12`. Runtime deps: biopython (+ stdlib HTTP for NCBI E-utilities).
- The library is standalone: **no arx-tools dependency at runtime**; arx-tools is used in the test suite as an acceptance check.

## Core rule
**Fail loudly.** Anything that is consumed but cannot be written out again, or that violates the rules below, raises an error.
Users of the library are trusted on sequence content (no re-translation check), but not on structure/naming.

## API sketch
```python
b = GenBankBuilder(
    identifier="FAM13496-i1-1.1",   # max 50 chars, no space : / \
    taxid=1308,                     # any valid NCBI taxid; name + lineage fetched via E-utilities (cached)
    strain="FAM13496",              # = Arx "organism"; optional, default = identifier; identifier must start with it; max 40 chars
    transl_table=11,                # default 11
    annotation_regex="annotation_regex.json",  # optional path; plain {TYPE: regex} OR Arx annotations.json format; default = packaged copy
    ncbi_email=None, ncbi_api_key=None,        # default: env vars NCBI_EMAIL / NCBI_API_KEY; email required
    restricted=False,               # organism.json: real bool
    tags=None,                      # organism.json: list of unique non-empty str; None = key omitted
    old_identifier=None,            # genome.json: previous name of the genome (non-empty line, <=50 chars, != identifier)
)
b.add_contig(name="FAM13496-i1-1.1_scf1",   # user supplies the name (re-use existing numbers); validated, raise otherwise sequence=..., topology="linear"|"circular",
             kind="scaffold"|"chromosome"|"plasmid")
b.add_gene(locus_tag="FAM13496-i1-1.1_000001", contig=..., location=...,   # Biopython-style, strand ±1,
           feature_type="CDS"|"tRNA"|"rRNA"|"ncRNA"|"tmRNA",               # join() and </> partial ends allowed
           product=..., gene=None, pseudo=False, translation=None,
           codon_start=1, **extra_qualifiers)   # only qualifiers that occur in the PGAP oracle
           # -> builder auto-creates the `gene` parent feature (same location, locus_tag, /pseudo if pseudo)
b.add_feature(contig=..., feature_type="repeat_region"|"regulatory", location=..., **qualifiers)  # no locus_tag
b.add_annotation(locus_tag, anno_type="KG", anno_id="K00001")
b.write(out_dir)   # validates everything, then writes gbk + organism.json + genome.json + sidecars
```

## Rules enforced
### Naming
- Contigs: `f'{identifier}_scf{n}'`, names supplied by the user (re-use their numbering); must be n = 1..N contiguous,
  **in order of `add_contig` calls**, no leading zeros. Raise otherwise.
- Locus tags: `f'{identifier}_{n:06d}'`, exactly 6 digits, unique, **increasing along contig order and position**.
- identifier starts with strain (= Arx organism; may be equal). Strain defaults to identifier.

### Features (mirror the PGAP oracle; reject anything else)
- Locus-tagged: `gene` (auto-created by the builder) + one child (`CDS`, `tRNA`, `rRNA`, `ncRNA`, `tmRNA`) sharing the locus_tag.
- Non-tagged: `repeat_region` (CRISPR), `regulatory`.
- Pseudo genes: `/pseudo` on gene + CDS, no translation.
- Every gene must have a nucleotide sequence (= extractable from its contig via its location).
- All protein (translation) and nucleotide sequences must end up in the gbk.
- Qualifiers on genes: `/product`, `/gene` (validated against the GC regex: 3–11 chars), and the PGAP ones
  (`codon_start`, `transl_table`, `translation`, `protein_id`, `inference`, `note`, `anticodon`,
  `ncRNA_class`, `regulatory_class`, `bound_moiety`, `ribosomal_slippage`, `db_xref` (RFAM), `rpt_*`).
- **No `/EC_number` or GO in the gbk.** EC/GO go only to sidecars.

### Header
- Required: taxid, organism/strain, topology per contig.
- `SOURCE` / `ORGANISM` = NCBI scientific name of taxid; ORGANISM lineage fetched from NCBI (current names,
  e.g. Bacillota — will differ from oracle's Firmicutes; accepted).
- `source` feature: `/organism`, `/mol_type="genomic DNA"`, `/strain` (**always**), `/db_xref="taxon:{taxid}"` (exactly one db_xref).
- Auto-generated: LOCUS (date = build date; handle Biopython's long-LOCUS-name limitation),
  DEFINITION (contains "chromosome"/"plasmid" per contig kind so Arx classifies contigs; chromosome requires circular),
  ACCESSION/VERSION = contig name,
  COMMENT = `Made by Abrinca GenBankBuilder {version}` (version via importlib.metadata; do **not** imitate PGAP)
  + a `##Genome-Annotation-Data##` block (Annotation Pipeline :: Abrinca GenBankBuilder, Annotation Software revision :: {version},
  Annotation Date) with the easy computed counts (Genes (total), CDSs (total), Genes (coding), CDSs (with protein),
  Genes (RNA), rRNAs, tRNAs, ncRNAs, Pseudo Genes (total), CDSs (without protein)).
  arx-tools `import_genome` regenerates ffn/faa with `strict=True` and raises if `Genes (total)` / `CDSs (with protein)`
  mismatch (`rename_genbank.py:329,332` -> `genbank_to_fasta.py:121,160-177`) — a free sanity check.
- Not emitted: KEYWORDS, REFERENCE, DBLINK.

## Custom annotation sidecars
- Format (Arx): `locustag\tanno1, anno2, anno3`, one file per type: `{identifier}.{TYPE}` (arx-tools globs `*.KG` etc.).
- Only annotated genes are written; sorted by locus tag; annotations sorted, `", "`-joined; no header, no blank lines
  (Arx stops reading at the first blank line!).
- Categories + regexes from a JSON file (default `annotation_regex.json`, same regexes as Arx's `annotations.json`).
  Users add categories by passing their own JSON path. Both formats accepted: `{TYPE: regex}` and Arx's
  `{TYPE: {name, color, regex, hyperlinks}}`.
- Duplicate annotation (same gene, type, id): silently ignored.
- Annotations allowed on CDS genes incl. pseudo; **rejected on RNA genes** (PGAP oracle: EC_number only on CDS, 319 normal + 2 pseudo, never on RNAs).
- Errors: unknown locus tag; regex mismatch; unknown category (message tells the user to add it to the JSON).
- **GC and GP sidecars are forbidden** (Arx derives them from `/gene` and `/product`; prefix handling is inconsistent).
- Example `{id}.KG`:
  ```
  EXAMPLE1234-2-1.1_000001	K00001
  EXAMPLE1234-2-1.1_000008	K00001, K00002
  ```

## Taxonomy helpers (NCBI E-utilities, cached on disk in `~/.cache/gbk-creator`)
Credentials: customer passes `ncbi_email` / `ncbi_api_key`, or sets env vars `NCBI_EMAIL` / `NCBI_API_KEY`.
The library never reads `.env` files itself. Locally: gitignored `.env` (template: `.env.example`),
used via `uv run --env-file .env pytest`. Never commit keys.
1. latin name → taxid (esearch; must be exactly one hit)
2. taxid → scientific name + rank + lineage (efetch)

## Test frontend + tests
Parser script: `.arx.gff` + `.fna` (+ `.faa`/`.ffn` for sequences) + `*.emapper.annotations` → builder → output.
- GFF: GFF3 can carry EC (`EC_number=`) and GO (`Ontology_term=GO:...`). Our gff has `EC_number` (bare, needs `EC:` prefix), no GO.
  These go to `.EC` / `.GO` sidecars.
- Known GFF loss: `join()` locations (ribosomal slippage) collapse to one span; partial `<`/`>` are in the coordinate columns.
- eggnog (emapper 1.0.3, 22 cols; mapping from Arx `load_eggnog_file`):
  GOs(6)→GO, EC(7)→`EC:`x, KEGG_ko(8)→KG (strip `ko:`), KEGG_Reaction(11)→KR,
  Preferred_name(5)→`EP:`x, eggNOG OGs(18)→`EO:`+part before first `@`, description(21)→`ED:`x.
  Query names have a `gnl|extdb|` prefix → strip.
- Tests (pytest):
  - **Semantic diff** vs oracle: parse both with Biopython; compare features/locations/qualifiers with an explicit
    allowlist of known differences (lineage, dates, COMMENT, joins, EC_number).
  - Unit tests for every rule violation above.
  - Acceptance: run arx-tools `check_files` / `import_genome` on output; final manual check by importing into a real Arx instance.

## Resolved decisions (formerly open)
- Contig names: user-supplied, validated (see Naming).
- `gene` parent: auto-created by `add_gene`.
- Annotation JSON: both formats.
- Duplicate annotations: ignored.
- Annotations on pseudo/RNA: copy PGAP (CDS incl. pseudo OK, RNA rejected).
- NCBI credentials: args or env vars; see Taxonomy helpers.
- Global ID clashes in Arx: fix in Arx (`website/models/GenomeContent.py:238,363,388` use `bulk_create(ignore_conflicts=True)`), not here.
