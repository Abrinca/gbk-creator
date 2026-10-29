# gbk-creator

Turn genome data (e.g. from your own database) into a **PGAP-canonical GenBank file plus custom-annotation
sidecars** that import into [Arx](https://github.com/Abrinca/arx-tools) with `import_genome` — no renaming needed.
See [Plan.md](Plan.md) for the full spec and rules.

```
out/FAM13496-i1-1.1/FAM13496-i1-1.1.gbk
out/FAM13496-i1-1.1/FAM13496-i1-1.1.KG     # locus_tag<TAB>anno1, anno2   (one file per annotation type)
out/FAM13496-i1-1.1/organism.json          # name (= strain), taxid, restricted, [tags]
out/FAM13496-i1-1.1/genome.json            # cds_tool, [old_identifier]
```

Nothing else is written; arx-tools derives fna/faa/ffn/gff/genome.json itself.

## Usage

```python
from gbk_creator import GenBankBuilder

b = GenBankBuilder(
    identifier="FAM13496-i1-1.1",   # max 50 chars, no space : / \ |
    taxid=1308,                     # name, rank and lineage come from NCBI (cached on disk)
    strain="FAM13496",              # Arx "organism"; identifier must start with it; default = identifier
    old_identifier="FAM13496-old",  # optional, for genome.json: links the new Arx genome to its old name
    restricted=False,               # optional, a real bool
    tags=["probiotic"],             # optional, unique non-empty strings; None = no "tags" key
)
b.add_contig("FAM13496-i1-1.1_scf1", sequence, topology="linear", kind="scaffold")  # names: {id}_scf1, _scf2, ... in order
b.add_gene("FAM13496-i1-1.1_000001", contig="FAM13496-i1-1.1_scf1",
           location=SimpleLocation(145, 586, strand=-1), feature_type="CDS",
           product="SprT family protein", translation="MNLTD...", protein_id="extdb:FAM13496-i1-1.1_000001")
b.add_feature("FAM13496-i1-1.1_scf1", "repeat_region", SimpleLocation(4827, 6116, 1), rpt_family="CRISPR", ...)
b.add_annotation("FAM13496-i1-1.1_000001", "KG", "K03095")
b.write("out/FAM13496-i1-1.1")      # validates everything first, then writes gbk + sidecars
```

Anything that cannot be written faithfully or violates the rules raises a `GbkCreatorError` (a `ValueError`) with
a message saying what to fix — e.g. an unknown annotation category tells you to add it to your annotation JSON.
`b.validate()` runs the global checks (e.g. locus tags increasing along contig order and position) without writing.

Annotation categories (`KG`, `EC`, `GO`, ...) and their regexes come from a JSON file, either `{"KG": "^K[0-9]{5}$"}`
or Arx's `annotations.json` format; pass it as `annotation_regex="my.json"` (default: the packaged
`annotation_regex.json`). Category names must be **exactly two upper-case letters** (Arx stores them in a
2-character field). `GC`/`GP` sidecars are forbidden (Arx derives them from `/gene` and `/product`);
EC numbers and GO terms go only to sidecars, never into the gbk.

#### Adding your own category (e.g. `ZZ`)

1. Add it to your annotation JSON (`{"ZZ": "^ZZ[0-9]{3}$"}`) and pass it as `annotation_regex=`; `add_annotation(lt, "ZZ", "ZZ001")` then writes `{identifier}.ZZ`.
2. The category must also exist in **Arx's `annotations.json`** (folder structure), or Arx rejects the file on import.
3. arx-tools' default `file_finder` only globs a fixed list of types (`GC GP EP ED EO EC KG KR GO SL OL`), so
   the new file would be copied but not registered. Point `ARX_IMPORT_SETTINGS` (or `import_genome --import_settings`)
   to a JSON file that repeats the default `custom_annotations` entries and adds yours:
   ```json
   {"file_finder": {"custom_annotations": [
       {"glob": "*.GC", "anno_type": "GC", "expected": false},
       "... one entry per default type (GP EP ED EO EC KG KR GO SL OL) ...",
       {"glob": "*.ZZ", "anno_type": "ZZ", "expected": false}
   ]}}
   ```
   (`tests/test_arx_acceptance.py::test_new_two_letter_category_needs_custom_import_settings` proves both the
   problem and the fix.)

Taxonomy helpers: `gbk_creator.name_to_taxid("Streptococcus thermophilus")` (exactly one NCBI hit required) and
`gbk_creator.fetch_taxon(1308)` (scientific name, rank, lineage).

### organism.json and genome.json

`write()` also writes two small metadata files that arx-tools merges in on `import_genome`:

* `organism.json`: `name` (= strain), `taxid`, `restricted` and `tags` (only if given); no `representative`, arx-tools
  sets `name` and `representative` itself. It is **merged over the existing organism's `organism.json`**: every
  non-null, non-empty value here wins (including `restricted: false` and `taxid`), whereas omitted `tags` and
  `tags=[]` leave existing values alone.
* `genome.json`: `cds_tool` (`"Abrinca GenBankBuilder v<version>"`) and, if `old_identifier=` is given (non-empty
  single line, max 50 characters, different from the identifier), `old_identifier` to link the new genome to its old
  name. arx-tools does not recognise our COMMENT block for its own `cds_tool` detection; should it ever, its
  `genome.json` merge runs afterwards, so our values win. No `contig_format` or file names are written.

Importing also copies both files into the genome folder (arx-tools then overwrites `genome.json` with the merged
result). The CLI exposes `--restricted`, `--tag` (repeatable) and `--old-identifier`.

### NCBI credentials

Pass `ncbi_email=` / `ncbi_api_key=` to `GenBankBuilder`, or set the environment variables `NCBI_EMAIL` (required) and
`NCBI_API_KEY` (optional, raises the rate limit). The library never reads `.env` files and never stores the key.
Results are cached in `$GBK_CREATOR_CACHE` (default `~/.cache/gbk-creator`); delete it to refresh.

## Test frontend (CLI)

`gbk-creator` converts an Arx-style sample (`.arx.gff` + `.fna`, optionally `.faa`, `.ffn`, `*.emapper.annotations`)
by driving the builder; it exists to exercise the library with real data:

```sh
uv run --env-file .env gbk-creator data/FAM13496-i1-1.1 --out out/FAM13496-i1-1.1
uv run --env-file .env gbk-creator --gff X.arx.gff --fna X.fna --faa X.faa --ffn X.ffn --eggnog X.emapper.annotations --out out/X
```

EC numbers in the GFF (`EC_number`) become `EC:`-prefixed `.EC` entries, `Ontology_term` becomes `.GO`; eggnog columns are
mapped like Arx's `load_eggnog_file`. With a `.ffn`, joins that the GFF collapsed (ribosomal slippage, tRNA
introns, features across the origin of circular contigs) are reconstructed.

The frontend fails loudly by default: leftover `.faa`/`.ffn` records, translations on non-CDS features, annotation
ids that do not match their category regex (PGAP's preliminary `EC:3.6.5.n1`), `gene` rows without CDS/RNA child
and unknown qualifiers all raise. `--lenient` (Python: `lenient=True`, optionally `dropped=[]` to collect what was
lost) is the single opt-in lossy mode: it drops leftover faa/ffn records, regex-invalid annotation ids and orphan
`gene` rows with a logged warning instead. Translations on non-CDS features always raise.

## Development

```sh
uv sync                                  # runtime deps + dev group (pytest, arx-tools)
cp .env.example .env                     # gitignored; set NCBI_EMAIL / NCBI_API_KEY
uv run --env-file .env pytest            # everything
uv run pytest -m "not network"           # offline only (unit tests never touch the network anyway)
```

* Unit tests (rules, format, taxonomy with a faked HTTP layer) run offline.
* Real-data tests (semantic diff against the PGAP oracle gbk, arx-tools acceptance incl. `import_genome`) read
  `data/{sample}/...` (gitignored) and skip when it is absent.
* `tests/test_taxonomy_transport.py` mocks `urlopen` and checks that credentials never leak (messages, tracebacks,
  logs, cache files) plus retry/backoff and rate limiting. Offline tests run with `NCBI_EMAIL`/`NCBI_API_KEY` removed.
* Tests marked `network` talk to NCBI E-utilities and skip without `NCBI_EMAIL`.
