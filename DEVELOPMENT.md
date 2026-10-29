# Development notes

User documentation lives in [README.md](README.md); the design spec and the rules the library enforces are in
[Plan.md](Plan.md).

## Setup and tests

```sh
uv sync                                  # runtime deps + dev group (pytest, arx-tools)
cp .env.example .env                     # gitignored; set NCBI_EMAIL / NCBI_API_KEY
uv run --env-file .env pytest            # everything
uv run pytest -m "not network"           # offline only (unit tests never touch the network anyway)
python examples/from_database.py         # the README example (offline)
```

* Unit tests (rules, format, taxonomy with a faked HTTP layer) run offline. Offline tests run with
  `NCBI_EMAIL`/`NCBI_API_KEY` removed from the environment, so real credentials can never show up in assertion output.
* Real-data tests (semantic diff against the PGAP oracle gbk, arx-tools acceptance incl. `import_genome`) read
  `data/{sample}/...` (gitignored) and skip when it is absent.
* `tests/test_taxonomy_transport.py` mocks `urlopen` and checks that credentials never leak (messages, tracebacks,
  logs, cache files) plus retry/backoff and rate limiting.
* `tests/test_examples.py` runs `examples/from_database.py`, and checks that the README's coordinate table, its embedded
  example excerpt and its feature/qualifier tables match the code (`ALLOWED_*` constants in `builder.py`).
* Tests marked `network` talk to NCBI E-utilities and skip without `NCBI_EMAIL`.
* The library never reads `.env` files and never stores the API key; never commit one.

## Test frontend (`gbk-creator` CLI)

`gbk-creator` converts an Arx-style sample (`.arx.gff` + `.fna`, optionally `.faa`, `.ffn`, `*.emapper.annotations`)
by driving the builder; it exists to exercise the library with real data (`src/gbk_creator/frontends/arx_gff.py`):

```sh
uv run --env-file .env gbk-creator data/FAM13496-i1-1.1 --out out/FAM13496-i1-1.1
uv run --env-file .env gbk-creator --gff X.arx.gff --fna X.fna --faa X.faa --ffn X.ffn --eggnog X.emapper.annotations --out out/X
```

Other flags: `--identifier`, `--strain`, `--taxid`, `--annotation-regex`, `--plasmid`, `--restricted`, `--tag` (repeatable),
`--old-identifier`, `--lenient`.

* EC numbers in the GFF (`EC_number`) become `EC:`-prefixed `.EC` entries, `Ontology_term` becomes `.GO`; eggnog columns are
  mapped like Arx's `load_eggnog_file`.
* With a `.ffn`, joins that the GFF collapsed (ribosomal slippage, tRNA introns, features across the origin of circular
  contigs) are reconstructed.
* It fails loudly by default: leftover `.faa`/`.ffn` records, translations on non-CDS features, annotation ids that do
  not match their category regex (PGAP's preliminary `EC:3.6.5.n1`), `gene` rows without CDS/RNA child and unknown
  qualifiers all raise. `--lenient` (Python: `lenient=True`, optionally `dropped=[]` to collect what was lost) is the
  single opt-in lossy mode: it drops leftover faa/ffn records, regex-invalid annotation ids and orphan `gene` rows with a
  logged warning. Translations on non-CDS features always raise. The oracle tests list the exact expected drops per
  sample (`EXPECTED_LENIENT_DROPS` in `tests/conftest.py`).

## What arx-tools does with our files (merge details)

* `organism.json` is merged over the organism's existing `organism.json`: every non-null, non-empty value wins
  (including `restricted: false` and `taxid`); omitted `tags`, `tags=[]` and null values leave existing ones alone.
* `genome.json` (`cds_tool`, optional `old_identifier`) is merged after arx-tools' own gbk-based detection, so ours wins;
  arx-tools does not recognise our COMMENT block (`Failed to discover annotation information` warning). We never write
  `contig_format` (arx-tools reads it from the import dir's `genome.json`).
* Both files are also copied into the genome folder; arx-tools then overwrites `genome.json` with the merged result.
* arx-tools only globs a fixed list of annotation types (`GC GP EP ED EO EC KG KR GO SL OL`): see the README section on
  custom categories; `test_new_two_letter_category_needs_custom_import_settings` covers it.
* `/old_locus_tag` is written to the `gene` and the child feature as the last qualifier (like arx-tools' rename mode);
  Arx reads the first one per locus (`Gene.old_locus_tag`, max 100 characters).
