# AGENTS.md — genetics.el

Emacs package: read, browse, annotate, compare and report on consumer
genetics raw-data exports (23andMe, AncestryDNA, MyHeritage/FTDNA CSV,
VCF / VCF.gz). Everything is local; README.md is the full reference.

## Layout

- `genetics.el` — entry point (`genetics-open`, summary buffer).
- `genetics-base.el` — customization, errors, JSON and output-file helpers.
- `genetics-core.el` — structs, chromosome/genotype helpers, kit access, registry.
- `genetics-detect.el` — format, build and reference-call detection, strand notes.
- `genetics-gzip.el` — gzip/BGZF decompression (gzip executable or zlib).
- `genetics-parse.el` — parsers, VCF offset index, cache, `genetics-parse-file`.
- `genetics-stats.el` — stats, sex inference, caveats text.
- `genetics-catalog.el` — annotation files (JSON/Org): loading, cache, lookup.
- `genetics-annotate.el` — curated-site resolution, risk-allele logic, APOE.
- `genetics-browse.el`, `-lookup.el`, `-report.el`, `-compare.el`, `-export.el`.
- `genetics-org.el` — Org dynamic blocks (genetics-summary, -hits, -apoe).
- `examples/` — sample Org report, report and screenshot regeneration scripts.
- `genetics-snpedia.el` — the ONLY file allowed to touch the network.

## Changing it

- `make test`, `make compile` (warnings are errors) and `make checkdoc`
  must pass (`make lint` runs the last two). Without make (Windows):
  `emacs -Q --batch -l test/run-tests.el [test|compile|checkdoc|all]`.
- Cross-platform (Linux, macOS, native Windows; CI matrix in
  `.github/workflows/test.yml`): find programs with `executable-find`, run
  them with argument lists (no shell), build paths with `expand-file-name`,
  write data files with `genetics--with-output-file` (UTF-8, LF).
- Fixtures in `test/fixtures/` are synthetic. Never read the real data
  directory (`genetics-data-directory`) from tests or code paths run in tests.
- Errors: `define-error` under `genetics-error`; never message-and-return-nil.
- Privacy: no network code outside genetics-snpedia.el (a test greps for it).
  SNPedia is opt-in and sends only the rsid.
- Genotype strings: uppercase letters, `--` for no-call; strand notes matter.
