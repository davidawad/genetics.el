# AGENTS.md — genetics.el

Emacs package: read, browse, annotate, compare and report on consumer
genetics raw-data exports (23andMe, AncestryDNA, MyHeritage/FTDNA CSV,
VCF / VCF.gz). Everything is local; README.md is the full reference.

## Layout

- `genetics.el` — entry point (`genetics-open`, summary buffer).
- `genetics-core.el` — customization, errors, structs, shared helpers.
- `genetics-parse.el` — format detection, parsers, VCF offset index, cache.
- `genetics-stats.el` — stats, sex inference, caveats text.
- `genetics-annotate.el` — annotation files (JSON/Org), risk-allele logic, APOE.
- `genetics-browse.el`, `-lookup.el`, `-report.el`, `-compare.el`, `-export.el`.
- `genetics-snpedia.el` — the ONLY file allowed to touch the network.

## Changing it

- `make test`, `make compile` (warnings are errors) and `make checkdoc`
  must pass (`make lint` runs the last two).
- Fixtures in `test/fixtures/` are synthetic. Never read the real data
  directory (`genetics-data-directory`) from tests or code paths run in tests.
- Errors: `define-error` under `genetics-error`; never message-and-return-nil.
- Privacy: no network code outside genetics-snpedia.el (a test greps for it).
  SNPedia is opt-in and sends only the rsid.
- Genotype strings: uppercase letters, `--` for no-call; strand notes matter.
