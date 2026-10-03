# genetics.el

Read, browse, annotate, compare and report on consumer genetics raw-data
exports inside Emacs. Everything runs locally; there is no network code
except an opt-in, rsid-only SNPedia lookup (see [Privacy](#privacy)).

> Informational only, not medical advice. Consumer genotyping arrays are not
> diagnostic. Confirm any finding with a clinical-grade test and a clinician.

## Overview

- Auto-detects and parses 23andMe, AncestryDNA, MyHeritage / FamilyTreeDNA
  CSV and VCF (plain and `.vcf.gz`).
- Multi-GB VCFs are never loaded: above a size limit the file is
  offset-indexed (rsid -> byte offset, per-chromosome ranges) and records are
  read lazily by seeking.
- Summary buffer with SNP count, no-call rate, heterozygous/homozygous counts,
  per-chromosome counts, inferred sex and strand/build caveats.
- Filterable tabulated browser, single-rsid lookup, annotation layer (JSON or
  Org), Org report with APOE haplotype, kit-to-kit comparison, CSV/JSON export.

## Install

Requires Emacs 29.1 or newer. For gzipped VCFs the `gzip` executable must be
on `exec-path`.

```elisp
(use-package genetics
  :load-path "~/src/genetics-el"
  :commands (genetics-open genetics-browse genetics-lookup
             genetics-report genetics-compare))
```

or, with `package-vc` (Emacs 29+):

```elisp
(package-vc-install "https://gitlab.com/davidawad/genetics-el")
```

## Quick start

```elisp
M-x genetics-open RET ~/dna/genome.txt RET   ; parse + summary buffer
M-x genetics-browse                          ; table of records
M-x genetics-lookup RET rs1801133 RET        ; one SNP in every loaded kit
M-x genetics-report                          ; Org report
M-x genetics-compare                         ; two loaded kits
```

`genetics-open` prompts starting in `genetics-data-directory`.

## Supported formats

Detection looks at header comments and columns, not the file extension.

| Format | Example | Notes |
|--------|---------|-------|
| 23andMe `.txt` | `rs4477212<TAB>1<TAB>82154<TAB>AA` | `#` comments; build read from "build 37"; chip v3/v4/v5 from a header hint or SNP count (about 960k / 570k / 630k), otherwise `unknown`; `--` is a no-call; single letters on X/Y/MT are hemizygous |
| AncestryDNA `.txt` | `rs4477212<TAB>1<TAB>82154<TAB>A<TAB>A` | chromosomes 23/24/25/26 become X/Y/XY/MT; `0` alleles are no-calls |
| MyHeritage CSV | `"rs4477212","1","82154","AA"` | `RSID,CHROMOSOME,POSITION,RESULT`, `#` comments mentioning MyHeritage |
| FamilyTreeDNA CSV | `"rs4477212","1","82154","AA"` | same columns, no comments (build 37 assumed, and flagged as assumed) |
| VCF / VCF.gz | `chr1 82154 rs4477212 A G . . . GT 0/1` | `GT` converted to letters (`0/1` -> `AG`, phased, multi-allelic, `./.` -> no-call, haploid); `chr` prefix removed; `.` ids stored as `chrom:pos`; build from `##reference` / `##contig` |

Genotypes are stored as upper-case letters (`AG`), `--` for a no-call. VCF
indels use `/` between alleles (`AT/A`). Comparison and filters are
order-insensitive (`AG` equals `GA`).

### Large VCFs

When a VCF is larger than `genetics-vcf-eager-limit` (200 MB by default) it is
indexed in chunks of `genetics-chunk-size` bytes without keeping records.
`.vcf.gz` above the limit is decompressed once into `genetics-cache-directory`
and indexed there. Lookups, browsing a chromosome and comparison work on such
kits; the no-call and heterozygosity statistics and sex inference are not
computed for them.

## Commands and keys

| Command | Purpose |
|---------|---------|
| `genetics-open` | Parse a file, register the kit, show the summary |
| `genetics-summary` | Show the summary buffer of a kit |
| `genetics-close` | Unload a kit |
| `genetics-browse` | Filterable table of records |
| `genetics-lookup` | Detail for one rsid across loaded kits |
| `genetics-report` | Org report (prefix argument: also write a file) |
| `genetics-compare` | Concordance between two kits |
| `genetics-export-csv`, `genetics-export-json` | Export all records, or the browser's filtered view |
| `genetics-reload-annotations` | Re-read annotation files |

Summary buffer (`genetics-summary-mode`): `b` browse, `r` report, `l` lookup,
`c` compare, `g` refresh; buttons for the same.

Browser (`genetics-browse-mode`, derived from `tabulated-list-mode`):

| Key | Command |
|-----|---------|
| `RET` | lookup the rsid at point |
| `c` | filter chromosome |
| `p` | filter position range |
| `s` | filter rsid by regexp |
| `g` | filter genotype (order-insensitive) |
| `n` / `h` / `o` | no-calls only / heterozygous only / homozygous only |
| `a` | toggle annotated only |
| `x` | clear all filters |
| `E` / `J` | export current view to CSV / JSON |
| `R` / `i` | report / summary |

Active filters are shown in the header line. At most `genetics-browse-limit`
rows are displayed; the header line says when the list is truncated. Exports
are not limited.

## Customization

| Variable | Default | Meaning |
|----------|---------|---------|
| `genetics-data-directory` | `~/Documents/Genetics/` | start directory when prompting for a file |
| `genetics-use-cache` | `t` | cache parsed kits as `.eld` files |
| `genetics-cache-directory` | `(locate-user-emacs-file "genetics-cache/")` | cache location (also holds SNPedia answers, decompressed VCFs) |
| `genetics-vcf-eager-limit` | 200 MB | VCFs above this are offset-indexed |
| `genetics-chunk-size` | 4 MB | read size when streaming |
| `genetics-browse-limit` | 5000 | maximum rows in the browser |
| `genetics-annotation-files` | the shipped curated JSON | annotation files, later ones override |
| `genetics-snpedia-enabled` | `nil` | allow SNPedia lookups |

The parse cache is an `.eld` file (printed Lisp, no sqlite) named from a hash
of the file's truename, size and modification time, so editing or replacing the
file invalidates it automatically.

## Annotation files

`genetics-annotation-files` lists JSON and/or Org files. The shipped
`annotations/genetics-curated.json` covers APOE (rs429358, rs7412), MTHFR
(rs1801133, rs1801131), Factor V Leiden (rs6025), HFE C282Y (rs1800562), LCT
(rs4988235) and CYP2C19*2 (rs4244285), each with its strand spelled out and
SNPedia / dbSNP links.

JSON: an array of objects.

```json
[{"rsid": "rs1801133", "gene": "MTHFR", "risk_allele": "A",
  "other_allele": "G", "effect": "...", "magnitude": 2,
  "strand": "+ strand G/A; gene strand C/T", "notes": "...",
  "url": "https://www.snpedia.com/index.php/Rs1801133",
  "genotypes": {"GG": "...", "AG": "...", "AA": "..."}}]
```

Only `rsid` is required. `other_allele` lets the package recognize palindromic
and strand-flipped genotypes. Genotype keys may be written in any allele order.

Org: one heading per SNP with a property drawer (see
`annotations/genetics-example.org`).

```org
* rs1805007 MC1R R151C
:PROPERTIES:
:GENE: MC1R
:RISK_ALLELE: T
:OTHER_ALLELE: C
:EFFECT: Associated with red hair and fair skin.
:MAGNITUDE: 2
:STRAND: + strand, GRCh37: C/T.
:URL: https://www.snpedia.com/index.php/Rs1805007
:GT_CT: One R151C allele.
:END:
Free text becomes the notes.
```

`RSID` may be omitted when the heading contains an `rsNNN` id. Files are
re-read automatically when they change.

### Risk alleles and strand

Risk-allele copies are counted on the strand the data uses. If neither allele
matches the risk allele but its complement does, the result is flagged
"possible strand flip" (with the copy count that would apply) and is never
flipped silently. Palindromic SNPs (A/T, C/G) are flagged as ambiguous.

APOE: rs429358 (T/C) and rs7412 (C/T) on the + strand give e2 = T+T, e3 = T+C,
e4 = C+C (e1 = C+T, very rare). A double heterozygote (CT + CT) is reported as
e2/e4 (most likely) or e1/e3 (very rare); phase cannot be resolved from
unphased genotypes.

## Caveats

- Strand: 23andMe reports the + (forward) strand of GRCh37; AncestryDNA,
  MyHeritage and FTDNA are assumed to do the same; VCF genotypes are relative
  to the file's REF on its own build. Many clinical sources quote gene-strand
  alleles (for example MTHFR C677T is `A` on the + strand and `T` on the gene
  strand). The curated annotations spell out both.
- Build: no liftover is performed. Comparing kits on different builds warns.
- Arrays measure a small fraction of the genome; no-calls and errors occur.

## Privacy

- All parsing, annotation, reporting and caching happens locally.
- Compressed VCFs are decompressed by the local `gzip` into a temp buffer or a
  file in your cache directory; nothing is uploaded.
- The only network code is `genetics-snpedia.el`. It is disabled unless
  `genetics-snpedia-enabled` is non-nil, asks for confirmation on first use per
  session, sends only the rsid (for example `.../api.php?action=parse&page=Rs429358&...`),
  never genotypes, positions or file names, and caches answers on disk.
- A test greps the package sources and fails if `url-retrieve`,
  `make-network-process` or `open-network-stream` appear anywhere else.
- Annotation links open in your browser only when you click them.
- The cache directory contains your parsed genotypes in plain text; protect it
  like the original file.

## Development

```sh
make test       # ERT suite (offline, synthetic fixtures only)
make compile    # byte-compile with byte-compile-error-on-warn
make checkdoc   # checkdoc on all non-test sources; fails on any warning
make lint       # compile + checkdoc
make clean
```

Files: `genetics.el` (entry point), `genetics-core.el`, `genetics-parse.el`,
`genetics-stats.el`, `genetics-annotate.el`, `genetics-browse.el`,
`genetics-lookup.el`, `genetics-report.el`, `genetics-compare.el`,
`genetics-export.el`, `genetics-snpedia.el`. Tests live in `test/*-test.el`,
fixtures (synthetic, fake genotypes) in `test/fixtures/`.

## License

MIT, Copyright (C) 2026 David Awad.
