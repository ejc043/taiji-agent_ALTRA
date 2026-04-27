# Bulk Taiji input sheet schema

Read this file when the user asks about the xlsx format in detail, or when you
need to validate/fix a sheet by hand rather than via the build script.

## Sheet 1: `Active`

One row per *assay track*, not per sample. A biological sample (`group`)
produces two or three rows: RNA-seq, ATAC-seq, and optionally HiC.

### Column order and semantics

The column order is fixed. Taiji parses by column name, but keeping the order
stable is helpful for humans reviewing the sheet.

```
type | id | group | rep | path | tags | format | cohort
```

**`type`** — one of `RNA-seq`, `ATAC-seq`, `HiC`. Anything else is an error.

**`id`** — unique across the whole sheet. The convention (which the build
script follows) is:
- RNA-seq → `RNA-{group}` (hyphen between prefix and group)
- ATAC-seq → `ATAC-{group}` (hyphen)
- HiC → `HiC_{group}` (underscore — yes, it's inconsistent; preserved to
  match the reference sheet format Taiji ships with)

**`group`** — the biological-sample key. This is the single most important
column because it's what ties RNA+ATAC(+HiC) rows together into one sample.
All rows that share a `group` are treated by Taiji as measurements of the
same sample.

**`rep`** — integer replicate number, typically `1`. If you have technical
replicates, give each its own row with the same `group` and incrementing
`rep`.

**`path`** — absolute path to the processed file on the filesystem that will
actually run Taiji. Relative paths are a bad idea: Taiji's working directory
isn't guaranteed to be where you think it is.

**`tags`** — semantic hint for Taiji about what's in the file:
- `GeneQuant` for RNA-seq (gene-level quantification, two-or-more column TSV
  with gene IDs in column 1)
- `ChromosomeLoop` for HiC
- *blank* for ATAC-seq (ATAC uses `format` instead of `tags`)

**`format`** — file format hint, currently only used for ATAC:
- `NarrowPeak` for standard ENCODE narrowPeak (BED6+4) ATAC peak files
- *blank* for RNA-seq and HiC

**`cohort`** — condition label used to construct Taiji contrasts
(e.g. `KO` vs `WT`, `naive` vs `effector`, `tumor` vs `normal`). All rows for
the same `group` must have the same `cohort`.

### Per-type row template

```
RNA-seq  | RNA-{group}  | {group} | 1 | /abs/path/{group}.rna.txt     |            |            | {cohort}
ATAC-seq | ATAC-{group} | {group} | 1 | /abs/path/{group}.narrowPeak  |            | NarrowPeak | {cohort}
HiC      | HiC_{group}  | {group} | 1 | /abs/path/loops.txt           | ChromosomeLoop |        | {cohort}
```

## Sheet 2: `active_metadata`

One row per `group`. All four columns are required.

```
Submitter_ID | Case_ID | vcf_Location | gtf_Location
```

- **`Submitter_ID`** and **`Case_ID`** — both set to the `group` value. Taiji
  treats them as separate fields but in practice they are redundant for most
  single-site analyses.
- **`vcf_Location`** — absolute path to the genome FASTA (despite the name).
  The column is historically named for variant annotation, but Taiji uses it
  as the reference sequence.
- **`gtf_Location`** — absolute path to the gene annotation GTF for the same
  build.

All rows in `active_metadata` should reference the same genome build within a
single run; mixing builds is a run-breaking error.

## Validation rules the build script enforces

The build script rejects a sheet that would violate these, unless
`--no-strict` is passed:

1. Every `group` in `Active` appears in `active_metadata` and vice versa.
2. Every `group` has at least one `RNA-seq` and one `ATAC-seq` row.
3. `id` values are unique.
4. `cohort` is consistent within `group`.
5. `tags` and `format` match the expected values for each `type`.
6. All `path`, `vcf_Location`, and `gtf_Location` files exist and are
   readable at the time the sheet is written.
7. All metadata rows use the same genome build (resolved from `--genome`).

Non-fatal warnings (printed but not errors even in strict mode):

- HiC row is missing for a group (HiC is optional).
- Genome config has a TODO path — the script will write the sheet but prints
  a prominent warning.
- A group has `rep > 1` for some types but not others — usually a sign of a
  dropped replicate file.

## Why the odd naming?

`vcf_Location` for a FASTA, `HiC_` with underscore while RNA/ATAC use hyphens,
the two redundant ID columns in `active_metadata` — these are inherited from
the reference sheet format used by the Wang lab's in-house Taiji runs. Keeping
them stable means the xlsx is drop-in compatible with existing configs. If
a future Taiji release drops the quirks, the build script can be updated in
one place.
