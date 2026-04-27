---
name: build-taiji-input
description: Construct the two-sheet xlsx input file that bulk Taiji consumes (`Active` + `active_metadata`) from a directory of processed RNA-seq / ATAC-seq / HiC files plus a small CSV of sample-to-cohort mappings. Use this skill whenever the user wants to prepare, build, generate, assemble, or create a Taiji input sheet, input xlsx, samplesheet, or config for a bulk Taiji run — including phrases like "make the Taiji input", "set up Taiji", "prep samples for Taiji", "convert my files for Taiji", or anything that implies turning a pile of per-sample files into the exact xlsx format bulk Taiji requires. Also trigger when the user mentions sample groups, GeneQuant, NarrowPeak, or ChromosomeLoop in the same breath as a run setup.
---

# Build Taiji input xlsx

## What this skill produces

A single `.xlsx` file with two sheets, matching exactly what bulk Taiji reads:

1. **`Active`** — one row per *assay track*. Each biological sample (a "group") contributes up to three rows: RNA-seq, ATAC-seq, and optionally HiC. Columns (in order):

   | Column   | Required | Content                                                                 |
   |----------|----------|-------------------------------------------------------------------------|
   | `type`   | yes      | `RNA-seq`, `ATAC-seq`, or `HiC`                                         |
   | `id`     | yes      | Unique across the sheet. Convention: `RNA-{group}`, `ATAC-{group}`, `HiC_{group}` |
   | `group`  | yes      | Biological sample key that pairs tracks together                        |
   | `rep`    | yes      | Replicate number (default `1`)                                          |
   | `path`   | yes      | Absolute path to the processed file                                     |
   | `tags`   | cond.    | `GeneQuant` for RNA, `ChromosomeLoop` for HiC, blank for ATAC           |
   | `format` | cond.    | `NarrowPeak` for ATAC, blank otherwise                                  |
   | `cohort` | yes      | Condition label (e.g. `KO`, `WT`)                                       |

2. **`active_metadata`** — one row per `group`, columns: `Submitter_ID`, `Case_ID`, `vcf_Location`, `gtf_Location`. The genome FASTA goes in `vcf_Location` (this is a Taiji naming quirk — it is a FASTA, not a VCF) and the gene annotation GTF in `gtf_Location`. Both are resolved from the chosen genome build.

## Inputs this skill accepts

Two things, plus a genome flag:

1. `--data-dir PATH` — a directory holding per-sample processed files. The skill globs for:
   - RNA: `{group}.rna.txt` (GeneQuant; two+ column TSV, gene IDs in column 1)
   - ATAC: `{group}.narrowPeak` (standard ENCODE narrowPeak / BED6+4)
   - HiC (optional): `{group}.hic.loops` or whatever the `--hic-pattern` says
2. `--samples PATH` — a CSV with at least these columns:
   - `group` (required): must match the filename stem in the data dir
   - `cohort` (required): condition label, e.g. `KO`, `WT`, `naive`, `effector`
   - `hic_path` (optional): absolute path to a per-sample HiC loops file, overrides the genome default
   - `rep` (optional): replicate number, default `1`
3. `--genome NAME` — one of the builds in `assets/genomes.yml` (e.g. `mm10`, `hg38`). Resolves FASTA, GTF, and the default HiC loops path for that species/build.

A template CSV lives at `assets/samples_template.csv`. The genome registry lives at `assets/genomes.yml` — **users of this skill at a new site must edit this file once** to point at their local reference paths. The example values are UCSD/TCF1-project paths and will not work on other clusters.

## How to run it

The skill bundles a single script. Invoke it with:

```bash
python scripts/build_taiji_input.py \
  --data-dir /stg3/data1/eunice/Projects/JohnChang/TCF1/data \
  --samples  ./samples.csv \
  --genome   mm10 \
  --out      ./taiji.input.xlsx
```

Optional flags worth knowing:

- `--rna-pattern "{group}.rna.txt"` — override the RNA filename pattern
- `--atac-pattern "{group}.narrowPeak"` — override the ATAC pattern
- `--hic-pattern "{group}.hic.loops"` — if set, the skill looks for per-sample HiC files; if absent, the per-genome default HiC path is used
- `--no-hic` — skip HiC rows entirely even if files exist
- `--genome-config PATH` — use a different genomes.yml (e.g. for a second site)
- `--strict` — fail on any missing RNA or ATAC file; otherwise warn and skip that group

The script:

1. Reads the samples CSV.
2. For each group, resolves RNA and ATAC file paths under `--data-dir` using the name patterns. Missing files become structured errors in `--strict` mode, warnings otherwise.
3. Resolves HiC: per-sample override > discovered file > genome default > omitted (if `--no-hic`).
4. Emits the `Active` sheet with the correct `tags`/`format` for each `type`, and IDs using the `RNA-{group}` / `ATAC-{group}` / `HiC_{group}` convention.
5. Emits the `active_metadata` sheet with FASTA and GTF from the chosen genome.
6. Prints a summary (N groups, N tracks, missing files, cohort breakdown).
7. Exits non-zero on validation errors (unless `--no-strict` was passed).

## When to reach for `references/sheet_schema.md`

If the user asks about the exact schema, column semantics, or wants to hand-write rows, open `references/sheet_schema.md` — it has the full annotated example and the rules about what is valid for each `type`. The SKILL.md body is kept short on purpose; detailed format rules live in the reference file.

## Interaction pattern

When a user brings messy inputs, work through this order:

1. **Confirm the genome build** before anything else. Mixed builds silently break Taiji's downstream network inference and the error messages are confusing, so make the build explicit up front and echo it back.
2. **Inspect what's in the data dir.** If the filename convention differs from the defaults (e.g. `sample.rna.bed` instead of `{group}.rna.txt`), propose a `--rna-pattern` override rather than renaming files.
3. **Ask for the samples CSV if none was provided.** The cohort column can't be inferred reliably — mis-labeling KO/WT ruins the analysis, so this is worth one round trip to confirm.
4. **Run the script with `--strict`** on the first pass so the user sees every missing file explicitly. Only drop to non-strict once the user has decided which samples to drop.
5. **Show the summary table** the script prints (groups × cohorts × tracks) so the user can eyeball that the counts match their expectations before handing the xlsx to Taiji.

## Why the defaults are what they are

The column set, tag conventions, and id format mirror an in-use Taiji bulk input sheet (mm10, TCF1 project) — the goal is that a sheet produced by this skill is drop-in compatible with existing Taiji configs, no manual fixups. Departures from these defaults (different patterns, extra columns) are supported but require explicit flags, because silent schema drift is the single most common reason Taiji runs fail on otherwise-fine data.

HiC is optional because many bulk Taiji runs proceed with RNA+ATAC alone; when HiC is included, most projects use a shared per-build reference loops file rather than per-sample chromatin conformation data, so the default path lives in the genome registry, not in the samples CSV. Per-sample HiC is supported via the `hic_path` column for the minority of projects that have it.
