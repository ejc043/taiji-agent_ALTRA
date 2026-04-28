---
name: pseudobulk-construct
description: Convert a single-cell object (.rds Seurat, .h5ad AnnData, or .h5mu MuData) into pseudobulk GeneQuant TSVs + per-cluster narrowPeak files that feed the bulk Taiji pipeline. Use this skill whenever the user wants to pseudobulk, aggregate, collapse, summarize, or bulk-ify a single-cell dataset — especially before handing off to build-taiji-input. Runs Seurat WNN clustering with scale-aware resolution tuning (targeting ~100-300 cells/cluster, adjusted for total cell count), filters small clusters (<20 RNA cells or <20 ATAC cells when assays are separate; <20 cells when paired), sums raw counts per (cluster × metadata) combination for RNA, and calls MACS2 peaks per cluster for ATAC. Trigger on phrases like "pseudobulk this", "construct pseudobulks", "cluster and aggregate", "collapse single cell to bulk", "make the scRNA bulk-compatible", "per-cluster peaks", "prep this for Taiji from single cell", or anywhere the user needs to bridge SC → bulk Taiji.
---

# Pseudobulk construction (scRNA + scATAC → bulk Taiji inputs)

## What this skill produces

Given a single-cell object and (for ATAC) a fragments file, emit per-cluster × per-metadata-group pseudobulk outputs in the exact layout `build-taiji-input` consumes:

```
output_dir/
├── rna/
│   ├── <cluster>__<metadata_col>__<value>.tsv         # gene_id<TAB>count
│   └── ...
├── atac/
│   ├── <cluster>__<metadata_col>__<value>_peaks.narrowPeak
│   └── ...
├── clusters.csv                 # per-cell: barcode, seurat_cluster, metadata cols
├── resolution_trace.json        # the resolution-search history (for audit)
└── manifest.tsv                 # samples sheet ready for build-taiji-input
```

The `manifest.tsv` is the bridge: each row is one (cluster × metadata_value) "sample" with paths to the RNA TSV, ATAC narrowPeak, cohort/group columns derived from the metadata, and genome tag. Feed it directly to `build-taiji-input --samples manifest.tsv`.

## What this skill refuses to do

- Bulk inputs. If `detect-dataset-type` classifies the directory as bulk, this skill errors and redirects to `build-taiji-input`.
- `sc-undetermined` single-cell inputs. If the skill can't tell whether RNA and ATAC are paired (multiome) or separate, it will not guess — it asks the user to clarify, or to provide a `.h5mu` / rename files with an explicit `multiome` / `rna` / `atac` token.
- `separate-assay` inputs that haven't been co-embedded yet. The skill checks for a transferred-label column in the object's `@meta.data` (default: `predicted.id`, configurable via `--transferred-label-col`) and refuses if it's missing — pointing the user at the Signac integrate_atac workflow first: <https://stuartlab.org/signac/articles/integrate_atac>.

## Inputs

| Flag                        | Purpose                                                                               |
|-----------------------------|---------------------------------------------------------------------------------------|
| `--input`                   | Path to the single-cell object (.rds / .h5ad / .h5mu). Required.                      |
| `--fragments`               | Path to `fragments.tsv.gz` (plus `.tbi`). Required if peaks are to be called.         |
| `--genome`                  | hg38 / hg19 / mm10 / mm39. Drives Signac::CallPeaks `effective.genome.size` and downstream-Taiji metadata. |
| `--metadata-cols`           | Comma-separated metadata columns to stratify on (e.g. `donor,condition`). Default: auto-detect categorical columns with 2-20 unique values, skipping cell-barcode-like columns. |
| `--output-dir`              | Where the `rna/`, `atac/`, `clusters.csv`, `manifest.tsv` land.                       |
| `--target-cluster-size`     | Target mean cluster size for the resolution search. Default: `200` (range `[100, 300]`). |
| `--min-cluster-cells`       | Drop clusters with fewer than this many cells. Default: `20`. Applied per modality for separate-assay, combined for multiome/multiome. |
| `--clustering-signal`       | `wnn` (default for multiome), `rna`, or `atac`. Auto-selected from detect-dataset-type's `sc_modality`. |
| `--transferred-label-col`   | Column in `@meta.data` that holds Signac-transferred labels for separate-assay. Default: `predicted.id`. |
| `--peak-caller`             | `macs2` or `macs3`. Auto-detected from PATH (prefers macs3); pseudobulk.py resolves the full binary path and passes it to Signac::CallPeaks via `macs2.path`. |
| `--skip-data-type-check`    | Bypass the detect-dataset-type pre-flight (use only if you know the input is SC).     |
| `--rna-only` / `--atac-only`| Restrict output to a single modality (useful for RNA-only scRNA objects).             |
| `--dry-run`                 | Plan the resolution search and list planned outputs without invoking peak calling.    |

## Dependencies the skill expects in the execution environment

This skill is not self-contained — it shells out to R and MACS2. The Python entry point checks these at startup and emits a clear error if they're missing:

- `Rscript` ≥ 4.2 on `PATH`.
- R packages: `Seurat` ≥ 5.0, `Signac` ≥ 1.12 (for ATAC), `Matrix`, `GenomicRanges`, `dplyr`. For `.h5ad` input: `SeuratDisk` or `anndata` + `sceasy`. For `.h5mu`: `MuDataSeurat`.
- `macs2` (or `macs3`) on `PATH`.
- `tabix` / `bgzip` available if fragments need to be subsetted.

On the user's SLURM environment these are typically loaded via `module load r/4.3 macs2/2.2.9.1` or an analogous Conda env. The skill does not install anything on its own.

## Interaction pattern

1. **Always run `detect-dataset-type` first** (the skill does this automatically unless `--skip-data-type-check` is passed). Use the returned `sc_modality` to pick the clustering signal: `multiome` → WNN; `separate-assay` → ATAC LSI clustering with prior label transfer required; `sc-undetermined` → refuse and ask for clarification.
2. **Confirm metadata columns with the user before spending R time**. The auto-detect step is intentionally conservative (2-20 unique categorical values, not a barcode, not continuous). Print the detected list and prompt for confirmation when running interactively. In non-interactive CI runs (`--yes`), the detected list is used as-is.
3. **Expose the resolution trace**. The binary search writes every (resolution, n_clusters, mean_size) tuple to `resolution_trace.json`. If the user is surprised by the chosen resolution, this is the audit artifact.
4. **Surface the cluster-size filter explicitly**. Clusters dropped by the 20-cell floor are logged with their sizes — small clusters often correspond to doublet pockets or lowly-represented cell types, and the user sometimes wants to know rather than have them silently disappear.
5. **On WNN failure, don't silently fall back**. If the object lacks a `ChromatinAssay` or `FindMultiModalNeighbors` errors, the skill surfaces the error and suggests running `--clustering-signal rna` explicitly rather than pretending WNN succeeded.

## Interoperability with the rest of the Taiji agent

- Input gated by `detect-dataset-type` (same soft-import pattern `build-taiji-input` uses).
- Output `manifest.tsv` has the exact columns `build-taiji-input --samples` expects (`name`, `rna_seq`, `atac_seq`, optional `hic`, `cohort`, `group`), so the two skills chain with zero glue code: `pseudobulk-construct` → `build-taiji-input` → bulk Taiji. `cohort` defaults to the metadata-column name, `group` to the metadata-value; both can be overridden.
- `--hic` is not produced (single-cell HiC is out of scope for this skill).

## When to reach for the references

- `references/clustering_strategy.md` — why WNN for multiome, how the scale-aware resolution search works, what to do when the search doesn't converge.
- `references/peak_calling.md` — MACS2 invocation per cluster, fragments-file requirements, what to do if fragments are missing.
- `references/output_format.md` — the exact on-disk layout, how the manifest maps to `build-taiji-input` columns, edge cases (metadata with missing values, clusters that survive one modality but not the other).

## Why the defaults are what they are

- **Scale-aware resolution**: the user pointed out that targeting a fixed ~200 cells/cluster with a naïve default-0.8 Louvain step won't work for a 1M-cell atlas (you'd get ~100 clusters that are each ~10k cells) or a 5k-cell pilot (you'd get 4 clusters of ~1200 each). The binary search seeds the initial resolution from N_cells / target_cluster_size and narrows from there; see `references/clustering_strategy.md` for the formula.
- **WNN over single-modality clustering for multiome**: this is the standard Seurat v5 recommendation for paired RNA+ATAC data. It gives noticeably better separation of rare cell types than clustering on either modality alone, at the cost of ~2× runtime for `FindMultiModalNeighbors`.
- **Hard 20-cell floor, no soft option**: peaks called from <20 cells are dominated by noise; summed counts from <20 cells are too sparse to be a meaningful "bulk" sample. Lowering the floor would silently degrade the quality of downstream Taiji.
- **Sum raw counts, not normalized**: bulk RNA-seq tools expect raw counts. Summing normalized values is never the right move.
- **Per-cluster Signac::CallPeaks, not global peaks + subset counting**: calling peaks on each pseudobulk individually catches cluster-specific regulatory elements that would be smoothed out by a single global call. It's the slower path but the biologically correct one. Going through Signac (instead of a direct gunzip|awk|macs pipe) buys us proper per-group `CreateFragmentObject` plumbing and per-group barcode reconciliation against `fragments.tsv.gz` — Seurat objects often carry suffixes (`AAACGCT-1`, `ATAC_AAACGCT`) that the raw 10x fragments file doesn't, and `call_peaks.R` strips them automatically (≥95% overlap required after correction, fail-loud otherwise).
