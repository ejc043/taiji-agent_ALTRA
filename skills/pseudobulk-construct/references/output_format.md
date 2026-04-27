# Output format

## Directory layout

After a successful run, `--output-dir` looks like:

```
output_dir/
├── clusters.csv                        # per-cell assignment + metadata
├── resolution_trace.json               # every (resolution, n_clusters, mean_size) tried
├── groups_plan.json                    # every (cluster x meta_col x value) group spec
├── per_cluster_barcodes/
│   ├── cluster0.txt                    # barcodes (one per line) used for peak calling
│   ├── cluster0__donor__D1.txt
│   └── ...
├── rna/
│   ├── cluster0__donor__D1.tsv         # gene_id<TAB>count, no header
│   ├── cluster0__donor__D2.tsv
│   └── _index.json                     # machine-readable map of group -> TSV path
├── atac/
│   ├── cluster0__donor__D1_peaks.narrowPeak
│   ├── cluster0__donor__D1_summits.bed
│   └── ...
└── manifest.tsv                        # feeds build-taiji-input --samples
```

## manifest.tsv — the handoff to build-taiji-input

The manifest has exactly the columns `build-taiji-input --samples` expects:

| Column    | Meaning                                                                         |
|-----------|---------------------------------------------------------------------------------|
| `name`    | Unique sample name. Defaults to `<cluster>__<meta_col>__<value>`.              |
| `cohort`  | Grouping axis. Defaults to the metadata column name.                            |
| `group`   | Grouping value. Defaults to `cluster<N>__<meta_value>`.                         |
| `rna_seq` | Absolute path to the GeneQuant TSV. Empty if `--atac-only`.                     |
| `atac_seq`| Absolute path to the narrowPeak file. Empty if `--rna-only` or MACS2 emitted no peaks. |
| `genome`  | Genome tag (hg38, mm10, etc.) — copied from `--genome`.                         |

After a successful pseudobulk run:

```bash
python build-taiji-input/scripts/build_taiji_input.py \
  --samples  <output_dir>/manifest.tsv \
  --data-dir <output_dir> \
  --genome   mm10 \
  --out      taiji_input.xlsx
```

No extra glue is required — all the TSV / narrowPeak paths in the manifest are absolute, so `build-taiji-input` doesn't need to know about the pseudobulk step at all.

## GeneQuant TSV format

```
ENSMUSG00000000001	42
ENSMUSG00000000003	0
ENSMUSG00000000028	117
...
```

Two columns, tab-separated, no header. Gene IDs are whatever was in the Seurat object's `rownames(counts)` — typically Ensembl IDs from cellranger. Counts are 32-bit integer sums of the raw (unnormalized) per-cell counts across the group's cells.

If your Seurat object uses gene symbols instead of Ensembl IDs, the TSV will have symbols — which will flow through `build-taiji-input` and into bulk Taiji as-is. Confirm the downstream genome annotation matches (Taiji's GENCODE gtf needs to use the same identifier space).

## narrowPeak format

Standard ENCODE narrowPeak (BED6+4):

```
chrom  start  end  name  score  strand  signalValue  pValue  qValue  peak
```

MACS2 emits `name` as `<group_name>_peak_N`. `signalValue` is the fold-enrichment, `pValue` and `qValue` are `-log10`. `peak` is the summit offset within the peak interval when `--call-summits` is passed.

## clusters.csv

Per-cell labels plus metadata columns used for stratification:

```
barcode,seurat_cluster,donor,condition
AAACGGGAGCGACGTA-1,0,D1,healthy
AAACGGGAGCGACGTA-2,0,D1,healthy
AAACGGGAGCTTACGT-1,3,D2,disease
...
```

Cells from dropped clusters (below `--min-cluster-cells`) are not included. The file is stable sort order (Seurat's `@meta.data` row order, which is insertion order), so diffs between runs only show what actually changed.

## resolution_trace.json

```json
{
  "signal": "wnn",
  "target_cluster_size": 200,
  "chosen_resolution": 0.825,
  "n_cells": 28143,
  "iterations": [
    {"iter": 1, "resolution": 0.35, "n_clusters": 42,  "mean_cluster_size": 670.1},
    {"iter": 2, "resolution": 0.85, "n_clusters": 104, "mean_cluster_size": 270.6},
    {"iter": 3, "resolution": 1.05, "n_clusters": 138, "mean_cluster_size": 203.9}
  ]
}
```

This file is the primary audit artifact for "why did the skill pick resolution X?" If a reviewer disagrees with the clustering, they can see the full search trajectory and rerun with an explicit `--target-cluster-size` pointing at the mean-size band they prefer.

## groups_plan.json

```json
{
  "signal": "wnn",
  "chosen_resolution": 0.825,
  "kept_clusters": ["0", "1", "3", "4", ...],
  "dropped_clusters": ["19"],
  "metadata_cols": ["donor", "condition"],
  "min_cluster_cells": 20,
  "groups": [
    {
      "name": "cluster0__donor__D1",
      "cluster": "0", "metadata_col": "donor", "metadata_value": "D1",
      "n_cells": 487, "cohort": "donor", "group": "cluster0__D1"
    },
    ...
  ]
}
```

This is what the Python orchestrator reads to drive MACS2 and to emit `manifest.tsv`. If you need to post-process the plan (merge two metadata values, rename a cohort, drop a specific group), edit this file and rerun `pseudobulk.py --resume` — though `--resume` is not implemented in v1, so for now the only way to post-process is to edit the manifest directly.

## Edge cases

- **No metadata columns detected**: the per-cluster-only path kicks in. `manifest.tsv` will have one row per cluster with `cohort=all` and `group=cluster<N>`. `build-taiji-input` handles this fine — the Active sheet just has a single cohort.
- **A metadata value appears in only one cluster**: the skill still emits a group for it. The resulting pseudobulk will be a single cluster's worth of cells for that metadata value; downstream Taiji treats it as one sample.
- **All cells in a cluster have missing (NA) values for a metadata column**: that cluster is skipped for that column but still contributes to other metadata-column stratifications. Warning logged per skipped group.
- **An RNA TSV is written but MACS2 produces no peaks for the same group**: common for rare, low-fragment clusters. The manifest row will have `rna_seq` set and `atac_seq` empty. `build-taiji-input` with `--no-hic --strict` will reject this row; without `--strict` it'll proceed with RNA-only.
- **An empty narrowPeak is written**: the skill treats a narrowPeak of size 0 as equivalent to "no peaks" and leaves `atac_seq` empty in the manifest rather than pointing to a useless file.
