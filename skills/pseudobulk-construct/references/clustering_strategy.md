# Clustering strategy

## Why WNN for multiome

Seurat's weighted nearest-neighbor (WNN) clustering learns per-cell weights for each modality and builds a joint graph that respects whichever modality is most informative for that particular cell. For multiome (RNA + ATAC from the same cells), this matters because some cell types are cleanly separated by expression but not accessibility (and vice versa). Clustering on RNA alone can collapse regulatory-distinct populations that share a transcriptomic program; clustering on ATAC alone can split populations that have the same accessibility landscape but different activation states.

The pipeline matches the Seurat v5 WNN vignette: RNA is normalized with `NormalizeData` → `ScaleData` → `RunPCA(50 PCs)`; ATAC goes through `RunTFIDF` → `FindTopFeatures(min.cutoff="q5")` → `RunSVD`; then `FindMultiModalNeighbors(reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:30))` builds the joint graph as `wsnn`, and `FindClusters(graph.name = "wsnn", resolution = r)` assigns cluster labels.

The first LSI component is dropped (`dims.list[[2]]` starts at 2) because in scATAC it typically correlates with read depth rather than biology — this is Signac's standard recommendation.

## Scale-aware resolution search

A naïve `FindClusters(resolution = 0.8)` call doesn't respect dataset size. For a 5k-cell pilot it produces a handful of large clusters; for a 1M-cell atlas it produces a couple thousand tiny ones. The user's spec targets a mean cluster size of ~100-300 cells, so the skill binary-searches `resolution` to land mean cluster size inside that band.

### Seeding

The initial resolution is seeded from the target number of clusters:

```r
target_n_clusters <- N_cells / target_cluster_size
r0 <- 0.15 * log2(N_cells / 1000) * sqrt(target_n_clusters / 20)
```

This heuristic is calibrated so that:

- At 5k cells targeting 200/cluster (~25 clusters), `r0 ≈ 0.3`.
- At 50k cells targeting 200/cluster (~250 clusters), `r0 ≈ 1.0`.
- At 1M cells targeting 200/cluster (~5000 clusters), `r0 ≈ 4.2` — clamped to 3.0.

### Iteration

After evaluating the seed, the skill adjusts in a bisection-style loop:

- If mean cluster size > target → cluster count is too low → increase resolution (move halfway toward `hi`, or multiply by 1.5 if `hi` is still the prior cap of 5.0).
- If mean cluster size < target → decrease resolution (move halfway toward `lo`, or multiply by 0.67).

Each iteration reuses the already-computed neighbor graph and just rerusns Louvain, so each step is cheap (seconds for a 100k-cell object vs. minutes for the initial `FindMultiModalNeighbors`).

Cap is 8 iterations. If the search hasn't converged by then, the last iteration is kept with a warning and the full trace is written to `resolution_trace.json` for audit.

### Acceptance band

When `--target-cluster-size` is the default 200, the accepted mean-size band is exactly the spec's `[100, 300]`. For any other target `t`, the band is `[0.5*t, 1.5*t]` — a symmetric 50% tolerance that makes sense as a general rule. If you need a tighter band for a specific dataset, override `--target-cluster-size` to a narrower target and let the skill re-center around it.

### When it doesn't converge

The usual cause is a dataset with a highly skewed cluster-size distribution: a couple of very large clusters and many tiny ones, where no single resolution puts the *mean* in the target band. `mean_size` isn't robust; `median_size` would be, but it shifts the interpretation of the spec in a way the user didn't ask for.

When the search doesn't converge the right move is usually to:

1. Inspect `resolution_trace.json` to see the size distribution at each step.
2. If one very large cluster is pulling the mean up, run with a slightly higher resolution than the last one tried — the big cluster will split, and the *new* mean will fit.
3. If the skew is extreme, consider subclustering the dominant cluster separately with Seurat's `FindSubCluster()` and running the pseudobulk step on the refined labels. The skill doesn't do this automatically because it's a manual curation decision.

## RNA-only and ATAC-only paths

- **RNA-only** (`--clustering-signal rna`): standard Seurat `FindNeighbors(reduction = "pca", dims = 1:30)` + `FindClusters`. Use this for scRNA-only objects, or when WNN fails and you want an explicit fallback.
- **ATAC-only** (`--clustering-signal atac`): Signac LSI pipeline + `FindNeighbors(reduction = "lsi", dims = 2:30)` + `FindClusters`. This is the canonical path for separate-assay data that has already been co-embedded via `integrate_atac` — clustering happens in the ATAC LSI space with a transferred-label column driving the RNA-side metadata axis.

The ATAC path enforces a safety check: it refuses to run unless the object carries a transferred-label column (default `predicted.id`), so separate-assay data that *hasn't* been integrated yet produces a loud error pointing at <https://stuartlab.org/signac/articles/integrate_atac> rather than silent bad clusters.

## Small-cluster filter

Any cluster with fewer than `--min-cluster-cells` cells (default 20) is dropped outright before pseudobulks are built. 20 is both the spec's floor and a pragmatic threshold: below ~20 cells, summed RNA counts are too sparse to distinguish from noise and MACS2 produces erratic peak calls. A user who wants to retain smaller clusters should reconsider whether they really want bulk Taiji on that group, or run a dedicated rare-cell analysis instead.

The filter runs separately at two levels:

1. At the cluster level: a cluster with <20 cells in total is dropped.
2. At the (cluster × metadata_col × metadata_value) level: a sub-group with <20 cells in that specific metadata value is dropped, even if the parent cluster passes.

This mirrors the spec. If `--metadata-cols donor,condition` is set and a cluster has 100 cells across 10 donors, the per-donor sub-groups will each have ~10 cells and will all be filtered out — the right behavior, because per-donor pseudobulks at that depth are worthless.

## Limitations worth calling out

- The resolution search uses **mean** cluster size as the acceptance criterion, not a percentile. Datasets with a single huge cluster may converge to a "passing" resolution where the median cluster is far from the target. Inspect `resolution_trace.json` if something looks off.
- The auto-detected metadata column list is heuristic (2-20 unique values, skip numerics with many levels, skip known technical columns). It can miss a meaningful column that happens to have 21 categories, or include a noisy one that happens to have 5. Always sanity-check the detected list before burning R time — the skill logs it clearly.
- The WNN graph is built on the default 30 RNA PCs + 29 ATAC LSI components (2:30). For very small datasets, 30 components is too many and the PCA/SVD may silently overfit; for very large atlases, 30 might miss biology. Tuning `dims.list` is currently not exposed as a flag — if this becomes a recurring issue, lift it out of the R script.
