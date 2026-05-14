# Skills — operational guidance for Claude

Supplement to the root CLAUDE.md. Covers recurring failure modes and mandatory checks discovered during live runs. Root CLAUDE.md takes precedence on architecture; this file governs skill invocation order and QC gates.

## Mandatory QC before co-embedding

**sc-qc must always run before coembed-construct.** Do not assume a `_prep.rds`, `_clean.rds`, or any other pre-named object is fully QC'd. Verify the distributions explicitly before coembedding.

### RNA — required checks

```r
rna <- readRDS("obj.rds")
cat("nFeature_RNA: min=", min(rna$nFeature_RNA), " max=", max(rna$nFeature_RNA), "\n")
cat("percent.mt: max=", round(max(rna$percent.mt), 2),
    " | n>=10%:", sum(rna$percent.mt >= 10), "\n")
```

Canonical filters (matches `JohnChang/TCF1/analysis/coembed_preprocess.R`):

| Filter | Operator | Threshold | Why |
|--------|----------|-----------|-----|
| `nFeature_RNA` | `> 200` and `< 5000` | strict exclusive | empty droplets (low) / doublets (high) |
| `percent.mt` | `< 10%` | strict exclusive | dead / lysed cells |

**`nCount_RNA` is not filtered.** The lab reference (`coembed_preprocess.R`) does not include a UMI count minimum, and we match that. Do not add a nCount filter without explicit instruction.

Common gap: `_prep` objects are often filtered for nFeature range only, leaving percent.mt unfiltered. That is the check most likely to be missing — always verify `sum(rna$percent.mt >= 10)` before coembedding.

### ATAC — required checks

```r
cat("nCount_peaks: min=", min(atac$nCount_peaks), "\n")
cat("nFeature_peaks: min=", min(atac$nFeature_peaks), "\n")
cat("TSS.enrichment: min=", min(atac$TSS.enrichment), "\n")
cat("nucleosome_signal: max=", max(atac$nucleosome_signal), "\n")
cat("blacklist_ratio: max=", max(atac$blacklist_ratio), "\n")
```

Standard filters:

| Filter | Threshold |
|--------|-----------|
| `nCount_peaks` | ≥ 1000 |
| `nFeature_peaks` | ≥ 500 |
| `TSS.enrichment` | ≥ 2 |
| `nucleosome_signal` | < 4 |
| `blacklist_ratio` | < 0.05 |

ATAC QC metrics (TSS, nucleosome_signal, blacklist_ratio) must be pre-computed upstream with Signac before sc-qc can apply them. If absent from meta.data, sc-qc will warn and skip — compute them first.

### Enforced pipeline order

```
sc-qc (RNA: --input rna.rds)
sc-qc (ATAC: --input atac.rds)          ← both required before coembedding
        |
coembed-construct --rna rna_filtered.rds --atac atac_filtered.rds
        |
chooseR_clustering.R                     ← resolution selection
        |
pseudobulk-construct --use-existing-clusters
```

Never pass unverified objects directly to coembed-construct. The downstream cost of filtering-after-coembedding (re-running 2+ hour jobs) far exceeds the cost of a 5-minute QC check upfront.

## chooseR resolution selection — interpretation guide

`chooseR_clustering.R` sweeps Louvain resolutions and selects the **plateau onset**, not the global silhouette maximum.

**Algorithm:**
1. For each resolution: `FindClusters` → subsample 2000 cells → Euclidean distances in 30-dim joint PCA → `cluster::silhouette()` → record mean silhouette width.
2. Detect a **leading spike**: if `sil[res=0.1] > 1.5 × sil[res=0.2]`, the first point is flagged as a spike and excluded from selection. The spike reflects the two dominant compartments (e.g. spleen vs siIEL) separating at coarse resolution — biologically real but too coarse for Taiji input.
3. Among remaining points, find the longest run of consecutive resolutions where |Δsil| < 0.05 (the stable plateau). Select the **midpoint** of that run — equidistant from the leading edge (where the curve may still be descending) and the trailing edge (where collapse begins).
4. Fallback if no plateau: highest silhouette among non-spike candidates.

**Reading the silhouette plot:** the spike (if flagged) is marked with a grey ×. The chosen resolution is marked "plateau onset" in red. The plateau region is where silhouette has stabilised before collapsing to ~0 at high resolutions.

**Typical curve shape for co-embedded spleen/siIEL data (81k cells):**
- res=0.10: spike (~0.32) — excluded
- res=0.60–1.20: stable plateau (~0.03–0.07) — selection zone
- res≥2.0: collapse to ~0

**Cluster×assay composition check** (always run after chooseR):

```r
obj <- readRDS("coembed_clustered.rds")
tab <- table(obj$seurat_clusters, obj$assay)
print(tab)
# Clusters with 0 ATAC cells are RNA-only satellite populations — expected
# when RNA N >> ATAC N. They contribute RNA TSVs only; no narrowPeak.
# Clusters with 1–19 ATAC cells will produce unreliable peaks — treat as RNA-only.
```

**Env vars to tune selection without editing the script:**
- `CHOOSER_SPIKE_FACTOR` (default 1.5): raise to be less aggressive about spike exclusion
- `CHOOSER_PLATEAU_TOL` (default 0.05): lower to require a tighter plateau before selecting

## _prep / _clean file naming — do not trust blindly

`_prep.rds`, `_clean.rds`, `_filtered.rds` indicate the file was processed upstream, not that it passed all required QC gates. Observed pattern: prep objects in this project had nFeature [200, 5000] applied but **percent.mt was not filtered**, leaving ~120 dead/damaged cells in 60k. Always run the RNA verification block above — specifically `sum(rna$percent.mt >= 10)` — before passing to coembed-construct.

## Gene activity for anchor transfer — use prep object's pre-built RNA assay when present

`_prep.rds` ATAC objects from this lab's pipeline typically carry a pre-built `"RNA"` assay (full-genome gene activity, ~21k genes). `coembed_preprocess.R` used this directly: `query.assay = "RNA"`.

**Do NOT recompute GeneActivity and store as "ACTIVITY" when a pre-built "RNA" assay exists.** Recomputing restricted to `VariableFeatures(rna)` (2000 genes) and only keeping the intersection drops ~18% of anchoring features (1634/2000 in D7). The lab reference used the full 21808-gene assay, giving FindTransferAnchors access to the complete variable feature set.

`coembed.R` auto-detects: if `"RNA" %in% Assays(atac)`, it normalizes/scales the existing assay and passes `query.assay = "RNA"`. If absent, it falls back to fresh GeneActivity.

## UMAP method — umap-learn with correlation metric

The lab reference used Python umap-learn (via reticulate) with `metric = "correlation"`. Seurat v5 defaults changed to uwot/cosine. The scripts enforce umap-learn when available; it must be installed in the Python env reticulate uses.

On this HPC: umap-learn lives in the taiji-agent conda env. Reticulate must point to it via `RETICULATE_PYTHON=/stg3/data1/eunice/.local/share/mamba/envs/taiji-agent/bin/python`. A stale virtualenv at `~/.virtualenvs/r-reticulate/` can shadow this — `coembed.R` sets `RETICULATE_PYTHON` explicitly before calling reticulate.

If umap-learn import still fails, the script falls back to uwot/cosine and logs a message. The UMAP shape will differ from the reference but clustering and downstream Taiji inputs are unaffected.

## Signac fragment path — always provide --fragments

ATAC `.rds` objects built on another machine embed an absolute fragment path that will not resolve at the new location. coembed.R will rebuild the Fragment handle automatically when `--fragments` is supplied, but it will fail with an opaque tabix error if the path is stale and no override is given. Always pass `--fragments` explicitly unless you built the object on the current machine.
