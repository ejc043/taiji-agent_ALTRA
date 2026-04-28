# Co-embedding strategy

## What this skill implements

End-to-end automation of the [Stuart 2019 / Signac integrate_atac vignette](https://stuartlab.org/signac/articles/integrate_atac), workflow B (full co-embedding), with one deliberate omission: cell-type label transfer.

## The vignette has two flows; we implement the second

| Flow | Steps | Output | Used here? |
|------|-------|--------|------------|
| **A: Label transfer only** | FindTransferAnchors + TransferData(refdata = pbmc.rna$seurat_annotations) | ATAC object gains `predicted.id`; two separate UMAPs | No — we deliberately skip the annotation projection |
| **B: Full co-embedding** | Flow A + impute RNA into ATAC + merge + joint PCA/UMAP | One Seurat object, both modalities in shared space, `assay` meta column tags origin | Yes — this is the entire skill |

## Why we skip the annotation transfer

`TransferData(refdata = <reference>$cell_type)` projects named labels from a reference atlas onto query cells. The output looks like a clean cluster annotation but is silently brittle:

- **Reference choice matters.** PBMC reference vs. tissue-matched reference vs. your own paired RNA → all three give different labels for the same query cell. There's no signal-to-noise ratio that tells you which is right.
- **Bad anchors are invisible.** Low `prediction.score.max` cells get labeled with the same confidence as high-score cells; downstream you can't tell which is which without inspecting the score distribution by hand.
- **Closely-related types confuse silently.** Naive vs. memory B cells, CD4 vs. CD8 effector subsets, monocyte subtypes — these get cross-classified at high confidence and produce misleading per-cell-type pseudobulks.

The mechanical co-embedding step (steps 1–8 in `coembed.R`) is deterministic given inputs and has clearer failure modes (anchor count too low, joint UMAP collapses, etc.). The cell-type *labeling* step is where biological judgment lives, and it stays user-side.

If you want named labels for downstream Taiji stratification (e.g. `--cohort-col cell_type` in `pseudobulk-construct`), assign them upstream:

```r
# After coembed-construct, before pseudobulk-construct:
coembed <- readRDS("output/coembed.rds")
# ... your reference-anchored TransferData call here ...
coembed$cell_type <- transfer_data_results$predicted.id
saveRDS(coembed, "output/coembed.rds")
```

Then pass `--metadata-cols cell_type,condition,tissue` to `pseudobulk-construct`.

## Why CCA reduction for FindTransferAnchors

The vignette uses `reduction = "cca"` for cross-modality anchoring. CCA finds linear combinations of features that maximize correlation across the two datasets — well-suited for the gene-expression vs. gene-activity comparison where the feature spaces are nominally the same (both indexed by gene symbol) but the values come from different molecular measurements.

For same-modality integration (e.g. cross-sample RNA-RNA), `RPCA` is usually faster and equally good. For cross-modality (RNA-ATAC), CCA is the established choice and what Signac documents.

## Why LSI as `weight.reduction` for TransferData

`TransferData` projects values from reference to query using anchor weights. The query-side weight reduction determines which cells are similar enough to share the projected value. For ATAC queries, the LSI reduction (computed before GeneActivity) is the canonical "neighborhood structure of ATAC cells" — using it preserves ATAC-side biology when imputing RNA.

We skip the first LSI component (`--lsi-skip-first 1`) because it typically correlates with sequencing depth rather than biological variation. This matches the vignette's `dims = 2:30`.

## Resolution binary search on the joint PCA

Identical algorithm to `pseudobulk-construct/scripts/load_and_cluster.R::binary_search_resolution`. Targets ~200 cells/cluster (range [100, 300]) by default; max 8 iterations. The seed resolution is `0.15 * log2(N_cells / 1000) * sqrt(target_n_clusters / 20)` clamped to [0.05, 3.0]. `--resolution <r>` overrides the search if the user wants a fixed value.

The cluster count this targets is calibrated for biological cell-type granularity. For a 10k-cell coembed, ~50 clusters; for 100k cells, ~500. You can adjust by passing `--target-cluster-size` (e.g. `--target-cluster-size 500` for fewer, larger clusters).

## When this skill refuses

- **Multiome input** (single .h5mu / cellranger-arc filtered_feature_bc_matrix + atac_fragments.tsv): doesn't need co-embedding — same cells, just use `pseudobulk-construct --clustering-signal wnn`.
- **Pre-co-embedded input** (one .rds with both modalities already merged + clustered): also doesn't need this skill — `pseudobulk-construct --use-existing-clusters` directly.
- **Bulk input**: the `detect-dataset-type` gate refuses upfront.

## Output schema details

`coembed_summary.json` keys:

```json
{
  "rna_input":  "/abs/path/rna.rds",
  "atac_input": "/abs/path/atac.rds",
  "genome":     "hg38",
  "n_cells_rna":   12453,
  "n_cells_atac":  11827,
  "n_cells_total": 24280,
  "n_anchors": 8743,
  "n_genes_used": 2000,
  "chosen_resolution": 0.83,
  "resolution_trace": [{"iter": 1, "resolution": 0.51, ...}, ...],
  "n_clusters_total": 47,
  "n_clusters_kept":  43,
  "dropped_clusters": ["44","45","46","47"],
  "metadata_cols": ["genotype","tissue"],
  "output_path": "/abs/path/coembed.rds"
}
```

## What needs first-run validation

The R script is syntax-clean but the sandbox can't run Seurat/Signac/EnsDb. First-run on Eunice's SLURM env should verify:

1. The `Annotation(atac) <- annotations` setter works against the user's actual ChromatinAssay (Signac API has been stable on this since v1.0).
2. The TransferData call's `weight.reduction = atac[["lsi"]]` doesn't throw a dim-mismatch when LSI dims aren't sequential.
3. EnsDb packages load (they should be in the conda env after rerunning `bin/install.sh --profile sc`; the new entries in `environment.sc.yml` add them).
4. The QC UMAP renders even on a headless cluster node (no DISPLAY) — `ggsave` does, but verify.
