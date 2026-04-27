# Classification rules

## Extension → class mapping

| Extension        | Class        | Notes                                                              |
|------------------|--------------|--------------------------------------------------------------------|
| `.tsv`           | bulk         | RNA-seq GeneQuant. Weak on its own — see "`.tsv` ambiguity".       |
| `.narrowPeak`    | bulk         | ATAC-seq peaks (ENCODE BED6+4 format).                             |
| `.bedpe`         | bulk         | HiC chromatin loops (6-column paired BED + optional score/metadata). |
| `.h5ad`          | single-cell  | AnnData / scanpy. HDF5-based, never used for bulk.                 |
| `.rds`           | single-cell  | R serialized object. In the Taiji context, effectively always Seurat / SingleCellExperiment. |
| `.h5mu`          | single-cell  | MuData. Multi-modal by construction, so always classified as `sc_modality=multiome`. |

All matching is **case-insensitive** and **.gz-aware**, so `.tsv`, `.TSV`, `.tsv.gz`, `.Narrowpeak.gz`, etc. all match.

## Filename-pattern overrides (SC takes priority over BULK)

Some 10x / cellranger outputs contain a `.tsv` or `.tsv.gz` suffix but are unambiguously single-cell. These filename patterns override the extension-level rule so that a cellranger-arc directory containing `atac_fragments.tsv.gz` is not miscounted as bulk RNA-seq:

| Filename substring                  | Class       | Context                                                   |
|-------------------------------------|-------------|-----------------------------------------------------------|
| `fragments.tsv`                     | single-cell | scATAC fragments (10x / Signac).                          |
| `barcodes.tsv`                      | single-cell | 10x GEX/ATAC cell barcodes.                               |
| `features.tsv`                      | single-cell | 10x GEX features.                                         |
| `genes.tsv`                         | single-cell | Older 10x GEX genes file (pre-CellRanger 3.0).           |
| `matrix.mtx`                        | single-cell | Sparse count matrix (10x / STARsolo).                    |
| `filtered_feature_bc_matrix.h5`     | single-cell | 10x GEX or Multiome combined output.                     |
| `filtered_peak_bc_matrix.h5`        | single-cell | 10x ATAC peak-barcode matrix.                            |
| `raw_feature_bc_matrix.h5`          | single-cell | 10x raw counts (unfiltered).                             |

Matching is a case-insensitive substring check on the basename.

## Single-cell modality sub-classification

When classification is `single-cell`, the result carries an additional `sc_modality` field that says whether RNA and ATAC come from the SAME cells (no co-embedding needed) or DIFFERENT cells (co-embedding required before any regulatory-network analysis).

Tier logic, first match wins:

1. **Any `.h5mu` file present** → `multiome`. MuData is designed to hold paired multi-modal data; there is no single-modality use case in this lab's pipelines.
2. **cellranger-arc signature**: both `filtered_feature_bc_matrix.h5` and `atac_fragments.tsv` present in the scanned tree → `multiome`. This is the canonical 10x Genomics Multiome output layout where RNA and ATAC are profiled from the same cells.
3. **Explicit multiome token** in any filename — one of `multiome`, `multi_omic`, `multi-omic`, `multimodal`, `arc` → `multiome`. Catches hand-named `.h5ad` / `.rds` files like `pbmc_multiome.h5ad`.
4. **Both RNA and ATAC hints** across filenames: any of `rna` / `gex` / `expression` / `scrna` / `gene_expression` AND any of `atac` / `scatac` / `peak` / `chromatin_accessibility` → `separate-assay`. The common case: `pbmc_rna.h5ad` + `pbmc_atac.h5ad` in the same directory.
5. Otherwise → `sc-undetermined`.

### Why tier order matters

The tiers are ordered from strongest evidence (file format or canonical directory layout) to weakest (filename tokens). A naïve implementation that checked RNA and ATAC hints first would flag a cellranger-arc directory as `separate-assay` (because the filenames contain both `feature` and `atac` substrings) — exactly the opposite of the truth. Structural signals always win.

### How separate-assay users should proceed

`separate-assay` is the case the [Signac integrate_atac vignette](https://stuartlab.org/signac/articles/integrate_atac) exists to solve: RNA and ATAC come from different cells (e.g. a scRNA-seq experiment and a separate scATAC-seq experiment on the same tissue type) and must be brought into a shared low-dimensional embedding before labels, expression, or gene activity can be transferred between them.

The workflow, at a high level:

1. Quantify the multiome/reference peaks in the scATAC-seq object (`FeatureMatrix` using the reference peak set) so both datasets share the same features.
2. On the scATAC-seq side, compute LSI (`FindTopFeatures` → `RunTFIDF` → `RunSVD`). The same is done for the multiome's ATAC assay.
3. Merge the two objects and recompute LSI on the merged data to get a baseline (uncorrected) embedding.
4. `FindIntegrationAnchors(object.list = list(pbmc.multi, pbmc.atac), reduction = "rlsi", dims = 2:30)` — reciprocal LSI projection finds anchors in each dataset's LSI space.
5. `IntegrateEmbeddings(anchorset = anchors, reductions = merged[["lsi"]], new.reduction.name = "integrated_lsi", dims.to.integrate = 1:30)` — integrates the LSI coordinates rather than the normalized count matrix, which is much better suited to the sparse scATAC matrix.
6. Optional follow-up: use `TransferData` / `GeneActivity` to push RNA labels and expression onto the query ATAC cells.

After integration, the user has a single object with a shared `integrated_lsi` reduction that can be fed into the single-cell Taiji workflow. `build-taiji-input` is bulk-only and will refuse the input regardless.

### Limitations of the filename-only heuristic

- Does not read `.h5ad`/`.rds` contents — a `.h5ad` whose filename gives no hint (e.g. `pbmc.h5ad`) will be `sc-undetermined` even if the file internally contains both RNA and ATAC assays.
- Cannot detect the case where a single `.h5ad` holds ONLY RNA but there's no sibling ATAC file — the skill would mark this `sc-undetermined` where a richer reader could say "pure scRNA-seq, route directly to SCENIC+".
- Tier 4 can be tripped by coincidental naming (e.g. a metadata file called `rna_and_atac_notes.md` next to a single `.h5ad`). Mitigated in practice because the file-iteration pool is restricted to recognized SC files plus matched 10x filenames, but a user with very loose naming could get a false `separate-assay` flag.

When in doubt, the right answer is: ask the user which layout they have. Don't guess silently.

## Classification logic

Let `B = {bulk extensions found in the directory}` and `S = {SC extensions found in the directory}`.

| `B`      | `S`      | Classification  | CLI exit      |
|----------|----------|-----------------|---------------|
| non-empty | empty    | `bulk`          | 0             |
| empty    | non-empty | `single-cell`   | 0 (with warning) |
| non-empty | non-empty | `mixed`         | 2 strict / 0 with `--allow-mixed` |
| empty    | empty    | `unknown`       | 0 (with top-extensions hint) |

Path errors (missing or unreadable root) short-circuit to exit `3` but still produce a best-effort result for whatever paths were readable.

## Why these exact extensions

- **`.tsv` (bulk RNA-seq).** The user's pipeline emits GeneQuant as TSVs (gene ID in column 1, counts/TPM in column 2+). Bulk Taiji reads them as `tags=GeneQuant`. A pure `.tsv`-only directory is flagged with a warning because nothing in the extension itself says "gene quant" — the same extension is used by every tabular tool on earth.
- **`.narrowPeak` (bulk ATAC-seq).** ENCODE-standard peak calls from MACS2 / MACS3. Taiji reads them with `format=NarrowPeak`. There is no SC equivalent — scATAC peak calls come out of ArchR / Signac and live inside `.rds` objects, not loose narrowPeak files.
- **`.bedpe` (bulk HiC).** Chromatin loops from HiC-Pro / Juicer / Cooler pipelines. Bulk Taiji reads them with `tags=ChromosomeLoop`. Single-cell HiC workflows are rare in this lab and are not in scope.
- **`.h5ad` (SC AnnData).** HDF5-backed, the scanpy/anndata ecosystem's default. Never used for bulk expression (bulk doesn't need the .obs/.var structure).
- **`.rds` (SC R object).** Technically `.rds` can contain any R object. In this lab's pipelines, it is effectively always a Seurat or SingleCellExperiment object. False-positive risk is low because bulk pipelines don't serialize R objects — they emit tabular files.

## Edge cases and caveats

### `.tsv` ambiguity

A directory of just `.tsv` files is flagged because:
- Bulk metadata sheets are often TSVs.
- Some scRNA tools export count matrices as TSVs (transposed, cells × genes).
- Old-style STAR / RSEM outputs use `.tsv` but so do tables of DEG results, sample metadata, QC reports, etc.

When this warning fires, confirm with the user that the TSVs are per-sample GeneQuant files (gene IDs in column 1, one row per gene, consistent row count across samples) before letting `build-taiji-input` run.

### `.rds` ambiguity

`.rds` is always classified as single-cell because in this lab's pipelines that's what it carries. It would be wrong in the general R ecosystem (a bulk DESeq2 result can be serialized as `.rds`), but the bias is protective: if the user has a bulk-only workflow they shouldn't be serializing to `.rds` anyway, so flagging it as SC nudges them to export a TSV instead.

### Mixed directories

The strict default exists because the common mixed-directory cases are user mistakes (copied the wrong folder, left an h5ad from a previous project, etc.). If the mixed state is intentional (a comparison project that stores both), pass `--allow-mixed` to downgrade to a warning.

### Path is a single file

Passing a single file works and classifies it as a one-element directory. Useful for quick one-liners like `python detect.py suspect.rds`.

### Empty or nearly-empty directories

`empty` and `unknown` are both classified as `unknown`. The CLI surfaces the top-5 non-matching extensions as a hint — this is usually enough to notice that the user is looking at a wrapper directory and the real data is one level down.

## Extending to new formats

To add a new recognized format, edit `BULK_EXTS` or `SC_EXTS` in `scripts/detect.py`. Keep two properties:

1. The extension must be **unambiguous in this lab's pipelines**. If a given extension is used for both bulk and SC workflows by the same lab, do not add it — force the user to disambiguate via a per-file check.
2. After adding an extension, add a fixture directory under `skills/detect-dataset-type-workspace/fixtures/` and an assertion in `evals/evals.json` so the eval loop covers it.

## What this skill does *not* do

- Does not open or parse file contents. No HDF5 reads, no tab-delimited sniffing, no R runtime.
- Does not move or rename files.
- Does not distinguish scRNA from scATAC within the SC class — both are `single-cell`.
- Does not validate that a bulk dataset is *complete* (i.e. that every sample has both RNA and ATAC). That's `build-taiji-input`'s job.
- Does not detect mis-labeled files (e.g. a `.tsv` containing single-cell data). The user's pipeline conventions are trusted at the extension level.
