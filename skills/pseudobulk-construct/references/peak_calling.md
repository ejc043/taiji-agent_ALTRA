# Peak calling (per-cluster MACS2)

## Why per-cluster, not global + subset

The two viable ATAC pseudobulk strategies are:

1. **Global peak call + per-cluster counting**: run MACS2 once on all fragments, then count reads per cluster in the resulting peak set.
2. **Per-cluster peak call** (what this skill does): subset fragments to the cluster's cells, call MACS2 on each subset, emit one narrowPeak per cluster.

(2) catches regulatory elements that are active in a minority cluster but get washed out in the global signal. For bulk Taiji the downstream tool wants per-sample peak sets anyway — it uses the narrowPeak files directly as the accessibility layer — so emitting them per cluster aligns naturally with the expected input contract. The tradeoff is runtime: a dataset with 40 clusters runs MACS2 40 times, which is fine on SLURM but noticeable locally.

If you do want strategy (1) for comparison, skip this skill and call MACS2 manually on the full fragments; the output will feed the same `build-taiji-input` manifest but with one ATAC file shared across all pseudobulks.

## MACS2 invocation

For each group, the skill runs:

```
gunzip -c fragments.tsv.gz \
  | awk -v BARCFILE=<cluster>.txt -F'\t' \
        'BEGIN {while ((getline b < BARCFILE) > 0) bc[b]=1; close(BARCFILE)}
         {if ($4 in bc) print $1"\t"$2"\t"$3}' \
  | macs2 callpeak \
      -t - -f BED \
      --nomodel --shift -100 --extsize 200 \
      -q 0.05 -g <hs|mm> \
      --call-summits \
      -n <group_name> --outdir <output_dir>/atac
```

Rationale for each flag:

- `-f BED`: fragments files are 3-column BED-like after the awk filter.
- `--nomodel --shift -100 --extsize 200`: the ATAC-seq standard. ATAC fragments are already Tn5 cut-site pairs; MACS2's default shifting model is built for ChIP-seq and doesn't apply. Shift −100 / extend 200 centers a 200 bp window on each cut site, which is what scATAC tutorials (10x, Signac, ArchR) all use.
- `-q 0.05`: 5% FDR. Matches ENCODE's per-replicate narrowPeak threshold. Can be tightened to 0.01 if downstream analysis is noise-sensitive.
- `-g hs|mm`: effective genome size. Mapped from `--genome` via a small lookup table in `pseudobulk.py` (hg19/hg38 → hs, mm9/mm10/mm39 → mm). If you need another organism, extend `macs2_genome_size()` — MACS2 takes any integer bp count.
- `--call-summits`: emit summit positions. Downstream Taiji doesn't strictly require them but they're cheap and occasionally useful for motif work.

## Fragments-file requirements

The skill needs a 10x / Signac-style `fragments.tsv.gz` with these columns (tab-separated):

```
chrom  start  end  barcode  duplicate_count
```

It only reads columns 1-4, so the duplicate-count column can be missing. The file must be gzipped (`gunzip -c` is the first step of the pipeline) and should ideally have a companion `.tbi` index if you plan to do anything with `tabix` later — the skill itself does a full scan, so the index isn't required for this step.

Barcodes in the fragments file must match barcodes in the Seurat object's `colnames()`. If your Seurat object was created from the same cellranger-arc output the fragments came from, they will match by construction. If the object was re-barcoded (prefixed with a sample name, for example), pass the same prefix handling through to the barcodes file or the match will silently fail — MACS2 will run on zero fragments and emit an empty narrowPeak.

## What the skill does when fragments are missing

If `--fragments` is not passed:

- RNA pseudobulks are still produced.
- The manifest's `atac_seq` column is empty for all rows.
- `build-taiji-input` will emit a warning about missing ATAC but produce the xlsx for RNA-only columns.

If `--fragments` is passed but the file is unreadable (bad path, truncated gzip, permission denied):

- The skill prints a warning per group, skips that group's peak call, and continues.
- Groups with missing narrowPeak files end up with empty `atac_seq` in the manifest.
- This is intentionally non-fatal so a batch run with one bad cluster doesn't lose the good outputs. For strict-mode behavior, check narrowPeak sizes before running `build-taiji-input`.

## Quality checks to run afterward

The skill does not automatically QC the peak calls. You should spot-check:

- **Peak count per cluster**: typical scATAC pseudobulks produce 20k-100k peaks. Clusters emitting <1k peaks probably had too few cells or too-shallow fragments; clusters emitting >200k probably picked up noise.
- **FRiP (fraction of reads in peaks)**: > 0.3 is healthy for scATAC pseudobulks. Lower suggests the peak call is dominated by background.
- **Peak overlap with a reference blacklist**: ENCODE blacklist regions shouldn't contribute more than a few percent of peaks.

None of these are gates — they're audit checks that catch the "something went wrong upstream" case.

## MACS3

The `--peak-caller macs3` flag swaps the binary name. MACS3 is a drop-in for the flags above; the scATAC recipe is unchanged. Output format is also narrowPeak. If you use MACS3, just verify that `macs3 --version` returns >= 3.0 on the SLURM node; older builds labeled "macs3" were pre-release and had compatibility bugs.

## Limitations

- The awk filter loads the entire barcodes list into memory, then streams fragments. For a per-cluster list of ~2000 barcodes and a 1M-line fragments file, this is fast. For a run with hundreds of clusters, the fragments file is scanned N times — there's no cached-seek optimization. If this becomes a bottleneck, the pragmatic fix is to sort fragments by barcode once and then iterate in sorted order; out of scope for v1 of this skill.
- Per-cluster peak calling on very small clusters (close to the 20-cell floor) produces noisy narrowPeak files even when MACS2 returns success. Consider raising `--min-cluster-cells` for better downstream peak quality.
- MACS2's effective-genome-size lookup is species-level only (`hs`/`mm`). For more exotic genomes, pass `--genome` as a bp integer and extend the mapping.
