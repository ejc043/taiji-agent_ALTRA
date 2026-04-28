# Peak calling (per-cluster Signac::CallPeaks)

## Why per-cluster, not global + subset

The two viable ATAC pseudobulk strategies are:

1. **Global peak call + per-cluster counting**: run MACS once on all fragments, then count reads per cluster in the resulting peak set.
2. **Per-cluster peak call** (what this skill does): subset fragments to the cluster's cells, call MACS on each subset, emit one narrowPeak per cluster.

(2) catches regulatory elements that are active in a minority cluster but get washed out in the global signal. For bulk Taiji the downstream tool wants per-sample peak sets anyway — it uses the narrowPeak files directly as the accessibility layer — so emitting them per cluster aligns naturally with the expected input contract. The tradeoff is runtime: a dataset with 40 clusters runs MACS 40 times, which is fine on SLURM but noticeable locally.

If you do want strategy (1) for comparison, skip this skill and call MACS manually on the full fragments; the output will feed the same `build-taiji-input` manifest but with one ATAC file shared across all pseudobulks.

## Implementation: `call_peaks.R`

Peak calling is driven by `skills/pseudobulk-construct/scripts/call_peaks.R`, which uses **`Signac::CallPeaks`** (Signac shells out to `macs2`/`macs3` under the hood). The R route — instead of a direct `gunzip | awk | macs2` pipe — is preferred because it gives us:

- **`CreateFragmentObject`-based plumbing**: per-group fragment scoping that stays consistent with the rest of the Signac/Seurat ecosystem.
- **Per-group barcode reconciliation**: handles the common case where the Seurat object's `colnames()` carry suffixes (`AAACGCT-1`, `ATAC_AAACGCT`) that don't appear in the raw 10x `fragments.tsv.gz`.
- **An `.rds` peak object alongside the narrowPeak**: useful if you want to reload peaks into Signac downstream.

### Per-group flow

For each entry in `groups_plan.json`:

1. Resolve cells from `clusters.csv` (cluster id [+ optional metadata col/value]).
2. Restrict to the ATAC subset:
   - If `meta.data` has an `assay` column (the integrate_atac convention), keep cells where `assay == 'ATAC'`.
   - Otherwise use `Cells(obj[[atac_assay]])`.
3. Reconcile barcodes against a sample of `fragments.tsv.gz` (default: first 5M rows, deduped). Try candidate strip patterns and pick the one that maximizes overlap. Patterns tried, in order:
   - identity (no correction)
   - split on `_` keep first piece (`ATAC_AAACGCT` → `AAACGCT`)
   - split on `-` keep first piece (`AAACGCT-1` → `AAACGCT`)
   - strip trailing `-<digits>` (`AAACGCT-1`, `AAACGCT-12`)
   - strip trailing `_<digits>` (`AAACGCT_1`, `AAACGCT_12`)
4. **Fail loud if best overlap < `--min-overlap` (default 95%)** — print sample barcodes from both sides so the user can diagnose. The most common cause is the wrong fragments file for the sample, not a fixable suffix.
5. `subset(obj, cells = ...)` → `RenameCells(new.names = corrected)` → `Fragments(obj_sub) <- CreateFragmentObject(path, cells = corrected)`. Renaming on the **subset** (not the parent) avoids duplicate-name collisions when both RNA and ATAC sides of a co-embedded object share bare barcodes.
6. `Signac::CallPeaks(object, macs2.path, effective.genome.size, combine.peaks = TRUE)`.
7. Write `<group>_peaks.rds` (saveRDS) and `<group>.narrowPeak` (ENCODE 10-col).

Idempotency: a group whose `<group>_peaks.rds` already exists is skipped, so partial reruns don't redo expensive MACS invocations.

## Fragments-file requirements

The skill needs a 10x / Signac-style `fragments.tsv.gz` with these columns (tab-separated):

```
chrom  start  end  barcode  duplicate_count
```

It only reads columns 1–4 for the barcode-overlap sample; for peak calling itself, Signac wants a properly bgzipped + tabix-indexed file (`.tbi` sibling required). If `tabix index` hasn't been run on your fragments file, do that once: `tabix -p bed fragments.tsv.gz`.

## What the skill does when fragments are missing

If `--fragments` is not passed:

- RNA pseudobulks are still produced.
- The manifest's `atac_seq` column is empty for all rows.
- `build-taiji-input` will emit a warning about missing ATAC but produce the xlsx for RNA-only columns.

If `--fragments` is passed but the file is unreadable (bad path, truncated gzip, missing `.tbi`):

- `call_peaks.R` fails loudly at startup (sample step or Signac fragment-object creation).
- Groups whose narrowPeak isn't written end up with empty `atac_seq` in the manifest.

## Effective genome sizes

`Signac::CallPeaks` takes a numeric `effective.genome.size`. The skill maps:

| Genome | Size       |
|--------|------------|
| hg19   | 2.7e9      |
| hg38   | 2.9e9      |
| mm9    | 1.87e9     |
| mm10   | 1.87e9     |
| mm39   | 1.87e9     |

Add new genomes to `effective_genome_size()` in `call_peaks.R`.

## Quality checks to run afterward

The skill does not automatically QC the peak calls. You should spot-check:

- **Peak count per cluster**: typical scATAC pseudobulks produce 20k–100k peaks. Clusters emitting <1k peaks probably had too few cells or too-shallow fragments; clusters emitting >200k probably picked up noise.
- **FRiP (fraction of reads in peaks)**: > 0.3 is healthy for scATAC pseudobulks. Lower suggests the peak call is dominated by background.
- **Peak overlap with a reference blacklist**: ENCODE blacklist regions shouldn't contribute more than a few percent of peaks.

None of these are gates — they're audit checks that catch the "something went wrong upstream" case.

## MACS3

`--peak-caller macs3` makes the orchestrator pick the `macs3` binary; pseudobulk.py resolves the full path via `shutil.which()` and passes it to `Signac::CallPeaks` as `macs2.path` (the kwarg name is misleading — Signac happily uses macs3 if you point it there). MACS3 is a drop-in for the macs2 invocation Signac builds; the scATAC defaults are unchanged. If you use MACS3, just verify that `macs3 --version` returns >= 3.0 on the SLURM node; older builds labeled "macs3" were pre-release and had compatibility bugs.

## Limitations

- The 5M-row fragment-barcode sample is enough to detect suffix patterns reliably for typical 10x libraries (millions of fragments, tens of thousands of cells). For unusually small fragments files where 5M rows would exceed the file, the sample shrinks to the file size — still fine. For unusually shallow libraries where 5M rows comes from <100 cells, you could get a false-negative reconciliation; raise `--frag-sample` if the heuristic looks suspicious.
- Per-cluster peak calling on very small clusters (close to the 20-cell floor) produces noisy narrowPeak files even when MACS returns success. Consider raising `--min-cluster-cells` for better downstream peak quality.
- The effective-genome-size table covers the common builds. For exotic genomes, extend `effective_genome_size()` in `call_peaks.R`.
- The Signac route requires `Signac::CallPeaks` (Signac ≥ 1.12). Older Signac versions had a slightly different `CallPeaks` API; pin to ≥ 1.12 to avoid surprise.
