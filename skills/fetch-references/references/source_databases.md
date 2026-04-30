# Source databases for reference data

## Primary sources (what the manifest uses)

### GENCODE FASTA + GTF — via EMBL-EBI mirror

URL pattern:

```
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_<species>/release_<N>/<filename>
```

Where `<species>` is `human` or `mouse` and `<N>` is the release tag (e.g. `45`, `M34`, `M25`, `19`).

- **Pros:** EBI mirror is reliably faster than the GENCODE primary site (UCSC), serves over HTTPS, supports range requests, and is content-addressable per release. Releases are immutable once published.
- **Cons:** EBI occasionally schedules brief downtime for maintenance (announced in advance). Fallback: same path under `https://ftp.ebi.ac.uk/...` or the GENCODE primary at `https://www.gencodegenes.org/`.
- **Version pins in the manifest:**
  - `hg38` → GENCODE v45 (released Mar 2024). Stable release; we're not chasing the rolling latest.
  - `hg19` → GENCODE v19 (released Feb 2014, frozen). v19 is the canonical "GRCh37 era" release; GENCODE doesn't do new releases on GRCh37.
  - `mm10` → GENCODE M25 (released Apr 2020, the last release on GRCm38). Newer releases (M26+) are on GRCm39 only.

### CIS-BP MEME — vendored at `cisbp_database/`

Files (all checked into the repo, no network needed):

- **Human:** `cisbp_database/cisbp_human_2.meme` — full CIS-BP, 4443 PWMs, gene-symbol primary IDs (e.g. `MOTIF POU3F2`). This is the **canonical** human motif file for the project — verified against a reference Nextflow run on the RA/OA hg38 chr22 dataset (Spearman 0.97, RMSE ~4e-5; identical TF set, top-50 100% match).
- **Mouse:** `cisbp_database/cisBP_mouse.meme`.
- **Legacy human:** `cisbp_database/cisBP_human.meme` (1882 PWMs) is also vendored for reproducing older runs but is not the default.

**Why CIS-BP and not HOCOMOCO or JASPAR?**

- **CIS-BP** is the project's only supported source. Primary IDs are TF gene symbols, so Taiji's `GeneRanks.tsv` is directly interpretable without any post-hoc motif-ID-to-gene mapping. CIS-BP also has the broadest TF coverage (~4400 vs HOCOMOCO's ~770 for human), which matters when downstream analyses ask whether a specific TF is in the network.
- **HOCOMOCO** was the prior default but was removed: switching motif sources changes both the TF set and PageRank values, so locking to one source keeps cross-run results comparable. Reproductions of the reference Nextflow run only match when CIS-BP is used.
- **JASPAR** is also human-curated and a reasonable alternative; if a future project needs JASPAR, copy `reference_manifest.yml` and add a `motif_jaspar:` block plus restore the `--motif-source` choices in `fetch.py`.

## Alternative sources (if primaries are unreachable)

| Asset            | Primary (in manifest)          | Mirror / alternative                                       |
|------------------|--------------------------------|------------------------------------------------------------|
| GENCODE FASTA    | `ftp.ebi.ac.uk/pub/databases/gencode/...` | `ftp.ebi.ac.uk/...` (same host, alt path) or `gencodegenes.org` direct |
| GENCODE GTF      | same                           | same                                                       |
| Genome FASTA     | (GENCODE primary assembly)     | UCSC: `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz` (chromosome ordering differs slightly) |
| Annotation GTF   | (GENCODE)                      | Ensembl: `https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.chr.gtf.gz` (ID conventions differ) |
| MEME motifs      | CIS-BP vendored at `cisbp_database/` | JASPAR 2024 (would need a manifest fork; not currently wired up) |

If you hot-swap a primary URL, also remember to update `approx_size_mb` in the manifest — UCSC's `hg38.fa.gz` is ~900 MB compressed but unzips to ~3.1 GB, which is close enough to GENCODE's ~3.0 GB to satisfy the ±5% tolerance, but Ensembl GTF is noticeably smaller than GENCODE's.

## Why specific GENCODE versions vs "always latest"

Pinning to specific releases instead of `release_latest`:

1. **Reproducibility.** A pipeline run today and a re-run in six months should consume the same reference. Latest pointers drift.
2. **Memory's worth of context.** Manuscript reviewers ask "which annotation version did you use?" — having `gencode.v45.annotation.gtf` baked into the path answers that for free.
3. **Bug surfaces.** Each new GENCODE release occasionally breaks downstream tools that hardcoded chromosome-name expectations or relied on a now-removed feature type. Pinning lets you bump deliberately.

To upgrade: edit the version number in `reference_manifest.yml` and re-run with `--force`. The skill will re-download into the same paths (overwriting the old files); the manifest itself is the change-tracked record of what was current.

## Why these particular sizes

`approx_size_mb` numbers in the manifest are the **uncompressed** sizes the files have on disk after gunzip. The size check is on what's at `<output_dir>/<genome>/<file>`, not on the .gz the URL serves. The numbers were determined by running the downloader once and recording the file sizes; they're rounded to the nearest 100 MB to leave room for minor re-cuts.

## What's NOT in the manifest

- **BWA / STAR / RSEM indexes.** Out of scope for this skill — the agent's design assumption is that all inputs to Taiji are pre-aligned/pre-quantified, so genome indexes are not needed.
- **Blacklist BEDs.** The ENCODE blacklist (hg38: `https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz`) is sometimes useful for filtering ATAC peaks but isn't required by Taiji directly. Add it to the manifest if your workflow needs it.
- **Chromosome size files.** Auto-generated from `genome.fa.fai` if you need them: `cut -f1,2 genome.fa.fai > chrom.sizes`.
- **Epitensor HiC loop files.** These are the per-genome fallback HiC loops that `build-taiji-input` uses when a sample doesn't supply its own bedpe. The Taiji-pipeline org distributes them at `https://github.com/Taiji-pipeline/Epitensor`; not auto-fetched here because most users either (a) have per-sample HiC bedpe, or (b) skip HiC entirely.

Add any of the above to the manifest if your workflow needs them.
