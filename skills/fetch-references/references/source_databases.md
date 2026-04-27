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
  - `mm39` → GENCODE M34 (released Mar 2024).

### HOCOMOCO v11 MEME — autosome.org

URL pattern:

```
https://hocomoco11.autosome.org/final_bundle/hocomoco11/full/<SPECIES>/mono/HOCOMOCOv11_full_<SPECIES>_mono_meme_format.meme
```

Where `<SPECIES>` is `HUMAN` or `MOUSE`.

- **Pros:** Stable URL since v11 release in 2018; small files (~5 MB); MEME format is the canonical input for TF-binding-motif scanning tools (FIMO, MEME-suite, Taiji's regulatory step).
- **Cons:** The autosome.org host occasionally hiccups; if that happens, the same files mirror at `https://opera.autosome.org/...` for a fallback path.
- **Why HOCOMOCO and not CIS-BP or JASPAR?** All three are valid choices:
  - **HOCOMOCO** (this skill's default): human-curated; high-confidence motifs only; smaller; widely used in regulatory-network work.
  - **JASPAR**: also human-curated; broader species coverage; MEME format available at `https://jaspar.genereg.net/download/data/<release>/JASPAR<release>_CORE_<species>_non-redundant_pfms_meme.txt`.
  - **CIS-BP**: largest motif catalog (includes inferred motifs); requires picking species + format from a UI rather than scripted download — awkward for an idempotent installer.

If you want JASPAR instead, copy `reference_manifest.yml` and replace the `motif:` URL. Reasonable JASPAR replacement for hg38:

```yaml
motif:
  url:    https://jaspar.genereg.net/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt
  target: hg38/JASPAR2024_vertebrates.meme
  gunzip: false
  approx_size_mb: 3
```

## Alternative sources (if primaries are unreachable)

| Asset            | Primary (in manifest)          | Mirror / alternative                                       |
|------------------|--------------------------------|------------------------------------------------------------|
| GENCODE FASTA    | `ftp.ebi.ac.uk/pub/databases/gencode/...` | `ftp.ebi.ac.uk/...` (same host, alt path) or `gencodegenes.org` direct |
| GENCODE GTF      | same                           | same                                                       |
| Genome FASTA     | (GENCODE primary assembly)     | UCSC: `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz` (chromosome ordering differs slightly) |
| Annotation GTF   | (GENCODE)                      | Ensembl: `https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.chr.gtf.gz` (ID conventions differ) |
| MEME motifs      | HOCOMOCO v11 autosome.org      | JASPAR 2024 (URL above); CIS-BP (manual download)         |

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
