---
name: fetch-references
description: Idempotently download the per-genome reference data Taiji needs — primary-assembly FASTA + GENCODE GTF + HOCOMOCO MEME motif file — from stable upstream hosts (EMBL-EBI for GENCODE, autosome.org for HOCOMOCO) into <work_dir>/dependencies_data/. Read URL specs from a single reference_manifest.yml, skip files that already exist with the right size, gunzip on the fly, optionally `samtools faidx` the FASTA, and optionally rewrite skills/build-taiji-input/assets/genomes.yml to point at the downloaded files. Trigger on phrases like "download genome", "fetch reference data", "get the FASTA", "stage the GTF", "where do I get HOCOMOCO", "download MEME motif file", "set up reference data for Taiji", or anywhere the user is missing genome.fa / genes.gtf / *.meme before running Taiji.
---

# Fetch references (per-genome reference-data downloader)

## What this skill produces

For each `--genome <name>`, downloads three files into `<output_dir>/<genome>/`:

| Kind     | Default name                       | Source                      | Approx size (uncompressed) |
|----------|------------------------------------|-----------------------------|----------------------------|
| `fasta`  | `genome.fa` (+ `.fai` if samtools) | GENCODE (EMBL-EBI mirror)   | ~3 GB (human), ~2.7 GB (mouse) |
| `gtf`    | `genes.gtf`                        | GENCODE (EMBL-EBI mirror)   | ~0.9-1.5 GB                |
| `motif`  | `HOCOMOCOv11_<species>.meme`       | autosome.org (HOCOMOCO v11) | ~5 MB                      |

Default output directory is `dependencies_data/` (under the current working directory, override with `--output`).

Available genomes are listed by `--list`: `hg38`, `hg19`, `mm10`, `mm39`. Add new ones by editing `scripts/reference_manifest.yml` — the schema is documented at the top of that file.

## Idempotency contract

- A target file present at the expected path **with size within ±5% of the manifest's `approx_size_mb`** is treated as already-downloaded and **skipped**. The skill never re-downloads files unnecessarily.
- `--check` reports what's present / missing without downloading anything.
- `--force` re-downloads everything for the chosen genome.
- `--dry-run` prints which URLs would be hit without touching the network.
- The skill never deletes files it didn't write itself. A `.part` temp file is used during download and atomically renamed on success; a partial file from a crashed run is cleaned up on the next attempt.

## Library API

```python
from fetch import fetch_genome, check_genome

# Download (or skip if present)
result = fetch_genome("hg38", output_dir="./dependencies_data")
# result.files: {"fasta": FileResult(target=..., status="downloaded"), ...}
# result.ok:    True if all files are present / downloaded

# Just check what's there
status = check_genome("hg38", output_dir="./dependencies_data")
```

## CLI

```bash
# Show what's available
python scripts/fetch.py --list

# Check what's already downloaded
python scripts/fetch.py --genome hg38 --output dependencies_data/ --check

# Download anything missing
python scripts/fetch.py --genome hg38 --output dependencies_data/

# Re-download even if present
python scripts/fetch.py --genome hg38 --output dependencies_data/ --force

# Dry-run: print what would be downloaded
python scripts/fetch.py --genome hg38 --output dependencies_data/ --dry-run

# Download + edit build-taiji-input/assets/genomes.yml in place
python scripts/fetch.py --genome hg38 --output dependencies_data/ --update-genomes-yml

# Machine-readable JSON report
python scripts/fetch.py --genome hg38 --output dependencies_data/ --json
```

## Source choices (why these URLs)

- **GENCODE FASTA + GTF via EMBL-EBI mirror** (`ftp.ebi.ac.uk/pub/databases/gencode/...`): the consortium URL pattern `release_N/...` has been stable for a decade; releases are immutable once published. EBI mirror is faster than the SciSchool primary.
- **HOCOMOCO v11 MEME from autosome.org**: hosted at the same `final_bundle/hocomoco11/full/<SPECIES>/mono/` path since 2018; versioned, no rolling drift. v11 is the most recent release with a stable URL.
- **Why not Ensembl?** Their URL paths bake the release number in differently and mouse genome filenames use a different convention than GENCODE. Picking one consortium and sticking with it makes URLs predictable across genomes.
- **Why not UCSC?** UCSC's `hg38.fa` is fine for many uses but the chromosome ordering and contig naming drift slightly across genome builds. GENCODE's "primary assembly" FASTA is the cleaner choice for downstream Taiji.

See `references/source_databases.md` for full justification + alternatives if a primary source is unreachable.

## Pre-flight: check before download

The skill is **always idempotent** — but in interactive use you usually want to know what you'd be downloading before kicking off ~5 GB of traffic. Recommended pattern:

```bash
# 1. See what's present and what's missing
python scripts/fetch.py --genome hg38 --output dependencies_data/ --check

# 2. If anything is missing, see what would be downloaded
python scripts/fetch.py --genome hg38 --output dependencies_data/ --dry-run

# 3. Download for real
python scripts/fetch.py --genome hg38 --output dependencies_data/
```

The skill prints a per-file report after every run with status (`present` / `downloaded` / `skipped` / `failed`), bytes on disk, and the target path.

## Integration with the rest of the agent

- **Auto-attaches to workflow-log.** When a workflow is active (the workflow-log skill's `active_run` pointer file exists), `fetch.py` writes a `fetch_references` stage entry with the per-file status, sizes, and target paths. When no run is active it logs nothing.
- **Optional `--update-genomes-yml`.** After a successful download, edits `skills/build-taiji-input/assets/genomes.yml` in place to point the chosen genome's `fasta:` and `gtf:` keys at the resolved paths. Existing `hic:` paths and other genome entries are preserved verbatim. Use this on first install so `build-taiji-input` doesn't trip the placeholder-paths warning.
- **The MEME motif file is NOT in `genomes.yml`.** Taiji's regulatory-network step reads the motif file from its YAML config, not from the xlsx. After this skill downloads the MEME file, point your Taiji config's `motif_file:` (or equivalent key — varies by Taiji version) at the resolved path.

## Interaction pattern

1. **Always run `--check` first** when a user is unsure what's already on disk. It's free (no network), tells you exactly what's missing, and surfaces files that exist but have wrong sizes (likely interrupted previous downloads).
2. **Default to GENCODE primary assembly + HOCOMOCO v11** unless the user has a specific reason for an alternative source. The manifest is editable; pinning to non-default sources is a deliberate choice that should live in a per-site fork of `reference_manifest.yml`.
3. **Watch out for the FASTA size.** ~3 GB uncompressed for human takes 3-10 minutes on a typical lab network. The downloader streams + decompresses on the fly so peak disk use is ~3 GB, not 6 GB.
4. **`samtools faidx` is best-effort.** If samtools is on PATH, the `.fai` index is built automatically after download. If not, the skill prints a NOTE and the user can run `samtools faidx genome.fa` manually later — the FASTA file is still usable without the index for some Taiji subcommands.

## When to reach for the reference doc

`references/source_databases.md` — full list of upstream hosts, alternatives if a primary is unreachable, and the rationale for the version pins (e.g. why GENCODE v45 specifically vs latest).

## Why the defaults are what they are

- **Manifest-driven, not hardcoded URLs.** Adding a new genome (CHM13, T2T-CHM13v2, drosophila, etc.) is a single YAML edit, not a code change.
- **`approx_size_mb` not exact byte counts.** Upstream files occasionally get re-released with cosmetic changes (whitespace, headers, etc.) that nudge the byte count by <1%. ±5% tolerance avoids false-positive "wrong size" failures while still catching obviously-truncated downloads.
- **No SHA256 by default.** Adding hashes would require maintenance every time GENCODE re-cuts a release. The size check catches the common failure modes (truncation, partial download) without that maintenance burden. If reproducibility matters more than convenience for a given site, add `sha256: <hex>` to the manifest entry — the downloader will verify when present.
- **gunzip on the fly, not "download .gz then unzip".** Saves ~3 GB of intermediate disk for the human FASTA.
- **`.part` then atomic rename.** A killed download leaves a `.part` file that the next run cleans up; the final target is only ever written atomically.

## Limitations

- HTTP only. No support for FTP, S3, GCS, or magnet links. The current sources all serve over HTTPS so this is fine; rip-out of `urllib.request` would be needed for protocol expansion.
- Single-stream downloads (no parallel chunks). For a 3 GB FASTA at ~50 MB/s this is fine; for slower links a multi-stream tool like `aria2c` would be faster.
- No checksum verification by default (see "approx_size_mb" rationale above). Add `sha256:` to the manifest entry to enable per-file verification.
- macOS-only edge case: on some macOS setups `urllib.request` doesn't trust the system root CA store and HTTPS hits silently fail. Workaround: `pip install certifi && export SSL_CERT_FILE=$(python -m certifi)`.
