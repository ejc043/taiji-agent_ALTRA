---
name: detect-dataset-type
description: Classify a dataset directory as bulk, single-cell, mixed, or unknown based purely on file extensions and filename patterns. For single-cell, additionally determine whether RNA and ATAC come from the SAME cells (multiome) or DIFFERENT cells (separate-assay, requiring Signac/Seurat co-embedding). Use this skill whenever the user wants to determine, detect, figure out, check, classify, identify, or tell the data type of a directory, folder, or project — especially before deciding which Taiji workflow to run. Bulk signatures are .tsv (RNA-seq), .narrowPeak (ATAC-seq), and .bedpe (HiC); single-cell signatures are .h5ad (AnnData), .rds (Seurat / SingleCellExperiment), .h5mu (MuData), plus 10x cellranger/cellranger-arc filename patterns (fragments.tsv, barcodes.tsv, matrix.mtx, filtered_feature_bc_matrix.h5, etc.). Trigger on phrases like "is this bulk or single-cell", "what kind of data is this", "detect the assay", "classify this dataset", "sc vs bulk", "h5ad or tsv", "multiome or separate", "do I need to co-embed", or anywhere the user needs a pre-flight check before a pipeline.
---

# Detect dataset type (bulk vs single-cell, with SC modality sub-classification)

## What this skill produces

A classification of the input path(s) into one of four labels, plus the evidence:

| Classification | Meaning                                                               |
|----------------|-----------------------------------------------------------------------|
| `bulk`         | Only bulk signatures present (.tsv, .narrowPeak, .bedpe).            |
| `single-cell`  | Only single-cell signatures present (see list below).                |
| `mixed`        | Both bulk and single-cell signatures present. Strict-mode error.     |
| `unknown`      | No recognized signatures. Returns the top extensions as a hint.      |

When classification is `single-cell`, the result also carries a `sc_modality`:

| sc_modality       | Meaning                                                                  |
|-------------------|--------------------------------------------------------------------------|
| `multiome`        | RNA + ATAC measured from the SAME cells (paired). No co-embedding needed. |
| `separate-assay`  | RNA and ATAC from DIFFERENT cells. Must be co-embedded via Signac/Seurat. |
| `sc-undetermined` | Single-cell data but modality layout couldn't be inferred from filenames. |

Downstream consumers (notably `build-taiji-input`) should proceed only on `bulk`. Any SC result — regardless of modality — is a redirect: bulk Taiji cannot consume single-cell inputs.

## Single-cell signatures recognized

**Extensions:** `.h5ad` (AnnData), `.rds` (Seurat / SingleCellExperiment), `.h5mu` (MuData — multi-modal by construction).

**Filename patterns** (matched case-insensitively as substrings; take priority over bulk extensions to prevent 10x cellranger-arc outputs like `atac_fragments.tsv.gz` from being misclassified as bulk):

- `fragments.tsv` — scATAC fragments (10x / Signac)
- `barcodes.tsv` / `features.tsv` / `genes.tsv` — 10x GEX/ATAC metadata
- `matrix.mtx` — sparse count matrix (10x / STARsolo)
- `filtered_feature_bc_matrix.h5` — 10x GEX or Multiome combined output
- `filtered_peak_bc_matrix.h5` — 10x ATAC peak-barcode matrix
- `raw_feature_bc_matrix.h5` — 10x raw counts

## Modality sub-classification: how it decides

Tier logic, first match wins:

1. **Any `.h5mu` present** → `multiome` (MuData is multi-modal by construction).
2. **cellranger-arc signature**: both `filtered_feature_bc_matrix.h5` and `atac_fragments.tsv` present → `multiome`.
3. **Explicit multiome token** in any filename (`multiome`, `multi_omic`, `multimodal`, `arc`) → `multiome`.
4. **Both RNA and ATAC hints** present across filenames (`rna`/`gex`/`expression`/`scrna` AND `atac`/`scatac`/`peak`/`chromatin_accessibility`) → `separate-assay`.
5. Otherwise → `sc-undetermined`.

This is a filename-only heuristic — it does not open `.h5ad`/`.rds` files to inspect assay slots. A generic `pbmc.h5ad` with no sibling ATAC file will end up `sc-undetermined`, which is the correct conservative behavior: the skill warns the user to clarify rather than guessing wrong.

## Inputs this skill accepts

One or more paths (files or directories). Directories are walked recursively by default, skipping `.git`, `__pycache__`, `.ipynb_checkpoints`, `node_modules`, `.venv`, hidden files, and `.DS_Store` / `Thumbs.db` noise.

Extension matching is **case-insensitive** and **.gz-aware** — so `.tsv`, `.TSV`, `.tsv.gz`, and `.NarrowPeak` all match. This is deliberate because ENCODE conventions are case-sensitive on paper (e.g. `.narrowPeak`) but in practice collaborators send a mix.

## How to run it

The skill bundles a single script. Library use (preferred when gating another pipeline):

```python
from pathlib import Path
import sys
sys.path.insert(0, str(Path("skills/detect-dataset-type/scripts")))
from detect import detect_dataset_type

result = detect_dataset_type(["/path/to/project"])
if not result.ok_for_bulk_taiji:
    raise SystemExit(f"build-taiji-input cannot proceed: {result.classification}")
```

CLI use:

```bash
python scripts/detect.py /path/to/project
python scripts/detect.py /path/to/project --format json
python scripts/detect.py /path/to/project --no-recursive
python scripts/detect.py /path/a /path/b          # multiple roots, combined report
python scripts/detect.py /path/to/project --allow-mixed  # mixed -> warning instead of error
```

Exit codes:
- `0` — classified (bulk, single-cell, or unknown). SC and unknown are non-fatal at the CLI level; a caller decides what to do with them.
- `2` — mixed dataset under strict mode (default). Intentionally loud.
- `3` — path error (root does not exist or is unreadable).

## Output shape (JSON)

```json
{
  "classification": "single-cell",
  "sc_modality":    "separate-assay",
  "paths_scanned": ["/abs/path"],
  "bulk_files":    {},
  "sc_files":      {".h5ad": ["/abs/path/pbmc_rna.h5ad", "/abs/path/pbmc_atac.h5ad"]},
  "other_exts":    {".md": 2, ".yml": 1},
  "sizes":         {"/abs/path/pbmc_rna.h5ad": 1234, ...},
  "warnings":      ["single-cell RNA-seq and scATAC-seq detected from DIFFERENT cell populations ... https://stuartlab.org/signac/articles/integrate_atac ..."],
  "errors":        []
}
```

`sc_modality` is `null` for anything that isn't `classification == "single-cell"`.

`sizes` is informational; the classifier itself is extension-only and does not use size as a tiebreaker. It's exposed so downstream tools can sanity-check for empty or suspiciously small files.

## When to reach for `references/classification_rules.md`

Open it when the user asks *why* a specific extension maps to a specific class, when handling an ambiguous case (e.g. "I have .tsv files that are actually scRNA count matrices"), or when deciding whether a new format should be added. The SKILL.md body is deliberately short; edge-case reasoning lives in the reference.

## Interaction pattern

When a user hands you a directory and asks "what is this":

1. **Run the classifier first, show the evidence.** Don't guess from the user's phrasing — even if they say "bulk", run the tool and confirm.
2. **On `mixed`, push back.** Default behavior is to fail; offer the user the three escape hatches in order of preference: (a) split directories, (b) remove the minority files, (c) pass `--allow-mixed` if they really meant to keep both.
3. **On `single-cell`, redirect based on `sc_modality`.** Don't pretend `build-taiji-input` will work.
   - **`multiome`**: RNA + ATAC are already paired at the cell level. Route to a single-cell Taiji workflow that consumes paired inputs (e.g. sc-Taiji with MuData/AnnData + fragments). No co-embedding step needed.
   - **`separate-assay`**: the user has RNA and ATAC from DIFFERENT cell populations. They MUST co-embed before any downstream regulatory-network work. Point them at Signac's integrate_atac workflow: <https://stuartlab.org/signac/articles/integrate_atac>. The recipe is: quantify the multiome peaks in the scATAC object with `FeatureMatrix`, compute LSI on the scATAC side (TF-IDF → SVD), merge the objects, then `FindIntegrationAnchors(reduction="rlsi")` and `IntegrateEmbeddings()` to produce a shared `integrated_lsi` reduction. Once the datasets share a low-dim embedding, RNA labels/expression can be transferred onto the ATAC cells via `TransferData`/`GeneActivity`, and the combined object is suitable for a single-cell Taiji run.
   - **`sc-undetermined`**: ask the user explicitly whether this is multiome or separate-assay — don't guess. Rename hint: adding `multiome`, `rna`, or `atac` tokens to filenames (or providing a `.h5mu`) is enough to disambiguate on the next run.
4. **On `unknown`, surface the top extensions.** Usually the directory is one level off — the real data is in a sibling folder, or the files are gzipped in a way the classifier didn't catch.
5. **On `bulk` with the `.tsv`-only warning, probe.** A directory of only `.tsv` files is a weak bulk signature — it could easily be metadata, logs, or a scRNA count matrix that someone exported to TSV. Ask the user to confirm the TSVs are GeneQuant files (gene IDs in column 1) before letting `build-taiji-input` run.

## Why the defaults are what they are

- **Extension-only, no content sniff.** The user's file-format conventions are unambiguous — .h5ad and .rds are never used for bulk, .narrowPeak and .bedpe are never used for single-cell, and .tsv is the acknowledged weak spot (handled with a warning). Buying deeper robustness would mean importing h5py and an R runtime, which is a big cost for a tiny win.
- **Strict on mixed by default.** Silently mixing SC into a bulk pipeline produces garbage with no error at Taiji-invocation time. Loud-fail-early is much cheaper than debug-at-analysis-time.
- **Gate, not filter.** The skill does not *remove* the unrecognized files; it classifies the whole directory. If the user really wants to split bulk from SC files in-place, that's a separate operation and should be explicit.
