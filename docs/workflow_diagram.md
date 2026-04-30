# Taiji-agent end-to-end workflow

## Where each step runs

```
╔════════════════════════════════════════════════════════════════════════════════════╗
║  SANDBOX (Linux aarch64)                  │   MAC / SLURM (Linux x86_64 or macOS)  ║
║  - skills + orchestration                 │   - actual Taiji binary execution      ║
║  - prep + validate                        │   - heavy compute                      ║
║  - workflow-log writer                    │   - reads same workspace folder        ║
╚════════════════════════════════════════════════════════════════════════════════════╝
                  │                                          │
                  └──────── shared workspace folder ─────────┘
                            (data/, dependencies_data/,
                             runs/<name>/Output/...)
```

## Setup phase (one-time per project)

```
                                       ╔══════════════════════╗
                                       ║  fetch-references    ║
                       (any machine)──▶║  --genome hg38       ║──▶ dependencies_data/hg38/
                                       ║  --output ...        ║       ├── genome.fa
                                       ║  --update-genomes-yml║       ├── genes.gtf
                                       ╚══════════════════════╝       └── cisbp_human_2.meme  (vendored)

                                       ╔══════════════════════╗
                                       ║  install-taiji.sh    ║──▶ binaries/
                       (any machine)──▶║  --system <os>       ║       ├── taiji-Ubuntu-x86_64
                                       ╚══════════════════════╝       ├── taiji-CentOS-x86_64
                                                                      └── taiji-macOS-Catalina-10.15
```

## Per-run pipeline

```
                                        user input
                                            │
                                            ▼
                              ┌──────────────────────────┐
                              │ runs/<name>/             │
                              │   ├── samples.csv        │   ◄── manual one-time edit
                              │   └── taiji_config       │
                              │       .template.yml      │   ◄── per-project formula
                              └────────────┬─────────────┘
                                           │
                                           ▼
              ╔════════════════════════════════════════════════════╗
              ║           bash bin/run-taiji.sh runs/<name>        ║
              ║                                                    ║
              ║  ┌────────────────────────────────────────────┐    ║
              ║  │ 1. pick binary by uname -sm                │    ║
              ║  │ 2. regenerate xlsx with REPO_ROOT paths    │    ║
              ║  │ 3. preflight (verify all paths exist)      │    ║
              ║  │ 4. workflow-log start                       │    ║
              ║  │ 5. delegate to taiji-runner ────────────┐  │    ║
              ║  │ 6. workflow-log finalize                │  │    ║
              ║  └────────────────────────────────────────────┘    ║
              ╚═══════════════════════════════════╪════════════════╝
                                                  │
        ┌─────────────────────────────────────────┴─────────────────────────────────┐
        │                                                                           │
        ▼                                                                           ▼
   ┌──────────┐  detect-dataset-type ────────▶  bulk?  ───yes──▶  build-taiji-input
   │ data/    │      │ ▲                                                  │
   │  RA_11_* │      │ │                                                  ▼
   │  OA_02_* │      ▼ │                                          taiji_input.xlsx
   └──────────┘   single-cell?                                    (Active +
                       │ │                                         active_metadata)
                       ▼ │                                                │
                pseudobulk-construct                                      │
                       │ │                                                │
                       └─┴───▶ manifest.tsv ────────────────▶ build-taiji-input
                                                                          │
                                                                          ▼
                                                              ┌──────────────────┐
                                                              │  taiji-runner    │
                                                              │  (per-sample)    │
                                                              └────────┬─────────┘
                                                                       │
                       ┌───────────────────────────────────────────────┼───────────────┐
                       │ runs/<name>/                                  │               │
                       │   Input/<sample>_input/                       │               │
                       │     ├── <sample>_input.tsv     ◄── filtered Active sheet    │
                       │     └── <sample>_config.yml    ◄── 3 placeholders filled    │
                       │   Output/Partial/<sample>_output/                            │
                       │     (taiji's cwd)                                            │
                       │   taiji_config_files.txt       ◄── for SLURM-array compat  │
                       └─────────────────────────────────────┬───────────────────────┘
                                                             │
                                          ┌──────────────────┴──────────────────┐
                                          │  for each sample (sequential by     │
                                          │  default; --parallel N optional):   │
                                          ▼                                     ▼
              ┌──────────────────────────────────┐            ┌──────────────────────────────────┐
              │ ./binaries/taiji-<os>            │            │ ./binaries/taiji-<os>            │
              │   run --config <sample>_config   │            │   run --config <sample>_config   │
              │     +RTS -N<threads>             │            │     +RTS -N<threads>             │
              │   cwd: <sample>_output/          │            │   cwd: <sample>_output/          │
              └──────────────┬───────────────────┘            └──────────────┬───────────────────┘
                             │                                               │
                             ▼                                               ▼
              ┌──────────────────────────────────┐            ┌──────────────────────────────────┐
              │ Output/Partial/RA_11_output/     │            │ Output/Partial/OA_02_output/     │
              │   ★ GeneRanks.tsv                │            │   ★ GeneRanks.tsv                │
              │   ★ GeneRanks_PValues.tsv        │            │   ★ GeneRanks_PValues.tsv        │
              │   Network/RA_11/                 │            │   Network/OA_02/                 │
              │     edges_combined.csv           │            │     edges_combined.csv           │
              │     edges_binding.csv            │            │     edges_binding.csv            │
              │     nodes.csv                    │            │     nodes.csv                    │
              │   ATACSeq/openChromatin.bed.gz   │            │   ATACSeq/...                    │
              │   RNASeq/expression_profile.tsv  │            │   RNASeq/...                     │
              │   GENOME/genome.index            │            │   GENOME/genome.index            │
              │   sciflow.db                     │            │   sciflow.db                     │
              └──────────────────────────────────┘            └──────────────────────────────────┘
```

## Side-channel: workflow-log

The `workflow-log` skill writes alongside the pipeline, capturing every stage.
Sibling skills auto-attach via `<run_dir>/log/active_run` pointer file.

```
              ┌─── runs/<name>/log/<run_id>.md       (human-readable narrative)
              ├─── runs/<name>/log/<run_id>.jsonl    (machine-parseable, per-stage)
              ├─── runs/<name>/log/index.jsonl       (cross-run start/finalize index)
              ├─── runs/<name>/log/active_run        (pointer for sibling auto-attach)
  workflow ──┤
   stages    ├─── runs/<name>/log/RA_11.taiji.stdout (per-sample)
              ├─── runs/<name>/log/RA_11.taiji.stderr
              ├─── runs/<name>/log/OA_02.taiji.stdout
              └─── runs/<name>/log/OA_02.taiji.stderr

  Stages emitted:
    init                       (workflow-log.start)
    detect                     (detect-dataset-type)
    pseudobulk                 (pseudobulk-construct, SC branch only)
    fetch_references           (fetch-references, optional)
    build_input                (build-taiji-input)
    taiji_run.<sample>         (taiji-runner, one per sample)
    final_summary              (workflow-log.finalize)
```

## Verified output ★ — what proof of success looks like

For sample RA_11 on the CHEM280 chr22 demo (4 min 44 s wall time on Mac):

```
GeneRanks.tsv   (765 TF rankings)
  ────────────────────────────────────────────
                       RA_11
  SP2_HUMAN.H11MO.0.A   5.26e-3  ◄── ubiquitous SP-family TF (top 1)
  SP1_HUMAN.H11MO.0.A   4.40e-3
  SP3_HUMAN.H11MO.0.B   3.94e-3
  PATZ1_HUMAN.H11MO.0.C 3.25e-3  ◄── chr22q11.21-resident TF (chr22 bias confirmed)
  WT1_HUMAN.H11MO.0.C   3.15e-3
  ...
  ────────────────────────────────────────────
  p < 0.01:    3 TFs
  p < 0.05:   71 TFs
  p < 0.10:  112 TFs
  p ≥ 0.10:  652 TFs
```

If `GeneRanks.tsv` exists with ~765 rows AND stderr's last line is `Program finishes successfully`, the run is healthy.

## Skill catalog

| # | Skill | What it does | Where it can run |
|---|-------|--------------|------------------|
| 1 | `detect-dataset-type` | Classify directory as bulk / SC / mixed / unknown; sub-classify SC modality (multiome / separate-assay / undetermined) | sandbox or Mac |
| 2 | `fetch-references` | Idempotent per-genome staging (FASTA + GTF from GENCODE/EMBL-EBI; MEME from vendored CIS-BP at `cisbp_database/`); size-tolerance check; optional samtools faidx | Mac (sandbox egress blocks the hosts) |
| 3 | `build-taiji-input` | Walk data dir + samples manifest → produce Active + active_metadata xlsx for bulk Taiji | sandbox or Mac |
| 4 | `pseudobulk-construct` | Bridge SC objects (.rds/.h5ad) to bulk-Taiji inputs via Seurat WNN clustering + per-cluster MACS3 peak calling | Mac (R + MACS3 not in sandbox) |
| 5 | `taiji-runner` | Per-sample Taiji 1.3.0 orchestrator: split xlsx → per-sample TSVs/configs → invoke binary per sample → validate `GeneRanks.tsv` presence | Mac (Taiji binary x86_64 only) |
| 6 | `workflow-log` | Per-run audit log (md + jsonl + index) auto-attached by sibling skills via `active_run` pointer | sandbox or Mac |

## Skill composition rules (load-bearing)

- **CLI is a thin wrapper over a pure library function.** Every skill exposes a callable Python (or R) function returning a dataclass; the CLI is just argparse + `print(result.summary())`.
- **No file content sniffing.** All classification routes on extensions and filenames.
- **Strict-on-mixed by default.** Silent SC contamination of a bulk pipeline produces garbage.
- **Soft-import siblings.** Cross-skill calls use `sys.path.insert(0, sibling/scripts)` + `try: from <m> import ...`. A standalone install of one skill always works.
- **Reference docs in `references/`, not in `SKILL.md`.** SKILL.md stays short.
- **Every recognized format/pattern ships with a fixture and an assertion** under `<skill>-workspace/fixtures/`.
- **Always cite limitations.** Eunice prefers explicit "this can fail when X" over silent confidence.
- **Output formats designed for downstream chaining.** When skill A's output feeds skill B, A emits exactly what B expects with no transformation step in between.
