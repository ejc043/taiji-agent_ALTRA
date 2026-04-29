# Per-sample directory layout

Mirrors the convention set by the existing UCSD wrapper at
`Taiji/taiji_wrapper-uniqueGen.py`. Drop-in compatible.

## Structure

For an xlsx with N samples (e.g. RA_11, OA_02), under `--run-dir`:

```
<run_dir>/
├── Input/
│   ├── <sample1>_input/
│   │   ├── <sample1>_input.tsv          # filtered Active sheet, sample1's rows only
│   │   └── <sample1>_config.yml         # template materialized for sample1
│   ├── <sample2>_input/
│   │   ├── <sample2>_input.tsv
│   │   └── <sample2>_config.yml
│   └── ...
├── Output/
│   └── Partial/
│       ├── <sample1>_output/            # taiji's cwd; all per-sample outputs
│       │   ├── GeneRanks.tsv               # the headline output (TF rankings)
│       │   ├── GeneRanks_PValues.tsv       # p-values for each TF rank
│       │   ├── Network/<sample1>/
│       │   │   ├── edges_combined.csv      # TF -> target edges
│       │   │   ├── edges_binding.csv       # raw binding edges
│       │   │   └── nodes.csv               # graph nodes
│       │   ├── ATACSeq/                    # ATAC intermediates (TFBS, peaks)
│       │   ├── RNASeq/                     # RNA intermediates
│       │   └── sciflow.db                  # Taiji workflow state DB
│       ├── <sample2>_output/
│       └── ...
├── taiji_config_files.txt               # newline-list of per-sample config paths
└── log/
    ├── <sample1>.taiji.stdout
    ├── <sample1>.taiji.stderr
    ├── <sample2>.taiji.stdout
    └── <sample2>.taiji.stderr
```

## Per-sample TSV format

Same 8 columns as the xlsx Active sheet:

```
type	id	group	rep	path	tags	format	cohort
RNA-seq	RNA-RA_11	RA_11	1	/abs/path/RA_11_RNA.tsv	GeneQuant		RA
ATAC-seq	ATAC-RA_11	RA_11	1	/abs/path/RA_11_REP1...narrowPeak		NarrowPeak	RA
HiC	HiC_RA_11	RA_11	1	/abs/path/RA_11_hicloops_chr22.bedpe	ChromosomeLoop		RA
```

One row per (sample × modality). For one sample with RNA + ATAC + HiC, that's three rows.

## Per-sample config (after substitution)

Starting from a template like:

```yaml
input:      [insert_input_filepath_here]
output_dir: [insert_output_directory_here]
genome:     [insert_genome_filepath_here]
tmp_dir:    /tmp/taiji_RA_OA_chr22
annotation: ${REPO_ROOT}/dependencies_data/hg38/genes.gtf
motif_file: ${REPO_ROOT}/cisbp_database/cisbp_human_2.meme
```

The runner produces (for sample RA_11 with REPO_ROOT=/path/to/agent):

```yaml
input:      /path/to/agent/runs/<run_dir>/Input/RA_11_input/RA_11_input.tsv
output_dir: /path/to/agent/runs/<run_dir>/Output/Partial/RA_11_output
genome:     /path/to/agent/dependencies_data/hg38/genome.fa  # from active_metadata.vcf_Location[RA_11]
tmp_dir:    /tmp/taiji_RA_OA_chr22
annotation: /path/to/agent/dependencies_data/hg38/genes.gtf
motif_file: /path/to/agent/cisbp_database/cisbp_human_2.meme
```

## taiji_config_files.txt

Newline-delimited list of per-sample config paths. Used by the existing
UCSD `Taiji_UniGen_array.sh` SLURM array script (`sed -n -e "$SLURM_ARRAY_TASK_ID p"`
on this file picks the right config for each array task).

This file isn't used by the runner itself when invoked locally — the runner
iterates over its own per-sample list — but it's emitted for compatibility
with the UCSD pipeline.

## Per-sample stdout/stderr

Each sample's output streams are captured separately:

- `log/<sample>.taiji.stdout` — Taiji's stdout (mostly progress messages)
- `log/<sample>.taiji.stderr` — Taiji's stderr (warnings, errors, RTS messages)

Persisted regardless of success/failure. The last 20 lines of stderr are
also embedded in the workflow-log entry as `stderr_tail` for quick scanning
without opening the file.

## Why this exact layout

- **Input/ and Output/Partial/ are separate** so you can `rm -rf Output/` to
  wipe results without losing your generated TSVs/configs. Useful when
  iterating on a config tweak.
- **Per-sample subdirs** so Taiji's intermediates (which can be hundreds of
  MB per run) don't collide between samples.
- **`taiji_config_files.txt` at the run-dir root** makes SLURM-array
  submission a one-line `sbatch -a 1-N%P <script> taiji_config_files.txt`
  (matching the existing UCSD pattern).
- **`log/` shared for all samples** so the workflow-log skill (which writes
  its `<run_id>.md` and `<run_id>.jsonl` there) lives alongside per-sample
  Taiji logs.
