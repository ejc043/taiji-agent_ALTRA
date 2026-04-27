# Per-stage field reference

Each stage emits a single JSONL object (and the corresponding Markdown section). The four canonical stages and their expected fields are below. New skills can call `TaijiLog.append_custom(stage, title, payload)` for non-canonical stages.

## init (the run header)

Written once by `TaijiLog.start(...)`.

| Field              | Source                                                 |
|--------------------|--------------------------------------------------------|
| `host`             | `socket.gethostname()`                                 |
| `user`             | `getpass.getuser()`                                    |
| `work_dir`         | absolute path passed to `start()`                      |
| `platform`         | `platform.platform()` (Linux kernel + distro hints)    |
| `python`           | major.minor.patch                                       |
| `genome`           | passed to `start()`                                    |
| `reference_paths`  | dict passed to `start()` (label -> absolute path)      |
| `conda_env`        | `$CONDA_DEFAULT_ENV` if set                            |
| `taiji_version`    | first line of `binaries/taiji --version` if installed  |
| `doctor_snapshot`  | `bin/doctor.sh --json --no-r` parsed; `{"ran": false, "reason": "..."}` if not run |
| `note`             | optional free-text                                      |

## detect (Stage A — detect-dataset-type)

Written by sibling skill via `log.append_detect(detect_result)`.

| Field                     | Type           | Notes                                       |
|---------------------------|----------------|---------------------------------------------|
| `classification`          | string         | `bulk` / `single-cell` / `mixed` / `unknown`|
| `sc_modality`             | string \| null | `multiome` / `separate-assay` / `sc-undetermined` / null |
| `sample_count_hint.bulk_total` | int      | Sum of all bulk-file counts                 |
| `sample_count_hint.sc_total`   | int      | Sum of all SC-file counts                   |
| `sample_count_hint.by_kind`    | dict     | `{ext_or_pattern: count}` for both bulk + SC |
| `warnings`                | list[str]      | from DetectResult.warnings                  |
| `errors`                  | list[str]      | from DetectResult.errors                    |
| `paths_scanned`           | list[str]      | what the detector walked                    |
| `evidence`                | dict           | full bulk_files + sc_files dicts (folded)   |
| `status`                  | enum           | `pass` if no errors, else `fail`            |

`sample_count_hint` is a hint, not a final tally — a 12-sample bulk dataset with RNA + ATAC + HiC reports `bulk_total: 36` because that's three files per sample. The actual sample count is locked in at `build_input`.

## pseudobulk (Stage B — pseudobulk-construct)

Written via `log.append_pseudobulk(summary)`. Recommended fields (all optional, but the more the better for audit):

| Field                     | Type           | Notes                                       |
|---------------------------|----------------|---------------------------------------------|
| `clustering_signal`       | string         | `wnn` / `rna` / `atac`                      |
| `n_cells_in`              | int            | Total cells loaded                          |
| `n_cells_kept`            | int            | After QC + small-cluster filter             |
| `chosen_resolution`       | float          | The resolution the binary search converged on |
| `resolution_trace`        | list[dict]     | Full search history from `resolution_trace.json` |
| `n_clusters_kept`         | int            | After `--min-cluster-cells` filter          |
| `n_clusters_dropped`      | int            | Below the floor                             |
| `dropped_cluster_sizes`   | dict           | `{cluster_id: n_cells}` for the dropped     |
| `mean_cluster_size`       | float          | Across kept clusters                        |
| `n_groups_planned`        | int            | Total `(cluster × metadata × value)` rows   |
| `manifest_path`           | string         | Path to the produced `manifest.tsv`         |
| `metadata_cols`           | list[str]      | Auto-detected or user-specified             |
| `transferred_label_col`   | string         | Only for separate-assay; usually `predicted.id` |
| `n_with_label`            | int            | Cells with a non-NA transferred label       |
| `coembedding_required`    | bool           | True if classification was separate-assay   |
| `status`                  | enum           | `pass` / `fail`                             |

For the SC case, this stage is the most information-dense entry in the log — it captures the answer to "what did clustering decide and why." Always include `chosen_resolution`, `resolution_trace`, and the `dropped_cluster_sizes` dict.

## build_input (Stage C — build-taiji-input)

Written via `log.append_build_input(summary)`.

| Field                     | Type           | Notes                                       |
|---------------------------|----------------|---------------------------------------------|
| `xlsx_path`               | string         | Absolute path to the produced xlsx          |
| `xlsx_sha256_short`       | string         | First 12 hex chars of SHA-256 (auto-added)  |
| `n_active_rows`           | int            | Rows in the `Active` sheet                  |
| `n_dropped_samples`       | int            | Samples filtered out (missing files etc.)   |
| `dropped_reasons`         | dict           | `{sample_name: reason}`                     |
| `samples_path`            | string         | Path to the input samples manifest          |
| `genome`                  | string         | Genome tag passed to the skill              |
| `gate_classification`     | string         | What detect-dataset-type returned at gate   |
| `gate_sc_modality`        | string \| null | sc_modality at gate, if any                 |
| `cohorts`                 | list[str]      | Distinct cohort values in the Active sheet  |
| `status`                  | enum           | `pass` / `fail`                             |

`xlsx_sha256_short` is appended automatically by `log.append_build_input` — caller doesn't need to provide it. Useful for spotting silent corruption between runs.

## taiji_run (Stage E — running the Taiji binary)

Written via `log.append_taiji_run(summary)`. The Taiji binary itself doesn't know about this skill, so the caller (a wrapper script or the future taiji submit skill) must populate the fields after the binary exits.

| Field                     | Type           | Notes                                       |
|---------------------------|----------------|---------------------------------------------|
| `command`                 | string         | The exact command line invoked              |
| `exit_code`               | int            | 0 = success                                 |
| `duration_s`              | float          | Wall-clock seconds                          |
| `output_dir`              | string         | Path to Taiji's output directory            |
| `output_dir_size_mb`      | float          | Auto-computed if `output_dir` exists        |
| `n_pagerank_files`        | int            | Count of `*_pagerank.tsv` in output_dir     |
| `n_edges_files`           | int            | Count of `*_edges_combined.csv`             |
| `peak_memory_kb`          | int            | If wrapped in `/usr/bin/time -v` or pulled from sacct |
| `stderr_tail`             | string         | Last ~20 lines, useful for failure diagnosis |
| `slurm_job_id`            | string         | If submitted via SLURM                      |
| `status`                  | enum           | `pass` if `exit_code == 0`, else `fail`     |

## <stage>.error

Written by any caller via `log.append_error(stage, message, stderr_tail, traceback)`. Used when a stage failed loudly enough to want a separate entry rather than just `status:"fail"` on the main entry.

| Field          | Type     | Notes                                          |
|----------------|----------|-----------------------------------------------|
| `error_message`| string   | One-line summary                               |
| `stderr_tail`  | string   | Last ~20 lines of stderr                       |
| `traceback`    | string   | Python traceback if applicable                 |
| `status`       | enum     | always `error`                                 |

## final_summary

Written once by `log.finalize(...)`. Always the last line in the JSONL.

| Field          | Type           | Notes                                              |
|----------------|----------------|----------------------------------------------------|
| `note`         | string \| null | Optional free-text                                 |
| `duration_s`   | float          | Wall-clock seconds for the whole workflow          |
| `stages_seen`  | list[str]      | Every distinct `stage` value observed in this run  |
| `status`       | enum           | `success` / `failure` / `partial`                  |

`stages_seen` is computed from the JSONL itself — useful for spotting "we said the run finished but never wrote a `taiji_run` entry, which means the actual binary call was skipped." A common diagnostic.

## Adding a new stage

For new skills, prefer `append_custom(stage, title, payload)` over a new top-level method. The `payload` dict gets folded into the JSONL entry verbatim (under `_to_jsonable` coercion). Stage names should be `snake_case` and unambiguous (`taiji_submit`, `network_post_filter`, etc.). Update this doc when you add one.
