---
name: workflow-log
description: Construct and update a per-run audit log for the Taiji-agent workflow. Each run gets a Markdown narrative + JSONL machine-readable record in <work_dir>/log/, and a one-line summary lands in log/index.jsonl for cross-run grep. Sibling skills (detect-dataset-type, build-taiji-input, pseudobulk-construct) auto-attach to the active run via a small `active_run` pointer file — when no run is active the skills run as before, when one is active they emit log entries automatically. Trigger on phrases like "log this run", "write a log file", "track the workflow", "audit the pipeline", "what did the workflow do", "summary of the run", or anywhere the user wants reproducibility / traceability across the detect → pseudobulk → build → taiji chain.
---

# Workflow log (per-run audit trail for the Taiji-agent pipeline)

## What this skill produces

For each workflow run, in `<work_dir>/log/`:

| File                   | Purpose                                                              |
|------------------------|----------------------------------------------------------------------|
| `<run_id>.md`          | Human-readable narrative — one section per stage, render in any editor |
| `<run_id>.jsonl`       | One JSON object per stage — machine-parseable, append-only           |
| `index.jsonl`          | One line per run (start + finalize) for cross-run grep / dashboards  |
| `active_run`           | Plain text holding the current `run_id`; deleted on `finalize`       |

`<run_id>` is `<YYYYMMDD_HHMMSS>_<4hex>` — sortable, unique enough that two runs in the same second still don't collide.

## What gets logged

**Run header (once at start):** UTC timestamp, hostname, user, working directory, platform, Python version, conda env, genome tag, Taiji binary version, full reference-data paths, optional free-text note, and a snapshot of `bin/doctor.sh --json --no-r` so you can see which deps were available when the run started.

**Per-stage entries** (each carries `ts`, `run_id`, `stage`, `status`):

- `detect` — classification, sc_modality, sample-count breakdown by file kind, warnings, errors, full evidence dict.
- `pseudobulk` — clustering signal, n_cells_in / n_cells_kept, chosen resolution, full resolution_trace, n_clusters_kept / n_clusters_dropped + dropped sizes, mean cluster size, n_groups_planned, manifest path.
- `build_input` — xlsx path (with SHA-256 short-hash), Active-sheet row count, dropped samples and reasons, samples manifest path, gate's classification + sc_modality.
- `taiji_run` — full command, exit code, wall-clock duration, output dir + total size, count of `*_pagerank.tsv` and `*_edges_combined.csv` files, optional peak memory (from `/usr/bin/time -v` or `sacct`), optional stderr tail.
- `<stage>.error` — error message, stderr tail, traceback. Logged when a stage explicitly throws.
- `final_summary` — overall status, optional duration_s, list of every stage seen.

## Library API (preferred — sibling skills use this)

```python
from log import TaijiLog

# At workflow start:
log = TaijiLog.start(
    work_dir=".",
    genome="hg38",
    reference_paths={"fasta": "/db/hg38.fa", "gtf": "/db/hg38.gtf",
                     "meme":  "/db/cisbp_human_2.meme"},
    note="CHEM280 WI26 demo dataset",
)

# Each skill, internally:
log = TaijiLog.attach_active(work_dir=".")     # returns None if no active run
if log:
    log.append_detect(detect_result)           # accepts DetectResult dataclass or dict
    log.append_pseudobulk({...summary...})
    log.append_build_input({...summary...})
    log.append_taiji_run({...summary...})
    log.append_error(stage="build_input", message="...", stderr_tail="...")

# At the end:
log.finalize(status="success", duration_s=1234.5)
```

`attach_active` is the bridge — it reads `<work_dir>/log/active_run` and binds to the current run. If no run is active it returns None, so calls degrade silently when the user runs a single skill outside a workflow.

## CLI

```bash
# Start a run (returns the run_id on stdout)
RUN_ID=$(python scripts/log.py start \
  --work-dir . \
  --genome hg38 \
  --reference-paths '{"fasta":"/db/hg38.fa","gtf":"/db/hg38.gtf"}' \
  --note "CHEM280 demo")

# Append a stage entry from JSON on disk or inline
python scripts/log.py append --stage detect --json-file detect_result.json
python scripts/log.py append --stage build_input \
  --json '{"xlsx_path":"taiji_input.xlsx","n_active_rows":12}'

# Check what's running
python scripts/log.py status

# Close the run
python scripts/log.py finalize --status success --duration-s 1234.5
```

## Interaction pattern

1. **Always start the log first.** A workflow that calls the skills before `TaijiLog.start(...)` runs will produce no log — which is fine, but means you can't reconstruct what happened. The recommended entry point for any non-trivial workflow is `python scripts/log.py start ...` immediately followed by the skill chain.

2. **Sibling skills auto-attach.** Once `active_run` exists, every call to `detect`, `pseudobulk`, or `build-taiji-input` writes its own log entry. No flags, no per-skill `--log-dir` arg.

3. **One run at a time per `<work_dir>`.** Calling `TaijiLog.start(...)` while another run is active overwrites the `active_run` pointer (the previous run's md/jsonl files stay intact, but new entries flow to the new run). Set a different `--work-dir` to run two pipelines in parallel.

4. **Finalize on every exit path.** Including failure. Use `--status failure` so the index reflects what really happened — searching `log/index.jsonl` for `"status":"failure"` is the canonical "show me the runs that broke" query.

5. **The `active_run` file is the source of truth.** If a workflow crashes mid-run, the file is left in place — meaning the next skill call will keep appending to the same log. Manually `rm log/active_run` to force a reset, or just call `finalize --status failure` to do it cleanly.

## When to reach for the references

- `references/log_format.md` — exact field schema for both Markdown and JSONL output. Use when writing a tool that consumes the log.
- `references/per_stage_fields.md` — what each stage records and why. Use when extending a skill to log additional fields.

## Why the defaults are what they are

- **Markdown + JSONL, not just one or the other.** Markdown is readable; JSONL is greppable. Building both costs nothing because every stage already structures its data — the same dict feeds both writers.
- **Append-only, no rewrites.** Logs that mutate are logs that lie. Every stage is a new line in JSONL and a new section in Markdown.
- **`active_run` as a file, not a global env var.** Files survive shell hops, env-var loss between Python and R, etc. Two skills running in two different processes both see the same `active_run` file.
- **`status` is enum-like (`pass`/`fail`/`error`/`info`)**, not free text. Lets you grep the JSONL for failure modes without regex acrobatics.
- **No silent failures from the log skill itself.** If `attach_active` raises, sibling skills get a `None` and continue without logging — they never crash the workflow because the log writer had a bug.

## Limitations

- File-locking is not used. Two processes appending to the same `<run_id>.jsonl` simultaneously can interleave bytes within a single line. Pragmatically rare because each stage is a single fast write, but worth knowing if you ever parallelize.
- Hostname / user / conda env are best-effort: containerized envs may report unhelpful values. They go into the header verbatim — don't trust them blindly for compliance evidence.
- The `doctor.sh` snapshot at start is `--no-r` to avoid a 30s R-package check on every workflow run. Run `doctor.sh` separately if you want the full R-side report.
