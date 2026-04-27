# Log file format

## Filenames

```
<work_dir>/log/
├── <run_id>.md            primary, human-readable
├── <run_id>.jsonl         primary, machine-readable
├── index.jsonl            cross-run summary index
└── active_run             plain text containing the current run_id
```

`<run_id>` = `YYYYMMDD_HHMMSS_<4hex>` — e.g., `20260426_142315_a3f9`. Sortable by name; collision-resistant within a single second.

## JSONL schema (one object per line)

Every line carries:

| Field    | Type     | Notes                                                    |
|----------|----------|----------------------------------------------------------|
| `ts`     | string   | ISO 8601 UTC timestamp (`...Z` suffix)                   |
| `run_id` | string   | The run identifier — same on every line in this file     |
| `stage`  | string   | `init` / `detect` / `pseudobulk` / `build_input` / `taiji_run` / `<stage>.error` / `final_summary` / custom |
| `status` | string   | `pass` / `fail` / `error` / `info`                       |

Stage-specific fields then follow at the top level. See `per_stage_fields.md` for the exact field set per stage.

### Example line (init / header)

```json
{"ts":"2026-04-26T14:23:15Z","run_id":"20260426_142315_a3f9","stage":"init",
 "host":"taiji-login01","user":"ejc043","work_dir":"/scratch/ejc043/run01",
 "platform":"Linux-4.18.0-553-x86_64","python":"3.11.6","genome":"hg38",
 "reference_paths":{"fasta":"/db/hg38.fa","gtf":"/db/hg38.gtf"},
 "conda_env":"taiji-agent","taiji_version":"taiji 1.3.0",
 "doctor_snapshot":{"ran":true,"ok":true,
   "skills":[{"skill":"detect-dataset-type","missing_required":[]}, ...]}}
```

### Example line (detect)

```json
{"ts":"2026-04-26T14:23:18Z","run_id":"20260426_142315_a3f9","stage":"detect",
 "status":"pass","classification":"bulk","sc_modality":null,
 "sample_count_hint":{"bulk_total":36,"sc_total":0,
   "by_kind":{".tsv":12,".narrowpeak":12,".bedpe":12}},
 "warnings":[],"errors":[],
 "paths_scanned":["/scratch/ejc043/run01/data"]}
```

### Example line (final_summary)

```json
{"ts":"2026-04-26T15:11:42Z","run_id":"20260426_142315_a3f9",
 "stage":"final_summary","status":"success","note":null,"duration_s":2907.3,
 "stages_seen":["init","detect","build_input","taiji_run","final_summary"]}
```

## Markdown layout

```
# Taiji workflow log — `<run_id>`

_Started: <UTC timestamp>_

## Run header

- **host:** ...
- **user:** ...
- ...

---

## [<UTC ts>] Stage A — detect-dataset-type

- **status:** `pass`
- **classification:** `bulk`
- **sc_modality:** _(none)_
- **sample_count_hint:** `{"bulk_total":36,"sc_total":0,...}`
- **evidence:**
<details><summary>show</summary>

```json
{ ... full evidence dict ... }
```
</details>

## [<UTC ts>] Stage B — pseudobulk-construct

...

## [<UTC ts>] Final summary

- **status:** `success`
- **duration_s:** `2907.3`
- **stages_seen:** `["init","detect","build_input","taiji_run","final_summary"]`
```

Fields in any single section come from the same dict that produced the JSONL line — they're rendered as bullet points if scalar / small, and folded into a `<details>` block if the JSON repr exceeds 200 chars (so the markdown stays scannable).

## index.jsonl schema

One object per event (`start` and `finalize`):

```json
{"event":"start","run_id":"20260426_142315_a3f9","ts":"2026-04-26T14:23:15Z",
 "work_dir":"/scratch/ejc043/run01","genome":"hg38","note":"CHEM280 demo"}
{"event":"finalize","run_id":"20260426_142315_a3f9","ts":"2026-04-26T15:11:42Z",
 "status":"success","note":null,"duration_s":2907.3}
```

Useful queries against `index.jsonl`:

```bash
# All failed runs
jq 'select(.event=="finalize" and .status!="success")' log/index.jsonl

# Run durations
jq 'select(.event=="finalize") | {run_id, duration_s}' log/index.jsonl

# Runs that started but never finalized (workflow crashed)
comm -23 \
  <(jq -r 'select(.event=="start")    | .run_id' log/index.jsonl | sort) \
  <(jq -r 'select(.event=="finalize") | .run_id' log/index.jsonl | sort)
```

## active_run

Plain text file containing exactly the active run's `run_id`, no newline necessary. Sibling skills `cat` it (or `Path.read_text()` it) to discover the run they should log into. Removed by `finalize`. If the file is missing or empty, `TaijiLog.attach_active(...)` returns `None` and skills skip logging.

## Stability guarantees

- **Field names are append-only.** New fields can be added to any stage. Existing fields are not renamed; if a field needs to change semantics, a new field is added and the old one is left in place for one release cycle.
- **`stage` and `status` enums are stable.** `init`, `detect`, `pseudobulk`, `build_input`, `taiji_run`, `<stage>.error`, `final_summary`. Status: `pass`, `fail`, `error`, `info`.
- **`ts` is always UTC, ISO 8601 with `Z` suffix.** Local-time stamps are an anti-pattern in audit trails because daylight-saving silently corrupts the ordering.
- **JSONL is append-only.** Existing lines are never rewritten. Markdown is also append-only at the section level (the header is written once at start and never edited).
