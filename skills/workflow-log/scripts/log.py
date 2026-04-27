"""TaijiLog — per-run audit log for the Taiji-agent workflow.

Two artifacts per run, both in <work_dir>/log/:

    <run_id>.md       human-readable narrative (one section per stage)
    <run_id>.jsonl    one JSON object per stage, machine-parseable

Plus two cumulative files:

    log/index.jsonl   one line per run (init + finalize), for cross-run grep
    log/active_run    text file containing the current run_id (deleted on finalize)

The `active_run` file is the bridge between the workflow-log skill and the
other three skills. Each sibling skill (detect-dataset-type, build-taiji-input,
pseudobulk-construct) calls `TaijiLog.attach_active(work_dir)` at runtime —
returns a TaijiLog instance bound to the current run if one exists, or None
otherwise. So:

    * Calling a skill outside a workflow run -> no log entries, silent.
    * Calling under a workflow that started a log -> auto-logs.

The skills don't need to know the run_id; they discover it via the file.

Library API
-----------

    from log import TaijiLog

    log = TaijiLog.start(work_dir=".", genome="hg38")
    log.append_detect(detect_result)
    log.append_pseudobulk(summary_dict)
    log.append_build_input(summary_dict)
    log.append_taiji_run(summary_dict)
    log.append_error(stage="detect", message="...", stderr_tail="...")
    log.finalize(status="success")

CLI
---

    python log.py start    --work-dir . --genome hg38 [--note "..."]
    python log.py append   --stage detect --json-file detect_result.json
    python log.py append   --stage build_input --json '{"xlsx": "...", "rows": 12}'
    python log.py finalize --status success [--note "..."]
    python log.py status                                      # show active run
"""

from __future__ import annotations

import argparse
import getpass
import hashlib
import json
import os
import platform
import socket
import subprocess
import sys
import uuid
from dataclasses import asdict, dataclass, field, is_dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _local_now_compact() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def _gen_run_id() -> str:
    """Sortable, unique-enough run id: <local-timestamp>_<4hex>."""
    return f"{_local_now_compact()}_{uuid.uuid4().hex[:4]}"


def _to_jsonable(obj: Any) -> Any:
    """Best-effort coerce dataclasses / Paths / sets / etc. to JSON-friendly."""
    if obj is None or isinstance(obj, (bool, int, float, str)):
        return obj
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, dict):
        return {k: _to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_jsonable(x) for x in obj]
    if isinstance(obj, set):
        return sorted(_to_jsonable(x) for x in obj)
    if is_dataclass(obj):
        return _to_jsonable(asdict(obj))
    # Fallback: stringify. Loses fidelity but never crashes.
    return repr(obj)


def _sha256_short(path: Path) -> str | None:
    """First 12 hex chars of the SHA256 of `path` (for output corruption checks)."""
    try:
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(65536), b""):
                h.update(chunk)
        return h.hexdigest()[:12]
    except OSError:
        return None


def _safe_call(cmd: list[str], timeout: int = 5) -> str | None:
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        if out.returncode == 0:
            return (out.stdout or "").strip().splitlines()[0] if out.stdout else None
    except Exception:
        pass
    return None


def _doctor_snapshot(work_dir: Path) -> dict:
    """Best-effort doctor.sh --json snapshot. Returns {'ran': False, ...} if
    doctor.sh isn't present or doesn't run. Never raises."""
    doctor = work_dir / "bin" / "doctor.sh"
    if not doctor.exists():
        return {"ran": False, "reason": "bin/doctor.sh not found"}
    try:
        out = subprocess.run(
            ["bash", str(doctor), "--json", "--no-r"],
            capture_output=True, text=True, timeout=30,
        )
        try:
            parsed = json.loads(out.stdout)
        except json.JSONDecodeError:
            return {"ran": False, "reason": "doctor output not JSON",
                    "stderr_tail": out.stderr[-400:]}
        return {"ran": True, "ok": parsed.get("ok"),
                "skills": [{"skill": r["skill"],
                            "missing_required": r.get("missing_required", [])}
                           for r in parsed.get("reports", [])]}
    except Exception as e:
        return {"ran": False, "reason": f"doctor.sh failed: {e}"}


# ---------------------------------------------------------------------------
# TaijiLog
# ---------------------------------------------------------------------------


@dataclass
class TaijiLog:
    work_dir: Path
    run_id: str
    log_dir: Path
    md_path: Path
    jsonl_path: Path
    active_path: Path = field(init=False)
    index_path: Path = field(init=False)

    def __post_init__(self):
        self.active_path = self.log_dir / "active_run"
        self.index_path = self.log_dir / "index.jsonl"

    # --- Lifecycle: start / attach / finalize ----------------------------

    @classmethod
    def start(
        cls,
        work_dir: str | Path = ".",
        genome: str | None = None,
        reference_paths: dict | None = None,
        note: str | None = None,
        run_id: str | None = None,
    ) -> "TaijiLog":
        """Begin a new logged workflow run. Writes the header section and
        marks <log_dir>/active_run so sibling skills auto-attach."""
        work_dir = Path(work_dir).resolve()
        log_dir = work_dir / "log"
        log_dir.mkdir(exist_ok=True)
        rid = run_id or _gen_run_id()
        log = cls(
            work_dir=work_dir,
            run_id=rid,
            log_dir=log_dir,
            md_path=log_dir / f"{rid}.md",
            jsonl_path=log_dir / f"{rid}.jsonl",
        )
        log._write_header(genome=genome,
                          reference_paths=reference_paths or {},
                          note=note)
        log.active_path.write_text(rid)
        log._append_index({"event": "start", "run_id": rid,
                           "ts": _utc_now(), "work_dir": str(work_dir),
                           "genome": genome, "note": note})
        return log

    @classmethod
    def attach_active(cls, work_dir: str | Path = ".") -> "TaijiLog | None":
        """Return a TaijiLog bound to the currently-active run, or None if
        no run is active. Sibling skills call this — None means "skip
        logging silently"."""
        work_dir = Path(work_dir).resolve()
        active = work_dir / "log" / "active_run"
        if not active.exists():
            return None
        try:
            rid = active.read_text().strip()
            if not rid:
                return None
        except OSError:
            return None
        log_dir = work_dir / "log"
        return cls(
            work_dir=work_dir,
            run_id=rid,
            log_dir=log_dir,
            md_path=log_dir / f"{rid}.md",
            jsonl_path=log_dir / f"{rid}.jsonl",
        )

    def finalize(self, status: str = "success", note: str | None = None,
                 duration_s: float | None = None) -> None:
        """Close the run. Writes a final_summary entry and removes the
        active_run pointer."""
        record = {
            "stage": "final_summary",
            "status": status,
            "note": note,
            "duration_s": duration_s,
            "stages_seen": self._stages_seen(),
        }
        self._write_stage("final_summary", "Final summary", record,
                          status=status)
        self._append_index({"event": "finalize", "run_id": self.run_id,
                            "ts": _utc_now(), "status": status, "note": note,
                            "duration_s": duration_s})
        try:
            self.active_path.unlink()
        except (FileNotFoundError, PermissionError, OSError):
            # Best-effort: a stale pointer from a previous session under
            # different uid won't be unlinkable, but the run is still finalized
            # (md/jsonl/index updates already happened above).
            pass

    # --- Per-stage append methods ----------------------------------------

    def append_detect(self, detect_result: Any) -> None:
        """Log the detect-dataset-type result. Accepts a DetectResult
        dataclass or a dict (so the sibling skill can pass either)."""
        d = _to_jsonable(detect_result)
        bulk_files = d.get("bulk_files", {}) or {}
        sc_files = d.get("sc_files", {}) or {}
        rec = {
            "stage": "detect",
            "classification": d.get("classification"),
            "sc_modality": d.get("sc_modality"),
            "sample_count_hint": {
                "bulk_total": sum(len(v) for v in bulk_files.values()),
                "sc_total":   sum(len(v) for v in sc_files.values()),
                "by_kind":    {**{k: len(v) for k, v in bulk_files.items()},
                               **{k: len(v) for k, v in sc_files.items()}},
            },
            "warnings": d.get("warnings", []),
            "errors":   d.get("errors", []),
            "paths_scanned": d.get("paths_scanned", []),
            "evidence":  {"bulk_files": bulk_files, "sc_files": sc_files},
        }
        status = "fail" if rec["errors"] else "pass"
        self._write_stage("detect", "Stage A — detect-dataset-type", rec,
                          status=status)

    def append_pseudobulk(self, summary: dict) -> None:
        """Log a pseudobulk-construct summary. `summary` is expected to
        contain at least: clustering_signal, n_cells_in, n_clusters_kept,
        n_clusters_dropped, dropped_cluster_sizes, chosen_resolution,
        resolution_trace, mean_cluster_size, n_groups_planned, manifest_path,
        and (for separate-assay) transferred_label_col + n_with_label."""
        s = _to_jsonable(summary)
        rec = {"stage": "pseudobulk", **s}
        status = s.get("status", "pass")
        self._write_stage("pseudobulk", "Stage B — pseudobulk-construct", rec,
                          status=status)

    def append_build_input(self, summary: dict) -> None:
        """Log a build-taiji-input summary. Expected fields: xlsx_path,
        n_active_rows, n_dropped_samples, dropped_reasons, samples_path,
        genome, gate_classification, gate_sc_modality."""
        s = _to_jsonable(summary)
        rec = {"stage": "build_input", **s}
        # Hash the xlsx for corruption-detection if it exists.
        xlsx = s.get("xlsx_path")
        if xlsx and Path(xlsx).exists():
            rec["xlsx_sha256_short"] = _sha256_short(Path(xlsx))
        status = s.get("status", "pass")
        self._write_stage("build_input", "Stage C — build-taiji-input", rec,
                          status=status)

    def append_taiji_run(self, summary: dict) -> None:
        """Log a Taiji binary run summary. Expected fields: command,
        exit_code, duration_s, output_dir, n_pagerank_files, n_edges_files,
        peak_memory_kb (optional), stderr_tail (optional)."""
        s = _to_jsonable(summary)
        rec = {"stage": "taiji_run", **s}
        out_dir = s.get("output_dir")
        if out_dir and Path(out_dir).exists():
            try:
                rec["output_dir_size_mb"] = round(
                    sum(p.stat().st_size for p in Path(out_dir).rglob("*")
                        if p.is_file()) / (1024 * 1024), 1)
            except Exception:
                pass
        status = "pass" if s.get("exit_code") == 0 else "fail"
        self._write_stage("taiji_run", "Stage E — Taiji binary run", rec,
                          status=status)

    def append_error(self, stage: str, message: str,
                     stderr_tail: str | None = None,
                     traceback: str | None = None) -> None:
        rec = {"stage": f"{stage}.error", "error_message": message,
               "stderr_tail": stderr_tail, "traceback": traceback}
        self._write_stage(f"{stage}.error", f"ERROR — {stage}", rec,
                          status="error")

    def append_custom(self, stage: str, title: str, payload: dict,
                      status: str = "info") -> None:
        """Escape hatch for stages outside the four canonical ones."""
        rec = {"stage": stage, **_to_jsonable(payload)}
        self._write_stage(stage, title, rec, status=status)

    # --- Internal: file writers ------------------------------------------

    def _write_header(self, genome: str | None,
                      reference_paths: dict, note: str | None) -> None:
        env_info = {
            "ts": _utc_now(),
            "run_id": self.run_id,
            "stage": "init",
            "host": socket.gethostname(),
            "user": getpass.getuser(),
            "work_dir": str(self.work_dir),
            "platform": platform.platform(),
            "python": sys.version.split()[0],
            "genome": genome,
            "reference_paths": _to_jsonable(reference_paths),
            "note": note,
            "conda_env": os.environ.get("CONDA_DEFAULT_ENV"),
            "taiji_version": _safe_call(
                [str(self.work_dir / "binaries" / "taiji"), "--version"]),
            "doctor_snapshot": _doctor_snapshot(self.work_dir),
        }
        # JSONL: header line first.
        with self.jsonl_path.open("a") as fh:
            fh.write(json.dumps(env_info) + "\n")
        # Markdown header.
        md = []
        md.append(f"# Taiji workflow log — `{self.run_id}`\n")
        md.append(f"_Started: {env_info['ts']}_\n")
        md.append(f"\n## Run header\n")
        md.append(f"- **host:** `{env_info['host']}`\n")
        md.append(f"- **user:** `{env_info['user']}`\n")
        md.append(f"- **work_dir:** `{env_info['work_dir']}`\n")
        md.append(f"- **platform:** `{env_info['platform']}`\n")
        md.append(f"- **python:** `{env_info['python']}`\n")
        md.append(f"- **conda_env:** `{env_info['conda_env'] or '(none)'}`\n")
        md.append(f"- **genome:** `{genome or '(unset)'}`\n")
        md.append(f"- **taiji_version:** `{env_info['taiji_version'] or '(not installed)'}`\n")
        if reference_paths:
            md.append("- **reference_paths:**\n")
            for k, v in reference_paths.items():
                md.append(f"  - {k}: `{v}`\n")
        if note:
            md.append(f"- **note:** {note}\n")
        snap = env_info["doctor_snapshot"]
        if snap.get("ran"):
            md.append(f"- **doctor:** {'OK' if snap.get('ok') else 'FAIL'}\n")
            for skill in snap.get("skills", []):
                missing = skill.get("missing_required", [])
                tag = "OK" if not missing else f"MISSING: {','.join(missing)}"
                md.append(f"  - {skill['skill']}: {tag}\n")
        else:
            md.append(f"- **doctor:** not run ({snap.get('reason')})\n")
        md.append("\n---\n")
        self.md_path.write_text("".join(md))

    def _write_stage(self, stage: str, title: str, record: dict,
                     status: str) -> None:
        # Always stamp + tag.
        record = {"ts": _utc_now(), "run_id": self.run_id,
                  "status": status, **record}
        # Append jsonl.
        with self.jsonl_path.open("a") as fh:
            fh.write(json.dumps(_to_jsonable(record)) + "\n")
        # Append markdown section.
        with self.md_path.open("a") as fh:
            fh.write(f"\n## [{record['ts']}] {title}\n\n")
            fh.write(f"- **status:** `{status}`\n")
            for key, value in record.items():
                if key in ("ts", "run_id", "status", "stage"):
                    continue
                fh.write(self._fmt_md_field(key, value))
            fh.write("\n")

    @staticmethod
    def _fmt_md_field(key: str, value: Any) -> str:
        # Big nested dicts/lists go in a collapsed json block; small scalars
        # render as inline bullets.
        if isinstance(value, (dict, list)) and len(json.dumps(value)) > 200:
            block = json.dumps(_to_jsonable(value), indent=2)
            return (f"- **{key}:**\n"
                    f"<details><summary>show</summary>\n\n"
                    f"```json\n{block}\n```\n\n</details>\n")
        if isinstance(value, (dict, list)):
            return f"- **{key}:** `{json.dumps(_to_jsonable(value))}`\n"
        if value is None:
            return f"- **{key}:** _(none)_\n"
        return f"- **{key}:** `{value}`\n"

    def _append_index(self, record: dict) -> None:
        with self.index_path.open("a") as fh:
            fh.write(json.dumps(_to_jsonable(record)) + "\n")

    def _stages_seen(self) -> list[str]:
        if not self.jsonl_path.exists():
            return []
        seen = []
        for line in self.jsonl_path.read_text().splitlines():
            try:
                seen.append(json.loads(line).get("stage"))
            except json.JSONDecodeError:
                continue
        return seen


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _cli_start(args):
    refs = json.loads(args.reference_paths) if args.reference_paths else {}
    log = TaijiLog.start(work_dir=args.work_dir, genome=args.genome,
                         reference_paths=refs, note=args.note,
                         run_id=args.run_id)
    print(log.run_id)


def _cli_append(args):
    log = TaijiLog.attach_active(work_dir=args.work_dir)
    if log is None:
        print(f"no active run in {args.work_dir}/log", file=sys.stderr)
        sys.exit(2)
    if args.json_file:
        payload = json.loads(Path(args.json_file).read_text())
    elif args.json:
        payload = json.loads(args.json)
    else:
        print("--json or --json-file required", file=sys.stderr)
        sys.exit(2)
    {
        "detect":      log.append_detect,
        "pseudobulk":  log.append_pseudobulk,
        "build_input": log.append_build_input,
        "taiji_run":   log.append_taiji_run,
    }.get(args.stage,
          lambda d: log.append_custom(args.stage, args.stage, d))(payload)


def _cli_finalize(args):
    log = TaijiLog.attach_active(work_dir=args.work_dir)
    if log is None:
        print(f"no active run in {args.work_dir}/log", file=sys.stderr)
        sys.exit(2)
    log.finalize(status=args.status, note=args.note,
                 duration_s=args.duration_s)
    print(f"finalized run {log.run_id} -> {args.status}")


def _cli_status(args):
    log = TaijiLog.attach_active(work_dir=args.work_dir)
    if log is None:
        print(f"no active run in {args.work_dir}/log")
        sys.exit(0)
    print(f"active run: {log.run_id}")
    print(f"  md: {log.md_path}")
    print(f"  jsonl: {log.jsonl_path}")
    print(f"  stages so far: {', '.join(log._stages_seen())}")


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("start", help="Begin a new run + write header.")
    s.add_argument("--work-dir", default=".")
    s.add_argument("--genome")
    s.add_argument("--reference-paths", help="JSON dict of label->path.")
    s.add_argument("--note")
    s.add_argument("--run-id", help="Override auto-generated run id.")
    s.set_defaults(func=_cli_start)

    a = sub.add_parser("append", help="Append a stage entry to the active run.")
    a.add_argument("--work-dir", default=".")
    a.add_argument("--stage", required=True,
                   help="detect|pseudobulk|build_input|taiji_run|<custom>")
    a.add_argument("--json",      help="Inline JSON payload.")
    a.add_argument("--json-file", help="Path to a JSON file with the payload.")
    a.set_defaults(func=_cli_append)

    f = sub.add_parser("finalize", help="Close the active run.")
    f.add_argument("--work-dir", default=".")
    f.add_argument("--status", default="success",
                   choices=["success", "failure", "partial"])
    f.add_argument("--note")
    f.add_argument("--duration-s", type=float)
    f.set_defaults(func=_cli_finalize)

    st = sub.add_parser("status", help="Print the active run id, if any.")
    st.add_argument("--work-dir", default=".")
    st.set_defaults(func=_cli_status)

    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
