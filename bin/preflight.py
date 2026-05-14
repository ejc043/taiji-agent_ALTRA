#!/usr/bin/env python3
"""preflight.py -- check what's installed and what's needed before running setup.

Reports env profile state, genome reference presence, and data classification.
Outputs JSON so callers (scripts, Claude) can decide what to install.

Works with system Python >= 3.6. Uses the taiji-agent conda env's Python for
subprocess calls to fetch.py and detect.py (both require Python >= 3.7).

Usage:
    python bin/preflight.py [--data-dir DIR] [--genome GENOME [...]]
                            [--refs-output DIR] [--env-name NAME]
                            [--system macos|centos|ubuntu]
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional

REPO_ROOT = Path(__file__).resolve().parent.parent

DETECT_PY = REPO_ROOT / "skills/detect-dataset-type/scripts/detect.py"
FETCH_PY  = REPO_ROOT / "skills/fetch-references/scripts/fetch.py"

PROFILE_REQUIRES = {
    "bulk":        "base",
    "single-cell": "sc",
    "mixed":       "sc",
    "unknown":     "base",
}

# sc requires base to already be installed
PROFILE_IMPLIES = {"sc": ["base", "sc"], "base": ["base"]}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _env_python(env_path):
    # type: (Optional[str]) -> Optional[str]
    """Return {env}/bin/python if the env exists, else None."""
    if not env_path:
        return None
    py = Path(env_path) / "bin" / "python"
    return str(py) if py.exists() else None


# ---------------------------------------------------------------------------
# Env detection  (mirrors install.sh get_env_path logic)
# ---------------------------------------------------------------------------

def _find_env(env_name):
    # type: (str) -> Optional[str]
    # 1. Currently activated
    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    if conda_prefix and Path(conda_prefix).name == env_name and Path(conda_prefix).is_dir():
        return conda_prefix

    # 2. Solver env list
    for solver in ("micromamba", "mamba", "conda"):
        if shutil.which(solver):
            try:
                out = subprocess.run(
                    [solver, "env", "list"],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                    universal_newlines=True, timeout=15,
                ).stdout
                for line in out.splitlines():
                    parts = line.split()
                    if parts and parts[0] == env_name:
                        p = parts[-1]
                        if Path(p).is_dir():
                            return p
            except Exception:
                pass

    # 3. Candidate paths
    user = os.environ.get("USER", "")
    mamba_root = os.environ.get("MAMBA_ROOT_PREFIX", "")
    home = Path.home()
    candidates = [
        mamba_root + "/envs/" + env_name if mamba_root else None,
        str(home / ".local/share/mamba/envs" / env_name),
        str(home / "micromamba/envs" / env_name),
        str(home / "miniforge3/envs" / env_name),
        str(home / "miniconda3/envs" / env_name),
        str(home / "anaconda3/envs" / env_name),
        str(home / "opt/miniconda3/envs" / env_name),
        "/opt/conda/envs/" + env_name,
        "/stg3/data1/" + user + "/.local/share/mamba/envs/" + env_name if user else None,
        "/stg3/data1/" + user + "/micromamba/envs/" + env_name if user else None,
        "/stg3/data1/" + user + "/miniconda3/envs/" + env_name if user else None,
    ]
    for c in candidates:
        if c and Path(c).is_dir():
            return c

    return None


def _read_installed_profiles(env_path):
    # type: (str) -> List[str]
    marker = Path(env_path) / ".taiji-agent-profiles"
    if not marker.exists():
        return []
    installed = []
    for line in marker.read_text().splitlines():
        parts = line.strip().split()
        if parts:
            installed.append(parts[0])
    return installed


# ---------------------------------------------------------------------------
# Check functions
# ---------------------------------------------------------------------------

def check_env(env_name="taiji-agent"):
    # type: (str) -> Dict
    path = _find_env(env_name)
    if not path:
        return {
            "path": None,
            "profiles_installed": [],
            "profiles_missing": ["base", "sc"],
            "ok": False,
        }
    installed = _read_installed_profiles(path)
    missing = [p for p in ["base", "sc"] if p not in installed]
    return {
        "path": path,
        "profiles_installed": installed,
        "profiles_missing": missing,
        "ok": len(missing) == 0,
    }


def check_data(data_dir, env_path=None):
    # type: (str, Optional[str]) -> Dict
    if not data_dir:
        return {
            "path": None,
            "classification": None,
            "sc_modality": None,
            "required_profile": None,
            "ok": None,
            "note": "no --data-dir provided",
        }
    p = Path(data_dir)
    if not p.exists():
        return {
            "path": data_dir,
            "classification": None,
            "sc_modality": None,
            "required_profile": None,
            "ok": False,
            "error": "directory does not exist: " + data_dir,
        }
    if not DETECT_PY.exists():
        return {
            "path": data_dir,
            "ok": False,
            "error": "detect.py not found at " + str(DETECT_PY),
        }

    env_py = _env_python(env_path)
    python_bin = env_py if env_py else sys.executable
    try:
        proc = subprocess.run(
            [python_bin, str(DETECT_PY), str(p), "--format", "json"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True, timeout=30,
        )
        result = json.loads(proc.stdout)
        classification = result.get("classification", "unknown")
        sc_modality = result.get("sc_modality")
        required_profile = PROFILE_REQUIRES.get(classification, "base")
        return {
            "path": data_dir,
            "classification": classification,
            "sc_modality": sc_modality,
            "required_profile": required_profile,
            "ok": True,
        }
    except Exception as e:
        return {"path": data_dir, "ok": False, "error": str(e)}


def _check_refs_direct(genome, refs_output):
    # type: (str, str) -> Dict
    """File-existence fallback when env Python is not available."""
    gdir = Path(refs_output) / genome
    present = []
    missing = []
    for kind, fname in [("fasta", "genome.fa"), ("gtf", "genes.gtf")]:
        f = gdir / fname
        if f.exists() and f.stat().st_size > 0:
            present.append(kind)
        else:
            missing.append(kind)
    meme_files = list(gdir.glob("*.meme")) if gdir.exists() else []
    if meme_files:
        present.append("motif")
    else:
        missing.append("motif")
    return {
        "output_dir": refs_output,
        "ok": len(missing) == 0,
        "files_present": present,
        "files_missing": missing,
        "note": "direct file check (no env Python available)",
    }


def check_refs(genomes, refs_output, env_path=None):
    # type: (List[str], str, Optional[str]) -> Dict
    if not genomes:
        return {}
    env_py = _env_python(env_path)
    python_bin = env_py if env_py else sys.executable
    can_use_fetch = FETCH_PY.exists() and (
        env_py is not None or sys.version_info >= (3, 7)
    )
    results = {}
    for genome in genomes:
        if not can_use_fetch:
            results[genome] = _check_refs_direct(genome, refs_output)
            continue
        try:
            proc = subprocess.run(
                [python_bin, str(FETCH_PY),
                 "--genome", genome, "--output", refs_output, "--check", "--json"],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                universal_newlines=True, timeout=30,
            )
            data = json.loads(proc.stdout)
            files = data.get("files", {})
            files_present = [k for k, v in files.items() if v.get("status") == "present"]
            files_missing = [k for k, v in files.items() if v.get("status") != "present"]
            results[genome] = {
                "output_dir": refs_output,
                "ok": len(files_missing) == 0,
                "files_present": files_present,
                "files_missing": files_missing,
            }
        except Exception:
            results[genome] = _check_refs_direct(genome, refs_output)
    return results


# ---------------------------------------------------------------------------
# Action planner
# ---------------------------------------------------------------------------

def build_actions(env, data, refs, system, refs_output):
    # type: (Dict, Dict, Dict, str, str) -> List[str]
    actions = []

    required = data.get("required_profile") or "base"
    installed = set(env.get("profiles_installed", []))
    need = set(PROFILE_IMPLIES.get(required, [required])) - installed
    if need:
        install_profile = "sc" if "sc" in need else "base"
        sys_flag = " --system " + system if system else ""
        actions.append("bash bin/install.sh --profile " + install_profile + sys_flag)

    for genome, ref_info in refs.items():
        if not ref_info.get("ok"):
            sys_flag = " --system " + system if system else ""
            out_flag = (" --fetch-output " + refs_output
                        if refs_output != "dependencies_data" else "")
            actions.append(
                "bash bin/install.sh --fetch-references --genome "
                + genome + out_flag + sys_flag
            )

    return actions


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Check what's installed and what's needed before setup."
    )
    ap.add_argument("--data-dir",    default="",                  help="data directory to classify")
    ap.add_argument("--genome",      nargs="*", default=[],       help="genome(s) to check refs for")
    ap.add_argument("--refs-output", default="dependencies_data", help="where genome refs live")
    ap.add_argument("--env-name",    default="taiji-agent")
    ap.add_argument("--system",      default="",                  help="macos|centos|ubuntu")
    args = ap.parse_args()

    env  = check_env(args.env_name)
    data = check_data(args.data_dir, env_path=env.get("path"))
    refs = check_refs(args.genome or [], args.refs_output, env_path=env.get("path"))
    actions = build_actions(env, data, refs, args.system, args.refs_output)

    report = {
        "env":            env,
        "data":           data,
        "references":     refs,
        "actions_needed": actions,
        "ready":          len(actions) == 0,
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
