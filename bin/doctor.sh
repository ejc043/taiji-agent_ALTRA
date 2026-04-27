#!/usr/bin/env bash
# doctor.sh — verify the runtime environment for every skill in this repo.
#
# Walks skills/*/dependencies.yml and checks:
#   - binaries: `which X` (with optional version probe)
#   - python_packages: `python -c "import X"`
#   - r_packages: `Rscript -e 'library(X)'` (only run if Rscript is on PATH)
#
# Emits a per-skill PASS/FAIL table and exits 0 only if all *required* checks
# pass. Optional deps and alternative-set entries are reported but never fail
# the build.
#
# Usage:
#   bash bin/doctor.sh                                # check all skills
#   bash bin/doctor.sh --skill pseudobulk-construct   # one skill only
#   bash bin/doctor.sh --profile base                 # only base-profile skills
#   bash bin/doctor.sh --profile sc                   # base + sc skills
#   bash bin/doctor.sh --profile full                 # base + sc + dev (everything)
#   bash bin/doctor.sh --no-r                         # skip R-package checks
#   bash bin/doctor.sh --quiet                        # only PASS/FAIL totals
#   bash bin/doctor.sh --json                         # machine-readable output

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SKILLS_DIR="${REPO_ROOT}/skills"

ONLY_SKILL=""
ONLY_PROFILE=""
SKIP_R=0
QUIET=0
JSON=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --skill)   ONLY_SKILL="$2"; shift 2 ;;
    --profile) ONLY_PROFILE="$2"; shift 2 ;;
    --no-r)    SKIP_R=1; shift ;;
    --quiet)   QUIET=1; shift ;;
    --json)    JSON=1; shift ;;
    -h|--help) sed -n '2,22p' "$0"; exit 0 ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done

# Heavy lifting in Python — bash can't parse YAML cleanly and the report
# logic is easier to read in one place.
PY=$(command -v python3 || command -v python || true)
if [[ -z "$PY" ]]; then
  echo "doctor: need python3 on PATH" >&2
  exit 2
fi

SKIP_R_ARG=$([[ $SKIP_R -eq 1 ]] && echo "1" || echo "0")
QUIET_ARG=$([[ $QUIET -eq 1 ]] && echo "1" || echo "0")
JSON_ARG=$([[ $JSON -eq 1 ]] && echo "1" || echo "0")

"$PY" - "$SKILLS_DIR" "$ONLY_SKILL" "$SKIP_R_ARG" "$QUIET_ARG" "$JSON_ARG" "$ONLY_PROFILE" <<'PY_EOF'
import json, os, shutil, subprocess, sys
from pathlib import Path

skills_dir, only_skill, skip_r_arg, quiet_arg, json_arg, only_profile = sys.argv[1:7]
SKIP_R = (skip_r_arg == "1")
QUIET = (quiet_arg == "1")
JSON = (json_arg == "1")
ONLY_PROFILE = only_profile.strip()
# Profile composition: passing --profile sc filters to {base, sc}; --profile dev
# filters to {base, dev}; --profile full means everything; empty means no filter.
if ONLY_PROFILE:
    PROFILE_FILTER = {
        "base": {"base"},
        "sc":   {"base", "sc"},
        "dev":  {"base", "dev"},
        "full": {"base", "sc", "dev"},
    }.get(ONLY_PROFILE)
    if PROFILE_FILTER is None:
        print(f"doctor: unknown --profile '{ONLY_PROFILE}' "
              "(expected: base|sc|dev|full)", file=sys.stderr)
        sys.exit(2)
else:
    PROFILE_FILTER = None  # no filter

# Minimal YAML loader without pyyaml dependency — fall back to pyyaml if
# available so users on the conda env get the real parser.
try:
    import yaml as _yaml
    def load_yaml(text): return _yaml.safe_load(text)
except ImportError:
    print("doctor: pyyaml not installed; install it for accurate parsing.",
          file=sys.stderr)
    print("doctor: temporarily using a naive parser — may misread complex YAML.",
          file=sys.stderr)
    def load_yaml(text):
        # Hard-fail rather than guess; we want a clear error.
        raise RuntimeError("pyyaml required: `pip install pyyaml` or `micromamba install pyyaml`")

def _ok(label):     return f"  [OK]   {label}"
def _miss(label):   return f"  [MISS] {label}"
def _opt(label):    return f"  [opt]  {label}"

def check_binary(spec):
    name = spec["name"]
    optional = bool(spec.get("optional", False))
    alts = spec.get("alternatives", []) or []
    found = shutil.which(name)
    if found:
        version = ""
        try:
            out = subprocess.run([name, "--version"], capture_output=True,
                                 text=True, timeout=5)
            version = (out.stdout or out.stderr).splitlines()[0][:80]
        except Exception:
            pass
        return True, optional, f"{name}: {found}" + (f"  ({version})" if version else "")
    if alts:
        for alt in alts:
            alt_path = shutil.which(alt)
            if alt_path:
                return True, optional, f"{name}|{alt}: {alt_path} (alternative)"
        return False, optional, f"{name} (alternatives tried: {','.join(alts)})"
    return False, optional, f"{name}"

def check_python_pkg(spec):
    name = spec["name"]
    # Some packages have different pip vs import names (e.g. pyyaml -> yaml).
    # Allow `import_name` to override the import target.
    import_name = spec.get("import_name", name)
    optional = bool(spec.get("optional", False))
    if name == "python":
        # Just verify the running interpreter satisfies any min_version.
        min_v = spec.get("min_version")
        if min_v:
            cur = ".".join(map(str, sys.version_info[:3]))
            return True, optional, f"python: {cur} (need >= {min_v})"
        return True, optional, f"python: {sys.version.split()[0]}"
    code = (
        f"import importlib, sys; "
        f"m = importlib.import_module('{import_name}'); "
        f"print(getattr(m, '__version__', 'unknown'))"
    )
    try:
        out = subprocess.run([sys.executable, "-c", code], capture_output=True,
                             text=True, timeout=10)
        if out.returncode == 0:
            label = f"{name}: {out.stdout.strip()}"
            if import_name != name:
                label += f"  (import {import_name})"
            return True, optional, label
        last_err = out.stderr.strip().splitlines()[-1] if out.stderr.strip() else "no message"
        return False, optional, f"{name} (import {import_name} failed: {last_err})"
    except Exception as e:
        return False, optional, f"{name} ({e})"

R_AVAILABLE = shutil.which("Rscript") is not None

def check_r_pkg(spec):
    name = spec["name"]
    optional = bool(spec.get("optional", False))
    if SKIP_R or not R_AVAILABLE:
        return None, optional, f"{name} (R not checked)"
    code = (
        f'suppressPackageStartupMessages(library({name})); '
        f'cat(as.character(packageVersion("{name}")))'
    )
    try:
        out = subprocess.run(["Rscript", "-e", code], capture_output=True,
                             text=True, timeout=30)
        if out.returncode == 0:
            return True, optional, f"{name}: {out.stdout.strip()}"
        last_err = (out.stderr.strip().splitlines() or [""])[-1]
        return False, optional, f"{name} (load failed: {last_err})"
    except Exception as e:
        return False, optional, f"{name} ({e})"

def run_for_skill(yml_path):
    spec = load_yaml(yml_path.read_text())
    if not spec:
        return {"skill": yml_path.parent.name, "results": [], "missing_required": []}
    skill = spec.get("skill", yml_path.parent.name)
    results = []
    missing_required = []

    for kind, checker in (("binaries", check_binary),
                          ("python_packages", check_python_pkg),
                          ("r_packages", check_r_pkg)):
        for s in spec.get(kind) or []:
            ok, optional, label = checker(s)
            results.append({"kind": kind, "ok": ok, "optional": optional,
                            "label": label})
            if ok is False and not optional:
                missing_required.append(f"{kind}:{s['name']}")
    return {"skill": skill, "results": results,
            "missing_required": missing_required}

# --- Walk skills ---
skills_root = Path(skills_dir)
yml_paths = sorted(skills_root.glob("*/dependencies.yml"))

if only_skill:
    yml_paths = [p for p in yml_paths if p.parent.name == only_skill]
    if not yml_paths:
        print(f"doctor: no dependencies.yml found for skill '{only_skill}'",
              file=sys.stderr)
        sys.exit(2)

# Apply profile filter: skip skills whose `profile:` field isn't in the
# requested profile set. Skills without a profile field are treated as `base`
# for backwards compat.
if PROFILE_FILTER is not None:
    filtered: list[Path] = []
    for p in yml_paths:
        try:
            spec = load_yaml(p.read_text()) or {}
            skill_profile = spec.get("profile", "base")
        except Exception:
            skill_profile = "base"
        if skill_profile in PROFILE_FILTER:
            filtered.append(p)
    yml_paths = filtered

reports = [run_for_skill(p) for p in yml_paths]
overall_ok = all(not r["missing_required"] for r in reports)

if JSON:
    print(json.dumps({"reports": reports, "ok": overall_ok}, indent=2))
    sys.exit(0 if overall_ok else 1)

# Human-readable output
sep = "=" * 72
for r in reports:
    print(sep)
    status = "PASS" if not r["missing_required"] else "FAIL"
    print(f"{status}  skill: {r['skill']}")
    print(sep)
    if not QUIET:
        for entry in r["results"]:
            tag = (_ok(entry["label"]) if entry["ok"] is True
                   else _opt(entry["label"]) if (entry["ok"] is None or entry["optional"])
                   else _miss(entry["label"]))
            print(tag)
    if r["missing_required"]:
        print(f"  -> missing required: {', '.join(r['missing_required'])}")
    print()

print(sep)
print(f"Overall: {'PASS' if overall_ok else 'FAIL'} "
      f"({sum(1 for r in reports if not r['missing_required'])}/{len(reports)} skills clean)")
print(sep)
if not R_AVAILABLE and not SKIP_R:
    print("note: Rscript not on PATH, R-package checks were skipped. "
          "Activate the conda env or pass --no-r to silence this notice.")

sys.exit(0 if overall_ok else 1)
PY_EOF
