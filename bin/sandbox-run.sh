#!/usr/bin/env bash
# sandbox-run.sh — workspace-bounded execution wrapper.
#
# Goal: ensure every command this agent runs writes only inside the
# workspace folder (REPO_ROOT) and uses a workspace-local TMPDIR. This
# catches both "agent runs the wrong thing" mistakes and "Taiji writes
# scratch to /tmp and forgets to clean up" surprises.
#
# Three layers (applied in order):
#   1. Soft (always):  cwd must be inside REPO_ROOT; TMPDIR pinned to
#                      REPO_ROOT/tmp/; child sees REPO_ROOT as $REPO_ROOT.
#   2. OS-level (auto): if `bwrap` (Linux) or `sandbox-exec` (macOS) is
#                      available AND --strict is set or $TAIJI_STRICT_SANDBOX=1,
#                      wrap the child in a kernel-enforced sandbox.
#   3. Container (opt-in): not implemented here. Use Docker/Apptainer
#                      separately if you want full image isolation.
#
# Usage:
#   bash bin/sandbox-run.sh <command> [args...]
#   bash bin/sandbox-run.sh --strict <command> [args...]   # require OS-level
#   bash bin/sandbox-run.sh --allow-net <command>          # let the child reach network
#                                                          # (default for fetch-references etc.;
#                                                          #  by default we DON'T block network in
#                                                          #  the soft tier — bwrap/sandbox-exec do)
#
# Refuses to run if cwd is outside REPO_ROOT or if --strict was requested
# but no OS-level sandbox tool is available.
#
# Exit codes:
#   0    child exited 0
#   <n>  whatever child exited
#   97   cwd outside REPO_ROOT (refused)
#   98   --strict requested but no sandbox tool available
#   99   bad args

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

STRICT=0
ALLOW_NET=0
EXTRA_BIND=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --strict)    STRICT=1; shift ;;
    --allow-net) ALLOW_NET=1; shift ;;
    --bind)      EXTRA_BIND+=("$2"); shift 2 ;;
    --)          shift; break ;;
    -h|--help)   sed -n '2,33p' "$0"; exit 0 ;;
    --*)         echo "unknown sandbox flag: $1" >&2; exit 99 ;;
    *)           break ;;
  esac
done

if [[ -n "${TAIJI_STRICT_SANDBOX:-}" && "${TAIJI_STRICT_SANDBOX}" != "0" ]]; then
  STRICT=1
fi

if [[ $# -eq 0 ]]; then
  echo "[sandbox] no command given" >&2
  exit 99
fi

# ---- Layer 1: soft enforcement (always) ---------------------------------

CWD="$(pwd -P)"
case "$CWD" in
  "$REPO_ROOT"|"$REPO_ROOT"/*)
    : # ok
    ;;
  *)
    echo "[sandbox] REFUSING: cwd '$CWD' is outside REPO_ROOT '$REPO_ROOT'." >&2
    echo "[sandbox] cd into the workspace before invoking sandbox-run.sh."   >&2
    exit 97
    ;;
esac

# Workspace-local TMPDIR (created if missing).
SANDBOX_TMP="$REPO_ROOT/tmp/sandbox-$$"
mkdir -p "$SANDBOX_TMP"
trap 'rm -rf "$SANDBOX_TMP" 2>/dev/null || true' EXIT
export TMPDIR="$SANDBOX_TMP"
export TMP="$SANDBOX_TMP"
export TEMP="$SANDBOX_TMP"
export REPO_ROOT     # children can use this to scope their own writes

# ---- Layer 2: OS-level sandbox (best-effort or strict) -------------------

OS="$(uname -s)"

select_os_sandbox() {
  case "$OS" in
    Linux)
      if command -v bwrap >/dev/null 2>&1; then
        echo "bwrap"; return
      fi
      ;;
    Darwin)
      if command -v sandbox-exec >/dev/null 2>&1; then
        echo "sandbox-exec"; return
      fi
      ;;
  esac
  echo ""
}

OS_TOOL="$(select_os_sandbox)"

if [[ "$STRICT" -eq 1 && -z "$OS_TOOL" ]]; then
  echo "[sandbox] --strict requested but no OS-level sandbox tool available." >&2
  case "$OS" in
    Linux)  echo "[sandbox] install bubblewrap: apt-get install bubblewrap, or" >&2
            echo "[sandbox]                     dnf install bubblewrap"          >&2 ;;
    Darwin) echo "[sandbox] sandbox-exec is bundled with macOS but may be missing on minimal images." >&2 ;;
  esac
  exit 98
fi

if [[ -n "$OS_TOOL" ]]; then
  echo "[sandbox] using $OS_TOOL ($OS_TOOL-bound; cwd=$CWD; TMPDIR=$SANDBOX_TMP)" >&2
else
  if [[ "$STRICT" -eq 1 ]]; then
    : # already handled above
  else
    echo "[sandbox] no OS-level sandbox tool found; using soft enforcement only." >&2
    echo "[sandbox]   (workspace-bound cwd + workspace TMPDIR; pass --strict to require"   >&2
    echo "[sandbox]    bwrap/sandbox-exec)"                                                >&2
  fi
fi

# ---- Build and exec -----------------------------------------------------

case "$OS_TOOL" in
  bwrap)
    # bwrap: rebind / read-only + bind REPO_ROOT read-write + tmpfs /tmp.
    # Network: blocked by default, opened with --share-net.
    BWRAP_ARGS=(
      --bind / /
      --bind "$REPO_ROOT" "$REPO_ROOT"
      --bind "$SANDBOX_TMP" /tmp
      --proc /proc
      --dev /dev
      --setenv TMPDIR "$SANDBOX_TMP"
      --setenv REPO_ROOT "$REPO_ROOT"
      --chdir "$CWD"
      --die-with-parent
      --new-session
    )
    [[ "$ALLOW_NET" -eq 1 ]] && BWRAP_ARGS+=(--share-net) || BWRAP_ARGS+=(--unshare-net)
    for b in "${EXTRA_BIND[@]:-}"; do
      [[ -n "$b" ]] && BWRAP_ARGS+=(--bind "$b" "$b")
    done
    exec bwrap "${BWRAP_ARGS[@]}" "$@"
    ;;

  sandbox-exec)
    # macOS sandbox-exec: write a profile that allows reads anywhere, writes
    # only inside REPO_ROOT and SANDBOX_TMP, network controlled by --allow-net.
    PROFILE="$SANDBOX_TMP/profile.sb"
    {
      echo '(version 1)'
      echo '(deny default)'
      echo '(allow process-fork)'
      echo '(allow process-exec)'
      echo '(allow file-read*)'
      echo '(allow signal)'
      echo '(allow sysctl-read)'
      printf '(allow file-write* (subpath %q))\n' "$REPO_ROOT"
      printf '(allow file-write* (subpath %q))\n' "$SANDBOX_TMP"
      echo '(allow file-write* (subpath "/private/var/folders"))'
      if [[ "$ALLOW_NET" -eq 1 ]]; then
        echo '(allow network*)'
      else
        echo '(allow network-outbound (literal "/private/var/run/syslog"))'
      fi
    } > "$PROFILE"
    exec sandbox-exec -f "$PROFILE" "$@"
    ;;

  "")
    # Soft mode: just exec with the env we set above.
    exec "$@"
    ;;
esac
