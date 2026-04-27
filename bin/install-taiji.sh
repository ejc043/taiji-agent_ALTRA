#!/usr/bin/env bash
# install-taiji.sh — download the Taiji binary for the chosen --system.
#
# Usage:
#   bash bin/install-taiji.sh                    # auto-detect OS
#   bash bin/install-taiji.sh --system centos    # explicit
#   bash bin/install-taiji.sh --system ubuntu
#   bash bin/install-taiji.sh --system macos
#   bash bin/install-taiji.sh --version v1.3.0   # pin a specific release
#   bash bin/install-taiji.sh --force            # re-download even if present
#
# Naming convention on disk: file is always saved as the bare GitHub asset
# name. No version suffix, no -latest. --version controls only which release
# tag the URL points at; the file overwrites in place on re-fetch.
#
#   --system centos  ->  binaries/taiji-CentOS-x86_64
#   --system ubuntu  ->  binaries/taiji-Ubuntu-x86_64
#   --system macos   ->  binaries/taiji-macOS-XX-XX  (manual download; varies by macOS version)
#
# Plus a convenience symlink at binaries/taiji -> the active binary, so
# bin/run-taiji.sh and other callers can find it via a stable path.
#
# By-name auto-detection: if a matching binary already exists in binaries/,
# the script symlinks it and skips the download (idempotent re-runs).
#
# Exit codes: 0 ok, 2 unknown --system, 3 download failed, 4 chmod failed.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BINDIR="${REPO_ROOT}/binaries"

SYSTEM=""
VERSION="latest"
FORCE=0

usage() {
  sed -n '2,15p' "$0"
  exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --system)  SYSTEM="$2"; shift 2 ;;
    --version) VERSION="$2"; shift 2 ;;
    --force)   FORCE=1;     shift   ;;
    -h|--help) usage 0 ;;
    *) echo "unknown arg: $1" >&2; usage 2 ;;
  esac
done

# Auto-detect OS if --system not given.
if [[ -z "$SYSTEM" ]]; then
  uname_s="$(uname -s)"
  case "$uname_s" in
    Linux)
      if [[ -f /etc/os-release ]]; then
        # shellcheck disable=SC1091
        . /etc/os-release
        case "${ID:-}" in
          ubuntu|debian)             SYSTEM="ubuntu" ;;
          centos|rhel|rocky|almalinux|fedora) SYSTEM="centos" ;;
          *) echo "unrecognized Linux distro '${ID:-unknown}'; pass --system explicitly" >&2; exit 2 ;;
        esac
      else
        echo "no /etc/os-release; pass --system explicitly" >&2; exit 2
      fi
      ;;
    Darwin) SYSTEM="macos" ;;
    *)      echo "unsupported OS '$uname_s'; pass --system explicitly" >&2; exit 2 ;;
  esac
  echo "[install-taiji] auto-detected --system $SYSTEM"
fi

# Map --system to the release asset name (or asset glob for macOS).
# ASSET_GLOB ends with `*` so the by-name detection picks up files named
# either `taiji-CentOS-x86_64` (no suffix), `taiji-CentOS-x86_64-v1.3.0`
# (explicit version), or `taiji-CentOS-x86_64-latest` (legacy).
case "$SYSTEM" in
  centos) ASSET="taiji-CentOS-x86_64";  ASSET_GLOB="taiji-CentOS-x86_64*" ;;
  ubuntu) ASSET="taiji-Ubuntu-x86_64";  ASSET_GLOB="taiji-Ubuntu-x86_64*" ;;
  macos)
    # macOS asset filenames vary by macOS major version (taiji-macOS-Catalina-10.15,
    # taiji-macOS-BigSur-11, etc.). We can't pick the right one to download
    # automatically, but we CAN auto-detect a matching file already in
    # binaries/ and symlink it without re-downloading.
    ASSET=""
    ASSET_GLOB="taiji-macOS-*"
    ;;
  *) echo "unknown --system '$SYSTEM'; expected centos|ubuntu|macos" >&2; exit 2 ;;
esac

mkdir -p "$BINDIR"

DEST="${BINDIR}/taiji"

# ---- Auto-detection: existing binary on disk -----------------------------
# Look for a pre-existing matching binary in BINDIR by name (not just the
# `taiji` symlink). Useful when the user has committed binaries directly,
# or when an earlier install already downloaded one.
EXISTING=""
if [[ "$FORCE" -eq 0 ]]; then
  # Search for any file matching the system's glob, prefer the most-recent.
  while IFS= read -r f; do
    [[ -f "$f" ]] && EXISTING="$f" && break
  done < <(ls -1t "$BINDIR/"$ASSET_GLOB 2>/dev/null | grep -vE '/taiji$' || true)
fi

if [[ -n "$EXISTING" ]]; then
  echo "[install-taiji] found existing matching binary: $EXISTING"
  if [[ -L "$DEST" || -e "$DEST" ]] && [[ "$(readlink "$DEST" 2>/dev/null || echo "")" == "$(basename "$EXISTING")" ]]; then
    echo "[install-taiji] symlink already correct: $DEST -> $(basename "$EXISTING")"
  else
    # `ln -sfn` atomically replaces an existing symlink without needing rm
    # first, and `-n` prevents dereferencing if $DEST happens to be a symlink
    # to a directory.
    if ln -sfn "$(basename "$EXISTING")" "$DEST" 2>/dev/null; then
      chmod +x "$EXISTING" 2>/dev/null || true
      echo "[install-taiji] symlinked: $DEST -> $(basename "$EXISTING")"
    else
      echo "[install-taiji] could not (re)create symlink at $DEST" >&2
      echo "[install-taiji]   if $DEST is a real file owned by another user, remove it manually" >&2
      echo "[install-taiji]   then re-run, OR call this binary directly: $EXISTING" >&2
      # Non-fatal — bin/run-taiji.sh's binary auto-pick logic doesn't
      # require the $DEST symlink (it picks by uname -sm). The symlink is
      # for convenience; failure to create it shouldn't block the install.
    fi
  fi
  echo "[install-taiji] skipping download (pass --force to re-fetch)"
  exit 0
fi

# ---- macOS download path: cannot pick automatically ----------------------
if [[ "$SYSTEM" == "macos" && -z "$ASSET" ]]; then
  echo "[install-taiji] no taiji-macOS-* binary found in $BINDIR." >&2
  echo "[install-taiji] macOS asset filename varies by macOS version, so we can't" >&2
  echo "[install-taiji] auto-download. Browse:" >&2
  echo "[install-taiji]   https://github.com/Taiji-pipeline/Taiji/releases/${VERSION}" >&2
  echo "[install-taiji] and grab the matching taiji-macOS-XX-XX asset, drop it at:" >&2
  echo "[install-taiji]   ${BINDIR}/taiji-macOS-<version>" >&2
  echo "[install-taiji] then re-run this script (or just create the symlink yourself):" >&2
  echo "[install-taiji]   ln -s taiji-macOS-<version> ${BINDIR}/taiji && chmod +x ${BINDIR}/taiji-macOS-<version>" >&2
  exit 0
fi

# ---- Linux download path -------------------------------------------------
# URL: GitHub serves /releases/latest/download/<asset> as a 302 to the
# real /releases/download/<tag>/<asset>; curl -L follows the redirect.
# `--version` controls which release tag the URL points at, but the file
# on disk is always saved as the bare asset name (e.g. taiji-CentOS-x86_64)
# — no version suffix, no -latest. Re-fetching overwrites the same file.
if [[ "$VERSION" == "latest" ]]; then
  URL="https://github.com/Taiji-pipeline/Taiji/releases/latest/download/${ASSET}"
else
  URL="https://github.com/Taiji-pipeline/Taiji/releases/download/${VERSION}/${ASSET}"
fi
TARGET="${BINDIR}/${ASSET}"

if [[ -e "$TARGET" && "$FORCE" -eq 0 ]]; then
  echo "[install-taiji] $TARGET already exists (pass --force to re-download)"
  ln -sfn "$ASSET" "$DEST" 2>/dev/null || true
  echo "[install-taiji] current binary: $("$TARGET" --version 2>/dev/null || echo unknown)"
  exit 0
fi

echo "[install-taiji] downloading $URL"
echo "[install-taiji]   -> $TARGET"
# `-fL`: fail on HTTP errors and follow redirects (releases/latest/download
#        always 302s to the actual release tag).
# `--retry`: handle transient 5xx/network blips.
# `-#`: progress bar to stderr.
if ! curl -fL --retry 3 --retry-delay 2 -# -o "$TARGET" "$URL"; then
  echo "" >&2
  echo "[install-taiji] download failed." >&2
  echo "[install-taiji]   url:    $URL" >&2
  echo "[install-taiji]   target: $TARGET" >&2
  echo "[install-taiji] common causes:" >&2
  echo "[install-taiji]   - sandbox / proxy blocking github.com (try from your Mac directly)" >&2
  echo "[install-taiji]   - asset filename changed in a Taiji release; verify at:" >&2
  echo "[install-taiji]       https://github.com/Taiji-pipeline/Taiji/releases/${VERSION}" >&2
  echo "[install-taiji]   - --version pin doesn't exist (pass a real tag like v1.3.0)" >&2
  exit 3
fi

chmod +x "$TARGET" || { echo "[install-taiji] chmod +x failed" >&2; exit 4; }

# Atomically (re)create the convenience symlink.
ln -sfn "$(basename "$TARGET")" "$DEST" 2>/dev/null \
  || echo "[install-taiji] WARN: could not create symlink at $DEST (binary still usable at $TARGET)" >&2

echo "[install-taiji] installed: $TARGET"
[[ -L "$DEST" ]] && echo "[install-taiji] symlink:   $DEST -> $(readlink "$DEST")"
echo "[install-taiji] smoke test: $TARGET --help"
"$TARGET" --help >/dev/null 2>&1 \
  && echo "[install-taiji] OK" \
  || echo "[install-taiji] WARNING: --help did not exit cleanly. The binary may still work; inspect manually."
