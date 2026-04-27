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
# The binary lands at <repo-root>/binaries/taiji and is symlinked to
# <repo-root>/binaries/taiji-<system>-<version> for traceability.
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

# Map --system to the release asset name.
case "$SYSTEM" in
  centos) ASSET="taiji-CentOS-x86_64" ;;
  ubuntu) ASSET="taiji-Ubuntu-x86_64" ;;
  macos)
    # The macOS asset filename varies by macOS major version (taiji-macOS-XX-XX).
    # We cannot reliably pick this from a script — let the user see the
    # release page and pass the exact suffix via --version <tag>.
    echo "[install-taiji] macOS asset filename varies by macOS version." >&2
    echo "[install-taiji] Browse https://github.com/Taiji-pipeline/Taiji/releases/${VERSION}" >&2
    echo "[install-taiji] and grab the matching taiji-macOS-XX-XX asset, then drop it at:" >&2
    echo "[install-taiji]   ${BINDIR}/taiji" >&2
    echo "[install-taiji] and \`chmod +x\` it. Skipping automatic download." >&2
    exit 0
    ;;
  *) echo "unknown --system '$SYSTEM'; expected centos|ubuntu|macos" >&2; exit 2 ;;
esac

# Build the URL.
if [[ "$VERSION" == "latest" ]]; then
  URL="https://github.com/Taiji-pipeline/Taiji/releases/latest/download/${ASSET}"
else
  URL="https://github.com/Taiji-pipeline/Taiji/releases/download/${VERSION}/${ASSET}"
fi

mkdir -p "$BINDIR"

DEST="${BINDIR}/taiji"
TRACEABLE="${BINDIR}/${ASSET}-${VERSION}"

if [[ -e "$DEST" && "$FORCE" -eq 0 ]]; then
  echo "[install-taiji] $DEST already exists (pass --force to re-download)"
  echo "[install-taiji] current version: $("$DEST" --version 2>/dev/null || echo unknown)"
  exit 0
fi

echo "[install-taiji] downloading $URL -> $TRACEABLE"
if ! curl -fL --retry 3 --retry-delay 2 -o "$TRACEABLE" "$URL"; then
  echo "[install-taiji] download failed" >&2
  exit 3
fi

chmod +x "$TRACEABLE" || { echo "[install-taiji] chmod +x failed" >&2; exit 4; }

# Replace any existing taiji symlink/file with one pointing at the new binary.
rm -f "$DEST"
ln -s "$(basename "$TRACEABLE")" "$DEST"

echo "[install-taiji] installed: $DEST -> $(readlink "$DEST")"
echo "[install-taiji] smoke test: $DEST --help"
"$DEST" --help >/dev/null 2>&1 \
  && echo "[install-taiji] OK" \
  || { echo "[install-taiji] WARNING: --help did not exit cleanly. Inspect manually."; }
