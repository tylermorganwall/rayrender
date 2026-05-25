#!/usr/bin/env bash
set -euo pipefail

repo_root="$(git rev-parse --show-toplevel)"
codex_tmp="${RAYRENDER_CODEX_TMP:-/private/tmp/codex-projects/rayrender}"

src_copy="$codex_tmp/src-copy"
r_lib="$codex_tmp/R-lib"
tmpdir="$codex_tmp/tmp"
ccache_dir="$codex_tmp/ccache"
makevars="$codex_tmp/Makevars"

mkdir -p "$src_copy" "$r_lib" "$tmpdir" "$ccache_dir"

rm -rf "$src_copy"

rsync -a \
  --delete \
  --exclude ".git" \
  --exclude ".Rproj.user" \
  --exclude "tools/benchmarks/results" \
  --exclude "*.o" \
  --exclude "*.so" \
  --exclude "*.dll" \
  --exclude "*.dylib" \
  "$repo_root/" "$src_copy/"

cat > "$makevars" <<'EOF'
CC = ccache clang
CXX = ccache clang++
CXX11 = ccache clang++
CXX14 = ccache clang++
CXX17 = ccache clang++
CXX20 = ccache clang++
EOF

export R_LIBS_USER="$r_lib"
export TMPDIR="$tmpdir"
export CCACHE_DIR="$ccache_dir"
export R_MAKEVARS_USER="$makevars"

if [ -d /opt/homebrew/opt/ccache/libexec ]; then
  export PATH="/opt/homebrew/opt/ccache/libexec:$PATH"
fi

if [ -d /usr/local/opt/ccache/libexec ]; then
  export PATH="/usr/local/opt/ccache/libexec:$PATH"
fi

R CMD INSTALL --preclean -l "$r_lib" "$src_copy"