#!/usr/bin/env bash
set -euo pipefail

repo_root="$(git rev-parse --show-toplevel)"
codex_tmp="${RAYRENDER_CODEX_TMP:-/private/tmp/codex-projects/rayrender}"

cd "$repo_root"

tools/codex/install-local.sh

export R_LIBS_USER="$codex_tmp/R-lib"
export TMPDIR="$codex_tmp/tmp"
export NOT_CRAN=true
export RAYRENDER_REPO_ROOT="$repo_root"

Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); options(device = function(...) grDevices::pdf(file = NULL)); library(rayrender); ns = asNamespace("rayrender"); for (name in ls(ns, all.names = TRUE)) assign(name, get(name, envir = ns), envir = .GlobalEnv); setwd(Sys.getenv("TMPDIR")); testthat::test_dir(file.path(Sys.getenv("RAYRENDER_REPO_ROOT"), "tests/testthat"), filter = "screen-(line|text)", reporter = "summary", load_package = "none")'

Rscript tools/spectral-tests/check-spectral-assets.R \
  --manifest docs/spectral/assets-manifest.csv \
  --build-source

Rscript tools/spectral-tests/capture-pr0-baselines.R \
  --compare docs/spectral/baselines/pr0-legacy-baselines.csv \
  --output "$codex_tmp/pr0-current-baselines.csv" \
  --environment-output "$codex_tmp/pr0-current-environment.md" \
  --artifact-dir "$codex_tmp/pr0-artifacts"
