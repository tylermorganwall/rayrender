#!/usr/bin/env bash
set -euo pipefail

# Release branches can omit the maintainer suite. Use a fixed revision from
# this repository, while measuring the package sources at the checked-out SHA.
if [[ ! -d tools/benchmarks ]]; then
  # Exclude the harness from the working-tree status recorded with results.
  printf '\n/tools/benchmarks/\n' >> "$(git rev-parse --git-path info/exclude)"
  git archive 98655fb48d372b1e8c4b221635608b96caf61db2 -- tools/benchmarks ':!tools/benchmarks/results' | tar -x
fi
