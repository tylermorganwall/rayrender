# Benchmark repair verification — 2026-09-18

## Confirmed against the checkout, deployment and available Actions logs

- The deployed [dashboard](https://www.rayrender.net/benchmarks/index.html) had
  212 rows: 64 scalar build failures, 22 historical SSE build failures, and 126
  successful SSE render rows. Its CSV exactly matched the fetched history branch.
- [Actions run 35189477771](https://github.com/tylermorganwall/rayrender/actions/runs/35189477771)
  reported success despite scalar `build_failed` rows. Its scalar install log
  shows `simd_sgn()` trying to pass `int[4]` storage to `__m128i` intrinsics. The
  checkout selects integer storage with `HAS_SSE`/`HAS_NEON`, but the function
  selected SSE using `__SSE2__`. A standalone x86 compilation reproduced this.
- The same standalone check with SSE4.1 exposed the analogous `simd_mul()` guard
  problem. Both guards now match the storage backend. The existing scalar sign
  implementation compares rather than negating/subtracting its argument, so
  `INT_MIN` and `INT_MAX` are safe.
- Actual CI compilation selected C++20 with `-g -O2`; the JSON set other
  standard-specific flags but omitted `CXX20FLAGS`. Corrected configurations
  have revision 2 and are separated from historical flag meanings.
- GNU time parsing split at the colon inside `(h:mm:ss or m:ss)`. Missing fields
  could also reach an `if (NA)` condition. Full-label/value splitting and strict
  duration parsing now cover minutes, hours, missing and malformed values.
- Both scenes explicitly assigned `NA` to BVH timing; there was no structured
  production timing path. Native constructor timing now reaches the worker JSON
  and CSV. The measurement is opt-in and its scope is documented in README.
- The dashboard was reproduced in Chromium both online and from generated local
  HTML. Neither produced an uncaught JavaScript error or failed network request.
  It did split every iteration by its timestamp, report misleading `n=1`
  summaries, hide build failures behind `NA`, and mix scene/settings series.
  Its raw script CSV also retained HTML entities in diagnostic strings.
- Collection accepted whatever artifacts existed and could succeed with no rows.
  New collection checks every expected configuration and retains explicit failure
  records. The workflow remains failed after publishing failed-attempt diagnostics.

## Confirmed design risk, not an observed stale deployment

The old pkgdown push workflow had no dependency on benchmark publication, so it
could read history before publication. The inspected live CSV was current;
therefore stale history was not demonstrated as a cause of that particular page.
The new completion chain and source/history freshness check address the race.

## Checks performed locally

- R/testthat benchmark regression suite: **136 passing assertions**, no warnings
  or skips at the last recorded run. Covers parsing, successful and failed
  validation, missing/duplicate iterations and artifacts, identities, settings,
  finite/range checks, compiler flags/backend verification, legacy data,
  idempotent CSV publication, reruns, latest valid comparison and site failures.
- Standalone C++ compile/execute tests: ARM scalar and NEON, plus x86 scalar,
  SSE2 and SSE4.1 under Rosetta. Integer signs cover zero, positive/negative values,
  `INT_MIN`, `INT_MAX` and their neighbors. All pass.
- Installed-package `test-benchmark-timing.R`: **7 passing assertions**. Ordinary
  renders lack timing attributes; benchmark renders have numeric intervals,
  multiple mesh/world trees are counted, and capture resets on subsequent renders.
- Playwright/Chromium fixtures: success, failed, partial, legacy, empty, latest
  attempt versus last complete comparison, filters, escaped diagnostic text,
  available metrics alongside absent BVH, desktop/mobile overflow and an explicit
  HTTP 503 load failure. No unexpected console errors or failed requests.
- Browser checks on regenerated deployed history: latest SSE medians aggregate
  three samples (0.978 s and 1.591 s); scalar failures are visible; absent BVH
  timing is explicit. Available render/memory metrics still render.
- Reduced **real scalar** run through the harness: two scenes × two iterations,
  16 × 16, two samples, 40 objects, zero warmups. It compiled and installed the
  package, validated four successful measurements, collected the expected local
  configuration, merged history twice (still four rows), generated the site and
  passed Chromium inspection. Production settings were not reduced.
  Actual Apple clang 16 commands selected `-std=gnu++20 -O3 -march=native` and
  preprocessing confirmed the scalar backend. Sphere/mixed BVH medians were
  0.000028167/0.000032271 seconds, with 1/8 trees respectively. Render medians
  were 0.018/0.018 seconds. These smoke measurements are not production comparisons.
- The local R Makeconf forces a compiler-cache wrapper that was blocked by the
  sandbox. Early attempts correctly failed validation and retained logs. The
  successful temporary config used `/usr/bin/clang++` and `MAKEFLAGS=-j8 CCACHE=`.
- `actionlint 1.7.7` passed both changed workflows; `yaml::read_yaml()` parsed
  both; `git diff --check` passed. Modified R files were formatted with `air`.

## Limits

No new GitHub Actions run, remote history push, or deployment was performed.
The full Linux scalar/SSE package matrix was not run on this ARM R installation;
SSE was verified by standalone x86 compilation/execution, with the production
matrix configured to exercise the full builds on Linux. GNU time parsing used
fixtures and existing Linux logs; this macOS smoke run does not report GNU RSS.
An optional attempt to run the entire pre-existing native Catch suite returned
no XML because the smoke package was built without `NOT_CRAN`; that broader
suite was not executed. The focused native sign and R timing tests did run.
