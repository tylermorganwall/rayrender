# Render benchmarks

The R harness installs each configuration into an isolated library, runs a fresh
R worker for every warmup/measured iteration, and writes CSV rows plus diagnostic
logs. Production settings remain in `configs/default.json`: 50 × 50 pixels,
32 samples, 3,000 objects, one render thread. CI compares
`o3_native_no_simd` and `o3_native_simd_sse` on Linux x86-64.

## Execution and validation

```sh
Rscript tools/benchmarks/run_render_benchmarks.R \
  --repo . --ref current --config tools/benchmarks/configs/default.json \
  --only-config o3_native_no_simd,o3_native_simd_sse \
  --benchmarks bvh_many_spheres,bvh_mixed_primitives \
  --iterations 3 --warmup 1 --extra-r-arg '--time-build true' \
  --output /tmp/render-benchmarks.csv --keep-workdirs --no-append
```

`current` includes working-tree changes; another `--ref` is exported with
`git archive`. Scene definitions come from the current harness. `--dry-run`
prints the plan. `--timeout-seconds`, `--workdir-root`, and `--extra-r-arg`
control execution. A reduced local check can use:

```sh
--iterations 2 --warmup 0 \
--extra-r-arg '--time-build true --width 16 --height 16 --samples 2 --sphere-count 40 --object-count 40'
```

On ARM, select scalar or the default auto-detected backend; an SSE configuration
must not be relabeled as a successful ARM/SSE measurement. If a local compiler
cache cannot run, use a temporary config with `MAKEFLAGS = "-j8 CCACHE="`.
Do not weaken the checked-in production configuration for smoke tests.

The harness writes `<output>.manifest.json` describing the request and validates
only this invocation's rows. It exits nonzero for failed builds/renders, missing
or duplicate scene/iteration pairs, wrong identity/settings/revision/backend,
or absent/nonfinite/out-of-range required metrics. Build, render, total and
process elapsed time and RSS must be positive; scene and BVH intervals may be
zero at timer resolution. A positive BVH construction count is required.
Linux runs require GNU time elapsed/RSS fields; unsupported platforms keep them
unavailable. Warmup failures also invalidate a run. CSVs and logs are written
before validation fails. GitHub summaries report statuses, medians, units and
sample counts. All validation and aggregation use R (`jsonlite`, `testthat` for
tests); Python is not used.

The worker enables the internal `rayrender.benchmark_timing` option. Native
`steady_clock` intervals cover `BVHAggregate` constructor bodies: primitive
bounds, building, primitive ordering, shadow classification, flattening and
BVH4 conversion. `bvh_build_seconds` sums these intervals on the scene-building
thread; `bvh_build_count` counts the world, mesh and volume-boundary trees.
Nested intervals are counted once. This excludes R scene creation, constructor
argument/member copies before the body, package compilation, ray tracing and
the distinct light-sampling hierarchy. No constructed tree means unavailable,
not zero. Ordinary renders do not acquire timing attributes. `render_seconds`
is the wall time of the public `render_scene()` call, including its R/native
preprocessing and postprocessing; `scene_build_seconds` times the R scene
description when requested; `total_seconds` combines those two intervals.

Configuration revision **2** adds `CXX20FLAGS` and corrected backend selection.
The scalar configuration disables configure's SIMD detection without enabling
`RAYSIMD`; the SSE configuration sets `RAYSIMD` and `HAS_SSE`. Both retain
`-O3 -march=native`, package macros, and package linker settings. Metadata records
the actual standard, compiler version, representative compile command, effective
flags, and backend obtained by preprocessing `simd.h` with that command's flags.
All emitted C++ compile commands are checked for consistent flags, final `-O3`,
and native architecture selection. Historical rows without a revision remain
legacy/unvalidated, separate from corrected builds.

## Collection and history

Data path:

1. `run_render_benchmarks.R` → isolated build → `run_one_render_benchmark.R` →
   scene → numeric worker JSON → configuration CSV.
2. Matrix jobs always upload available CSV/manifest, worker JSON and log artifacts.
3. `combine_benchmark_csvs.R` independently reconstructs the workflow request and
   requires exactly one artifact for each expected configuration. Missing,
   malformed, duplicated or mismatched artifacts create explicit failure records.
   Original artifacts remain available; invalid identities cannot overwrite an
   unrelated historical run.
4. The collector writes `attempt_status=complete` only when the entire comparison
   passes. Otherwise it writes diagnostics and exits nonzero. The publisher
   retains that failure result after publishing diagnostics.
5. `append_benchmark_history.R` merges `data/render_benchmarks.csv`, the RDS copy,
   and `latest.json` on `benchmark-history`. Identity includes run, attempt,
   commit, branch, configuration/revision, scene, settings and iteration. Repeated
   publication is idempotent; rerun attempts are distinct. Established historical
   rows (including failures, legacy duplicates and missing newer fields) remain.
   `latest.json` separates the latest attempt from the latest complete comparison.

The collector requires `--input-dir`, `--output`, `--config`, `--configs`,
`--benchmarks`, `--iterations` and `--warmup`. Pass the same `--extra-r-arg` as
execution. Outside CI, supply `--run-id`, `--run-attempt`, `--commit-sha` and
`--branch-name` from the execution manifest/CSV. Input filenames are
`<configuration>.csv`, including within artifact subdirectories.

```sh
Rscript tools/benchmarks/append_benchmark_history.R \
  --history-dir /tmp/benchmark-history --new-results /tmp/combined-results.csv
Rscript tools/benchmarks/build_benchmark_site.R \
  --history-csv /tmp/benchmark-history/data/render_benchmarks.csv \
  --site-dir /tmp/benchmark-site
```

Open the generated `index.html` directly or serve it over HTTP. The generator
writes R summaries as embedded JSON and `data/dashboard.json`, with a raw CSV
link. Failed and nonfinite samples are excluded independently for each metric;
missing BVH timing never hides available render measurements. The page defaults
to the latest complete comparison on the newest attempt's branch (or the latest
attempt if none is validated) and keeps latest-attempt diagnostics visible.
Filters select branch, run/attempt, scene and configuration. Medians group by
run, attempt, commit, configuration/revision, scene, full settings, compiler and
platform—not iteration timestamps. Each metric has its own sample count.
Package installation time is counted once per build, not per scene row.
Historical points are not connected across incompatible settings. Hosted-runner
hardware/load still adds noise; these results indicate trends, not precise
microbenchmark differences. Output hashes remain diagnostic rather than a gate.

## Deployment ordering and trust

`render-benchmarks.yml` still benchmarks pushes, PRs and manual runs. Only trusted
default-branch push/manual runs publish history. `gh-pages` and
`benchmark-history` pushes are excluded. The history job is serialized and fails
if established history cannot be fetched; it does not initialize over a fetch
failure or erase the branch.

`pkgdown.yaml` preserves push, PR, release and manual documentation builds and
also listens for completion of **Render Benchmarks**. For that event it requires
the same repository, default branch and a push/manual source event. It checks
the exact workflow attempt's successful **Commit benchmark-history branch** step,
so failed measurements can deploy diagnostics after successful publication.
It does not depend on a `GITHUB_TOKEN` push triggering another workflow.

Non-PR site builds check out current trusted default-branch code; no benchmark
artifact is executed. Builds have read-only repository permission. A separate
non-PR job receives deployment permission and only publishes generated site
files. All non-PR documentation deployments share workflow concurrency, and the
final deployment checks both source and history SHAs to reject superseded builds.
Missing/corrupt established history fails generation/deployment instead of
substituting an empty CSV.

## Regression checks

```sh
Rscript -e 'testthat::test_dir("tools/benchmarks/tests", stop_on_failure = TRUE)'
npm ci --prefix tools/benchmarks/tests/browser
npx --prefix tools/benchmarks/tests/browser playwright install chromium
Rscript tools/benchmarks/tests/build_fixtures.R /tmp/benchmark-fixtures
node tools/benchmarks/tests/browser/dashboard.cjs /tmp/benchmark-fixtures
actionlint .github/workflows/render-benchmarks.yml .github/workflows/pkgdown.yaml
```

Playwright is a pinned development-only dependency; the dashboard has no runtime
JavaScript dependencies. The suite tests successful, failed, partial, legacy and
empty datasets; filters, sample counts, failure text, safe JSON embedding and
HTTP load errors; compilation/sign correctness for scalar and supported SSE/NEON
backends; validation, missing artifacts, rerun identity and idempotent history.
`tests/testthat/test-benchmark-timing.R` exercises the native-to-R measurement
path against the installed package, including multiple trees and capture reset.

See [repair verification](REPAIR_VERIFICATION.md) for the investigation evidence
and local checks performed for this repair.
