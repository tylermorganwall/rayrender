# Render Benchmark Harness

This directory contains an R-only benchmark harness for measuring rayrender
render time across source refs, build settings, and R-defined scenes.

The orchestrator copies or checks out a source tree, installs rayrender into an
isolated temporary R library with `R CMD INSTALL -l`, then calls the R
worker once per warmup or measured iteration. Measured iterations are appended
to a CSV file. Each data row is one render iteration.

## Current Working Tree

Run the default smoke benchmark against the current checkout, including
uncommitted changes:

```sh
Rscript tools/benchmarks/run_render_benchmarks.R \
  --repo . \
  --ref current \
  --config tools/benchmarks/configs/default.json \
  --benchmarks bvh_many_spheres \
  --iterations 1 \
  --warmup 0 \
  --output /tmp/rayrender_benchmark_smoke.csv
```

Run the scalar optimized config and record scene construction timing when the
scene supports it:

```sh
Rscript tools/benchmarks/run_render_benchmarks.R \
  --repo . \
  --ref current \
  --config tools/benchmarks/configs/default.json \
  --only-config o3_native_no_simd \
  --benchmarks bvh_many_spheres \
  --iterations 1 \
  --warmup 0 \
  --extra-r-arg "--time-build true" \
  --output /tmp/rayrender_benchmark_smoke.csv
```

Run both bundled scenes with more iterations:

```sh
Rscript tools/benchmarks/run_render_benchmarks.R \
  --repo . \
  --ref current \
  --config tools/benchmarks/configs/default.json \
  --benchmarks bvh_many_spheres,bvh_mixed_primitives \
  --iterations 5 \
  --warmup 1 \
  --output tools/benchmarks/results/render_benchmark_times.csv
```

## Git Refs

Use any commit, tag, or branch in `--ref`:

```sh
Rscript tools/benchmarks/run_render_benchmarks.R \
  --repo . \
  --ref HEAD~1 \
  --config tools/benchmarks/configs/default.json \
  --benchmarks bvh_many_spheres \
  --iterations 5 \
  --warmup 1 \
  --output tools/benchmarks/results/render_benchmark_times.csv
```

For non-current refs, the harness exports the requested tree with `git archive`
into a temporary source directory. Benchmark scene files are read from the
current repository, so older commits do not need to contain this benchmark
harness.

## Release Commits

Run the bundled benchmark across every release commit in the current
first-parent history whose commit message contains a `v0.x.y` version:

```sh
Rscript tools/benchmarks/run_release_render_benchmarks.R \
  --repo . \
  --config tools/benchmarks/configs/default.json \
  --benchmarks bvh_many_spheres,bvh_mixed_primitives \
  --iterations 5 \
  --warmup 1 \
  --output tools/benchmarks/results/release_render_benchmark_times.csv
```

The release wrapper writes the benchmark rows to `--output` and writes a
sidecar commit manifest next to it with a `_commits.csv` suffix. Use
`--max-releases N --dry-run` for a quick plan check before starting a long run,
or `--all-refs` to search branches, tags, and remotes instead of just the
current first-parent history. By default, the wrapper uses `--cxx-std auto`:
it reads each release commit's `DESCRIPTION` and sets `CXX_STD` to the declared
C++20 when required and otherwise uses `CXX17` for older releases so current R
toolchains and Rcpp headers do not compile them as pre-C++11.

## Build Configurations

Build configs live in JSON under `tools/benchmarks/configs`. The default config
contains one baseline build and smoke-sized render settings so the smoke test
writes one data row quickly. Increase `width`, `height`, `samples`,
`sphere_count`, or `object_count` when collecting comparison data.

Additional configs can be added to `build_configs`:

```json
{
  "name": "o3_native_simd",
  "env": {},
  "makevars": {
    "CXX17FLAGS": "-O3 -march=native -DRAYSIMD",
    "CXX14FLAGS": "-O3 -march=native -DRAYSIMD",
    "CXXFLAGS": "-O3 -march=native -DRAYSIMD"
  }
}
```

Supported Makevars keys are `CXX_STD`, `CXX`, `CXX11`, `CXX14`, `CXX17`,
`CXX20`, `CC`, `CFLAGS`, `CXXFLAGS`, `CXX11FLAGS`, `CXX14FLAGS`,
`CXX17FLAGS`, `CXX20FLAGS`, `PKG_CXXFLAGS`, `PKG_LIBS`, and `MAKEFLAGS`.

The `o3_native_no_simd` config sets `RAYRENDER_DISABLE_SIMD=true`, which tells
the package configure script to skip auto-detected SSE/NEON flags. The
`o3_native_simd_sse` config defines the four-wide SIMD path with `-DRAYSIMD`
and `-DHAS_SSE`.

Filter configs with:

```sh
Rscript tools/benchmarks/run_render_benchmarks.R ... --only-config default
```

## Adding Scenes

Add a file at `tools/benchmarks/scenes/<name>.R`. It must define:

```r
run_benchmark = function(settings) {
  scene = rayrender::generate_ground()
  elapsed = system.time({
    image = rayrender::render_scene(
      scene,
      width = settings$width,
      height = settings$height,
      samples = settings$samples,
      preview = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      denoise = FALSE
    )
  })[["elapsed"]]
  attr(image, "render_seconds") = as.numeric(elapsed)
  image
}
```

Use `settings$time_build` if the scene should support timing both scene
construction and rendering. The bundled scenes build scene descriptions outside
the timed render by default, but `render_scene()` still includes user-facing
preprocessing and BVH construction.

Scenes may return either the rendered object directly or a list/data frame with
optional metric fields. Supported optional fields are `render_seconds`,
`bvh_build_seconds`, `scene_build_seconds`, `total_seconds`, `output_hash`, and
`artifact_path`. If BVH construction cannot be separated through the public R
API, set `bvh_build_seconds` to `NA` and report scene setup or total timing
instead.

## Useful Options

- `--keep-workdirs`: keep temporary source, logs, artifacts, and R library.
- `--workdir-root`: choose where temporary workdirs are created.
- `--timeout-seconds`: stop an individual render after a timeout.
- `--extra-r-arg`: pass an extra option to the R worker, for example
  `--extra-r-arg "--time-build true"`.
- `--no-append`: overwrite the output CSV.
- `--dry-run`: print the planned matrix and build commands.

## CSV Schema

The CSV includes:

- Run identity: `timestamp_utc`, `run_id`, `run_attempt`, `workflow`,
  `event_name`, `source_ref`, `commit_label`, `commit_sha`, `git_dirty`,
  `branch_name`, `pull_request_number`, `pull_request_head_ref`,
  `pull_request_base_ref`.
- Benchmark identity: `benchmark_name`, `benchmark_source`,
  `benchmark_settings_json`, `iteration_index`, `iterations`,
  `warmup_iterations`.
- Render settings: `seed`, `width`, `height`, `samples`, `threads`.
- Build metadata: `build_config_name`, `compiler_id`, `cc`, `cxx`,
  `cxxflags`, `cxx14flags`, `cxx17flags`, `pkg_cxxflags`, `makeflags`,
  `extra_env_json`, `simd_enabled`, `raysimd_defined`, `has_sse_defined`,
  `has_avx_defined`.
- Platform metadata: `runner_os`, `platform_system`, `platform_release`,
  `platform_machine`, `platform_processor`, `cpu_count_logical`,
  `python_version`, `r_version`. The compatibility field `python_version` is
  recorded as `NA` by this R harness.
- Package metadata: `package_name`, `package_version`.
- Timing and memory: `package_build_seconds`, `package_install_seconds`,
  `scene_build_seconds`, `bvh_build_seconds`, `render_seconds`,
  `total_seconds`, `process_elapsed_seconds`, `process_user_seconds`,
  `process_system_seconds`, `max_rss_kb`, `max_rss_mb`. Compatibility aliases
  `build_seconds` and `install_seconds` are also retained.
- Status and artifacts: `status`, `error`, `output_hash`, `artifact_path`,
  `build_log_path`, `benchmark_log_path`, `time_log_path`, `workdir_path`.

Rows with `status = build_failed` represent an intended benchmark that could not
run because the package build failed. Rows with `status = render_failed`
represent failed measured render iterations.

On Linux, each measured worker process is wrapped with `/usr/bin/time -v`.
`render_seconds` still comes from the R worker; `/usr/bin/time` supplies max RSS
and optional outer process timing. If GNU time is unavailable, those fields are
recorded as `NA`.

## GitHub Actions Dashboard

The workflow `.github/workflows/render-benchmarks.yml` runs on:

- `workflow_dispatch`;
- pushes to `main` and `deferred`;
- pull requests.

Pull requests run the benchmark matrix and upload CSV/log artifacts, but do not
push history or deploy Pages. Pushes to `main`/`deferred` and manual runs with
`publish=true` append successful artifacts to the `benchmark-history` branch and
deploy a static site through GitHub Actions Pages.

To run manually:

1. Open Actions -> Render Benchmarks.
2. Choose benchmark names, iterations, warmups, and whether to publish.
3. Run the workflow.

Set GitHub Pages to use GitHub Actions:

1. Open repository Settings -> Pages.
2. Set Source to `GitHub Actions`.
3. Run the workflow with `publish=true`.

Persistent files live on the `benchmark-history` branch:

- `data/render_benchmarks.csv`;
- `data/render_benchmarks.rds`;
- `data/latest.json`;
- `site/index.html`;
- `site/data/render_benchmarks.csv`;
- `site/assets/dashboard.js`;
- `site/assets/style.css`.

The generated dashboard is static HTML, CSS, and JavaScript. It shows latest
median rows, render time trends, BVH build time trends when available, max RSS,
and package build time. Filters are available for branch, benchmark, and build
configuration. The raw CSV is linked from the page.

You can test the history append and site generation locally:

```sh
mkdir -p /tmp/benchmark-history-test
Rscript tools/benchmarks/append_benchmark_history.R \
  --history-dir /tmp/benchmark-history-test \
  --new-results /tmp/rayrender_benchmark_smoke.csv

Rscript tools/benchmarks/build_benchmark_site.R \
  --history-csv /tmp/benchmark-history-test/data/render_benchmarks.csv \
  --site-dir /tmp/benchmark-history-test/site
```

Open `/tmp/benchmark-history-test/site/index.html` in a browser to inspect the
static page.

## Interpreting Results

Prefer medians over single timings. Use several iterations and at least one
warmup when comparing commits. Benchmarking renderers is sensitive to thermal
throttling, CPU frequency scaling, battery state, background processes, and
other users of the same machine. Keep the machine plugged in, close unrelated
CPU-heavy work, and compare refs in repeated alternating runs when changes are
small.

The worker records an output hash for the returned render object, but hash
changes do not fail the benchmark. Use them as a signal that a performance
comparison may also include a rendering behavior change.

GitHub-hosted runners are noisy. Treat the dashboard as trend detection, not
precise microbenchmarking. Compare medians across multiple iterations and avoid
over-interpreting small changes. Memory is Linux max RSS for the worker process,
so it is most comparable within the same workflow environment. Pull requests
from forks should not push benchmark history.
