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

Supported Makevars keys are `CXX`, `CXX11`, `CXX14`, `CXX17`, `CC`, `CFLAGS`,
`CXXFLAGS`, `CXX11FLAGS`, `CXX14FLAGS`, `CXX17FLAGS`, `PKG_CXXFLAGS`,
`PKG_LIBS`, and `MAKEFLAGS`.

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

- Benchmark identity: `benchmark_name`, `benchmark_source`,
  `benchmark_settings_json`, `iteration_index`, `iterations`,
  `warmup_iterations`.
- Render settings: `seed`, `width`, `height`, `samples`, `threads`.
- Platform metadata: `platform_system`, `platform_release`,
  `platform_machine`, `platform_processor`, `cpu_count_logical`,
  `python_version`, `r_version`. The compatibility field `python_version` is
  recorded as `NA` by this R harness.
- Package and git metadata: `package_name`, `package_version`, `source_ref`,
  `commit_label`, `commit_sha`, `git_dirty`.
- Build metadata: `build_config_name`, `compiler_id`, `cc`, `cxx`,
  `cxxflags`, `cxx17flags`, `pkg_cxxflags`, `makeflags`, `extra_env_json`,
  `build_seconds`, `install_seconds`, `build_log_path`.
- Render result: `render_seconds`, `status`, `error`, `output_hash`,
  `artifact_path`, `benchmark_log_path`, `workdir_path`.

Rows with `status = build_failed` represent an intended benchmark that could not
run because the package build failed. Rows with `status = render_failed`
represent failed measured render iterations.

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
