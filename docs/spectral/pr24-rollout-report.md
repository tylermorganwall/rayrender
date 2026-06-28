# PR24 Rollout Report

PR24 defines the user-facing spectral rollout boundary. It does not make the
spectral renderer the default and does not connect schema-v2 compiler input to a
runtime pixel renderer.

## Default Policy

The default remains:

```r
render_mode = "rgb_legacy"
```

Spectral mode remains opt-in. The current public contract is:

```r
render_scene(scene, render_mode = "spectral", return_result = TRUE)
render_animation(scene, camera_motion, render_mode = "spectral", return_result = TRUE)
spectral_render_capabilities()
```

Calling still or animation spectral mode with `return_result = FALSE` fails
clearly because the direct spectral renderer bridge is not connected.

The spectral renderer can become the default only after the PR24 acceptance
gates have rendered data, package checks pass on supported platforms, one stable
release exposes spectral mode as opt-in, and maintainers approve the appearance
change policy. `render_mode = "rgb_legacy"` remains the escape hatch and legacy
renderer code must not be removed in the default-switch PR.

## Still and Animation

Still and animation now share the R-side spectral compiler-input contract:

- scene schema/compiler;
- spectrum and image cache descriptors;
- RGB-to-spectrum table policy;
- Film/PixelSensor descriptors;
- Material/Light/Medium registries;
- integrator descriptor;
- output color conversion descriptor;
- region initialization metadata.

`render_animation(..., render_mode = "spectral", return_result = TRUE)` returns
frame-zero compiler input using the first requested camera-motion row. Runtime
animation is still disabled until the direct C++ spectral bridge exists.

## AOVs, Adaptive Sampling, and Denoising

PR24 records these contracts in `spectral_render_capabilities()` and in returned
compiler input under `render_defaults$spectral_policy`.

- Adaptive decisions must use Film/sensor linear channels, never wavelength
  packet components.
- Beauty and albedo denoiser inputs are output-linear RGB from the same
  Film/PixelSensor basis. Normals are geometric vectors.
- Spectral AOVs are defined for beauty, alpha, albedo, normal, depth, emission,
  and variance, but runtime emission of these AOVs remains disabled.
- Alpha is a stochastic cutout policy before Material closure and region
  transition. A transmissive dielectric is not alpha.

## Documentation and Migration Examples

New opt-in spectral scene construction:

```r
scene = ray_scene_v2(
  sphere(material = diffuse(reflectance = spectrum_rgb("red", role = "albedo")))
) |>
  add_light(point_light(c(0, 4, -3), power = spectrum_rgb("white", role = "illuminant")))

compiler_input = render_scene(
  scene,
  render_mode = "spectral",
  return_result = TRUE,
  integrator = path_integrator(max_depth = 8),
  sampler = sobol_sampler(pixel_samples = 16),
  film = rgb_film(width = 64, height = 64)
)
```

Legacy scene migration remains explicit:

```r
legacy_scene = sphere(material = glossy(color = "#884422"))
spectral_scene = legacy_scene_to_schema_v2(legacy_scene)
```

Legacy materials with ambiguous physical meaning report conversion warnings once
per scene. The conversion report is preserved on the schema-v2 scene.

## Intentional Deviations and Missing Runtime Features

The pbrt conformance matrix records per-subsystem deviations. PR24 keeps these
runtime features disabled:

- spectral still rendering;
- spectral animation rendering;
- preview/progress/cancellation through spectral Film;
- adaptive sampling through spectral Film;
- denoising and AOV file output;
- rendered pbrt reference-image comparison;
- default switch to spectral.

These are blocked on the direct schema-v2 compiler-to-C++ render bridge, not on
R descriptor serialization.

## Gate Evidence

Fast PR24 evidence:

- `air format R/spectral_rollout.R R/render_scene.R R/render_animation.R R/schema_v2_descriptors.R tests/testthat/test-schema-v2-descriptors.R`
- `Rscript -e "devtools::load_all('.', quiet=TRUE); testthat::test_file('tests/testthat/test-schema-v2-descriptors.R')"`: 158 passes.
- `R CMD INSTALL .`
- `R CMD build --no-build-vignettes --no-manual /Users/tyler/Desktop/R/rayrender`
- `R CMD check --no-manual --no-vignettes /private/tmp/rayrender_0.41.3.tar.gz`: Status OK; installed size 17.4 MiB with 9.5 MiB `extdata` and 6.6 MiB `libs`.

The package check no longer reports the Rd warnings present before PR24.
