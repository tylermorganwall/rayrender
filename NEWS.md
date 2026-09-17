# rayrender 0.42.0

* Allow `sky_light()` and `sky_light_image()` to position the Sun with `elevation`
  and `azimuth` instead of latitude, longitude, and date/time. Add rendered light
  management, standalone Sun/Moon, and medium examples, and install Prague sky
  datasets before building the pkgdown site.

* Move EXR white-balance baking controls from the render functions to
  `sky_light_image()`. Cache each adapted sky for use in stills and animations,
  reusing the generated sky when only the target white point changes.

* Fix Windows `near`/`far` macro collisions in geometry and volume code. Remove
  unused preview state and signedness warnings in BVH construction and spline
  interpolation, and make haze-cut sorting bounds explicit. Correct RGB-to-HSV
  conversion to include blue when finding the maximum channel and initialize
  hue on every path.

* Automatically select the NEE integrator for scenes containing native atmospheric
  skies or attached media, including clouds and instanced media, instead of
  requiring an explicit `integrator_type = "nee"` selection. Recognize emitting
  media when deciding whether to add automatic ambient illumination.

* Fix Windows macro collisions in atmospheric lighting and preview test builds
  using X11. Emit the transform error-bound overloads needed by volume boundaries
  explicitly so optimized GCC builds load reliably.

* Select finite lights with a spatial light BVH in the NEE integrator, using
  distance, estimated power and emission direction. Retain the previous
  scattering context for matching MIS probabilities and mix in the original
  distribution to preserve support for textured and animated emitters.

* Accelerate eligible opaque NEE shadow connections with early-exit BVH traversal
  and reduced triangle intersection work. Select explicit emitters with matching
  MIS weights, avoiding full light-mixture PDF scans and supporting multiple and
  instanced lights. Preserve ordered visibility for media and alpha masks, and
  retain the directional-mixture estimator for atmospheric transport.

* Store BVH4 leaf primitive ranges in compact records, reducing BVH memory use.

* Reduce BVH traversal bookkeeping while preserving intersection order, and
  defer unused medium transforms in surface hit records.

* Split the interactive preview status bar into camera information and rendering
  flags on separate rows.

* Reuse rendering workers across samples and wake on task completion, removing
  fixed polling delays from beauty and denoising feature passes while retaining
  preview cancellation and R interrupt handling.

* Reduce the default Moon disk resolution to 256 pixels in `moon_light()` and
  automatic sky lights, retaining the Moon renderer's 2x antialiasing.

* Skip grid emission lookups for media with zero RGB emission and no temperature
  field, or with emission disabled by a zero scale.

* Keep interpolated camera-up vectors at keyframes to avoid brief roll wobbles
  when damping camera motion. Keyed positions and lookat targets are retained.

* Allow `sky_light()` below a Sun elevation of -4.2 degrees. Prague contributes
  black sky and no solar haze, while enabled Moon, stars, and planets retain
  their light, atmospheric transmission, and Earth occlusion.

* Animate `cloud()` by translating sampling coordinates through fixed Perlin
  fields with `t`; `animation_seed` selects the direction. Broad and fine noise
  move together while the envelope and boundary stay fixed. This changes the
  seeded density fields from the previous simplex-perturbation implementation.
  Whole-cloud placement remains controlled by x/y/z.

* Expose sky model settings directly in `sky_light()` and `sky_light_image()`.
  Both accept `altitude`, `visibility`, and other named settings instead of a
  `sky_args` list. Image resolution and Hosek settings belong to
  `sky_light_image()`; native atmospheric sampling uses `sampling_resolution`.

* Add `cloud()`, a procedural cumulus/stratus volume centered at the origin, with
  standard position, rotation, and scale controls. The density field and its
  boundary transform together. Sky examples now use the exported constructor.
* Correct the Z rotation matrix used by R-side group and animation transforms
  to match the positive-angle convention of direct objects and instances.

* Enable altitude-dependent lighting, finite haze, and deferred haze sampling by
  default in `sky_light()`. Expose `haze` and `query_altitude` controls, and reduce
  horizon bands with `haze_filter`, sampling one nearby atmospheric path per haze
  query. Reuse exact spectral queries and precompute transmission reconstruction
  with `cache_spectra`, `transmission_table`, and `transmission_table_max_mb`.
* Exclude haze from volume interiors by default with `haze_in_volumes = FALSE`.
  Individual media also support `haze` and `haze_density_threshold`; clouds
  default to haze off, with a density cutoff of 0.05 when enabled.
* Add capsule terrain, trees, and cloud examples comparing altitude-dependent
  lighting and haze, with additional atmospheric controls described in the help.

* Fix Moon disk and legacy sky preparation in skymodelr so the upper limb
  remains visible after the center sets. Normalize the complete disk before
  horizon clipping, with continuous horizon attenuation and color.

* Add `sun_light()` and `moon_light()` scene elements using skymodelr positions,
  apparent sizes, Prague solar radiance, and detailed phase-dependent lunar
  textures. Dedicated disk sampling preserves detail independently of sky-map
  resolution, with consistent solid-angle PDFs and light-selection weights.
  Both lights use the public skymodelr disk generators for radiance and geometry.
* Add image-based `infinite_light()` scene elements, named light management,
  and additive sampling of multiple environments. Existing environment render
  arguments remain supported.
* Add `sky_light()` for native Prague location/time skies, with automatically
  sampled Sun and Moon disks. Add `sky_light_image()` for cached skymodelr EXRs
  shared across still renders, animation, and preview.
* Make animation environment rotation use the same direction as still renders,
  and isolate its preview environment transform from geometry transforms.
* Preserve preview RGB colors when exporting snapshots through rayimage.
* `nee` now traces surfaces and RGB participating media with null scattering and
  path-space MIS. Its images change; `rtiow` remains the default and legacy fog
  remains available in `rtiow` and `basic`.
* Add reusable homogeneous, dense-grid and uncompressed float NanoVDB medium
  descriptions, independent medium attachments, nested containment, anisotropic
  scattering, and RGB or normalized blackbody emission.
* Add volume transmittance alpha, alpha-aware adaptive sampling, color-only
  denoising for participating media, and consistent RGBA animation and snapshots.
* Correct NEE light PDFs, sampler dimension reuse, segment-based dielectric
  absorption, deterministic volume boundary intersections, and realistic-camera
  shutter-time preservation.
