# rayrender 0.42.3.9000

* Diffuse materials now always use the analytic normal-mapping model of
  Schuessler et al. (2017), including imported diffuse face materials and vertex
  colors. Removed the legacy diffuse bump implementation and the beta
  `normal_mapping` argument; ordinary `diffuse(bump_texture = ...)` uses the
  updated model. Zero-height bump maps preserve smooth normals and UV scale.
* Replaced the approximate Oren–Nayar BRDF with energy-preserving Oren–Nayar
  (EON). `sigma` maps to roughness as `min(sigma / 90, 1)`; zero remains
  Lambertian. EON includes color-dependent multiple scattering. Nonfinite,
  negative, or nonscalar `sigma` values fail validation.
* Fixed the Sampler triangle bump path to use the computed UV derivatives,
  and made vertex-color coordinates consistent between both hit overloads.
* Fixed the RNG triangle consistent-normal heuristic's dependence on ray length
  by matching the Sampler overload's normalized local direction.

## New features

- `render_scene()` Adds experimental `integrator_type = "metal"` for GPU wavefront
  rendering with diffuse, metal, microfacet reflection/transmission, glossy,
  dielectric, translucent, hair, and layered OpenPBR materials. Homogeneous
  `subsurface()` walks and `subsurface_diffusion()` respect priority overlaps.
  It supports image roughness maps on supported macOS GPUs. Analytic spheres,
  ellipsoids, cylinders, disks, and cylindrical/ribbon curves are triangulated
  automatically. Camera movement,
  environment rotation, progressive denoising, and cancellation remain interactive.
  Opaque scenes use specialized pipelines without interior-boundary handling.
  Unsupported scenes warn and use CPU NEE, which remains the default integrator.

- `subsurface()` Supports fully capped cylinders as closed scattering boundaries.

- `render_scene()` Adds `bvh_type = "hlbvh"` for parallel CPU BVH construction
  and `bvh_type = "metal"` for optional Metal GPU construction on macOS. Both
  retain CPU ray traversal and are available to animation and AO renders.
  SAH remains the default; HLBVH trades traversal quality for faster building.
  Metal uses cooperative sorting and parallel singleton-leaf construction, with
  shared node buffers retained through CPU packing to avoid a full arena copy.
  Both HLBVH builders prepare bounds in one pass and pack large SIMD trees in
  parallel; batched CPU Morton generation reduces scheduling overhead.
  Further construction optimizations reduce primitive-reordering and leaf-bound
  overhead. Metal sorting, grouping, and lower-tree construction now share one
  GPU submission, eliminating an intermediate CPU synchronization.
  Large triangle meshes also prepare bounds and shadow classifications in
  parallel. Both HLBVH backends validate compact bounds and compute Morton
  extents in parallel for large inputs. HLBVH packing initializes traversal
  storage on its workers, avoiding a serial initialization pass.

- `vdb_medium()` Loads native OpenVDB and uncompressed, ZIP or BLOSC NanoVDB
  volumes through `openvdbr`'s registered runtime interface. Sparse storage and
  native affine transforms are preserved, with direct grid sampling during
  rendering. `nanovdb_medium()` remains available as an alias, and PBRT imports
  can now render compressed NanoVDB smoke and temperature fields.

- `texture_ptex()` Adds filtered per-face textures through the required `ptexr`
  runtime package, using its installed headers through `LinkingTo` without
  bundling the Ptex implementation in rayrender. `read_pbrt()` preserves Ptex
  face IDs in memory and lazy PLY assets, imports color/roughness and material
  bump maps, and matches PBRT's texture encoding. A shared bounded cache avoids
  loading all texels into R.

- `read_pbrt()` Adds optional `mesh_storage = "ply"` to keep compatible inline
  meshes in file-backed assets instead of retaining their geometry arrays in R.

- `read_pbrt()` Supports serial PBRT `Import`, forward and cross-import object
  references, and exact conversion of uniform cubic B-splines to Bezier curves.

- `read_pbrt()` Reduces memory use and conversion overhead for large instanced
  scenes through shared prototypes, bounded parsed-file caching, and compact
  storage for compatible placements, meshes, and curves. Geometry attributes
  are preserved, and transform values remain numeric throughout parsing.

- `read_pbrt()` Retains illumination from constant-zero-alpha area emitters
  using invisible lights instead of discarding them. Visibility after diffuse
  bounces differs from PBRT and is explicitly diagnosed with `strict = FALSE`.

- `texture_mix()`, `texture_direction_mix()`, `texture_noise()`, `texture_checker()`,
  `texture_gradient()`, and `texture_image()` build composable scalar/color textures
  with independent UV, object, or world mappings. `diffuse()` and `microfacet()`
  accept color graphs; `microfacet()` and `openpbr()` accept roughness graphs, and
  `openpbr()` accepts base-color graphs. Images support explicit linear/sRGB decoding.
- `read_pbrt()` Preserves `directionmix` textures on supported color and roughness
  inputs, including nested image/constant inputs and anisotropic glass roughness.
  PBRT roughness textures are remapped after evaluation instead of replaced by a constant.

- `translucent()` Adds thin-sheet diffuse reflection and transmission with
  independent colors and textures, without a volumetric walk. `read_pbrt()`
  imports diffuse-transmission materials instead of replacing them with opaque
  surfaces.

- `read_pbrt()` Converts legacy PBRT Disney materials to OpenPBR, retaining
  base-color, bump, and roughness textures and mapping metallic, anisotropic
  specular, clearcoat, sheen, transmission, and subsurface controls. Model
  differences and unsupported texture controls are reported with `strict = FALSE`.

- `read_pbrt()` Improves opaque Disney-to-OpenPBR conversion with fitted diffuse,
  partial-metal, tint, sheen, and clearcoat mappings. Reduces excessive fuzz rims
  and coat highlights while preserving GGX widths. Thin, transmissive, and
  subsurface materials retain their physical interface IOR and prior conversion.

- `read_pbrt()` Imports PBRT scene files in R, returning an editable scene and
  camera/render arguments. Supports includes, transforms, instances, common
  geometry and materials, textures, lights, and homogeneous media, with explicit
  diagnostics for approximated or unsupported features.

- `point_light()` and `spot_light()` Add geometry-free lights with inverse-square
  falloff, PBRT v4 spot cones, shadows, and medium attenuation. Attach and manage
  them with `add_light()`, `get_light()`, `list_lights()`, and `remove_light()`.
  PBRT point and spot lights now use these native emitters.

- `read_pbrt()` Maps coated conductors to OpenPBR, imports hair, grid/NanoVDB
  media, realistic lens cameras, and animated object/instance transforms. Maps
  material displacement to bump mapping and PLY displacement through rayvertex,
  with diagnostics for model and tessellation approximations. Animated camera
  imports render a single exposure using both shutter endpoint poses.

- `hair()` Accepts RGB absorption and color inputs; scalar absorption is expanded
  to all three channels. Still images from camera motion now respect the next
  pose when camera motion blur is enabled. Realistic lens files accept spaces
  as well as tabs.

- Camera rays now respect `shutteropen` and `shutterclose` for animated geometry,
  while camera pose interpolation retains normalized shutter fractions. PBRT
  imports sample their full shutter interval.

- `openpbr()` Adds a separate OpenPBR Surface 1.1.1 uber shader with layered
  diffuse, metal, transmission, subsurface, coat, fuzz, thin-film, and emission
  controls. Solid interiors use random walks and dielectric priorities; thin
  surfaces use the thin-walled model. Supports base-color, roughness, and bump
  textures and automatically selects the NEE integrator.

- `subsurface()` Adds opt-in `accelerate = TRUE` to limit geometry queries to
  sampled free flights and simplify well-conditioned collision weights. Both
  random and guided walks retain their physical
  scattering model, embedded objects, and dielectric-priority interfaces.

- `subsurface_diffusion()` Adds fast normalized diffusion for optically thick
  solids, with a dielectric surface and priority-aware glass/liquid contacts.
  Its profile radius controls the scattering spread; use `subsurface()` when
  physical random-walk transport is needed. Importance-sampling surface exits
  reduces bright outliers during repeated glass/liquid interactions. Exit
  directions importance-sample Fresnel transmission for the adjacent material,
  avoiding directions outside the transmitting cone without changing the lobe.

* Add `disk_light()` for a uniform infinite disk with color, intensity, angular
  diameter, and direction controls, without image files or sky datasets.

## Bugfixes

- `render_scene()` Avoids scanning finite triangle emitters when a Metal ray
  escapes the scene, removing a major slowdown with triangulated area lights.
  OpenPBR surfaces retain their layered BSDF instead of being mistaken for
  ordinary dielectrics during GPU export.

- `subsurface_diffusion()` Resolves coincident glass/liquid contacts consistently
  for camera and shadow rays, including containment, refraction, and depth
  accounting. Distinct closed bodies retain both memberships when their first
  containment intersections have identical distances.

- `diffuse()` Preserves smooth mesh shading by applying analytic bump mapping
  relative to the interpolated surface normal. Geometric sidedness remains
  enforced, and bump textures retain the energy-conserving model.

- `render_scene()` Correctly selects equal-count BVH partitioning for
  `bvh_type = "equal"` instead of treating the selector as a leaf-size value.
- Overall: Makes the scalar float-to-integer lane bit cast agree with SSE/NEON
  in builds with SIMD disabled.

- Overall: Preserves shipped configuration scripts during Windows cleanup and
  removes generated configuration for both platforms.

- Overall: Retains generated vignette files in source packages so installation
  can create the vignette index successfully.

- `infinite_light()` Wraps latitude-longitude images across the longitude seam
  using texel centers. Bright sources straddling that seam now retain their
  energy, also when supplied through `render_scene(environment_light = ...)`.
  Image-light sampling accounts for the interpolated footprint of bright texels
  to reduce rare, high-weight samples around small HDR sources.

- Overall: Requires rayimage 0.28.0 for correct HBD display encoding and limits
  renderer threads to two when R package checks request a core limit.

- `read_pbrt()` Preserves camera handedness, converts EXR environment primaries
  to linear sRGB, and corrects the film normalization of blackbody volume
  emission. Constant spectral absorption/scattering coefficients and custom
  NanoVDB density/temperature grid names are now retained. `camera()` and
  `render_scene()` support `camera_flip_x` for matching imported camera frames.

- Overall: Resolves package-check warnings and notes for missing C++ override
  declarations, R utility calls, and example lines truncated in the PDF manual.

- `ply_model()` Reads supplied vertex normals from PLY files.

- Overall: Corrects PLY vertex, normal, and UV element counts, preventing excess
  allocation and phantom vertices in refinement. Mesh construction transfers
  temporary triangle storage into the BVH and releases parser/build buffers
  when they are no longer needed. OBJ consistent-normal arrays use face counts.

- Overall: Fixes transform-cache lookup and growth so collisions do not lose
  existing entries or allocate duplicate transforms. Signed zeros hash
  consistently with matrix equality.

- Overall: Reuses shared instance geometry and its BVH within each render while
  retaining independent placement transforms, emitter identities and medium
  boundary offsets. Large static PBRT populations use packed affine transforms.
  Instances reuse the inverse matrix already stored in their forward transform
  instead of allocating a second transform. Subsurface preparation preserves
  unchanged nested prototypes so ordinary instanced geometry remains shared.

- `dielectric()` Transmits direct illumination through ordinary glass when the
  priority-selected IORs match, preserving colored absorption and surface tint.
  Point, spot, area, and environment lights share this behavior, and matched
  boundaries preserve camera-path depth and light-sampling weights.

- Overall: Corrects bump-map slopes to use UV units and signed texture-repeat
  scaling with PBRT-style filtered height differences and curved-surface normal
  derivatives. Camera ray differentials propagate through transforms, reflection,
  refraction, and priority-skipped interfaces to select bump footprints; diffuse
  and volume scattering clear them. Height maps preserve floating-point precision
  and signed values, including independently scaled PBRT bump maps. Existing bump
  intensities may need adjustment.

- Overall: Mirrored spheres, ellipsoids, cylinders, and disks combine transform
  handedness with explicit normal reversal, correcting refraction through imported
  PBRT glass. Cylinder cap normals are transformed once and texture normals follow
  the same orientation.
- `ellipsoid()` Applies the axis scaling to normals in scalar-vector builds.
  `cylinder()` Correctly intersects cap-only axial rays and rays leaving a cap
  from inside.
- `microfacet()` Corrects rough-glass radiance scaling, Fresnel probabilities,
  visible-normal PDFs, and reflection/transmission evaluation. Index-matched
  transmission passes straight through without an undefined half-vector.

- `texture_noise()` Keeps its noise field continuous across zero and negative
  coordinate boundaries, removing seams from 3D procedural patterns.

- `read_pbrt()` Preserves shape alpha textures and fractional opacity, including
  independent UV scale and offset, so plant and flower meshes retain their cutout
  silhouettes instead of rendering the whole polygon.
  Rays passing through a cutout also continue through the remaining geometry
  inside an object instance.

- `hair()` Corrects directional sampling, PDF evaluation and projected scattering,
  preventing extreme noise and brightness errors. Hair now provides a bounded
  albedo guide to the denoiser and handles exactly grazing directions.

- `bezier_curve()` Caches curve segments and avoids matrix inversion for ray
  intersections. Correctly retains the nearest hit and respects ray intervals.
  Adds `split_depth = 3L`, matching PBRT CPU subdivision; use zero to save memory
  and setup time at the cost of more intersection work.

- `read_pbrt()` Preserves curve subdivision and ribbon-chain endpoint normals.
  Imported scenes start Russian roulette at bounce 5 to avoid tracing negligible
  contributions through high-depth scenes.

- `read_pbrt()` Preserves Film ISO exposure and uses PBRT's smooth default for
  coated diffuse materials, correcting underexposure and overly dull coatings.

- Meshes honor reversed orientation and mirrored transforms in geometric and
  shading normals and area-light selection. PBRT `ReverseOrientation` now
  points triangle emitters toward the intended side, fixing dark city windows.

- Texture-capable materials accept `image_offset = c(0, 0)` for native UV
  translation of their image maps. `read_pbrt()` carries `udelta` and `vdelta`
  into texture lookups, including independent bump and roughness offsets,
  without resampling image pixels. Zero offsets no longer produce diagnostics.

- `read_pbrt()` Uses a PEGTL grammar through `piton` to parse PBRT files,
  retaining large numeric arrays without duplicate character-token buffers.
  Reuses decoded parameter declarations and avoids revalidating meshes already
  checked during conversion, reducing setup costs in large scenes.
  Point, spot, and infinite lights are collected in bounded lists and attached
  once, avoiding repeated validation of growing light lists.
- `read_pbrt()` Loads PFM environment and texture images, preserving their
  floating-point values, scale, byte order, and row orientation.
- Fix texture selection after nested instances, which could attach another
  object's bump map or dereference a null bump buffer during rendering.
- `render_scene()` Avoids generating unused temporary texture paths for
  untextured objects, reducing preparation time and memory in large scenes.
- `grid_medium()` Accepts spatial emission multipliers without expanding them
  into RGB grids. `read_pbrt()` preserves uniform-grid `Lescale` arrays.
- `read_pbrt()` Accumulates large scene imports in bounded lists and combines
  them once, avoiding repeated copying of scene rows and diagnostics. Object
  definitions are combined once and reused by their instances.
- Overall: Avoids repeated list-column copying during scene preparation and
  skips unnecessary subsurface updates for ordinary surfaces. Resolves
  subsurface interiors once per preparation pass while preserving nested
  instances, explicit media, and material replacement.

- `read_pbrt()` Converts Disney subsurface color and diffusion distance together
  using Hyperion's dielectric-aware fit, reducing excessive absorption in the
  OpenPBR conversion. Per-channel extinction mean free paths are preserved;
  the fit's IOR assumption and random-walk depth requirements are diagnosed.
  Specular tint now correctly has no effect at pure dielectric/metallic endpoints.

- `read_pbrt()` Accepts PBRT execution options such as `Option "wavefront" true`,
  recording unsupported options in the conversion diagnostics instead of failing
  to parse otherwise renderable scenes. Parameter lists can continue after an
  `Include`, including nested and repeated includes.

- `read_pbrt()` Imports one-pixel environment maps without crashing the EXR writer.

- `read_pbrt()` Accepts identical repeated parameters. With `strict = FALSE`,
  unresolved material names use a diagnosed diffuse fallback and repeated
  texture names are redefined with a scope warning, allowing more previews to
  finish. Strict imports still reject unresolved or repeated resource names.

- Image textures: JPEG, HDR, and EXR files no longer fail during PNG-only alpha
  inspection before rendering.

- `subsurface_diffusion()` Keeps interactive rendering responsive when the camera
  enters an active diffusion region. Affected rays terminate with recorded
  diagnostics and a console warning instead of aborting the render.

- Overall: Denoising follows perfect glass and mirror paths to the first
  non-specular interaction for both albedo and normal guides. Nested glass no
  longer forces white guides after the first interface. Volume and subsurface
  renders now use guides from their scattering paths, with auxiliary prefiltering
  for final images.

- Overall: Parallel rays no longer produce invalid rectangle intersections.
  This prevents runaway memory use when diffusion probes lie in a box face;
  invalid ordered medium crossings terminate the affected path with diagnostics.

* Recover from classified medium and subsurface transport failures by terminating
  only the affected ray path instead of aborting the render. Retain accumulated
  light and keep the sample in the image average without retrying it. Scene
  validation errors and unclassified exceptions still stop rendering.
* Keep transport diagnostics silent by default, attached to the returned image
  as `path_warnings`, without console warnings or log files. Set
  `RAYRENDER_DEBUG_PATHS=true` (or `1`) to enable a summarized warning after image
  output and a diagnostic log, whose filename is attached as `path_warning_log`.
  Stills and animation frames retain failure counts and up to eight detailed
  examples per category. Collection runs only on failed paths, with no diagnostic
  locking, counting, or formatting on successful paths.
* Fix medium membership at glass/liquid contact edges and corners, and correct
  grazing subsurface collisions that round onto a boundary plane. These changes
  prevent spurious repeated-entry failures while preserving sampled flight
  distances and scattering weights. Added precision and path-recovery regressions.
* Fix a sign error in `microfacet(transmission = TRUE)` that produced negative
  transmitted-light contributions and could make rough transparent surfaces
  render black. Added a transmitted-radiance regression test.
* Add reproducible scenes in `tools/subsurface/` for milk and chocolate milk,
  iced coffee with dielectric-priority ice, droplet and microfacet condensation,
  and Monterey Bay terrain styled as a cookie in a bowl of milk.

* Mesh construction for `extruded_path()` and `extruded_polygon()` now lives in
  rayvertex's `extruded_path_mesh()` and `extruded_polygon_mesh()`. Rayrender's
  public wrappers preserve their arguments, materials, transforms, and SSS
  boundary identities. Geometry tests moved with the builders; rayrender retains
  integration and rendering tests. Requires rayvertex >= 0.16.0.
  
* Fixed `extruded_polygon()` cap winding in every plane, reflected scales and
  reversed heights, multipart holes, direct multiple-hole indices, plain
  `SpatialPolygons`, per-feature heights, positive x offsets, and uppercase plane
  names. Shared cap/wall vertices and local-coordinate winding calculations
  produce consistent closed meshes, including at large coordinate offsets.
  Redundant vertices are removed; invalid rings, holes, heights, and scales now
  produce explicit errors. Added geometry tests, SSS rendering regressions, and
  before/after examples in `tools/polygon/`.
* Preserve precise triangle distances when ordering medium intersections. Entry
  and exit hits at grazing mesh edges no longer swap when their single-precision
  distances round to the same value, preventing rare SSS boundary-state errors.
  Containment probes also stay on their original ray instead of offsetting past
  neighbouring faces at pointed mesh corners.

* Fixed `extruded_path()` closed seams, exact trim endpoints, terminal twist and
  morphing, taper normals, polygon winding, wrapped cap choices, and shared
  material identities. Zero-width tips now use nondegenerate triangle fans.
* Extracted an internal indexed sweep builder with adaptive arc-length sampling,
  added `initial_normal`, `smooth_angle`, and `arc_tolerance`, and added geometry
  and rendering regressions. Invalid widths, degenerate profiles, and incompatible
  closed seams now produce explicit errors.

* Added `subsurface()` for homogeneous RGB scattering inside a neutral dielectric
  boundary, with artist reflectance/extinction-distance controls or direct physical
  absorption and scattering coefficients. Interiors are prepared automatically,
  including inside instances, with explicit-medium conflict diagnostics.
* Added `subsurface(priority = ...)` for overlaps with glass and other SSS bodies.
  Lower values win; hidden interfaces pass through without extra reflection,
  refraction, or scattering. Glass/liquid contact can be modeled by extending the
  lower-priority liquid into the glass wall instead of leaving an air gap.
* Added separate internal SSS depth and compensated roulette, plus selectable
  ordinary and Dwivedi direction/distance proposals with joint probability and
  direct-light MIS correction. `random_walk` remains the default pending broader
  production benchmarks. The current guide is restricted to isotropic scattering;
  other anisotropy values retain their physical HG phase in an unguided fallback.
  Rough boundaries use a single-scattering GGX model.
* Added deterministic proposal/boundary tests, independent homogeneous slab
  references, and an opt-in raw-radiance benchmark in `tools/subsurface/validate.R`.

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

## Bugfixes

- `render_scene()` Opaque renders containing media or subsurface materials no
  longer trace unused transparency or wait for the 64-sample alpha convergence
  minimum. Transparent renders retain their coverage estimator and minimum.

## Documentation

- Overall: Adds a distinct rendered example directly to each new function's
  roxygen documentation, including texture composition, OpenPBR, diffusion,
  translucent materials, light management, and PBRT import. Examples use 16
  samples with denoising, directional lighting, and before/after renders for
  light editing and removal. Texture details explain each operation's purpose,
  practical uses, and the distinction between scaling values and coordinates.

- Overall: Adds a composable-texture vignette with 26 rendered examples,
  reproducible recipes, coordinate and image-mapping comparisons, and current
  material-input and filtering limitations.

## Other

- Overall: Speeds up scene-row assembly, including large PBRT imports, while
  preserving descriptor data, names, metadata, and type compatibility checks.

- Overall: Obtains GLM and the OpenPBR BSDF headers from the separate
  `glmheaders` and `openpbr` packages through `LinkingTo`, replacing bundled
  header copies while retaining third-party attribution.

- `render_scene()` Speeds up triangle coverage tests in medium and subsurface
  paths by reducing conditional branches, while preserving double-precision
  intersections and inclusive boundary handling.
- `render_scene()` Speeds up medium and subsurface membership checks by pruning
  distant closed objects and classifying individual boundaries from their first
  crossing. Retains full crossing replay for instances and exact contacts.
- `render_scene()` Distributes expensive regions across smaller rendering jobs,
  improving core utilization in scenes with localized glass or subsurface
  scattering. Preserves pixel sampling and adaptive convergence regions.
