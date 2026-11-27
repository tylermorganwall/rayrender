# rayrender 0.43.0

rayrender can now path trace on supported Mac GPUs, animate rigged glTF models,
import PBRT scenes, and save editable scenes in OpenUSD. New materials include
OpenPBR, thin translucent sheets, and two ways to render subsurface scattering.

## New features

### Rendering on the GPU and CPU

- `render_scene()` and `render_animation()` now default to
  `integrator_type = "auto"`. They use the new Metal wavefront path tracer on a
  ray-tracing-capable macOS GPU and CPU next-event estimation (NEE) elsewhere.
  Unsupported scene features or sampling modes trigger a warning and a CPU
  render of the whole scene. Set `integrator_type = "nee"` to choose the CPU
  explicitly. The preview identifies the active renderer, and the returned
  image's `wavefront` attribute records fallback and timing information.
- Metal renders instanced meshes, analytic spheres and ellipsoids, Bézier
  curves, constructive solid geometry (CSG), and layered materials. It supports
  open water surfaces, transparent backgrounds, Prague skies and haze, and
  homogeneous, dense, and sparse volumes with scattering and thermal emission.
  Perspective, orthographic, panoramic, and physical lens cameras support
  motion blur and interactive preview.
- `bvh_type = "hlbvh"` builds the bounding volume hierarchy in parallel on
  the CPU. On supported Macs, `bvh_type = "metal"` builds it on the GPU.
  Both options retain CPU ray traversal and can speed up setup for large scenes
  in still, animation, and ambient occlusion renders. The default builder for
  CPU tracing remains `"sah"`.

### Animation and cameras

- `animation_keyframes()` animates position, rotation, and scale on independent
  schedules built with `animation_keyframe_track()` and
  `animation_transform_track()`. Tracks support stepped, linear, and cubic
  interpolation, including quaternion rotations.
  `animation_clip_player()` can hold a pose, repeat a
  clip, play it backward, or ping-pong. In-place playback removes selected
  root motion, so a walk cycle can follow a path you define.
  `animation_loop_seam()` reports jumps between a loop's endpoints.
- `animation_read_gltf()` (also `read_gltf()`) imports local glTF and GLB
  triangle models with rigid node animation, skeletal skinning, and morph
  targets. Supply a rayrender material when importing. Set
  `recalculate_normals = TRUE` to recompute smooth normals on deformed poses.
- `animation_model_asset()` combines a hierarchy of parts authored with
  `animation_model_node()` and named clips built with `animation_model_clip()`.
  Multiple `animation_model_instance()` placements share geometry and can
  each play a different clip or start at a different time.
  `animation_list_clips()` lists clip IDs and targets, and
  `animation_model_pose()` inspects joint transforms at an exact time.
- `render_scene(time = ...)` renders a still at a specified scene time in
  seconds. `render_animation(times = ...)` accepts times in any order,
  including repeated times. For explicit exposures, the default
  `animation_motion_blur = "full"` samples moving parts, skin, and morphs
  throughout the shutter interval on both CPU and Metal. Use `"placement"`
  to freeze the internal pose while blurring the model's movement through
  the scene.
- `render_scene(timeline = ...)` adds scrubbing, frame stepping, reverse
  playback, and speed controls. `camera(keyframes = ...)` retains an editable
  camera path. The preview follows camera and projection changes at the
  authored times, including physical lens cameras with fixed optics during
  an exposure. Interactive navigation scales to the visible scene. Lowering
  the field of view to zero switches to orthographic projection.
- `render_animation()` draws screen-space text and lines with the camera for
  each frame, preserving pixel size and occlusion in previews and files.
  File output now overlaps rendering with post-processing and writing one
  frame at a time. Use `async_write = FALSE` for synchronous writes or
  `return_frames = TRUE` to retain the image arrays as well as the filenames.

### Scene import and export

- `read_pbrt()` imports PBRT v4 scenes into an editable R scene with camera
  and render settings. It handles includes, serial imports, shared instances,
  animated transforms, realistic cameras, meshes, curves and hair, lights,
  textures, and volumes. `strict = FALSE` returns diagnostics for approximated
  or unsupported features, including conversions from Disney and coated
  materials to OpenPBR. For large scenes, `mesh_storage = "ply"` keeps
  compatible inline meshes on disk. Choose a persistent `asset_dir` when
  saving an import for a later R session.
- `write_usd()` exports geometry, cameras, instancers, materials, animation,
  and volume references to binary or text OpenUSD stages. By default it
  collects referenced files beside the stage. `read_usd()` restores stages
  written by rayrender, including renderer-specific settings stored as typed
  USD metadata. Both functions use the optional `libopenusd` package
  (>= 26.8.0-0), which can be installed after rayrender. Export also uses
  `digest`.
- `set_render_settings()` stores render defaults with a scene and carries
  them through a USD round trip. Arguments supplied to a render call take
  precedence.

### Materials, textures, geometry, and lights

- `openpbr()` implements OpenPBR Surface 1.1.1, with layered diffuse and
  metallic reflection, transmission, subsurface scattering, coat, fuzz,
  thin-film interference, and emission. `translucent()` models diffuse
  reflection and transmission through a thin sheet with independently
  textured colors on the two lobes.
- `subsurface()` traces scattering inside closed solids. Specify reflectance
  and an extinction distance, or supply absorption and scattering coefficients
  directly. Material priorities handle contacts such as ice in a drink.
  Random and guided walks support optional `accelerate = TRUE`.
  `subsurface_diffusion()` offers a faster approximation for optically thick
  objects, with a radius that sets how far light spreads.
- Compose textures from images, noise, checkerboards, gradients, and
  constants. Mix by weight or surface direction, extract a color channel, and
  scale values independently of UV, object, or world coordinates. Use the
  result for color or roughness on supported materials. `texture_image()`
  accepts an explicit color encoding, and materials with image maps accept
  `image_offset` for UV translation.
  Color inputs also accept `"rgb(r, g, b)"` strings with channels from 0 to 255.
- `texture_ptex()` reads per-face textures with filtering through `ptexr`.
  PBRT imports retain the original face IDs, including for file-backed
  meshes. Metal approximates native Ptex filtering with a bounded,
  mipmapped texture atlas and reports invalid face lookups.
- `vdb_medium()` reads sparse OpenVDB and NanoVDB grids, including ZIP- and
  BLOSC-compressed NanoVDB, and preserves grid transforms, density, and
  temperature without expanding the volume into an R array.
  `nanovdb_medium()` remains an alias. `grid_medium()` also accepts spatial
  emission multipliers.
- `point_light()` and `spot_light()` add inverse-square lights with shadows
  and attenuation through media. `projection_light()` projects an image onto
  surfaces or through a volume. `disk_light()` places a uniform disk at
  infinity with a specified direction and angular diameter.
  `add_light()`, `get_light()`, `list_lights()`, and `remove_light()` manage
  scene lights by name.
- `csg_curve()` makes closed, tapered tubes along 3D Bézier paths.
  `csg_object()` supports a material for each component, and `csg_combine()`
  blends materials across smooth joins. Set `blend_materials = FALSE` for
  a sharp transition. `bezier_curve()` exposes `split_depth` to trade setup
  time and memory for intersection speed.

### Updating existing scenes

- Mesh constructors now count actual Loop subdivision passes.
  `subdivision_levels = 0` is the default and leaves the mesh as supplied.
  `1` performs one pass. Change an explicit `1` to `0` if it was used to
  disable subdivision.
- `diffuse()` replaces its legacy `fog` and `fogdensity` arguments with
  explicit media. Attach a `homogeneous_medium()`, `grid_medium()`, or
  `vdb_medium()` to the object with `set_medium()`.
- `diffuse()` uses energy-preserving Oren–Nayar for rough surfaces.
  `sigma = 0` is Lambertian, and roughness reaches its maximum at
  `sigma = 90`. Bump textures now use energy-conserving normal mapping
  throughout, including on imported diffuse faces and vertex colors. Drop
  the development-only `normal_mapping` argument. Existing bump strengths
  may need adjustment because slopes now use the texture's UV scale.

## Bugfixes

- `glossy()` fixes excessive brightness and reflections through bumped
  surfaces. Diffuse and glossy bumps retain direct light and smooth normals.
  Texture filtering follows reflected and refracted rays, and `render_ao()`
  uses the viewing direction for colored shading.
- `microfacet()` fixes rough glass that rendered black or produced negative
  radiance. Index-matched glass transmits direct light with the correct
  absorption and tint. Glass and subsurface boundaries now resolve coincident
  or grazing contacts consistently, including rough interfaces on GCC builds.
- `hair()` fixes bright outliers at grazing angles and provides an albedo
  guide for denoising. It accepts RGB color and absorption. Denoising now
  follows mirror and glass paths to the first non-specular interaction,
  and adaptive sampling correctly averages guides for volumes and surfaces.
- `extruded_polygon()` fixes cap orientation, holes, per-feature heights, and
  reflected scales. `extruded_path()` fixes seams on closed paths, trimmed
  endpoints, twist, taper normals, and zero-width tips. It adds
  `initial_normal`, `smooth_angle`, and `arc_tolerance` for sweep geometry.
  Both functions use the mesh builders in rayvertex.
- `ply_model()` reads supplied vertex normals and counts vertices and UVs
  correctly. Mirrored meshes and primitives retain their normal and emitter
  orientation. Cylinders handle axial cap hits, Bézier curves retain the
  nearest valid hit, and distant or thin CSG surfaces remain visible.
- `read_pbrt()` preserves film ISO, camera handedness (`camera_flip_x`),
  curve subdivision, ribbon normals, cutout opacity, texture offsets, and
  emitter orientation. It preserves HDR values from PFM images, converts
  EXR environment primaries to linear sRGB, and retains named volume grids
  and blackbody emission scaling. Nested instances keep their own textures.
- `infinite_light()` preserves the energy of bright texels at the longitude
  seam. Metal fixes dark skies seen through empty high-altitude clouds,
  bright outliers at the Sun's edge, and spotlights' angular falloff.
  Procedural noise is continuous across negative coordinates, and JPEG, HDR,
  and EXR texture loading handles images without PNG alpha.
- `generate_camera_motion()` keeps fixed look-at targets stable during orbits
  and omits the duplicate closing frame on a loop. The final frame blurs
  toward the first. Preview fixes also address camera jumps, reversed controls
  on reflected cameras, physical-lens keyframes, and brightness after a move.
- Camera rays use the authored shutter interval for animated geometry. Metal
  fixes ghost images in transform motion blur and keeps shadows and moving
  lights at the same time as the surface seen by each ray.
- `render_scene()` isolates recoverable medium and subsurface failures to the
  affected ray. The returned image records counts in `path_warnings`. Set
  `RAYRENDER_DEBUG_PATHS=true` for a summary and diagnostic log.
  Opaque volume renders can finish without waiting for alpha convergence.

## Documentation

- The rendered vignettes *Composing surface textures* and *Authoring and
  reusing animated models* cover texture mapping and filtering, independent
  animation playback, pose inspection, and full-pose motion blur. The
  installed `animation/` examples include glTF fixtures and Blender scripts
  for a flapping bird, a pulsing jellyfish, and a walking character.

## Other

- Large scenes take less time and memory to assemble with `add_object()`,
  `group_objects()`, and `rbind()`. Instances share geometry and acceleration
  structures. Metal animation reuses buffers across frames and refits changing
  geometry. CPU and GPU rendering also reduce work for image lights, many
  finite emitters, heterogeneous volumes, and localized subsurface scattering.
  `bvh_type = "equal"` now selects equal-count partitioning correctly.
- Source builds use `stbimageheaders` for image I/O, `cgltf` for glTF,
  `piton` for PBRT parsing, and the `glmheaders` and `openpbr` packages for
  their C++ headers. Volume support requires `openvdbr` >= 13.1.0-5.
  Windows cleanup preserves shipped configuration files, and Linux links
  the atomic runtime where required. WebAssembly builds link
  the target OpenEXR and Imath libraries. For offline builds,
  `RAYRENDER_WASM_RLIB` can point to an R library containing those Wasm binaries.
- Release checks exercise installation without OpenUSD. The benchmark workflow
  measures scalar and SIMD builds of the triggering commit.

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
