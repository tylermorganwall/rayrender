# rayrender 0.42.3.9000

* Diffuse materials now always use the normal-mapping model of
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

- `animation_model_instance()` Freezes complete internal poses at frame time
  during nonzero exposures, while outer placement and other scene objects retain
  motion blur. Rendering reports this once as a message instead of rejecting
  skinned or morph-animated models. Offset point exposures remain exact.
- Overall: New animation authoring, import, and inspection functions use the
  `animation_` prefix. The five unprefixed `model_*()` development aliases are
  removed; use the corresponding `animation_model_*()` functions.

- `read_gltf()` Adds `recalculate_normals = TRUE` to derive smooth normals from
  each deformed pose, preventing skin-weight gradients from leaving shading
  normals inconsistent with the surface. Authored normals remain the default.

- Blender example generators now create a feathered bird flap, a pulsing
  jellyfish, and a walking cycle from a supplied Rigify character, with editable
  scenes, GLB exports, preview movies, and export/loop validation. Walk checks
  cover forward knee alignment, foot spacing, extended heel-first steps,
  upright posture, reduced vertical bob, stance contact, and a smooth pelvic
  crease. `inst/animation/blender-models.R`
  renders all three GLBs natively in rayrender, with two independently phased
  instances following separate multi-key placement paths.

- `read_gltf()` Imports local glTF/GLB triangle models and rigid node clips into
  shared model assets, including sparse/interleaved accessors, independent
  property schedules, and authored cubic quaternion tracks. The parser headers
  are supplied by the separate `gltfheader` package through `LinkingTo`.
  Native morph and skeletal animation support all supplied joint-weight sets,
  inverse binds, rest weights, and STEP/linear/cubic morph channels. Morphs
  precede skinning; normals and tangents follow the blended deformation.
  Exact-time snapshots update CPU bounds and Metal geometry, including during
  preview seeks. Mesh light sampling uses the current triangle areas and
  matching probability densities. Deforming models require explicit scene times. Nonzero exposures freeze
  internal poses while retaining outer placement motion blur.

- `animation_model_asset()` Shares rigid model geometry and named clips across independently
  clocked `animation_model_instance()` placements. `animation_model_node()` and
  `animation_model_clip()` author
  parent/child parts with local pivots and independent property schedules;
  `animation_model_pose()` evaluates exact local, global, and world matrices without a
  renderer. `list_animations()` exposes stable clip IDs and target metadata.
  `clip_player(root_motion = "in_place")` removes explicitly selected translation
  axes and optional yaw while preserving vertical bob and child articulation.
  Saved scenes retain shared native geometry after reload, and copied instances
  keep independent playback settings through nested scene composition.

- `render_scene(timeline = ...)` Adds a preview slider, reverse playback, frame
  stepping, speed controls, and progressive convergence when paused. Authored
  camera rows and saved camera playback use scene time, including shutter
  exposure, with perspective, orthographic, or panoramic projections.
  Camera exposures visit interior keys on CPU and Metal and remain independent
  of seek history in preview, explicit stills, and offline output.
  Projection changes cut at the authored key, including within an exposure.
  Physical-camera tracks support exact rigid exposure and preview playback with
  constant aperture and focus, preserving initialized lens optics while seeking.
  Zero-width physical-camera playback also supports changing aperture and focus;
  complete lens preparation preserves the same optics after reverse seeks.
  Frozen explicit-time lens-camera snapshots retain unit exposure weighting
  instead of becoming black; legacy lens exposure weighting is unchanged.
- `clip_player()` Maps scene time to independent hold, repeat, or pingpong
  clocks, including reverse and frozen playback. `animate_keyframes(player = ...)`
  applies the clock to each shutter ray on CPU and Metal; `animation_loop_seam()`
  diagnoses endpoint pose jumps without altering authored keys.
- `animate_keyframes()` Adds multi-key outer object transforms with independently
  timed position, quaternion rotation, and scale channels. `keyframe_track()`
  and `transform_track()` support STEP, linear/SLERP, cubic, and smooth rotation.
  `render_scene(time = ...)` evaluates snapshots in scene seconds, with optional
  relative shutter offsets and preserved legacy two-endpoint motion.
- `render_animation(times = ...)` Renders explicit, possibly nonmonotonic scene
  times while retaining compiled geometry. Camera tracks use explicit time
  columns or timestamps derived from `fps`. Wavefront motion packs independent
  property tracks and rejects key schedules that lose device time precision.
  Exact snapshots preserve the authored side of imported repeat seams when
  the clip period is not exactly representable on the device.
- `projection_light()` Projects an image onto surfaces and through participating
  media on CPU and Metal, with image orientation, field of view, and intensity
  controls. `read_pbrt()` imports projection lights.
- `csg_object()` Supports per-component surface materials, with `NULL` inheriting
  the enclosing object's material. `csg_combine()` blends these materials across
  smooth joins on CPU and Metal; `blend_materials = FALSE` retains sharp material
  transitions. Dielectric and subsurface materials cannot be mixed.
- `csg_curve()` Creates closed, tapered tubes along 3D Bezier paths for tentacles,
  stems, and other curved solids that can participate in CSG operations.

- `dielectric()` Supports open interfaces on Metal, including water surfaces
  bounded by opaque terrain, with refraction, absorption, and material priorities.
  Open mesh instances retain shared geometry; closed subsurface and volume
  boundaries keep their containment behavior.

- `render_scene()` and `render_animation()` Default to `integrator_type = "auto"`,
  selecting Metal when a supported device is available and CPU NEE otherwise.
  Explicit integrator choices are preserved. PBRT imports use the same default.

- `render_scene()` Sizes interactive navigation from visible geometry, with
  finite scene bounds as a fallback. Held movement uses elapsed time; orbit,
  viewport-scaled pan, and proportional dolly controls retain E/C speed overrides.
  Picking a new target updates navigation scale without changing lens focus on
  right clicks.

- `render_scene()` Supports legacy Perlin noise, nested checker textures, and
  RGB/HSV gradients in UV or world coordinates with the Metal integrator,
  preserving the CPU texture definitions instead of falling back to CPU rendering.

- `texture_ptex()` Supports Metal through a bounded, mipmapped GPU atlas.
  Filtering approximates native Ptex; atlas resolution and memory limits are
  configurable without changing CPU rendering.

- `render_scene()` Intersects flat, cylindrical, and ribbon Bézier curves directly
  on Metal, preserving the CPU curve model and avoiding tube/strip tessellation.
  Curves support hair shading, tapered widths, textures, and animated transforms.

- `render_scene()` Uses a spatial light BVH for Metal scenes with many finite
  emitters, reducing noise from uneven lighting. Metal also honors guided
  direction and free-flight sampling in `subsurface(method = "guided")`.

- `render_scene()` Supports invisible area lights and transparent backgrounds on
  Metal. Hidden emitters still illuminate the scene; camera coverage includes
  alpha cutouts, volume extinction, and atmospheric haze.

- `csg_object()` Renders directly on Metal, including smooth blends, shells,
  procedural textures, emissive shapes, and animated transforms. Closed CSG
  solids support dielectric priorities and subsurface or volume interiors.

- `render_scene()` Supports physical lens cameras, panoramic cameras, camera motion
  blur, and animated geometry on Metal. Camera and geometry motion use independent
  per-pixel shutter times. Density-threshold haze runs on the GPU for dense and
  sparse volumes, retaining every crossing of the interpolated density cutoff.

- `render_scene()` Supports homogeneous, dense-grid, and sparse OpenVDB/NanoVDB
  volumes on Metal, including transformed density fields, scattering, and RGB or
  thermal emission. Prague atmospheres evaluate skymodelr's spectral model and
  finite haze directly on Metal, including atmospheric light before shadow blockers.
  Multiple Prague skies retain batched CPU atmosphere queries.
  Interactive haze and altitude controls remain available; `wavefront$atmosphere_seconds`
  reports the CPU atmosphere portion of rendering time.

- `render_scene()` Samples Metal image environments with luminance-weighted alias
  tables, reducing noise from compact bright sources and preserving interactive
  rotation. Mesh instances share GPU geometry and acceleration structures while
  retaining per-placement texture mapping, emission, and medium priorities.

- `render_scene()` Intersects spheres and ellipsoids analytically on Metal,
  retaining smooth silhouettes and accurate texture coordinates without
  tessellation. Transformed instances, area lights, and priority-aware SSS
  use the same analytic geometry.

- `render_scene()` Evaluates composable color and roughness textures on Metal,
  including UV/object/world mappings, seeded noise, direction mixes, and image
  decoding/addressing. Metal also supports image alpha coverage and filtered bump
  maps, including imported bump graphs and energy-conserving diffuse normal mapping.

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

- Overall: Animated skin and morph poses update resident Metal geometry and
  refit existing acceleration structures, preserving BVH topology and static
  scene buffers across frames and timeline seeks.

- Overall: Adaptive sampling now preserves correctly averaged denoising guides,
  fixing residual noise and inconsistent smoothing in Metal renders and previews.

- `render_scene()` Saved physical-camera keyframes retain their lens projection,
  aperture, and viewing direction. Timeline picking uses the displayed lens
  pose, and rigid saved-camera navigation restores the authored position.

- `spot_light()` Uses the same smooth angular falloff on Metal as on the CPU.

- `render_scene()` Preserves camera elevation when orbiting horizontally with
  A/D, including repeated taps and held keys, in CPU and Metal previews.

- `render_scene()` Keeps interactive preview brightness consistent immediately
  after camera, lighting, or resolution changes. CPU and Metal restart sampling
  from the cleared film's actual sample count instead of dividing the first
  new frame by two.

- `generate_camera_motion()` Samples closed paths without a duplicated closing
  frame. With `closed = TRUE`, `render_animation()` automatically blurs the last
  frame toward the first, including when rendering a subset of the animation.
  Open paths retain their stationary final shutter endpoint.

- `render_scene()` Prevents occasional large camera jumps when tapping movement
  keys during slow preview updates. Separate presses remain separate taps, and
  only uninterrupted holds accumulate elapsed movement time.

- `render_animation()` Retains Metal scene buffers and acceleration structures
  across compatible frames, avoiding repeated uploads and construction. Frame
  accumulation and diagnostics reset independently; returned frames report
  resource reuse and per-frame setup and sampling times.

- `read_pbrt()` Preserves persistent alpha-mask paths in instance prototypes,
  allowing imported scenes saved with `saveRDS()` to render in later R sessions.

- `render_scene()` Corrects reversed horizontal navigation and roll controls
  for imported PBRT cameras with horizontally reflected film coordinates.

- `render_scene()` Fixes coherent ghost poses in Metal geometry motion blur by
  sampling time independently for each camera ray. Shading, shadows, moving
  emitters, textures, and volume transport retain that time throughout the path.
  Nested rotation and scale/shear use native transform interpolation without
  per-sample acceleration-structure refits.

- `render_scene()` Places the interactive status bar below the image and identifies
  the active integrator, including CPU fallback when Metal is unavailable.

- `render_scene()` Updates the interactive `F` preview after each low-resolution
  sample instead of waiting for four samples, allowing more frequent camera
  interaction while retaining progressive accumulation and denoising.

- Overall: Diffuse and glossy bump mapping includes up to 16 local facet
  interactions on CPU and Metal, reducing darkening from truncated scattering.
  This applies to both Lambertian and EON diffuse reflection and retains material
  color and absorption. Local interactions do not consume `max_depth`; very steep
  bumps can still lose energy through truncation or geometric clipping.

- Overall: Higher-order normal mapping retains direct illumination from point,
  spot, area, and environment lights. Light sampling and surface scattering use
  matching evaluation and directional probabilities for multiple importance
  sampling (MIS). CPU basic and guided rendering also use the higher-order
  scattering evaluation.

- `glossy()` Corrects excessive brightness on CPU and Metal by applying the
  outgoing surface cosine once and matching sampling probabilities to the
  diffuse/microfacet mixture. This affects unbumped surfaces as well as bump maps,
  with both GGX and Beckmann distributions. CPU samples use the local shading
  hemisphere, fixing incorrect directions on rotated surfaces.

- `glossy()` Uses energy-conserving bump mapping relative to the unperturbed
  smooth surface normal on CPU and Metal. Evaluation and sampling share that
  frame, preserving smooth mesh shading and preventing back-facing bumps from
  amplifying light or reflecting through the surface. Texture reflectance is
  bounded to [0, 1], with nonfinite values treated as absorbing.

- `render_scene()` Preserves camera ray differentials when normalizing rays in
  the CPU NEE integrator, restoring the pixel footprints used to filter bump maps.

- Overall: Correctly transforms normal derivatives without normalizing them as
  surface normals. Flat surfaces retain finite zero derivatives instead of NaNs,
  and curved surfaces retain the correct bump response under nonuniform scaling.

- `render_ao()` Uses the incoming viewing ray when preparing colored surface
  scattering. Colored AO and bounce-count diagnostics use the same higher-order
  normal-mapping evaluation as the CPU renderers.

- `render_scene()` Applies consistent mesh normals on Metal to match CPU glossy
  and metallic shading, including nonuniform mesh placements. Prague sunlight
  uses matching disk support for radiance and sampling density, preventing bright
  outliers and distorted denoised images near the solar disk's edge.

- `render_scene()` Preserves the sky observer position across invisible volume
  boundaries and null collisions on Metal when finite haze is disabled. Empty
  parts of high-altitude clouds no longer darken the sky or incident lighting.
  Spectral shadow rays share geometry intersections across their GPU lanes.
- `csg_object()` Measures marching steps from the bounding-box entry so small
  objects remain intersectable from distant ray origins on CPU and Metal.

- `csg_object()` Preserves later crossings on clipped rays and resolves thin
  shells with scale-aware tolerances. Successive queries no longer count the
  same grazing surface repeatedly, preventing excess diffusion brightness.

- `render_scene()` Preserves large ground spheres and their checker-pattern
  alignment with analytic Metal intersections.

- `render_scene()` Avoids scanning finite triangle emitters when a Metal ray
  escapes the scene, removing a major slowdown with triangulated area lights.
  OpenPBR surfaces retain their layered BSDF instead of being mistaken for
  ordinary dielectrics during GPU export.

- `subsurface_diffusion()` Resolves coincident glass/liquid contacts consistently
  for camera and shadow rays, including containment, refraction, and depth
  accounting. Distinct closed bodies retain both memberships when their first
  containment intersections have identical distances.

- `diffuse()` Preserves smooth mesh shading by applying bump mapping
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
  path-space MIS. Its images change; `rtiow` remains the default.
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

- Overall: Adds a rendered model-animation vignette covering nodes, clips, shared
  assets, independently placed instances, and exact-time pose inspection.

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

- `render_scene()` Reuses identical motion groups throughout Metal scenes,
  including matching objects separated by scene traversal order. This avoids
  redundant ray traversal and preserves hardware-mask acceleration for repeated
  animations without changing their motion or shutter timing.

- `render_scene()` Uses less memory and renders repeated primitive instances
  faster with Metal. Animated scenes with up to seven motion groups also use
  hardware masks to avoid unnecessary intersection filtering; larger motion
  hierarchies retain the existing traversal path.

- `csg_group()` Accelerates Metal rendering of grouped curves with a shared
  distance-query BVH while preserving geometry, smooth joins, and material
  selection. Other CSG expressions retain their existing evaluation.

- `diffuse()` Removes the `fog` and `fogdensity` arguments. Create volumes with
  `set_medium()` and an explicit `homogeneous_medium()`, `grid_medium()`, or
  `vdb_medium()` instead. Removes the legacy material-based fog implementations
  from both CPU and Metal scene preparation.

- Overall: Tracks material, normal-mapping, and sampling header dependencies in
  Unix and Windows builds so incremental compilation rebuilds affected C++ files
  instead of mixing incompatible scattering-record layouts.

- `render_scene()` Further reduces Metal render times for heterogeneous volumes
  and Prague skies by reusing shadow grid traversal and skipping finite-segment
  atmosphere work when haze is disabled.

- `render_scene()` Accelerates Metal lighting through volumes and alpha-mapped
  surfaces with Prague skies when finite haze is disabled, separating shadow
  visibility and RGB attenuation from spectral light evaluation.

- `render_scene()` Reduces Metal density-threshold haze work by sharing density
  lookups and reusing cutoff polynomials. Animated scenes upload and refit only
  changing geometry; scenes with one shared motion transform rays into an
  immutable reference acceleration structure.

- `render_scene()` Uses constant-time alias sampling for CPU image environments,
  with finer blue-noise samples for image-light proposals. Large maps trade
  additional table setup and storage for faster directional sampling.

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
