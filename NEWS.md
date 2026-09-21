# rayrender 0.42.0

* Clamp interactive orbit and pitch at the camera up axis by default. Add
  quaternion-based free rotation with consistent local pitch/roll controls,
  available in the Camera pane and via `camera_rotation = "free"`.

* Align realistic-camera navigation and horizontal image orientation with the
  perspective camera, including viewport picking and animation.

* Fix spurious repeated-entry errors at translated fog-box boundaries, including
  Fast preview after volume edits near a surface.

* Native editor volumes now expose density, scattering, absorption, anisotropy,
  emission and haze controls, with live validation, undo and exported overrides.

* Organize the left pane into Camera and Sky tabs. Move Start final render and
  Export R code together above Undo/Redo, keeping shared render settings visible
  in either tab. Requires rayimgui 0.0.17 (ABI 1.12).

* Stop animation playback when the viewport is clicked or receives manual input.
  Clicks restore the starting camera; navigation restores it before moving.
  Keep M for pause/resume and F for preview quality, including during playback.

* Handle F in every focused editor pane and during animation playback. Share
  the Fast preview checkbox state, toggle once per press, and preserve text-field
  keyboard ownership and undo/redo.

* Add Loop playback to the Animation panel to repeat saved frames continuously.
  Toggle it while playing or paused; turning it off finishes the current pass.
  Preserve pause/resume, Stop restoration, and undo/redo of the loop setting.

* Right-click a keyframe to replace its saved camera and thumbnail with the
  current viewport view. Preserve its order and timing, with undo/redo support.
  Requires rayimgui 0.0.16 (ABI 1.11).

* Keep M available for animation playback after clicking a keyframe thumbnail.
  The Animation panel handles its own shortcut while active inputs retain
  keyboard focus. Add prominent Play/Pause/Resume and Stop buttons; pausing
  holds the current frame and resuming continues the same path. Requires
  rayimgui 0.0.15 (ABI 1.10).

* Add editable frame intervals between camera keyframe snapshots, including the
  closing transition, with undo/redo and synchronized position, lens, and
  orientation timing. `generate_camera_motion()` accepts `segment_frames` for
  spline, linear, quad, cubic, and exp interpolation.

* Preserve physical ray geometry above Prague's 15 km coefficient range.
  Downward haze no longer disappears in a growing circle, and celestial
  visibility follows the observer's actual horizon. Out-of-range coefficients
  still use their nearest available altitude. Remote rays advance along their
  physical direction to atmospheric entry, preserving vacuum distances and
  avoiding invalid native queries beyond the outer atmosphere.

* Sample Hosek's Sun as a celestial disk alongside an atmosphere-only image,
  using skymodelr 0.6.5. Initial renders, live Fast Sun drags, undo and exports
  keep the same solar size and illumination, including the 89.9-degree limit.
  `sun_light()` now supports Hosek turbidity and manual Sun angles.

* Update Sun direction continuously during drags with temporary Fast quality.
  Add a live Camera pane with validated position, target, up direction and lens
  controls, synchronized with navigation, keyframes, undo/redo and export.

* Keep the selection outline visible with its matching viewport image during
  object drags and render restarts, preventing flicker between completed samples.

* Add validated Prague atmosphere scale (meters per scene unit) and base altitude
  (meters above sea level) inputs to the editor. Edits rebuild lighting after
  release, persist across sky-model changes, and support undo/redo and export.

* Restart preview sample counts correctly after camera or scene edits and quality
  changes. This prevents auto-exposure from doubling after movement and avoids
  brightness drift with manual exposure.

* Hovering the editor export status shows the complete saved path. Repeated live
  object transforms retain affine matrices and consistent inverses, preventing
  accumulated shading errors and rejected translation edits.

* Preview object transforms during drags at Fast quality, restoring the chosen
  quality on release. Valid material inputs update immediately; invalid drafts
  show red fields with hover diagnostics. Texture paths include Browse dialogs.
  Requires rayimgui 0.0.13 (ABI 1.9).

* Escape closes the native editor from every panel, including active text inputs,
  popups, transform drags and animation playback. Requires rayimgui 0.0.12 (ABI 1.8).

* Add native editor undo/redo for inputs, scene edits, camera/render/sky settings
  and keyframes. Control/Command-Z undoes; add Shift to redo. Drags coalesce into
  one edit, text fields retain local history, and failed restores preserve state.
  Requires rayimgui 0.0.11 (ABI 1.7).

* Add a bottom native animation panel with camera keyframe navigation, saving,
  deletion and playback, camera motion blur/shutter and open/closed paths.
  Saved keyframes have clickable immutable thumbnails. Requires rayimgui 0.0.10.

- Add Export R code to the native editor. Save cumulative object and material
  changes, nested instance overrides, sky, camera and viewport appearance as a
  runnable script with companion scene data. Preserve temporary assets and
  replay exports in headless renders through `apply_scene_edits()`.

- Cap sky Sun elevation at 89.9 degrees before image generation or native Prague
  sampling. Keep the editor angle and separate solar disks synchronized with
  the clamped direction, including typed values and model changes.

- Preserve fractional roughness-map values instead of truncating them in the
  shared byte cache. Each material samples its own range/flip settings, and the
  editor retains the original map path and mapping when applying other edits.

- Use native Prague atmosphere transport when selected in the editor, including
  scenes that start with a Hosek image. Expose haze and altitude-query controls,
  preserve their settings across model changes, and reserve volume transport
  before interactive sampling begins.

- Add visual Sun controls: a quarter-circle elevation handle and circular
  azimuth dial, each with a numeric input inside. Requires rayimgui 0.0.9 (ABI 1.5).

- Restore the native editor sky panel for image and atmospheric skies, with a
  Hosek/Prague selector, Sun direction, and location/UTC controls. The demo starts
  with Hosek; failed sky loads leave the current lighting intact.

- Expand the native material inspector with type-specific surface parameters,
  procedural textures, and color/alpha/bump/roughness file inputs. Edits preserve
  per-instance isolation and apply atomically; alpha changes refresh mesh shadow
  classification and texture paths are included in exported scene edits.

- Show preview FPS beside the sample count in the left render panel. The counter
  averages completed preview frames over half a second, including denoising.

- Add a Denoise checkbox to the native editor's left render panel. It updates
  normal and fast previews plus the final render, preserves accumulated samples,
  and is disabled when denoising support is unavailable.

- Restore Shift-click viewport selection. Repeated Shift-clicks descend through
  groups and instance contents, with child edits isolated to the clicked copy.
  Nested objects appear in the scene tree and use matching visibility masks.
  Idle transform handles yield to Shift-click; active drags retain ownership.

- Editing a material on an Instances hierarchy node now updates the matching
  source slot on every placement, while preserving untouched per-instance values.
- Add an occlusion-aware selection outline for viewport and hierarchy
  selections, including groups and instance collections. A continuous dark/light
  border comes from cached binary coverage, with a transparent interior that
  preserves the object's rendered colors. Render noise cannot flip outline
  pixels, and preview-only resets reuse the mask. Requires rayimgui 0.0.8 (ABI 1.4).

- Add native Shift-click object selection, transform gizmos, and a material-specific
  panel. Selection stops at outer groups, instance placements, and mesh roots.
  Edits rebuild geometry, volume boundaries, and light sampling after workers
  drain; returned images carry committed `scene_edits` metadata. Requires rayimgui
  0.0.6 or later.

- Native preview: fix macOS keyboard focus with rayimgui 0.0.5, add a Fast preview
  checkbox synchronized with F, and make UTC date/time components draggable and
  directly editable before applying location/time.

* Connect standard camera movement, lens, picking and keyframe controls to the
  native preview. Add Sun elevation/azimuth and latitude/longitude/UTC controls
  for atmospheric skies, applying scene changes between completed samples.

* Use the selected tone map and the sRGB display transfer in live previews and
  preview snapshots, matching final image color processing instead of using
  a square-root gamma approximation. This applies to stills and animations,
  including the legacy and native preview windows.

* Add an optional native preview through rayimgui. Choose `gui = "auto"`,
  `"imgui"`, `"legacy"`, or `"none"` in `render_scene()`. The native window
  displays progressive CPU renders with exposure and render/atmosphere controls;
  rayrender still builds and renders without rayimgui.

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

* Add continuous local evolution to `cloud()` with `t` and an independent
  `animation_seed`. Broad and fine noise evolve smoothly while position remains
  controlled by x/y/z. Time zero preserves existing clouds exactly.

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
