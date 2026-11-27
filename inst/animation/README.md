# Animation examples

For visual demonstrations with modeled geometry, use the
[Blender cycle generators](blender/README.md): a feathered bird, a pulsing
jellyfish, and a walking cycle using a supplied character rig. They produce
editable Blender scenes, animated GLBs, and rendered loop movies.
`blender-models.R` imports those GLBs through native node, morph, and skeletal
animation. It places two independently phased copies of each asset on separate
five-key paths. No Blender process or baked geometry sequence is used during
rayrender rendering. The supplied character is not distributed with the package.

## Render the authored models

After generating the GLBs as described in [the Blender README](blender/README.md):

```r
library(rayrender)
source(system.file("animation", "blender-models.R", package = "rayrender"))
assets = path.expand("~/Desktop/3dmodels/rayrender-cycles")

# A native morph snapshot at an absolute scene time, with an authored camera.
jellyfish = blender_jellyfish_scene(assets)
render_scene(jellyfish$scene, time = .75, shutter_offsets = c(0, 0),
  width = 640, height = 480, samples = 64, preview = FALSE)

# Generate three four-second movies and a combined movie (requires ffmpeg).
# Use integrator_type = "nee" for CPU. Metal requires a supported macOS device.
render_blender_model_examples(assets, file.path(assets, "rayrender-examples"),
  samples = 32, integrator_type = "metal")

walkers = blender_walk_scene(assets)
# The repeated time returns the same pose after a backward seek.
frames = render_animation(walkers$scene, times = c(1, 3, .25, 1),
  shutter_offsets = c(0, 0), width = 320, height = 240, samples = 32,
  preview = FALSE, progress = FALSE, plot_scene = FALSE)
if (interactive()) {
  render_scene(walkers$scene, time = 0, timeline = c(0, 4),
    shutter_offsets = c(0, 0), samples = 128, interactive = TRUE)
}
```

The movie renderer preflights every requested Metal time at low resolution,
then records explicit scene times and actual backend use in `render-report.json`.
It saves PNG frames and per-frame backend reports before rejecting any fallback
in the full-size render. A requested Metal render that falls back to CPU is an
error in this example. Each frame is evaluated from immutable rest geometry.
Walker placement advances at approximately 0.75 units/second while the two
internal walk clips loop with different phases; changing placement does not
change the gait. The bird and jellyfish also use different clip speeds.

Imported glTF materials are replaced by explicit rayrender materials. The
examples assign colors to authored parts and use rayrender studio lights.
Deforming assets require explicit `time`/`times`. Nonzero shutter intervals
freeze all internal node transforms, skins, and morphs at frame time while
sampling outer placement at each ray's shutter time. This is whole-model motion
blur, not deformation blur. A message is emitted once per render/preview
session. Zero-width exposures sample the exact frame time plus shutter offset.

The canonical animation API uses the `animation_` prefix; the earlier unprefixed
`model_*()` aliases have been removed; use `animation_model_*()` instead. Other
legacy names remain compatibility aliases. Existing render entry points are unchanged.

For a complete authoring workflow with rendered figures, see
`vignette("model-animation", package = "rayrender")`. Its reproducible script is
`vignettes/scripts/model-animation.R` in the source checkout.

The small examples below are implementation and regression fixtures.
`primitive-path.R` runs a four-key primitive path, exposure spanning an interior
key, and independently timed rotation/position authoring. `offline-times.R`
renders nonmonotonic requested times with a static camera and checks that
repeating a time reproduces the same frame. Both run without external assets:

```r
source(system.file("animation", "primitive-path.R", package = "rayrender"))
source(system.file("animation", "offline-times.R", package = "rayrender"))
```

## Original rigid animation example

`rigid-birds.R` builds one shared primitive model with a body and two animated
wings. Its two instances use different five-key 3D paths, clip phases, and speeds.
The geometry and motion are original; no external model or texture is required.

```r
library(rayrender)
source(system.file("animation", "rigid-birds.R", package = "rayrender"))
birds = rigid_birds_scene()
animation_list_clips(birds$asset)
animation_model_pose(birds$first, c(2, .4, 2))

# Exact requested-time snapshots with a static scene camera.
frames = render_animation(birds$scene, times = (0:47) / 12, fps = 12,
  width = 320, height = 240, samples = 64, preview = FALSE,
  progress = FALSE, plot_scene = FALSE, denoise = FALSE, tonemap = "hbd")
rayimage::ray_write_image(frames[[1]], filename = "birds-first.png")

# Play, reverse, change speed, or drag to any scene time in the existing preview.
if (interactive()) {
  render_scene(birds$scene, time = 0, timeline = c(0, 4), fps = 30,
    preview = TRUE, interactive = TRUE, samples = 128)
}
```

Positions and clip clocks use seconds. Changing preview speed does not change
the authored clip clocks or shutter duration. Zero-width exposures are snapshots;
`shutter_offsets = c(-.05, .05)` blurs outer placement with a frozen internal pose.

`imported-models.R` uses the original local glTF/GLB fixtures through animation_read_gltf():

```r
library(rayrender)
source(system.file("animation", "imported-models.R", package = "rayrender"))
scene = imported_models_scene()
frames = render_animation(scene, times = c(1, .25, 2, .25),
  width = 320, height = 240, samples = 32, preview = FALSE,
  progress = FALSE, plot_scene = FALSE, denoise = FALSE)
if (interactive()) {
  render_scene(scene, time = 0, timeline = c(0, 4), interactive = TRUE)
}
```

animation_read_gltf() supports rigid node, skeletal, and morph animation with explicit
rayrender material overrides. Morph deltas are applied before linear skinning,
using every supplied paired JOINTS/WEIGHTS set after total-weight normalization.
The skin palette is inverse(mesh global) * joint global * inverse bind; outer
placement is applied afterward. Missing inverse binds default to identity.
By default, normals use the inverse transpose of the blended skin map; tangents use its
linear part, are orthogonalized, and retain reflection handedness. Singular
blends and malformed influences fail explicitly. CPU and Metal acceleration
structures are refitted when a pose changes. Metal updates the affected mesh
buffers while retaining triangle ordering, BVH topology, static geometry and
GPU work buffers. Placement alone and repeated identical poses need no reskinning
or geometry refit. Frame reports expose `acceleration_builds`,
`acceleration_refits`, and `geometry_bytes_updated`; after initial construction,
pose changes report zero builds. Refits preserve correctness for large pose
changes, although the original BVH partition can become less efficient.
Use `animation_read_gltf(..., recalculate_normals = TRUE)` to replace authored normals
with smooth, corner-angle-weighted normals from each posed triangle surface.
This includes the effect of varying skin weights and position-only morphs on
shading, while retaining imported vertex splits. The walking example enables
this option to keep foot and body shading consistent with the deformed mesh.
Source glTF/GLB buffers are decoded into owned data; reloaded assets do not need
those source files. Copies share immutable native geometry and retain independent
clip settings. See fixtures/README.md for the original test-data provenance.

## Frozen internal poses with placement motion blur

`frozen-model-blur.R` uses the authored walking GLB from the Blender generator.
It renders two independently phased walkers alongside a fast moving orange
sphere. The sphere and walking paths blur over the shutter while each internal
walking pose stays fixed within that frame. It writes 12 PNG frames and a
sharp/blurred comparison without modifying the source model.

```sh
Rscript /path/to/rayrender/inst/animation/frozen-model-blur.R /path/to/rayrender-cycles output-directory
```

For interactive seeking with the same nonzero exposure:

```r
source(system.file("animation", "frozen-model-blur.R", package = "rayrender"))
scene = frozen_model_blur_scene("/path/to/rayrender-cycles")
if (interactive()) {
  render_scene(scene, time = 0, timeline = c(0, 2),
    shutter_offsets = c(-.1, .1), interactive = TRUE)
}
```
