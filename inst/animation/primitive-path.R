library(rayrender)

# Run with Rscript inst/animation/primitive-path.R. No external assets or display
# are required. Change preview/plot_scene below to inspect a selected snapshot.
path = animation_transform_track(
  time = c(0, .7, 2.3, 4),
  position = rbind(c(-2, 0, 0), c(-1, 1, 0), c(1, 0, 0), c(2, 1, 0)),
  angle = rbind(c(0, 0, 0), c(0, 40, 0), c(0, 80, 0), c(0, 120, 0)),
  position_interpolation = "cubic",
  rotation_interpolation = "smooth"
)
scene = animation_keyframes(
  cube(x = .2, width = .6, material = diffuse(color = "red")),
  path
) |>
  add_object(xz_rect(y = -.6, xwidth = 12, zwidth = 12)) |>
  add_object(sphere(y = 5, radius = 2, material = light(intensity = 10)))

# An exact snapshot uses a zero-width exposure by default.
snapshot = render_scene(
  scene,
  time = 1.4,
  integrator_type = "nee",
  width = 120,
  height = 90,
  samples = 16,
  lookfrom = c(0, 3, 9),
  lookat = c(0, 0, 0),
  aperture = 0,
  fov = 45,
  preview = FALSE,
  plot_scene = FALSE,
  progress = FALSE,
  denoise = FALSE
)
stopifnot(all(is.finite(snapshot)))

# The exposure crosses the interior .7-second key; bounds enclose the full curve.
exposure = render_scene(
  scene,
  time = .7,
  shutter_offsets = c(-.3, .3),
  integrator_type = "nee",
  width = 120,
  height = 90,
  samples = 32,
  lookfrom = c(0, 3, 9),
  lookat = c(0, 0, 0),
  aperture = 0,
  fov = 45,
  preview = FALSE,
  plot_scene = FALSE,
  progress = FALSE,
  denoise = FALSE
)
stopifnot(all(is.finite(exposure)))

# Two full turns need intermediate orientations; identical endpoints alone do
# not encode winding. A separate position channel keeps its own key schedule.
turns = animation_transform_track(
  time = 0:8,
  angle = cbind(0, (0:8) * 90, 0),
  rotation_interpolation = "smooth",
  position = animation_keyframe_track(
    c(0, 3, 8),
    rbind(c(0, 0, 0), c(1, 2, 0), c(3, 0, 0)),
    interpolation = "cubic"
  )
)
stopifnot(inherits(turns, "ray_transform_track"))
