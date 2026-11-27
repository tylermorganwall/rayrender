library(rayrender)

# Run without a display or external model. Output frame order is independent
# of scene time; frames 2 and 4 request the same pose after different seeks.
path = animation_transform_track(
  c(0, .7, 2.3, 4),
  position = rbind(c(-1, 0, 0), c(0, 1, 0), c(1, 0, 0), c(-1, 0, 0)),
  position_interpolation = "cubic"
)
scene = animation_keyframes(
  sphere(radius = .4, material = light(intensity = 1)),
  path
) |>
  add_camera(camera(
    lookfrom = c(0, 0, 6),
    aperture = 0,
    lookat = c(0, 0, 0),
    fov = 0,
    ortho_dimensions = c(4, 4)
  ))
frames = render_animation(
  scene,
  times = c(2, .2, -1, .2, 4),
  width = 40,
  height = 40,
  samples = 4,
  integrator_type = "nee",
  preview = FALSE,
  plot_scene = FALSE,
  denoise = FALSE,
  progress = FALSE,
  ambient_light = FALSE,
  backgroundhigh = "black",
  backgroundlow = "black",
  bloom = FALSE
)
stopifnot(length(frames) == 5, identical(frames[[2]], frames[[4]]))

# To write files, add filename = "flight"; output files use sequence indices,
# not the requested seconds. start_frame/end_frame subset that output sequence.
