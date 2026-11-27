## ---- model-setup ----
library(rayrender)

## ---- model-geometry ----
# Build each part in its own coordinates, with its joint at the origin.
blue = diffuse("#287F94")
orange = diffuse("#E6A045")
dark = diffuse("#263745")
base_geometry = cylinder(
  y = .10,
  radius = .43,
  length = .20,
  material = dark
) |>
  add_object(cylinder(y = .23, radius = .29, length = .10, material = blue))
upper_geometry = sphere(radius = .15, material = orange) |>
  add_object(cube(
    y = .45,
    xwidth = .17,
    ywidth = .90,
    zwidth = .22,
    material = blue
  ))
forearm_geometry = sphere(radius = .14, material = dark) |>
  add_object(cube(
    y = .40,
    xwidth = .14,
    ywidth = .80,
    zwidth = .19,
    material = orange
  ))
gripper_geometry = sphere(radius = .12, material = dark) |>
  add_object(cube(
    y = .17,
    xwidth = .46,
    ywidth = .12,
    zwidth = .19,
    material = blue
  )) |>
  add_object(cube(
    x = -.20,
    y = .33,
    xwidth = .08,
    ywidth = .30,
    zwidth = .13,
    material = dark
  )) |>
  add_object(cube(
    x = .20,
    y = .33,
    xwidth = .08,
    ywidth = .30,
    zwidth = .13,
    material = dark
  ))

## ---- model-nodes ----
base = animation_model_node("base", base_geometry)
upper = animation_model_node(
  "upper",
  upper_geometry,
  parent = "base",
  transform = animation_transform_track(
    0,
    position = c(0, .32, 0),
    angle = c(0, 0, -25)
  )
)
forearm = animation_model_node(
  "forearm",
  forearm_geometry,
  parent = "upper",
  transform = animation_transform_track(
    0,
    position = c(0, .90, 0),
    angle = c(0, 0, 65)
  )
)
gripper = animation_model_node(
  "gripper",
  gripper_geometry,
  parent = "forearm",
  transform = animation_transform_track(
    0,
    position = c(0, .80, 0),
    angle = c(0, 0, -30)
  )
)

## ---- model-clip ----
# Channels are named by node ID. These rotations replace the rest rotations;
# the absent position channels retain each joint's rest offset.
reach = animation_model_clip(
  "reach",
  channels = list(
    upper = animation_transform_track(
      c(0, .5, 1, 1.5, 2),
      angle = rbind(
        c(0, 0, -25),
        c(0, 0, -55),
        c(0, 0, -25),
        c(0, 0, 15),
        c(0, 0, -25)
      ),
      rotation_interpolation = "smooth"
    ),
    forearm = animation_transform_track(
      c(0, .7, 1.3, 2),
      angle = rbind(c(0, 0, 65), c(0, 0, 105), c(0, 0, 30), c(0, 0, 65)),
      rotation_interpolation = "smooth"
    ),
    gripper = animation_transform_track(
      c(0, .4, 1, 1.6, 2),
      angle = rbind(
        c(0, 0, -30),
        c(0, 45, -10),
        c(0, 0, -30),
        c(0, -45, -50),
        c(0, 0, -30)
      ),
      rotation_interpolation = "smooth"
    )
  ),
  name = "Reach and turn"
)

## ---- model-asset ----
robot = animation_model_asset(
  nodes = list(base, upper, forearm, gripper),
  clips = list(reach),
  name = "Inspection robot"
)
animation_list_clips(robot)

## ---- model-rest ----
# No selected clip means the authored rest pose, even when clips are available.
rest = animation_model_instance(robot)
studio = xz_rect(xwidth = 200, zwidth = 200, material = diffuse("#DDD7CC")) |>
  add_object(sphere(
    x = -3,
    y = 6,
    z = 4,
    radius = 1.4,
    material = light("#FFF3DF", intensity = 5)
  ))
rest_scene = add_object(studio, rest) |>
  add_camera(
    camera(
      lookfrom = c(3.7, 2.9, 5),
      lookat = c(0, 1.1, 0),
      fov = 32,
      aperture = 0
    ),
    name = "rest"
  )

## ---- model-paths ----
# Both bases stay on the floor. Paths use scene seconds, independently of clips.
left_path = animation_transform_track(
  c(0, 1.2, 2.8, 4),
  position = rbind(
    c(-1.4, 0, .2),
    c(-1.7, 0, -.5),
    c(-.9, 0, -.7),
    c(-1.4, 0, .2)
  ),
  angle = rbind(c(0, -15, 0), c(0, 25, 0), c(0, 50, 0), c(0, -15, 0)),
  position_interpolation = "cubic",
  rotation_interpolation = "smooth"
)
right_path = animation_transform_track(
  c(0, .8, 2.4, 4),
  position = rbind(c(1.4, 0, -.3), c(1.7, 0, .3), c(.9, 0, .5), c(1.4, 0, -.3)),
  angle = rbind(c(0, -35, 0), c(0, -5, 0), c(0, 25, 0), c(0, -35, 0)),
  position_interpolation = "cubic",
  rotation_interpolation = "smooth"
)

## ---- model-instances ----
left_robot = animation_model_instance(
  robot,
  transform = left_path,
  animation = animation_clip_player("reach", loop = "repeat")
)
right_robot = animation_model_instance(
  robot,
  transform = right_path,
  animation = animation_clip_player(
    "reach",
    offset = .7,
    speed = .5,
    loop = "repeat"
  )
)
scene = add_object(studio, left_robot) |>
  add_object(right_robot) |>
  add_camera(
    camera(
      lookfrom = c(3.7, 3.2, 7),
      lookat = c(0, 1, 0),
      fov = 31,
      aperture = 0
    ),
    name = "pair"
  )

## ---- model-inspect ----
# Query forward, backward, and again. Compare pose values, not change flags.
requested_times = c(1.5, .25, 1.5)
poses = animation_model_pose(left_robot, requested_times)
stopifnot(identical(poses[[1]]$world, poses[[3]]$world))
left_pose = animation_model_pose(left_robot, 1)[[1]]
right_pose = animation_model_pose(right_robot, 1)[[1]]
data.frame(
  instance = c("left", "right"),
  scene_time = 1,
  clip_time = c(left_pose$clip_time, right_pose$clip_time)
)
# The translation column gives the gripper origin in world coordinates.
left_pose$world$gripper[1:3, 4]

## ---- model-render-settings ----
# Fixed settings make the supplied figures reproducible on the CPU renderer.
render_settings = list(
  width = 720,
  height = 480,
  samples = 128,
  sample_method = "sobol",
  min_variance = 0,
  max_depth = 12,
  ambient_light = TRUE,
  backgroundhigh = "#C7D9E8",
  backgroundlow = "#E5E1D9",
  preview = FALSE,
  progress = FALSE,
  plot_scene = FALSE,
  denoise = FALSE,
  bloom = FALSE,
  tonemap = "raw",
  integrator_type = "nee",
  shutter_offsets = c(0, 0)
)

## ---- model-render-rest ----
set.seed(20261018)
rest_image = do.call(
  render_scene,
  c(
    list(scene = rest_scene, time = 0),
    render_settings
  )
)
rayimage::ray_write_image(rest_image, "figures/model-animation/rest.png")

## ---- model-render-frames ----
# One scene, four requested times, one static scene camera.
frame_times = c(0, .5, 1, 1.5)
set.seed(20261018)
frames = do.call(
  render_animation,
  c(
    list(scene = scene, times = frame_times),
    render_settings
  )
)
for (i in seq_along(frames)) {
  rayimage::ray_write_image(
    frames[[i]],
    sprintf("figures/model-animation/pair-%02d.png", i)
  )
}
