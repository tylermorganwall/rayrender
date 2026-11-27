# Original local glTF example. Source this file after library(rayrender).
imported_models_scene = function() {
  file = system.file(
    "animation",
    "fixtures",
    "rigid.glb",
    package = "rayrender"
  )
  stopifnot(nzchar(file))
  asset = animation_read_gltf(file, material = diffuse("#e9a23b"))
  first = animation_model_instance(
    asset,
    animation_transform_track(
      c(0, 1, 2, 3),
      position = rbind(c(-3, 1, 0), c(-2, 1.5, 1), c(-3, 2, -.5), c(-3, 1, 0)),
      angle = rbind(c(0, 0, 0), c(15, 45, 0), c(-10, 90, 15), c(0, 0, 0)),
      position_interpolation = "cubic",
      rotation_interpolation = "smooth"
    ),
    animation_clip_player("clip-0", loop = "repeat")
  )
  second = animation_model_instance(
    asset,
    animation_transform_track(
      c(0, .8, 2, 3),
      position = rbind(c(2, 1, 0), c(3, 2, -.5), c(2, 1, 1), c(2, 1, 0)),
      angle = rbind(
        c(0, -30, 0),
        c(-15, -60, -10),
        c(15, -90, 0),
        c(0, -30, 0)
      ),
      position_interpolation = "cubic",
      rotation_interpolation = "smooth"
    ),
    animation_clip_player("clip-0", offset = .5, speed = -1.2, loop = "repeat")
  )
  scene = add_object(first, second) |>
    add_object(xz_rect(
      y = -.05,
      xwidth = 18,
      zwidth = 14,
      material = diffuse(
        checkercolor = "grey25",
        color = "grey80",
        checkerperiod = 1
      )
    )) |>
    add_object(sphere(
      x = -4,
      y = 9,
      z = 5,
      radius = 1.5,
      material = light(intensity = 40)
    )) |>
    add_camera(camera(lookfrom = c(11, 9, 17), lookat = c(1, 2, 0), fov = 40))
  # Reloading preserves model phases, hierarchy, cameras, and shared native meshes.
  scene = unserialize(serialize(scene, NULL))
  scene
}
