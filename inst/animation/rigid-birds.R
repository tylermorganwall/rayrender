# Original primitive-only fixture; no downloaded models or textures.
# Source this file after library(rayrender), then call rigid_birds_scene().
rigid_birds_scene = function() {
  body_geometry = ellipsoid(
    a = .25,
    b = .23,
    c = .55,
    material = diffuse("#176B80")
  ) |>
    add_object(sphere(
      y = .13,
      z = -.43,
      radius = .23,
      material = diffuse("#ECA338")
    )) |>
    add_object(cube(
      y = .10,
      z = -.68,
      xwidth = .16,
      ywidth = .09,
      zwidth = .25,
      angle = c(12, 0, 0),
      material = diffuse("#C84D28")
    )) |>
    add_object(sphere(
      x = -.18,
      y = .20,
      z = -.54,
      radius = .035,
      material = diffuse("#172434")
    )) |>
    add_object(sphere(
      x = .18,
      y = .20,
      z = -.54,
      radius = .035,
      material = diffuse("#172434")
    )) |>
    add_object(cube(
      z = .56,
      xwidth = .38,
      ywidth = .06,
      zwidth = .35,
      angle = c(-15, 0, 0),
      material = diffuse("#176B80")
    ))
  body = animation_model_node("body", body_geometry)
  left = animation_model_node(
    "left-wing",
    cube(
      x = -.55,
      xwidth = 1.1,
      ywidth = .075,
      zwidth = .4,
      material = diffuse("#D5822D")
    ),
    parent = "body",
    transform = animation_transform_track(0, position = c(-.19, .02, 0))
  )
  right = animation_model_node(
    "right-wing",
    cube(
      x = .55,
      xwidth = 1.1,
      ywidth = .075,
      zwidth = .4,
      material = diffuse("#D5822D")
    ),
    parent = "body",
    transform = animation_transform_track(0, position = c(.19, .02, 0))
  )
  # The body advances locally and bobs. In-place removes just x/z locomotion;
  # the multi-key placement paths below remain completely independent.
  clip = animation_model_clip(
    "flight",
    list(
      body = animation_transform_track(
        c(0, .5, 1),
        position = rbind(c(0, 0, 0), c(0, .10, -.2), c(0, 0, -.4)),
        position_interpolation = "cubic"
      ),
      `left-wing` = animation_transform_track(
        c(0, .35, .65, 1),
        angle = rbind(c(0, 0, -32), c(0, 0, 42), c(0, 0, 15), c(0, 0, -32)),
        rotation_interpolation = "smooth"
      ),
      `right-wing` = animation_transform_track(
        c(0, .4, .7, 1),
        angle = rbind(c(0, 0, 32), c(0, 0, -42), c(0, 0, -15), c(0, 0, 32)),
        rotation_interpolation = "smooth"
      )
    ),
    name = "Flapping flight"
  )
  asset = animation_model_asset(
    list(body, left, right),
    list(clip),
    name = "Solid bird"
  )
  first_path = animation_transform_track(
    c(0, .9, 2, 3.1, 4),
    position = rbind(
      c(-1.8, 1.1, .7),
      c(-1.5, 2, -.5),
      c(-.7, 2.2, -1),
      c(-1.1, 1.4, .1),
      c(-1.8, 1.1, .7)
    ),
    angle = rbind(
      c(0, -25, -10),
      c(-12, 20, 18),
      c(5, 60, 8),
      c(12, -10, -18),
      c(0, -25, -10)
    ),
    position_interpolation = "cubic",
    rotation_interpolation = "smooth"
  )
  second_path = animation_transform_track(
    c(0, .7, 1.8, 3, 4),
    position = rbind(
      c(1.6, 1.7, -.8),
      c(1.1, 1.1, .6),
      c(1.4, 1.8, 1),
      c(1.7, 2.4, -.4),
      c(1.6, 1.7, -.8)
    ),
    angle = rbind(
      c(-8, 35, 12),
      c(8, -15, -20),
      c(-10, -35, 5),
      c(0, 30, 20),
      c(-8, 35, 12)
    ),
    position_interpolation = "cubic",
    rotation_interpolation = "smooth"
  )
  first = animation_model_instance(
    asset,
    first_path,
    animation_clip_player(
      "flight",
      loop = "repeat",
      root_motion = "in_place",
      locomotion_root = "body"
    )
  )
  second = animation_model_instance(
    asset,
    second_path,
    animation_clip_player(
      "flight",
      offset = .37,
      speed = 1.4,
      loop = "repeat",
      root_motion = "in_place",
      locomotion_root = "body"
    )
  )
  scene = add_object(first, second) |>
    add_object(cube(
      y = -.18,
      xwidth = 8,
      ywidth = .3,
      zwidth = 6,
      material = diffuse(
        "#DDD5C4",
        checkercolor = "#A5B5B8",
        checkerperiod = .5
      )
    )) |>
    add_object(sphere(
      x = -3,
      y = 6,
      z = 2,
      radius = 1.2,
      material = light("#FFF2DC", intensity = 30)
    )) |>
    add_camera(
      camera(
        lookfrom = c(6.8, 4.5, 8),
        lookat = c(0, 1, 0),
        fov = 38,
        aperture = 0
      ),
      name = "flight"
    )
  list(asset = asset, first = first, second = second, scene = scene)
}
