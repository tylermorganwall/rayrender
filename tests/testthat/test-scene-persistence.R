test_that("saved defaults respect explicit values including NULL", {
  scene = structure(
    sphere(),
    ray_render_settings = list(
      width = 10,
      height = 20,
      environment_light = "sky.exr"
    )
  )
  probe = function(scene, width = 1, height = 2, environment_light = NULL) {
    call = saved_scene_render_call(
      scene,
      names(as.list(match.call())[-1]),
      environment(),
      probe
    )
    if (!is.null(call)) {
      return(do.call(probe, call))
    }
    list(width = width, height = height, environment_light = environment_light)
  }
  expect_identical(
    probe(scene, width = 50, environment_light = NULL),
    list(width = 50, height = 20, environment_light = NULL)
  )
  expect_identical(
    attr(
      structure(scene, ray_render_settings = list(environment_light = NULL)),
      "ray_render_settings"
    )$environment_light,
    NULL
  )
})

test_that("camera keyframes preserve editorial fields and reject invalid inputs", {
  keys = camera()$motion
  keys$exposure = 1.4
  keys$env_rotation = 0.7
  original = get("keyframes", envir = ray_environment)
  cam = camera(
    keyframes = keys,
    keyframe_motion_args = list(frames = 120, closed = TRUE)
  )
  expect_equal(cam$keyframes, as.data.frame(keys))
  expect_equal(cam$keyframe_motion_args$frames, 120)
  expect_true(cam$keyframe_motion_args$closed)
  expect_identical(get("keyframes", envir = ray_environment), original)
  keys$x = NA_real_
  expect_error(camera(keyframes = keys), "finite numeric")
})

test_that("saved keyframes reach the native preview buffer", {
  original = get("keyframes", envir = ray_environment)
  on.exit(assign("keyframes", original, envir = ray_environment))
  keys = as.data.frame(camera()$motion)
  keys$exposure = 1.3
  keys$env_rotation = 0.6
  scene = add_camera(sphere(), camera(keyframes = keys))
  expect_message(
    render_scene(
      scene,
      width = 8,
      height = 8,
      samples = 1,
      preview = FALSE,
      parallel = FALSE,
      denoise = FALSE,
      plot_scene = FALSE,
      progress = FALSE
    ),
    "Saving camera keyframes"
  )
  expect_equal(get_saved_keyframes()[names(keys)], keys)
})

test_that("scene assembly and grouping preserve render settings", {
  a = structure(
    sphere(),
    ray_render_settings = list(width = 10, environment_light = NULL)
  )
  b = structure(cube(), ray_render_settings = list(width = 20, height = 30))
  scene = add_object(a, b)
  expect_identical(
    attr(scene, "ray_render_settings"),
    list(width = 10, environment_light = NULL, height = 30)
  )
  expect_identical(
    attr(group_objects(scene), "ray_render_settings"),
    attr(scene, "ray_render_settings")
  )
})
