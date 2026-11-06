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

test_that("set_render_settings merges defaults and preserves explicit NULL", {
  original = add_camera(sphere(), camera(iso = 2.5)) |>
    set_render_settings(
      samples = 128,
      width = 1920,
      height = 1080,
      environment_light = "sky.exr"
    )
  updated = original |>
    set_render_settings(samples = 512, environment_light = NULL)
  expect_identical(
    attr(updated, "ray_render_settings"),
    list(samples = 512, width = 1920, height = 1080, environment_light = NULL)
  )
  expect_identical(attr(original, "ray_render_settings")$samples, 128)
  expect_identical(
    attr(original, "ray_render_settings")$environment_light,
    "sky.exr"
  )
  expect_identical(get_camera(updated), get_camera(original))
  expect_identical(set_render_settings(updated), updated)
  expect_error(set_render_settings(list(), samples = 1), "ray_scene")
  expect_error(set_render_settings(sphere(), 512), "must be named")
  expect_error(
    set_render_settings(sphere(), samples = 1, samples = 2),
    "unique"
  )
  expect_error(
    set_render_settings(sphere(), sampels = 1),
    "Unknown render settings: sampels"
  )
})

test_that("updated USD defaults render with explicit overrides", {
  skip_if_not_installed("libopenusd", "26.8.3")
  original_path = tempfile(fileext = ".usdc")
  updated_path = tempfile(fileext = ".usdc")
  on.exit(unlink(c(original_path, updated_path)))
  original = add_camera(sphere(), camera(lookfrom = c(0, 0, 5))) |>
    set_render_settings(
      samples = 1,
      width = 10,
      height = 8,
      preview = FALSE,
      parallel = FALSE,
      denoise = FALSE,
      plot_scene = FALSE,
      progress = FALSE
    )
  write_usd(original, original_path, copy_assets = FALSE)
  updated = read_usd(original_path) |>
    set_render_settings(width = 12, height = 6)
  expect_equal(dim(render_scene(updated))[1:2], c(6, 12))
  expect_equal(dim(render_scene(updated, width = 8))[1:2], c(6, 8))
  write_usd(updated, updated_path, copy_assets = FALSE)
  expect_identical(read_usd(updated_path), updated)
  expect_identical(read_usd(original_path), original)
})
