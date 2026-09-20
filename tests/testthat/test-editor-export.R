export_test_state = function() {
  list(
    version = 1L,
    objects = list(),
    sky = NULL,
    camera = list(
      x = 4,
      y = 2,
      z = -7,
      dx = 0,
      dy = 1,
      dz = 0,
      upx = 0,
      upy = 1,
      upz = 0,
      fov = 55,
      aperture = 0,
      focal = 8,
      orthox = 3,
      orthoy = 2,
      env_rotation = 15
    ),
    denoise = FALSE,
    exposure = 0.4,
    camera_motion_blur = FALSE,
    shutter_speed = 2,
    integrator_type = "nee"
  )
}

export_test_args = function() {
  list(
    width = 12,
    height = 8,
    samples = 2,
    min_variance = 0,
    rotate_env = 20,
    tonemap = "raw",
    sample_method = "random",
    parallel = FALSE,
    progress = FALSE,
    ambient_light = TRUE
  )
}

test_that("exports freeze the camera, exposure and output settings", {
  state = export_test_state()
  args = native_export_render_args(export_test_args(), state)
  expect_equal(args$lookfrom, c(4, 2, -7))
  expect_equal(args$lookat, c(0, 1, 0))
  expect_equal(args$camera_up, c(0, 1, 0))
  expect_equal(args$rotate_env, 35)
  expect_equal(native_export_render_args(list(), state)$rotate_env, 15)
  expect_equal(args$exposure, 0.4)
  expect_false(args$auto_exposure)
  expect_false(args$bloom)
  expect_false(args$denoise)
  expect_identical(args$gui, "none")
  expect_identical(args$integrator_type, "nee")
})

test_that("exported scripts run without the original workspace and preserve temporary maps", {
  directory = tempfile("scene export '")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  texture = tempfile(fileext = ".png")
  png::writePNG(array(.4, c(4, 4, 3)), texture)
  original = group_objects(sphere(
    material = microfacet(roughness_texture = texture)
  ))
  state = export_test_state()
  filename = file.path(directory, "scene 'view.R")
  script = write_native_scene_export(
    original,
    export_test_args(),
    state,
    filename,
    getwd()
  )
  expect_true(file.exists(script))
  expect_length(list.files(directory, pattern = "rds$"), 1L)
  saved = readRDS(file.path(directory, "scene 'view-scene.rds"))
  expect_false(identical(
    saved$material[[1]]$roughness_texture,
    original$material[[1]]$roughness_texture
  ))
  unlink(texture)
  # Evaluate the exact generated script with a render spy to inspect its complete
  # inputs, then render those inputs normally after leaving the export directory.
  captured = NULL
  local_mocked_bindings(render_scene = function(scene, ...) {
    captured <<- list(scene = scene, args = list(...), directory = getwd())
    invisible(NULL)
  })
  previous = getwd()
  expect_no_error(source(script, local = new.env(parent = baseenv())))
  expect_identical(getwd(), previous)
  expect_equal(captured$args$lookfrom, c(4, 2, -7))
  expect_equal(attr(captured$scene, "rayrender_scene_edits"), list())
  expect_true(file.exists(file.path(
    captured$directory,
    saved$material[[1]]$roughness_texture
  )))
  # A fresh R process sources the export after the original texture is gone.
  check = file.path(directory, "check.R")
  writeLines(
    c(
      paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
      "pdf(NULL)",
      paste0("result = source(", encodeString(script, quote = '"'), ")$value"),
      "stopifnot(all(is.finite(result)), max(result[,,1:3]) > 0)"
    ),
    check
  )
  output = system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(check)),
    stdout = TRUE,
    stderr = TRUE
  )
  expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))
  second = write_native_scene_export(
    sphere(),
    export_test_args(),
    state,
    filename,
    getwd()
  )
  expect_false(identical(script, second))
  expect_true(file.exists(script))
})

test_that("headless rendering replays transform edits and rejects invalid paths", {
  transform = diag(4)
  transform[1, 4] = 1.5
  edits = list(list(
    row = 1L,
    instance = 1L,
    transform = transform,
    materials = list()
  ))
  original = sphere(material = diffuse(color = "red"))
  changed = apply_scene_edits(original, edits)
  expect_null(attr(original, "rayrender_scene_edits"))
  args = list(
    width = 12,
    height = 8,
    samples = 2,
    gui = "none",
    plot_scene = FALSE,
    lookfrom = c(3, 2, -6),
    lookat = c(0, 0, 0),
    fov = 55,
    aperture = 0,
    parallel = FALSE,
    progress = FALSE,
    bloom = FALSE,
    denoise = FALSE,
    ambient_light = TRUE,
    min_variance = 0,
    sample_method = "random"
  )
  set.seed(234)
  replay = do.call(render_scene, c(list(scene = changed), args))
  set.seed(234)
  expected = do.call(
    render_scene,
    c(list(scene = sphere(x = 1.5, material = diffuse(color = "red"))), args)
  )
  expect_equal(as.numeric(replay), as.numeric(expected), tolerance = 1e-6)
  expect_equal(attr(replay, "scene_edits")[[1]]$transform, transform)
  edits[[1]]$row = 2L
  expect_error(
    do.call(
      render_scene,
      c(list(scene = apply_scene_edits(original, edits)), args)
    ),
    "Invalid editor index"
  )
  edits[[1]]$row = 1L
  edits[[1]]$transform[1, 1] = 0
  expect_error(
    do.call(
      render_scene,
      c(list(scene = apply_scene_edits(original, edits)), args)
    ),
    "singular|nonzero"
  )
})

test_that("saved manual sky direction replays without an editor", {
  sky = sky_light_image(
    40,
    -74,
    as.POSIXct("2026-06-21 16:00:00", tz = "UTC"),
    resolution = 16,
    moon = FALSE
  )
  state = list(
    model = 0L,
    latitude = 35,
    longitude = -80,
    datetime = "2026-06-22 15:00:00",
    manual = TRUE,
    elevation = 30,
    azimuth = 135,
    haze = TRUE,
    altitude = TRUE
  )
  restored = native_sky_restore(list(sky), state)
  expect_equal(restored$controls$elevation, 30)
  expect_true(restored$controls$manual)
  expect_identical(restored$lights[[1]]$type, "image")
  pixels = rayimage::ray_read_image(restored$lights[[1]]$filename)
  expect_true(all(is.finite(pixels)))
  expect_gt(max(pixels[,, 1:3]), 0)
})

test_that("temporary OBJ exports preserve material and texture dependencies", {
  source = tempfile("model-")
  destination = tempfile("assets-")
  dir.create(source)
  on.exit(unlink(c(source, destination), recursive = TRUE))
  texture = file.path(source, "color.png")
  png::writePNG(array(.5, c(3, 3, 3)), texture)
  writeLines(
    c("newmtl surface", "map_Kd -s 1 1 1 color.png"),
    file.path(source, "surface.mtl")
  )
  model = file.path(source, "model.obj")
  writeLines(
    c(
      "mtllib surface.mtl",
      "v 0 0 0",
      "v 1 0 0",
      "v 0 1 0",
      "usemtl surface",
      "f 1 2 3"
    ),
    model
  )
  paths = native_export_assets(
    list(model = model),
    destination,
    basename(destination),
    getwd()
  )
  unlink(source, recursive = TRUE)
  copied_model = file.path(dirname(destination), paths$model)
  material = sub("^mtllib ", "", readLines(copied_model)[1L])
  lines = readLines(file.path(destination, material))
  expect_match(lines[2L], "map_Kd -s 1 1 1")
  copied_texture = sub("^map_Kd -s 1 1 1 ", "", lines[2L])
  expect_true(file.exists(file.path(destination, copied_texture)))
})
