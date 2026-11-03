test_that("typed USD preserves renderer descriptors and camera animation", {
  skip_if_not_installed("libopenusd", "26.8.2")
  skip_if_not_installed("digest")
  directory = tempfile()
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  keys = as.data.frame(camera()$motion)
  keys$exposure = 1.4
  keys$env_rotation = 0.3
  scene = structure(
    add_camera(
      add_object(
        sphere(material = diffuse(color = "red")),
        cube(x = 2, xwidth = 2)
      ),
      camera(keyframes = keys, iso = 5)
    ),
    ray_render_settings = list(
      width = 24,
      height = 16,
      samples = 2,
      environment_light = NULL
    )
  )
  for (extension in c("usda", "usdc")) {
    file = file.path(directory, paste0("scene.", extension))
    write_usd(scene, file)
    expect_identical(read_usd(file), scene)
    expect_error(write_usd(scene, file), "Destination exists")
  }
  value = list(
    numbers = c(NA_real_, NaN, Inf, -Inf),
    integers = c(1L, NA_integer_),
    text = c("a", NA_character_),
    logical = c(TRUE, NA),
    empty = character(),
    raw = as.raw(0:255),
    matrix = matrix(1:6, 2),
    date = as.POSIXct("2025-01-01", tz = "UTC")
  )
  metadata_file = file.path(directory, "metadata.usdc")
  metadata_source = file.path(directory, "metadata.usda")
  usd_stream_native(
    list(
      "#usda 1.0\ndef Xform \"World\" (customData = { dictionary probe = ",
      usd_metadata(value),
      "} ) {}"
    ),
    metadata_source
  )
  usd_export_file_native(metadata_source, metadata_file)
  metadata = usd_inspect_native(metadata_file, "/World", character())[[
    1
  ]]$custom_data$probe
  expect_identical(usd_unpack_value(metadata), value)
  expect_error(
    usd_stream_native(usd_metadata(new.env()), metadata_source),
    "cannot preserve"
  )
})

test_that("USD collects images and NanoVDB files without embedding them in R", {
  skip_if_not_installed("libopenusd", "26.8.2")
  skip_if_not_installed("digest")
  source = tempfile()
  dir.create(source)
  bundle = tempfile()
  dir.create(bundle)
  on.exit(unlink(c(source, bundle), recursive = TRUE))
  image = file.path(source, "texture.png")
  png::writePNG(array(0.5, c(4, 4, 3)), image)
  volume = file.path(source, "density.nvdb")
  file.copy(test_path("fixtures/volumes/transformed-blosc.nvdb"), volume)
  scene = add_camera(
    add_object(
      sphere(material = diffuse(image_texture = image)),
      set_medium(cube(x = 3), vdb_medium(volume))
    ),
    camera()
  )
  file = file.path(bundle, "scene.usdc")
  original_scene = scene
  expect_warning(write_usd(scene, file), "NanoVDB")
  expect_identical(scene, original_scene)
  unlink(source, recursive = TRUE)
  restored = read_usd(file)
  expect_true(file.exists(restored$material[[1]]$image))
  expect_true(file.exists(restored$shape_info[[2]]$medium$filename))
  expect_true(startsWith(
    restored$shape_info[[2]]$medium$filename,
    normalizePath(bundle)
  ))
})

test_that("OBJ instances collect relative material dependencies", {
  skip_if_not_installed("libopenusd", "26.8.2")
  skip_if_not_installed("digest")
  source = tempfile()
  dir.create(source)
  bundle = tempfile()
  dir.create(bundle)
  on.exit(unlink(c(source, bundle), recursive = TRUE))
  png::writePNG(array(0.4, c(2, 2, 3)), file.path(source, "paint.png"))
  writeLines(
    c("newmtl paint", "Kd 1 0 0", "map_Kd paint.png"),
    file.path(source, "model.mtl")
  )
  obj = file.path(source, "model.obj")
  writeLines(
    c(
      "mtllib model.mtl",
      "v 0 0 0",
      "v 1 0 0",
      "v 0 1 0",
      "usemtl paint",
      "f 1 2 3"
    ),
    obj
  )
  scene = structure(
    add_camera(
      create_instances(obj_model(obj), x = c(-1, 1)),
      camera(lookfrom = c(0, 0, 5), aperture = 0)
    ),
    ray_render_settings = list(
      width = 12,
      height = 8,
      samples = 1,
      preview = FALSE,
      denoise = FALSE,
      parallel = FALSE,
      plot_scene = FALSE,
      progress = FALSE
    )
  )
  file = file.path(bundle, "scene.usdc")
  expect_warning(write_usd(scene, file), "OBJ material")
  unlink(source, recursive = TRUE)
  restored = read_usd(file)
  prototype = restored$shape_info[[1]]$shape_properties$original_scene[[1]]
  expect_true(file.exists(prototype$shape_info[[1]]$fileinfo))
  expect_identical(
    prototype$shape_info[[1]]$shape_properties$basename,
    dirname(prototype$shape_info[[1]]$fileinfo)
  )
  expect_equal(length(list.files(file.path(bundle, "scene_assets"))), 3)
  expect_equal(dim(render_scene(restored))[1:2], c(8, 12))
})

test_that("USD render round trips preserve images and explicit settings", {
  skip_if_not_installed("libopenusd", "26.8.2")
  skip_if_not_installed("digest")
  directory = tempfile()
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  scene = structure(
    add_camera(
      sphere(material = diffuse(color = "purple")),
      camera(lookfrom = c(0, 1, 5), aperture = 0)
    ),
    ray_render_settings = list(
      width = 16,
      height = 12,
      samples = 1,
      preview = FALSE,
      denoise = FALSE,
      parallel = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      debug_channel = "normals"
    )
  )
  file = file.path(directory, "scene.usdc")
  write_usd(scene, file)
  restored = read_usd(file)
  set.seed(1)
  before = render_scene(scene)
  set.seed(1)
  after = render_scene(restored)
  expect_equal(after, before, tolerance = 1e-07)
  expect_equal(dim(render_scene(restored, width = 20))[1:2], c(12, 20))
})

test_that("USD import rejects foreign stages instead of inventing renderer data", {
  skip_if_not_installed("libopenusd", "26.8.2")
  skip_if_not_installed("digest")
  file = tempfile(fileext = ".usda")
  on.exit(unlink(file))
  writeLines(
    c("#usda 1.0", "(defaultPrim = \"World\")", "def Xform \"World\" {}"),
    file
  )
  expect_error(read_usd(file), "rayrender USD stage")
})

test_that("native export validates assets and preserves an existing stage on failure", {
  skip_if_not_installed("libopenusd", "26.8.2")
  skip_if_not_installed("digest")
  directory = tempfile()
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  file = file.path(directory, "scene.usdc")
  scene = sphere()
  write_usd(scene, file)
  before = unname(tools::md5sum(file))
  invalid = scene
  attr(invalid, "unsupported") = new.env()
  expect_error(write_usd(invalid, file, overwrite = TRUE), "cannot preserve")
  expect_identical(unname(tools::md5sum(file)), before)
  write_usd(cube(), file, overwrite = TRUE)
  expect_identical(read_usd(file), cube())
  source = file.path(directory, "texture \"quoted\".png")
  png::writePNG(array(0.4, c(2, 2, 3)), source)
  textured = sphere(material = diffuse(image_texture = source))
  write_usd(textured, file, copy_assets = FALSE, overwrite = TRUE)
  expect_identical(
    read_usd(file)$material[[1]]$image,
    normalizePath(source, winslash = "/")
  )
  unlink(source)
  expect_error(read_usd(file), "Missing USD asset")
})

test_that("native metadata preserves precision, strings, attributes and empty columns", {
  skip_if_not_installed("libopenusd", "26.8.2")
  skip_if_not_installed("digest")
  file = tempfile(fileext = ".usdc")
  on.exit(unlink(file))
  scene = sphere()
  attr(scene, "edge_cases") = list(
    precise = c(pi, 1 + .Machine$double.eps, .Machine$double.xmin),
    text = c("quote\" and slash\\", "line\nnext\ttab", "café"),
    empty = list(logical(), integer(), numeric(), character(), raw(), list()),
    named = list(`odd key : with space` = matrix(1:6, 2)),
    values = c(NA_real_, NaN, Inf, -Inf)
  )
  write_usd(scene, file)
  expect_identical(read_usd(file), scene)
  empty = scene[0, ]
  write_usd(empty, file, overwrite = TRUE)
  expect_identical(read_usd(file), empty)
})

test_that("USD preserves column names and custom attributes", {
  skip_if_not_installed("libopenusd", "26.8.2")
  scene = add_object(sphere(), cube(x = 2))
  names(scene$x) = c("first", "second")
  names(scene$material) = c("one", "two")
  attr(scene$x, "units_note") = "scene units"
  path = tempfile(fileext = ".usdc")
  on.exit(unlink(path))
  write_usd(scene, path, copy_assets = FALSE)
  expect_identical(read_usd(path), scene)
})

test_that("USD streams nested geometry and metadata without string assembly", {
  skip_if_not_installed("libopenusd", "26.8.2")
  skip_if_not_installed("digest")
  directory = tempfile()
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  scene = add_object(sphere(), cube(x = 2))
  scene$shape_info[[1]]$probe = list(
    numbers = rep(c(pi, NA_real_, NaN, Inf, -Inf), 1000),
    matrix = matrix(seq_len(3000), ncol = 3)
  )
  for (extension in c("usda", "usdc")) {
    path = file.path(directory, paste0("streamed.", extension))
    write_usd(scene, path, copy_assets = FALSE)
    expect_identical(read_usd(path), scene)
  }
  expect_setequal(
    list.files(directory, all.files = TRUE, no.. = TRUE),
    c("streamed.usda", "streamed.usdc")
  )
})

test_that("a failed streamed export leaves the previous stage and no scratch files", {
  skip_if_not_installed("libopenusd", "26.8.2")
  skip_if_not_installed("digest")
  directory = tempfile()
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  path = file.path(directory, "scene.usdc")
  scene = sphere()
  write_usd(scene, path, copy_assets = FALSE)
  invalid = scene
  invalid$shape_info[[1]]$invalid = new.env()
  expect_error(
    write_usd(invalid, path, copy_assets = FALSE, overwrite = TRUE),
    "cannot preserve"
  )
  expect_identical(read_usd(path), scene)
  expect_identical(
    list.files(directory, all.files = TRUE, no.. = TRUE),
    "scene.usdc"
  )
})

test_that("streamed numeric properties reject missing values", {
  path = tempfile(fileext = ".usda")
  on.exit(unlink(path))
  for (value in list(NA_integer_, NA, NA_real_, Inf)) {
    expect_error(
      usd_stream_native(usd_property("int[]", "probe", value), path),
      "Missing|Nonfinite"
    )
  }
})

test_that("USD preserves instance, texture-graph, and environment assets", {
  skip_if_not_installed("libopenusd", "26.8.2")
  skip_if_not_installed("digest")
  source = tempfile(fileext = ".png")
  bundle = tempfile()
  dir.create(bundle)
  on.exit(unlink(c(source, bundle), recursive = TRUE))
  png::writePNG(array(.5, c(4, 8, 3)), source)
  prototype = sphere(material = diffuse(image_texture = source))
  scene = create_instances(prototype, x = c(-1, 1)) |>
    add_object(sphere(
      y = -2,
      material = diffuse(color = texture_image(source))
    )) |>
    add_infinite_light(infinite_light(source)) |>
    add_camera(camera(lookfrom = c(0, 1, 6), aperture = 0))
  attr(scene, "ray_render_settings") = list(
    environment_light = source,
    width = 12,
    height = 8,
    samples = 1,
    preview = FALSE,
    parallel = FALSE,
    denoise = FALSE,
    plot_scene = FALSE,
    progress = FALSE
  )
  path = file.path(bundle, "scene.usdc")
  expect_warning(write_usd(scene, path), "Procedural materials")
  expect_length(list.files(file.path(bundle, "scene_assets")), 1L)
  unlink(source)
  restored = read_usd(path)
  expect_true(file.exists(
    attr(restored, "ray_render_settings")$environment_light
  ))
  expect_true(file.exists(attr(restored, "ray_infinite_lights")[[1]]$filename))
  expect_equal(dim(render_scene(restored, width = 16))[1:2], c(8, 16))
  testthat::local_mocked_bindings(
    render_animation_camera = function(...) list(...),
    .package = "rayrender"
  )
  args = render_animation(restored, width = 20)
  expect_equal(args$width, 20)
  expect_equal(args$height, 8)
  expect_true(file.exists(args$environment_light))
})
