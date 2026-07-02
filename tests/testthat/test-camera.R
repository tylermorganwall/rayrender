test_that("camera creates one-frame static motion", {
  cam = camera(
    name = "main",
    lookfrom = c(0, 1, -10),
    lookat = c(0, 0, 0),
    fov = 35
  )

  expect_s3_class(cam, "ray_camera")
  expect_s3_class(cam$motion, "ray_camera_motion")
  expect_s3_class(cam$motion, "data.frame")
  expect_equal(nrow(cam$motion), 1)
  expect_equal(cam$motion$fov, 35)
  expect_equal(
    cam$motion$focal,
    sqrt(sum((c(0, 1, -10) - c(0, 0, 0))^2))
  )
})

test_that("camera preserves generated motion frame count", {
  motion = generate_camera_motion(
    positions = list(c(0, 1, -10), c(1, 2, -8), c(0, 1, -6)),
    lookats = list(c(0, 0, 0), c(0, 0, 0), c(0, 0, 0)),
    frames = 5,
    type = "linear",
    progress = FALSE
  )

  cam = camera(name = "orbit", motion = motion)

  expect_equal(nrow(cam$motion), 5)
})

test_that("generate_camera_motion returns camera motion data frame", {
  motion = generate_camera_motion(
    positions = list(c(0, 1, -10), c(1, 2, -8), c(0, 1, -6)),
    frames = 4,
    type = "linear",
    progress = FALSE
  )

  expect_s3_class(motion, "ray_camera_motion")
  expect_s3_class(motion, "data.frame")
})

test_that("add_camera attaches named active camera", {
  scene = generate_ground() |>
    add_camera(camera(name = "main"))

  expect_named(attr(scene, "ray_cameras"), "main")
  expect_equal(attr(scene, "active_camera"), "main")
})

test_that("add_object preserves cameras and active camera", {
  scene = generate_ground() |>
    add_camera(camera(name = "main"))

  scene = add_object(scene, sphere())

  expect_named(attr(scene, "ray_cameras"), "main")
  expect_equal(attr(scene, "active_camera"), "main")
})

test_that("duplicate camera names error unless replace is TRUE", {
  scene = generate_ground() |>
    add_camera(camera(name = "main"))

  expect_error(
    add_camera(scene, camera(name = "main")),
    "already exists"
  )

  replaced = add_camera(
    scene,
    camera(name = "main", fov = 55),
    replace = TRUE
  )
  expect_equal(get_camera(replaced, "main")$motion$fov, 55)
})

test_that("get_camera returns active camera", {
  scene = generate_ground() |>
    add_camera(camera(name = "wide"), active = FALSE) |>
    add_camera(camera(name = "main"), active = TRUE)

  expect_equal(get_camera(scene)$name, "main")
})

test_that("get_camera errors with multiple cameras and no active camera", {
  scene = generate_ground() |>
    add_camera(camera(name = "wide"), active = FALSE) |>
    add_camera(camera(name = "main"), active = FALSE)

  expect_error(get_camera(scene), "Multiple cameras")
})

test_that("list_cameras reports camera metadata", {
  motion = generate_camera_motion(
    positions = list(c(0, 1, -10), c(1, 2, -8)),
    frames = 3,
    type = "linear",
    progress = FALSE
  )
  scene = generate_ground() |>
    add_camera(camera(name = "wide", filename = "wide.png"), active = FALSE) |>
    add_camera(
      camera(name = "flythrough", motion = motion, filename = "fly_%04d.png"),
      active = TRUE
    )

  cameras = list_cameras(scene)

  expect_equal(cameras$name, c("wide", "flythrough"))
  expect_equal(cameras$frames, c(1L, 3L))
  expect_equal(cameras$animated, c(FALSE, TRUE))
  expect_equal(cameras$filename, c("wide.png", "fly_%04d.png"))
  expect_equal(cameras$active, c(FALSE, TRUE))
})

test_that("camera_frame_filenames expands animation filenames", {
  empty = camera_frame_filenames(NA_character_, 2)
  expect_equal(empty$filenames, c("", ""))
  expect_false(empty$write_image)

  expect_equal(
    camera_frame_filenames("foo", 2)$filenames,
    c("foo1.png", "foo2.png")
  )
  expect_equal(
    camera_frame_filenames("foo.png", 2)$filenames,
    c("foo1.png", "foo2.png")
  )
  expect_equal(
    camera_frame_filenames("foo_%04d.png", 2)$filenames,
    c("foo_0001.png", "foo_0002.png")
  )
  expect_equal(
    camera_frame_filenames("foo_%04d", 2)$filenames,
    c("foo_0001.png", "foo_0002.png")
  )
  expect_equal(
    camera_frame_filenames(c("one.png", "two"), 2)$filenames,
    c("one.png", "two.png")
  )
})

test_that("duplicate multi-camera output filenames are detected", {
  cameras = list(
    one = camera(name = "one", filename = "same.png"),
    two = camera(name = "two", filename = "same.png")
  )

  expect_error(
    validate_camera_output_filenames(cameras, mode = "image"),
    "Duplicate camera output filenames"
  )
})

test_that("camera_batch_plan combines cameras and supports mixed writing", {
  cameras = list(
    wide = camera(name = "wide", filename = "wide.png"),
    preview = camera(name = "preview", filename = NA_character_)
  )

  plan = camera_batch_plan(cameras, mode = "auto")

  expect_s3_class(plan$motion, "ray_camera_motion")
  expect_equal(nrow(plan$motion), 2)
  expect_equal(plan$filenames, c("wide.png", ""))
  expect_true(plan$write_image)
  expect_equal(plan$camera_index, c(1L, 2L))
  expect_equal(plan$camera_names, c("wide", "preview"))
})

test_that("camera_batch_plan combines static and animated cameras", {
  motion = generate_camera_motion(
    positions = list(c(0, 1, -10), c(1, 2, -8), c(0, 1, -6)),
    frames = 3,
    type = "linear",
    progress = FALSE
  )
  cameras = list(
    wide = camera(name = "wide", filename = "wide"),
    orbit = camera(name = "orbit", motion = motion, filename = "orbit_%04d")
  )

  plan = camera_batch_plan(cameras, mode = "auto")

  expect_equal(nrow(plan$motion), 4)
  expect_equal(
    plan$filenames,
    c("wide.png", "orbit_0001.png", "orbit_0002.png", "orbit_0003.png")
  )
  expect_equal(plan$camera_index, c(1L, 2L, 2L, 2L))
  expect_equal(plan$frame_numbers, c(1L, 1L, 2L, 3L))
})

test_that("camera_batch_metadata_compatible detects differing metadata", {
  expect_true(camera_batch_metadata_compatible(list(
    camera(name = "one"),
    camera(name = "two")
  )))
  expect_false(camera_batch_metadata_compatible(list(
    camera(name = "one"),
    camera(name = "two", iso = 200)
  )))
})

test_that("resolver validates filename vectors against attached camera frames", {
  motion = generate_camera_motion(
    positions = list(c(0, 1, -10), c(1, 2, -8)),
    frames = 2,
    type = "linear",
    progress = FALSE
  )
  scene = generate_ground() |>
    add_camera(camera(name = "move", motion = motion), active = TRUE)

  resolved = resolve_scene_camera(scene)
  resolved = lapply(
    resolved,
    apply_camera_overrides,
    overrides = list(filename = c("one.png", "two.png")),
    supplied = c(filename = TRUE)
  )

  expect_equal(
    camera_written_filenames(resolved[[1]], mode = "animation"),
    c("one.png", "two.png")
  )
})

test_that("resolve_scene_camera returns attached active camera", {
  scene = generate_ground() |>
    add_camera(camera(name = "main"))

  resolved = resolve_scene_camera(scene)

  expect_equal(resolved[[1]]$name, "main")
})

test_that("print.ray_camera summarizes static cameras", {
  cam = camera(
    name = "main",
    lookfrom = c(0, 1, -10),
    lookat = c(0, 0, 0),
    fov = 35,
    filename = "render.png"
  )

  output = capture.output(print(cam))

  expect_match(output[1], "ray_camera <main>", fixed = TRUE)
  expect_true(any(grepl("type: perspective static", output, fixed = TRUE)))
  expect_true(any(grepl("frames: 1", output, fixed = TRUE)))
  expect_true(any(grepl("lookfrom: c(0, 1, -10)", output, fixed = TRUE)))
  expect_true(any(grepl("output: render.png", output, fixed = TRUE)))
})

test_that("print.ray_camera summarizes animated cameras", {
  motion = generate_camera_motion(
    positions = list(c(0, 1, -10), c(1, 2, -8), c(0, 1, -6)),
    lookats = list(c(0, 0, 0), c(0, 0, 0), c(0, 0, 0)),
    frames = 5,
    type = "linear",
    progress = FALSE
  )
  cam = camera(name = "orbit", motion = motion, filename = "orbit_%04d.png")

  output = capture.output(print(cam))

  expect_match(output[1], "ray_camera <orbit>", fixed = TRUE)
  expect_true(any(grepl("type: perspective animation", output, fixed = TRUE)))
  expect_true(any(grepl("frames: 5", output, fixed = TRUE)))
  expect_true(any(grepl("last lookfrom:", output, fixed = TRUE)))
  expect_true(any(grepl("output: orbit_%04d.png", output, fixed = TRUE)))
})
