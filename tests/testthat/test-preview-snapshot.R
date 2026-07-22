test_that("preview snapshots use the first available default number", {
  snapshot_dir = withr::local_tempdir()
  withr::local_dir(snapshot_dir)

  expect_equal(
    next_preview_snapshot_filename(),
    "rayrender_snapshot1.png"
  )

  file.create("rayrender_snapshot1.png")
  file.create("rayrender_snapshot3.png")

  expect_equal(
    next_preview_snapshot_filename(),
    "rayrender_snapshot2.png"
  )
})

test_that("preview snapshots insert a number before the source extension", {
  snapshot_dir = withr::local_tempdir()
  source_filename = file.path(snapshot_dir, "render.final.jpg")
  file.create(file.path(snapshot_dir, "render.final1.jpg"))
  file.create(file.path(snapshot_dir, "render.final2.jpg"))

  expect_equal(
    next_preview_snapshot_filename(source_filename),
    file.path(snapshot_dir, "render.final3.jpg")
  )
})

test_that("preview snapshots use the local default without a source extension", {
  snapshot_dir = withr::local_tempdir()
  withr::local_dir(snapshot_dir)

  expect_equal(
    next_preview_snapshot_filename("render"),
    "rayrender_snapshot1.png"
  )
})

test_that("preview snapshots write the current preview image", {
  snapshot_dir = withr::local_tempdir()
  source_filename = file.path(snapshot_dir, "render.png")
  preview_image = array(0.5, dim = c(2, 3, 3))

  expect_message(
    {
      saved_filename = save_preview_snapshot(preview_image, source_filename)
    },
    "Saved preview snapshot:"
  )

  expect_equal(saved_filename, file.path(snapshot_dir, "render1.png"))
  expect_true(file.exists(saved_filename))
  expect_equal(dim(png::readPNG(saved_filename)), c(2, 3, 3))
})
