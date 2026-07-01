test_that("animation frame post-processing matches scene orientation", {
  rgb_mat = list(
    r = matrix(1:6, nrow = 2, ncol = 3),
    g = matrix(11:16, nrow = 2, ncol = 3),
    b = matrix(21:26, nrow = 2, ncol = 3)
  )

  frame = post_process_frame(
    rgb_mat = rgb_mat,
    debug_channel = 2,
    filename = "",
    tonemap = "raw",
    bloom = FALSE,
    write_file = FALSE,
    plot_scene = FALSE
  )

  expect_equal(frame[,, 1], fliplr(flipud(t(rgb_mat$r))))
  expect_equal(frame[,, 2], fliplr(flipud(t(rgb_mat$g))))
  expect_equal(frame[,, 3], fliplr(flipud(t(rgb_mat$b))))
})

test_that("render_animation respects denoise = FALSE", {
  skip_on_cran()
  skip_if_not(has_denoiser())

  camera_motion = data.frame(
    x = 278,
    y = 278,
    z = -800,
    dx = 278,
    dy = 278,
    dz = 277.5,
    aperture = 0,
    fov = 40,
    focal = 1077.5,
    orthox = 1,
    orthoy = 1,
    upx = 0,
    upy = 1,
    upz = 0
  )
  scene = generate_cornell()
  raw_prefix = tempfile("animation-raw-")
  denoised_prefix = tempfile("animation-denoised-")

  set.seed(1)
  render_animation(
    scene,
    camera_motion = camera_motion,
    filename = raw_prefix,
    samples = 4,
    width = 8,
    height = 8,
    denoise = FALSE,
    preview = FALSE,
    plot_scene = FALSE,
    parallel = FALSE,
    progress = FALSE,
    bloom = FALSE
  )
  set.seed(1)
  render_animation(
    scene,
    camera_motion = camera_motion,
    filename = denoised_prefix,
    samples = 4,
    width = 8,
    height = 8,
    denoise = TRUE,
    preview = FALSE,
    plot_scene = FALSE,
    parallel = FALSE,
    progress = FALSE,
    bloom = FALSE
  )

  raw_frame = png::readPNG(paste0(raw_prefix, "1.png"))
  denoised_frame = png::readPNG(paste0(denoised_prefix, "1.png"))

  expect_false(identical(raw_frame, denoised_frame))
})
