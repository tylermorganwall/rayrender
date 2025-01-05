

test_that("Test rendering with and without OIDN", {
  set.seed(1)
  rgb_output = generate_cornell() |> 
    render_scene(samples=16, denoise = FALSE, preview = FALSE,
      print_debug_info = TRUE, return_raw_array = TRUE)
  if(abs(sum(rgb_output) - 232160.35) > 100) {
    warning(sprintf("Sum of rgb_output is %f", sum(rgb_output)))
  }
  expect_true(abs(sum(rgb_output) - 232160.35) < 100)
  if(rayrender:::cppdef_HAS_OIDN()) {
    set.seed(1)
    rgb_output_oidn = generate_cornell() |> 
      render_scene(samples=16, preview = FALSE, return_raw_array = TRUE)
    if(abs(sum(rgb_output_oidn) - 235740.581) > 100) {
      warning(sprintf("Sum of rgb_output_oidn is %f", sum(rgb_output_oidn)))
    }
    expect_true(abs(sum(rgb_output_oidn) - 235740.581) < 100)
  }
})

test_that("Test rendering with and without OIDN, snapshot", {
  set.seed(1)
  path = tempfile(fileext = ".png")
  generate_cornell() |> 
    render_scene(samples=16, denoise = FALSE, 
      preview = FALSE, filename=path)
  expect_snapshot_file(path, variant = Sys.info()[["sysname"]])
  if(rayrender:::cppdef_HAS_OIDN()) {
    set.seed(1)
    path_oidn = tempfile(pattern="oidn", fileext = ".png")

    generate_cornell() |> 
      render_scene(samples=16, preview = FALSE, filename=path_oidn)
    expect_snapshot_file(path_oidn, variant = Sys.info()[["sysname"]])
  }
})
