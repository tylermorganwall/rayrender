test_that("BVH timing is opt-in, numeric, and reset for each render", {
  render = function(scene) {
    render_scene(
      scene,
      width = 8,
      height = 8,
      samples = 1,
      parallel = FALSE,
      preview = FALSE,
      progress = FALSE,
      plot_scene = FALSE,
      denoise = FALSE
    )
  }
  withr::local_options(rayrender.benchmark_timing = FALSE)
  expect_null(attr(render(sphere()), "bvh_build_seconds"))
  options(rayrender.benchmark_timing = TRUE)
  start = proc.time()[["elapsed"]]
  first = render(sphere())
  elapsed = proc.time()[["elapsed"]] - start
  expect_true(is.finite(attr(first, "bvh_build_seconds")))
  expect_gte(attr(first, "bvh_build_seconds"), 0)
  expect_lte(attr(first, "bvh_build_seconds"), elapsed)
  expect_equal(attr(first, "bvh_build_count"), 1)
  # A triangle mesh builds an internal BVH in addition to the world tree.
  multiple = render(triangle())
  expect_gt(attr(multiple, "bvh_build_count"), 1)
  expect_equal(attr(render(sphere()), "bvh_build_count"), 1)
})
