

test_that("Test rendering with and without OIDN", {
  set.seed(1)
  generate_cornell() |> 
    render_scene(samples=16, denoise=F, print_debug_info = TRUE) ->
  rgb_output
  if(abs(sum(rgb_output) - 232160.35) > 0.01) {
    message(sprintf("Sum of rgb_output is %f", sum(rgb_output)))
  }
  expect_true(abs(sum(rgb_output) - 232160.35) < 0.01)
  if(cppdef_HAS_OIDN()()) {
    set.seed(1)
    generate_cornell() |> 
      render_scene(samples=16) ->
    rgb_output_oidn
    generate_cornell() |> 
      render_scene(samples=16, preview = FALSE) ->
    rgb_output_oidn
    if(abs(sum(rgb_output_oidn) - 235740.581) > 0.01) {
      message(sprintf("Sum of rgb_output_oidn is %f", sum(rgb_output_oidn)))
    }
    expect_true(abs(sum(rgb_output_oidn) - 235740.581) < 0.01)
  }
})
