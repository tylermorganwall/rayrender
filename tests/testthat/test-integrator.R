test_that("Test different rendering integrator types, snapshot", {
  types = c("rtiow", "nee", "basic")
  for(i in seq_along(types)) {
    set.seed(i)
    expect_snapshot_output({generate_cornell() |> 
      render_scene(samples=4, denoise = FALSE, preview = FALSE,  width=6, height=6,
                  filename = tempfile(fileext = ".png"), parallel = FALSE,
                  integrator_type = types[i])}, variant = Sys.info()[["sysname"]])
  }
  for(i in seq_along(types)) {
    set.seed(i)
    expect_snapshot_output({generate_cornell() |> 
      render_scene(samples=4, denoise = TRUE, preview = FALSE,  width=6, height=6,
                  filename = tempfile(fileext = ".png"), parallel = FALSE,
                  integrator_type = types[i])}, variant = Sys.info()[["sysname"]])
  }
})