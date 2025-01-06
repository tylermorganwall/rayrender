test_that("Test different rendering integrator types, snapshot", {
  set.seed(1)
  types = c("rtiow", "nee", "basic")
  for(i in seq_along(types)) {
    set.seed(i)
    expect_snapshot_output({generate_cornell() |> 
      render_scene(samples=1, denoise = FALSE, preview = FALSE,  width=6, height=6,
                  filename = tempfile(fileext = ".png"), parallel = FALSE,
                  integrator_type = types[i])}, variant = Sys.info()[["sysname"]])
  }
})