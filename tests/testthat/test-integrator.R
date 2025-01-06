test_that("Test different rendering integrator types, snapshot", {
  set.seed(1)
  types = c("rtiow", "nee", "basic")
  for(i in seq_along(types)) {
    set.seed(i)
    expect_snapshot({generate_cornell() |> 
      render_scene(samples=4, denoise = FALSE, preview = FALSE,  width=6, height=6,
                  filename = tempfile(fileext = ".png"), 
                  integrator_type = types[i])}, variant = Sys.info()[["sysname"]])
  }
  if(rayrender:::cppdef_HAS_OIDN()) {
    types_oidn = c("rtiow", "nee", "basic")

    for(i in seq_along(types)) {
      set.seed(i+3)
       expect_snapshot({generate_cornell() |> 
        render_scene(samples=4, denoise = TRUE, preview = FALSE,  width=6, height=6,
                     filename = tempfile(fileext = ".png"),
                     integrator_type = types_oidn[i])}, variant = Sys.info()[["sysname"]])
    }
  }
})