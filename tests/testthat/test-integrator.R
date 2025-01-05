test_that("Test different rendering integrator types, snapshot", {
  set.seed(1)
  types = c("rtiow", "nee", "basic")
  # paths = tempfile(pattern=types, fileext = ".png")
  for(i in seq_along(types)) {
    set.seed(1)
    generate_cornell() |> 
      render_scene(samples=16, denoise = FALSE, preview = FALSE, 
                  integrator_type = types[i])
    expect_snapshot(paths[i], variant = Sys.info()[["sysname"]])
  }
  if(rayrender:::cppdef_HAS_OIDN()) {
    types_oidn = c("rtiow", "nee", "basic")
    # paths_oidn = tempfile(pattern=paste0("denoise_",types_oidn), fileext = ".png")

    for(i in seq_along(types)) {
      set.seed(1)
      generate_cornell() |> 
        render_scene(samples=16, denoise = TRUE, preview = FALSE, 
                     integrator_type = types_oidn[i])
       expect_snapshot(paths_oidn[i], variant = Sys.info()[["sysname"]])
    }
  }
})