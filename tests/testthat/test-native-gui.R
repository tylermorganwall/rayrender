test_that("non-native choices never probe the provider", {
  resolve = resolve_native_gui
  environment(resolve) = list2env(
    list(requireNamespace = function(...) stop("Unexpected provider probe")),
    parent = environment(resolve)
  )
  expect_identical(resolve("none", TRUE)$mode, "none")
  expect_identical(resolve("auto", FALSE)$mode, "none")
  expect_identical(resolve("legacy", TRUE)$mode, "legacy")
  expect_identical(resolve("legacy", FALSE)$mode, "legacy")
})

test_that("provider absence falls back only in auto mode", {
  resolve = resolve_native_gui
  environment(resolve) = list2env(
    list(requireNamespace = function(...) FALSE),
    parent = environment(resolve)
  )
  expect_message(result <- resolve("auto", TRUE), "provider_absent")
  expect_identical(result, list(mode = "legacy", api = NULL, fallback = FALSE))
  expect_error(resolve("imgui", TRUE), "provider_absent")
  expect_error(resolve("imgui", FALSE), "provider_absent")
})

test_that("none disables even an explicitly requested interactive preview", {
  calls = list()
  local_mocked_bindings(
    resolve_native_gui = function(gui, preview) {
      calls[[length(calls) + 1L]] <<- list(gui = gui, preview = preview)
      list(mode = "none", api = NULL, fallback = FALSE)
    }
  )
  scene = sphere(material = diffuse(color = "orange"))
  render = function(gui, preview, interactive, deferred_render) {
    set.seed(42)
    render_scene(
      scene,
      width = 8,
      height = 8,
      samples = 2,
      sample_method = "random",
      ambient_light = TRUE,
      lookfrom = c(0, 0, 5),
      lookat = c(0, 0, 0),
      parallel = FALSE,
      progress = FALSE,
      plot_scene = FALSE,
      denoise = FALSE,
      gui = gui,
      preview = preview,
      interactive = interactive,
      deferred_render = deferred_render
    )
  }
  none = render("none", TRUE, TRUE, TRUE)
  batch = render("auto", FALSE, FALSE, FALSE)
  expect_identical(none, batch)
  expect_true(all(is.finite(none)))
  expect_identical(calls[[1]], list(gui = "none", preview = FALSE))
  expect_identical(calls[[2]], list(gui = "auto", preview = FALSE))
})

test_that("native hierarchy preserves nested groups without merging independent groups", {
  inner = group_objects(rbind(sphere(), cube(x = 2)))
  outer = group_objects(rbind(inner, sphere(x = 4)))
  separate = group_objects(sphere(x = 8))
  scene = native_editor_scene(process_scene(rbind(outer, separate))$scene)
  expect_identical(scene$preview_paths[[1]], scene$preview_paths[[2]])
  expect_length(scene$preview_paths[[1]], 2L)
  expect_identical(scene$preview_paths[[3]], scene$preview_paths[[1]][1])
  expect_false(scene$preview_paths[[4]][1] %in% scene$preview_paths[[1]])
})
