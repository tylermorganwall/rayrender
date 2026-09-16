test_that("outer groups retain identity without changing the random stream", {
  set.seed(14)
  before = .Random.seed
  one = group_objects(add_object(sphere(x = -1), sphere(x = 1)))
  two = group_objects(add_object(sphere(x = -1), sphere(x = 1)))
  joined = add_object(one, two)
  prepared = native_editor_scene(process_scene(joined)$scene)
  expect_identical(prepared$preview_groups, c(1L, 1L, 2L, 2L))
  outer = group_objects(joined)
  expect_identical(
    native_editor_scene(process_scene(outer)$scene)$preview_groups,
    rep(1L, 4)
  )
  expect_identical(.Random.seed, before)
  expect_null(joined$preview_groups)
})

test_that("editor preparation recursively isolates materials and preserves input", {
  child = sphere(material = diffuse(color = "red"))
  child$shape_info[[1]]$material_id = 4L
  input = process_scene(create_instances(child, x = c(-2, 2)))$scene
  before = input$shape_info[[1]]$shape_properties$original_scene[[1]]
  prepared = native_editor_scene(input)
  after = prepared$shape_info[[1]]$shape_properties$original_scene[[1]]
  expect_true(is.na(after$shape_info[[1]]$material_id))
  expect_identical(before$shape_info[[1]]$material_id, 4L)
  expect_identical(
    input$shape_info[[1]]$shape_properties$original_scene[[1]],
    before
  )
  expect_identical(after$preview_groups, 0L)
})
