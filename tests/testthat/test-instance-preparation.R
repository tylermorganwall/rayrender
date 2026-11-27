test_that("preparation preserves unchanged nested prototype identity", {
  skip_if_not(capabilities("profmem"))
  for (leaf in list(
    sphere(),
    set_medium(sphere(), homogeneous_medium(sigma_a = .1))
  )) {
    nested = create_instances(create_instances(leaf))
    scene = add_object(nested, nested)
    before = scene$shape_info[[1]]$shape_properties$original_scene[[1]]
    identity = tracemem(before)
    prepared = prepare_subsurface(scene)
    expect_identical(prepared, scene)
    for (i in 1:2) {
      after = prepared$shape_info[[i]]$shape_properties$original_scene[[1]]
      expect_identical(tracemem(after), identity)
      untracemem(after)
    }
    untracemem(before)
  }
})

test_that("nested changed subsurface interiors are still prepared", {
  scene = create_instances(create_instances(sphere()))
  child = scene$shape_info[[1]]$shape_properties$original_scene[[1]]
  leaf = child$shape_info[[1]]$shape_properties$original_scene[[1]]
  leaf$material[[1]] = subsurface(sigma_a = .1, sigma_s = 2)[[1]]
  child$shape_info[[1]]$shape_properties$original_scene[[1]] = leaf
  scene$shape_info[[1]]$shape_properties$original_scene[[1]] = child
  prepared = prepare_subsurface(scene)
  prepared_child = prepared$shape_info[[1]]$shape_properties$original_scene[[1]]
  prepared_leaf = prepared_child$shape_info[[
    1
  ]]$shape_properties$original_scene[[1]]
  expect_identical(prepared_leaf$shape_info[[1]]$medium_owner, "subsurface")
  expect_null(leaf$shape_info[[1]]$medium)
  expect_identical(prepare_subsurface(prepared), prepared)
  leaf$shape_info[[1]]$medium = homogeneous_medium()
  child$shape_info[[1]]$shape_properties$original_scene[[1]] = leaf
  scene$shape_info[[1]]$shape_properties$original_scene[[1]] = child
  expect_error(prepare_subsurface(scene), "explicit medium")
})
