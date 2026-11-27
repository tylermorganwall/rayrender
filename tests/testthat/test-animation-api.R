test_that("remaining animation aliases retain their implementations", {
  names = c(
    keyframe_track = "animation_keyframe_track",
    transform_track = "animation_transform_track",
    animate_keyframes = "animation_keyframes",
    clip_player = "animation_clip_player",
    read_gltf = "animation_read_gltf",
    list_animations = "animation_list_clips"
  )
  for (old in names(names)) {
    expect_identical(
      getExportedValue("rayrender", old),
      getExportedValue("rayrender", names[[old]])
    )
  }
  asset = animation_model_asset(list(animation_model_node("body", sphere())))
  instance = animation_model_instance(asset)
  expect_equal(animation_model_pose(instance, 1)[[1]]$world$body, diag(4))
  expect_equal(nrow(animation_list_clips(asset)), 0)
})

test_that("model authoring and inspection expose only prefixed names", {
  old = paste0("model_", c("asset", "clip", "instance", "node", "pose"))
  exports = getNamespaceExports("rayrender")
  expect_true(all(paste0("animation_", old) %in% exports))
  expect_false(any(old %in% exports))
  for (name in old) {
    expect_false(exists(
      name,
      envir = asNamespace("rayrender"),
      inherits = FALSE
    ))
  }
})
