rigid_model_fixture = function(emissive = FALSE) {
  material = if (emissive) light("orange", intensity = 1) else diffuse("orange")
  body = animation_model_node(
    "body",
    ellipsoid(a = .25, b = .2, c = .45, material = material),
    transform = transform_track(0, position = c(0, .4, 0))
  )
  left = animation_model_node(
    "left",
    cube(x = -.45, xwidth = .9, ywidth = .06, zwidth = .3, material = material),
    parent = "body",
    transform = transform_track(0, position = c(-.2, 0, 0))
  )
  right = animation_model_node(
    "right",
    cube(x = .45, xwidth = .9, ywidth = .06, zwidth = .3, material = material),
    parent = "body",
    transform = transform_track(0, position = c(.2, 0, 0))
  )
  clip = animation_model_clip(
    "flap",
    list(
      body = transform_track(
        c(2, 3, 4),
        position = rbind(c(0, .4, 0), c(.5, .6, .2), c(1, .4, .4)),
        angle = rbind(c(0, 0, 0), c(0, 30, 0), c(0, 60, 0))
      ),
      left = transform_track(
        c(2, 2.5, 4),
        angle = rbind(c(0, 0, -35), c(0, 0, 35), c(0, 0, -35)),
        rotation_interpolation = "smooth"
      ),
      right = transform_track(
        c(2, 3.25, 4),
        angle = rbind(c(0, 0, 35), c(0, 0, -35), c(0, 0, 35)),
        rotation_interpolation = "smooth"
      )
    ),
    name = "Flight"
  )
  animation_model_asset(
    list(right, body, left),
    list(clip),
    name = "Rigid bird"
  )
}

flatten_model_pose = function(instance, time) {
  pose = animation_model_pose(instance, time)[[1]]
  asset = instance$animation_info[[1]]$model_instance$asset
  parts = lapply(asset$nodes, function(node) {
    scene = node$scene
    if (is.null(scene)) {
      return(NULL)
    }
    for (i in seq_len(nrow(scene))) {
      group = scene$transforms[[i]]$group_transform[[1]]
      scene$transforms[[i]]$group_transform[[1]] = pose$world[[node$id]] %*%
        if (nrow(group) == 4L) group else diag(4)
    }
    scene
  })
  Reduce(add_object, Filter(Negate(is.null), parts))
}
