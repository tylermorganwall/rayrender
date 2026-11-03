test_that("batch binding preserves ordinary rows and scene metadata", {
  a = add_camera(sphere(), camera(lookfrom = c(0, 1, 5)), name = "first")
  b = add_camera(cube(), camera(lookfrom = c(1, 1, 5)), name = "second")
  attr(b, "cornell") = TRUE
  batch = list(a, NULL, b, sphere(x = 2))
  expected = Reduce(add_object, Filter(Negate(is.null), batch))
  expect_identical(do.call(rbind, batch), expected)
  expect_identical(rbind(a), a)
  expect_identical(rbind(a, NULL), a)
  expect_identical(rbind(a[0, ], b), add_object(a[0, ], b))
  expect_error(rbind(a, a), "Duplicate camera")

  # Metadata-only fragments must participate even when they have no rows.
  lights = sphere()[0, ]
  lights = add_infinite_light(
    lights,
    sun_light(
      lat = 0,
      long = 0,
      datetime = as.POSIXct("2025-01-01", tz = "UTC"),
      name = "key"
    )
  )
  points = sphere()[0, ]
  points = add_light(points, point_light(name = "fill"))
  expect_identical(
    rbind(lights, a, points),
    Reduce(add_object, list(lights, a, points))
  )
  expect_error(rbind(a, lights, lights), "Duplicate infinite light")
  expect_error(rbind(a, points, points), "Duplicate point/spot light")
})

test_that("batch binding remaps nested instances in order without mutating inputs", {
  asset = rigid_model_fixture()
  instance = animation_model_instance(asset, animation = clip_player("flap"))
  original_id = model_instance_ids(instance)
  nested = create_instances(instance)
  result = rbind(instance, sphere(), nested, instance, nested)
  ids = model_instance_ids(result)
  expect_length(ids, 4)
  expect_identical(ids[[1]], original_id[[1]])
  expect_equal(anyDuplicated(ids), 0L)
  expect_identical(model_instance_ids(instance), original_id)
  expect_identical(model_instance_ids(nested), original_id)

  # References in animated child nodes must follow the remapped identity.
  check_nodes = function(scene, binding = NULL) {
    for (i in seq_len(nrow(scene))) {
      row_binding = scene$animation_info[[i]]$model_instance
      active = if (is.null(row_binding)) binding else row_binding
      node = scene$animation_info[[i]]$model_node
      if (!is.null(node)) {
        expect_identical(node$instance$id, active$id)
      }
      if (!is.null(row_binding)) {
        expect_identical(row_binding$asset, asset)
      }
      prototype = scene$shape_info[[i]]$shape_properties$original_scene
      if (!is.null(prototype)) check_nodes(prototype[[1]], active)
    }
  }
  check_nodes(result)
})

test_that("batched group transforms preserve composition and input ownership", {
  original = rbind(sphere(), group_objects(cube(), translate = c(1, 2, 3)))
  before = unserialize(serialize(original, NULL))
  pivot = c(1, -2, 3)
  translation = c(-3, 2, 1)
  angle = c(15, 25, 35)
  scale = c(2, .5, 1)
  matrix = generate_translation_matrix(translation) %*%
    generate_translation_matrix(pivot) %*%
    generate_rotation_matrix(angle, c(1, 2, 3)) %*%
    diag(c(scale, 1)) %*%
    generate_translation_matrix(-pivot)
  result = group_objects(
    original,
    pivot_point = pivot,
    translate = translation,
    angle = angle,
    scale = scale
  )
  expect_identical(result$transforms[[1]]$group_transform[[1]], matrix)
  expect_identical(
    result$transforms[[2]]$group_transform[[1]],
    matrix %*% original$transforms[[2]]$group_transform[[1]]
  )
  expect_identical(original, before)
  expect_identical(class(result$transforms), class(original$transforms))
  expect_identical(group_objects(original[0, ]), original[0, ])
})


test_that("scene column binding matches base binding without mutating records", {
  rows = list(sphere(x = 1), cube(y = 2), create_instances(sphere(), x = 1:3))
  before = unserialize(serialize(rows, NULL))
  expect_identical(bind_ray_scene_rows(rows), do.call(rbind.data.frame, rows))
  expect_identical(rows, before)
  expect_identical(
    bind_ray_scene_rows(list(rows[[1]][0, ], rows[[2]])),
    rbind.data.frame(rows[[1]][0, ], rows[[2]])
  )
  expect_identical(bind_ray_scene_rows(list(NULL, rows[[1]])), rows[[1]])
  # Non-scene vector classes retain base coercion and factor handling.
  ordinary = list(data.frame(x = factor("a")), data.frame(x = factor("b")))
  expect_identical(
    bind_ray_scene_rows(ordinary),
    do.call(rbind.data.frame, ordinary)
  )
})
