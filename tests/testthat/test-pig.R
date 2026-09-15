pig_test_triangles = function(scene) {
  sum(vapply(
    scene$shape_info[scene$shape == "mesh3d"],
    function(s) nrow(s$mesh_info[[1]]$indices),
    integer(1)
  ))
}

test_that("pig preserves the legacy positional arguments and transform pivot", {
  expect_identical(
    names(formals(pig))[1:9],
    c(
      "x",
      "y",
      "z",
      "emotion",
      "spider",
      "angle",
      "order_rotation",
      "scale",
      "diffuse_sigma"
    )
  )
  scene = pig(
    2,
    3,
    -4,
    "worried",
    FALSE,
    c(20, 40, -10),
    c(3, 1, 2),
    c(.5, 2, 3),
    15
  )
  reference = group_objects(
    sphere(),
    translate = c(2, 3, -4),
    angle = c(20, 40, -10),
    order_rotation = c(3, 1, 2),
    scale = c(.5, 2, 3),
    pivot_point = c(0, 1, 0)
  )
  expect_true(all(vapply(
    scene$transforms,
    function(t) identical(t, reference$transforms[[1]]),
    logical(1)
  )))
  expect_equal(pig(scale = 2), pig(scale = c(2, 2, 2)))
})

test_that("pig assembles the compact classic, ski and spider assets", {
  expect_equal(pig_test_triangles(pig()), 3299)
  expect_equal(pig_test_triangles(pig(ski = TRUE)), 5212)
  expect_equal(pig_test_triangles(pig(spider = TRUE, hair_count = 0)), 4275)
  expect_equal(pig(ski = TRUE), pig(ski = TRUE, emotion = "excited"))
  expect_equal(pig(emotion = "neutral"), pig(emotion = "cheerful"))
  expect_equal(
    pig(ski = TRUE, emotion = "neutral")$shape,
    pig(ski = TRUE, emotion = "cheerful")$shape
  )
  for (ski in c(FALSE, TRUE)) {
    for (emotion in c(
      "neutral",
      "skeptical",
      "worried",
      "angry",
      "surprised",
      "excited"
    )) {
      meshes = lapply(
        pig(ski = ski, emotion = emotion)$shape_info,
        function(s) s$mesh_info[[1]]
      )
      expect_true(all(vapply(
        meshes,
        function(m) all(is.finite(m$vertices)),
        logical(1)
      )))
      expect_true(all(vapply(
        meshes,
        function(m) nrow(m$normals) == nrow(m$vertices),
        logical(1)
      )))
      expect_true(all(vapply(
        meshes,
        function(m) all(m$indices >= 0 & m$indices < nrow(m$vertices)),
        logical(1)
      )))
    }
  }
})

test_that("spiders can't ski warns once and returns the unskied spider", {
  messages = character()
  scene = withCallingHandlers(
    pig(ski = TRUE, spider = TRUE, hair_count = 5),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(messages, "spiders can't ski")
  expect_equal(scene, pig(spider = TRUE, hair_count = 5))
})

test_that("spider hair is deterministic native hair and shares the mesh transform", {
  set.seed(271)
  before = .Random.seed
  scene = pig(
    spider = TRUE,
    hair_count = 24,
    x = 2,
    angle = c(0, 30, 0),
    scale = .5
  )
  expect_identical(.Random.seed, before)
  expect_equal(sum(scene$shape == "curve"), 24)
  again = pig(
    spider = TRUE,
    hair_count = 24,
    x = 2,
    angle = c(0, 30, 0),
    scale = .5
  )
  expect_identical(scene, again)
  different = pig(spider = TRUE, hair_count = 24, hair_seed = 9)
  expect_false(identical(scene$shape_info, different$shape_info))
  expect_true(all(vapply(
    scene$transforms,
    function(t) identical(t, scene$transforms[[1]]),
    logical(1)
  )))
  curves = scene$shape == "curve"
  expect_true(all(vapply(
    scene$material[curves],
    function(m) identical(m$type, hair()[[1]]$type),
    logical(1)
  )))
  expect_true(all(vapply(
    scene$shape_info[curves],
    function(s) s$shape_properties$curvetype == 1,
    logical(1)
  )))
  custom = hair(pigment = .3, red_pigment = 1.3)
  override = pig(spider = TRUE, hair_count = 3, hair_material = custom)
  expect_true(all(vapply(
    override$material[override$shape == "curve"],
    function(m) identical(m, custom[[1]]),
    logical(1)
  )))
  for (emotion in c("skeptical", "worried", "angry", "surprised", "excited")) {
    expect_equal(
      sum(
        pig(spider = TRUE, emotion = emotion, hair_count = 3)$shape == "curve"
      ),
      3
    )
  }
})

test_that("diffuse_sigma reaches matte surfaces without changing geometry", {
  plain = pig()
  rough = pig(diffuse_sigma = 35)
  geometry = function(scene) {
    lapply(scene$shape_info, function(s) {
      s$mesh_info[[1]][c("vertices", "indices", "normals")]
    })
  }
  expect_identical(geometry(plain), geometry(rough))
  diffuse_type = diffuse(sigma = 35)[[1]]$type
  matte = vapply(
    rough$material,
    function(m) identical(m$type, diffuse_type),
    logical(1)
  )
  expect_true(any(matte))
  expect_true(all(vapply(
    rough$material[matte],
    function(m) isTRUE(all.equal(m$sigma, 35 * pi / 180)),
    logical(1)
  )))
  expect_identical(plain$material[!matte], rough$material[!matte])
})

test_that("pig rejects invalid inputs before generating geometry", {
  expect_error(pig(ski = NA), "TRUE or FALSE")
  expect_error(pig(spider = 1), "TRUE or FALSE")
  expect_error(pig(emotion = "unknown"), "arg")
  expect_error(pig(x = Inf), "finite")
  expect_error(pig(scale = c(1, 2)), "scale")
  expect_error(pig(scale = 0), "nonzero")
  expect_error(pig(angle = NA_real_), "angle")
  expect_error(pig(order_rotation = c(1, 1, 3)), "permutation")
  expect_error(pig(diffuse_sigma = -1), "nonnegative")
  expect_error(pig(spider = TRUE, hair_count = -1), "hair_count")
  expect_error(pig(spider = TRUE, hair_count = .5), "hair_count")
  expect_error(pig(spider = TRUE, hair_length = 0), "hair_length")
  expect_error(pig(spider = TRUE, hair_seed = NA_real_), "hair_seed")
  expect_error(pig(spider = TRUE, hair_material = diffuse()), "hair_material")
})

test_that("all pig variants render from embedded meshes without OBJ files", {
  withr::local_options(cores = 1)
  for (variant in c("classic", "ski", "spider")) {
    scene = pig(
      ski = variant == "ski",
      spider = variant == "spider",
      hair_count = 8
    )
    path = withr::local_tempfile(fileext = ".png")
    expect_no_error(render_scene(
      scene,
      width = 32,
      height = 32,
      samples = 2,
      lookfrom = c(9, 6, 9),
      lookat = c(0, 1, 0),
      max_depth = 4,
      ambient_light = TRUE,
      preview = FALSE,
      interactive = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      bloom = FALSE,
      filename = path
    ))
    expect_true(file.exists(path))
  }
})
