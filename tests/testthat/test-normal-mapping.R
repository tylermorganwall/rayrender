test_that("diffuse has one normal-mapping model and validates roughness", {
  expect_false("normal_mapping" %in% names(formals(diffuse)))
  expect_null(diffuse()[[1]]$normal_mapping)
  expect_error(diffuse(normal_mapping = "legacy"), "unused argument")
  for (sigma in list(-1, Inf, NA_real_, NaN, numeric(), c(1, 2), "rough")) {
    expect_error(diffuse(sigma = sigma), "sigma")
  }
  rough = diffuse(sigma = 20)
  expect_identical(rough[[1]]$type, 4L)
  expect_equal(rough[[1]]$sigma, 20 * pi / 180)
  expect_equal(diffuse(sigma = 900)[[1]]$sigma, pi / 2)
  expect_identical(unserialize(serialize(rough, NULL)), rough)
  expect_no_error(prepare_scene_list(sphere(material = rough)))
  expect_no_error(prepare_scene_list(sphere(material = diffuse(fog = TRUE))))
})

test_that("imported diffuse face materials and vertex colors render with the default model", {
  folder = tempfile()
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  obj = file.path(folder, "triangle.obj")
  writeLines(
    c("newmtl matte", "Kd 0.3 0.6 0.2", "illum 1"),
    file.path(folder, "triangle.mtl")
  )
  writeLines(
    c(
      "mtllib triangle.mtl",
      "v -2 -2 0 1 0 0",
      "v 2 -2 0 0 1 0",
      "v 0 2 0 0 0 1",
      "usemtl matte",
      "f 1 2 3"
    ),
    obj
  )
  mesh = rayvertex::sphere_mesh()
  scenes = list(
    obj_model(obj),
    obj_model(
      obj,
      load_material = FALSE,
      material = diffuse(color = "grey60", sigma = 45)
    ),
    obj_model(obj, vertex_colors = TRUE),
    raymesh_model(mesh, override_material = FALSE),
    raymesh_model(
      mesh,
      override_material = TRUE,
      material = diffuse(color = "grey60", sigma = 45)
    ),
    bezier_curve(
      c(-1, -1, 0),
      c(-1, 1, 0),
      c(1, -1, 0),
      c(1, 1, 0),
      width = .2,
      material = diffuse(color = "grey60")
    )
  )
  for (scene in scenes) {
    set.seed(551)
    image = render_scene(
      scene,
      width = 12,
      height = 12,
      samples = 16,
      lookfrom = c(0, 0, 4),
      lookat = c(0, 0, 0),
      fov = 50,
      integrator_type = "basic",
      ambient_light = TRUE,
      backgroundhigh = "white",
      backgroundlow = "white",
      tonemap = "raw",
      denoise = FALSE,
      bloom = FALSE,
      min_variance = 0,
      max_depth = 3,
      parallel = FALSE,
      preview = FALSE,
      progress = FALSE,
      plot_scene = FALSE
    )
    expect_true(all(is.finite(image)))
    expect_gt(mean(image[,, 1:3]), 0)
    expect_lt(min(image[,, 1:3]), .99)
  }
})

test_that("physical smooth diffuse agrees across basic, mixture and NEE estimators", {
  obj = tempfile(fileext = ".obj")
  on.exit(unlink(obj))
  writeLines(
    c(
      "v -20 -20 0",
      "v 20 -20 0",
      "v 20 20 0",
      "v -20 20 0",
      "vn 0.6 0 0.8",
      "f 1//1 2//1 3//1",
      "f 1//1 3//1 4//1"
    ),
    obj
  )
  for (sigma in c(0, 11.5, 90)) {
    scene = obj_model(
      obj,
      load_material = FALSE,
      calculate_consistent_normals = TRUE,
      material = diffuse(
        color = rep(.7, 3),
        sigma = sigma
      )
    ) |>
      add_object(xy_rect(
        x = 2,
        y = 1,
        z = 2,
        xwidth = 2,
        ywidth = 2,
        flipped = TRUE,
        material = light(intensity = 5)
      ))
    means = matrix(0, 4, 3)
    for (j in seq_along(c("basic", "rtiow", "nee"))) {
      for (seed in 1:4) {
        set.seed(seed + 860)
        image = render_scene(
          scene,
          width = 12,
          height = 12,
          samples = 128,
          sample_method = "random",
          lookfrom = c(0, 0, 3),
          lookat = c(0, 0, 0),
          fov = 20,
          integrator_type = c("basic", "rtiow", "nee")[j],
          ambient_light = FALSE,
          backgroundhigh = "black",
          backgroundlow = "black",
          tonemap = "raw",
          denoise = FALSE,
          bloom = FALSE,
          min_variance = 0,
          max_depth = 3,
          parallel = FALSE,
          preview = FALSE,
          progress = FALSE,
          plot_scene = FALSE
        )
        expect_true(all(is.finite(image)))
        means[seed, j] = mean(image[,, 1:3])
      }
    }
    for (j in 1:2) {
      difference = means[, j] - means[, 3]
      # Independent seeds, conservative four-standard-error comparison.
      tolerance = 4 * sd(difference) / sqrt(nrow(means)) + .004
      expect_lt(abs(mean(difference)), tolerance)
    }
    expect_gt(min(means), .01)
  }
})

test_that("physical diffuse coexists with unchanged glass, SSS and a medium", {
  scene = sphere(
    x = -.8,
    radius = .5,
    material = diffuse()
  ) |>
    add_object(sphere(
      x = .7,
      radius = .5,
      material = subsurface(sigma_s = 2, sigma_a = .1)
    )) |>
    add_object(sphere(x = 0, y = -.6, radius = .3, material = dielectric())) |>
    add_object(set_medium(
      cube(width = 5),
      homogeneous_medium(sigma_s = .01)
    )) |>
    add_object(xy_rect(
      z = 2,
      y = 3,
      flipped = TRUE,
      material = light(intensity = 8)
    ))
  set.seed(991)
  image = render_scene(
    scene,
    width = 16,
    height = 16,
    samples = 16,
    lookfrom = c(0, 0, 4),
    lookat = c(0, 0, 0),
    integrator_type = "nee",
    tonemap = "raw",
    denoise = FALSE,
    bloom = FALSE,
    min_variance = 0,
    parallel = TRUE,
    preview = FALSE,
    progress = FALSE,
    plot_scene = FALSE
  )
  expect_true(all(is.finite(image)))
  expect_gt(mean(image[,, 1:3]), 0)
})

test_that("overridden and instanced physical meshes render under structured environment light", {
  image_file = tempfile(fileext = ".png")
  on.exit(unlink(image_file))
  environment = array(.15, c(8, 16, 3))
  environment[2:4, 3:7, ] = .9
  png::writePNG(environment, image_file)
  mesh = raymesh_model(
    rayvertex::sphere_mesh(),
    override_material = TRUE,
    scale = c(1.2, .8, .7),
    material = diffuse()
  )
  scene = create_instances(mesh, x = c(-.8, .8)) |>
    add_infinite_light(infinite_light(image_file))
  for (integrator in c("basic", "rtiow", "nee")) {
    set.seed(107)
    result = render_scene(
      scene,
      width = 16,
      height = 12,
      samples = 32,
      lookfrom = c(0, 0, 5),
      lookat = c(0, 0, 0),
      fov = 35,
      integrator_type = integrator,
      tonemap = "raw",
      denoise = FALSE,
      min_variance = 0,
      bloom = FALSE,
      parallel = TRUE,
      preview = FALSE,
      progress = FALSE,
      plot_scene = FALSE
    )
    expect_true(all(is.finite(result)))
    expect_gt(mean(result[,, 1:3]), .01)
  }
})
