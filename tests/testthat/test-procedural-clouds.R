cloud_example_helpers = new.env()
cloud_example_file = testthat::test_path(
  "..",
  "..",
  "inst",
  "examples",
  "procedural-clouds.R"
)
if (!file.exists(cloud_example_file)) {
  cloud_example_file = system.file(
    "examples",
    "procedural-clouds.R",
    package = "rayrender"
  )
}
sys.source(cloud_example_file, envir = cloud_example_helpers)

local_cloud_example_scene = function(spatial = TRUE, ...) {
  skip_if_not_installed("rayshader")
  skip_if_not_installed("ambient")
  skip_if_not("get_scene_metadata" %in% getNamespaceExports("rayshader"))
  withr::local_options(list(rgl.useNULL = TRUE), .local_envir = parent.frame())
  withr::defer(rgl::close3d(), envir = parent.frame())
  terrain = matrix(seq(100, 179), 10, 8)
  args = list(
    hillshade = rayshader::height_shade(terrain),
    heightmap = terrain,
    zscale = 100,
    shadow = FALSE,
    windowsize = c(240, 200),
    ...
  )
  if (spatial) {
    args$extent = c(xmin = 500000, xmax = 500900, ymin = 0, ymax = 1400)
    args$crs = 32631
  }
  do.call(rayshader::plot_3d, args)
  testthat::local_mocked_bindings(
    render_highquality = function(...) list(...),
    .package = "rayshader",
    .env = parent.frame()
  )
  rayshader::get_scene_metadata()
}

test_that("physical cloud sizes follow each horizontal axis and exaggerated elevations", {
  info = local_cloud_example_scene(
    geographic_aspect = FALSE,
    vertical_exaggeration = 2
  )
  args = cloud_example_helpers$render_perlin_cloud_scene(
    altitude = 1000,
    cloud_width = 400,
    cloud_depth = 600,
    cloud_height = 500,
    resolution = 24,
    denoise = FALSE,
    samples = 7,
    lat = -36.87593,
    long = 174.7647,
    datetime = as.POSIXct("2026-01-15 13:00:00", tz = "Pacific/Auckland"),
    sky_args = list(hosek = FALSE)
  )
  cloud = args$scene_elements
  size = cloud$shape_info[[1]]$shape_properties$boxinfo
  horizontal_scale = info$geographic_aspect$cell_meters
  expect_equal(size[1] * horizontal_scale[["x"]], 400)
  expect_equal(size[3] * horizontal_scale[["z"]], 600)
  expect_equal(size[2], 500 / info$effective_zscale)
  origin_y = mean(rgl::par3d("bbox")[3:4])
  expect_equal(
    cloud$y - size[2] / 2 + origin_y,
    1000 / info$effective_zscale
  )
  terrain_relative = cloud_example_helpers$render_perlin_cloud_scene(
    altitude = 1000,
    altitude_reference = "terrain",
    cloud_width = 400,
    cloud_depth = 600,
    cloud_height = 500,
    resolution = 24
  )
  expect_equal(
    terrain_relative$scene_elements$y - cloud$y,
    info$scene_bounds["y", "max"]
  )
  expect_identical(args$integrator_type, "nee")
  expect_identical(args$samples, 7)
  expect_false(args$denoise)
  expect_null(args$camera_location)
  expect_identical(args$sky_args, list(hosek = FALSE))
  expect_equal(c(args$lat, args$long), c(-36.87593, 174.7647))
  expect_equal(
    args$datetime,
    as.POSIXct("2026-01-15 13:00:00", tz = "Pacific/Auckland")
  )
})

test_that("cloud defaults cover the footprint and physical unit choices agree", {
  info = local_cloud_example_scene()
  args = list(altitude = 200, cloud_height = 100, resolution = 24)
  metric = do.call(cloud_example_helpers$render_perlin_cloud_scene, args)
  expect_false(any(
    c(
      "sample_method",
      "max_depth",
      "min_variance",
      "clamp_value",
      "denoise",
      "bloom",
      "tonemap",
      "preview",
      "plot",
      "lightsize",
      "lightintensity",
      "lightcolor",
      "ambient_light",
      "backgroundhigh",
      "backgroundlow",
      "camera_location",
      "camera_lookat",
      "ortho_dimensions"
    ) %in%
      names(metric)
  ))
  kilometers = cloud_example_helpers$render_perlin_cloud_scene(
    altitude = 0.2,
    cloud_height = 0.1,
    units = "kilometers",
    resolution = 24
  )
  feet = cloud_example_helpers$render_perlin_cloud_scene(
    altitude = 200 / 0.3048,
    cloud_height = 100 / 0.3048,
    units = "feet",
    resolution = 24
  )
  expect_equal(metric$scene_elements, kilometers$scene_elements)
  expect_equal(metric$scene_elements, feet$scene_elements)
  properties = metric$scene_elements$shape_info[[1]]$shape_properties$boxinfo
  expect_equal(properties[1], info$scene_dimensions[["x"]])
  expect_equal(properties[3], info$scene_dimensions[["z"]])
})

test_that("absolute altitude, feet elevations, and spatial centers map into the scene", {
  info = local_cloud_example_scene()
  extent = info$extent_3d
  args = cloud_example_helpers$render_perlin_cloud_scene(
    altitude = 1000 * 0.3048,
    altitude_reference = "absolute",
    elevation_units = "feet",
    center = c(extent[["xmax"]], extent[["ymax"]]),
    cloud_height = 200 * 0.3048,
    resolution = 24
  )
  cloud = args$scene_elements
  size = cloud$shape_info[[1]]$shape_properties$boxinfo
  expect_equal(cloud$x, info$scene_bounds["x", "max"], ignore_attr = TRUE)
  expect_equal(cloud$z, info$scene_bounds["z", "min"], ignore_attr = TRUE)
  expect_equal(size[2], 200 / info$effective_zscale)
  expect_equal(
    cloud$y - size[2] / 2 + mean(rgl::par3d("bbox")[3:4]),
    1000 / info$effective_zscale
  )
})

test_that("panning and display scaling preserve physical attachment and the rgl view", {
  info = local_cloud_example_scene()
  baseline = cloud_example_helpers$render_perlin_cloud_scene(
    resolution = 24
  )
  view = rgl::par3d()
  pan = c(2, 3, -1)
  matrix = view$userMatrix
  matrix[1:3, 4] = matrix[1:3, 1:3] %*% pan
  rgl::par3d(userMatrix = matrix, scale = c(2, 1.5, 0.75))
  shifted = cloud_example_helpers$render_perlin_cloud_scene(
    resolution = 24
  )
  a = baseline$scene_elements
  b = shifted$scene_elements
  expect_equal(c(b$x, b$y, b$z), (c(a$x, a$y, a$z) + pan) * c(2, 1.5, 0.75))
  expect_equal(rgl::par3d("userMatrix"), matrix)
  expect_equal(rgl::par3d("scale"), c(2, 1.5, 0.75))
  expect_equal(rayshader::get_scene_metadata()$scene_bounds, info$scene_bounds)
})

test_that("plain matrices need an explicit physical scale and render overrides work", {
  local_cloud_example_scene(spatial = FALSE)
  expect_error(
    cloud_example_helpers$render_perlin_cloud_scene(resolution = 24),
    "Physical cloud sizes require spatial distance metadata"
  )
  args = cloud_example_helpers$render_perlin_cloud_scene(
    resolution = 24,
    meters_per_scene_unit = c(100, 200),
    camera_location = c(10, 20, 30),
    camera_lookat = c(1, 2, 3),
    scene_elements = sphere()
  )
  expect_equal(nrow(args$scene_elements), 2)
  expect_equal(args$camera_location, c(10, 20, 30))
  expect_equal(args$camera_lookat, c(1, 2, 3))
  expect_null(args$ortho_dimensions)
  expect_error(
    cloud_example_helpers$render_perlin_cloud_scene(
      resolution = 24,
      meters_per_scene_unit = 100,
      center = c(0, 0, 0)
    ),
    "center must be a finite"
  )
  expect_error(
    cloud_example_helpers$render_perlin_cloud_scene(
      resolution = 24,
      meters_per_scene_unit = 100,
      cloud_height = -1
    ),
    "cloud_height must"
  )
  expect_error(
    cloud_example_helpers$render_perlin_cloud_scene(
      resolution = 24,
      meters_per_scene_unit = 100,
      integrator_type = "basic"
    ),
    "require integrator_type"
  )
})

test_that("clouds preserve the rgl bounds and camera for both projections", {
  for (fov in c(0, 45)) {
    local_cloud_example_scene(fov = fov)
    view = rgl::par3d()
    fields = c(
      "bbox",
      "userMatrix",
      "modelMatrix",
      "projMatrix",
      "observer",
      "zoom",
      "FOV",
      "scale"
    )
    for (crop in c(FALSE, TRUE)) {
      args = cloud_example_helpers$render_perlin_cloud_scene(
        altitude = 2000,
        cloud_width = 5000,
        cloud_depth = 6000,
        crop_to_extent = crop,
        resolution = 24
      )
      expect_false(any(
        c("camera_location", "camera_lookat", "ortho_dimensions") %in%
          names(args)
      ))
      expect_equal(rgl::par3d()[fields], view[fields])
    }
  }
})

test_that("extent cropping preserves the field while clipping transformed boundaries", {
  info = local_cloud_example_scene(geographic_aspect = FALSE)
  view = rgl::par3d()
  pan = c(2, 3, -1)
  matrix = view$userMatrix
  matrix[1:3, 4] = matrix[1:3, 1:3] %*% pan
  scale = c(2, 1.5, 0.75)
  rgl::par3d(userMatrix = matrix, scale = scale)
  args = list(
    altitude = 200,
    cloud_height = 300,
    cloud_width = 3000,
    cloud_depth = 4000,
    center = c(500700, 400),
    resolution = 24
  )
  original = do.call(
    cloud_example_helpers$render_perlin_cloud_scene,
    args
  )$scene_elements
  cropped = do.call(
    cloud_example_helpers$render_perlin_cloud_scene,
    c(args, list(crop_to_extent = TRUE))
  )$scene_elements
  original_position = c(original$x, original$y, original$z)
  cropped_position = c(cropped$x, cropped$y, cropped$z)
  size = cropped$shape_info[[1]]$shape_properties$boxinfo
  cropped_bounds = sweep(rbind(-size / 2, size / 2), 2, cropped_position, "+")
  offset = rowMeans(matrix(rgl::par3d("bbox"), nrow = 3, byrow = TRUE)) - pan
  expected = sweep(sweep(t(info$scene_bounds), 2, offset), 2, scale, "*")
  expect_equal(unname(cropped_bounds[, c(1, 3)]), unname(expected[, c(1, 3)]))
  expect_equal(cropped$y, original$y)
  expect_equal(size[2], original$shape_info[[1]]$shape_properties$boxinfo[2])
  a = original$shape_info[[1]]$medium
  b = cropped$shape_info[[1]]$medium
  expect_identical(a$density, b$density)
  expect_identical(
    a[c("sigma_a", "sigma_s", "g")],
    b[c("sigma_a", "sigma_s", "g")]
  )
  expect_equal(
    sweep(a$bounds, 2, original_position, "+"),
    sweep(b$bounds, 2, cropped_position, "+")
  )
  expect_false(
    "crop_to_extent" %in%
      names(do.call(
        cloud_example_helpers$render_perlin_cloud_scene,
        c(args, list(crop_to_extent = TRUE))
      ))
  )
})

test_that("extent cropping leaves contained clouds alone and clips a partial overlap", {
  info = local_cloud_example_scene()
  args = list(
    cloud_width = 400,
    cloud_depth = 500,
    resolution = 24
  )
  original = do.call(
    cloud_example_helpers$render_perlin_cloud_scene,
    args
  )$scene_elements
  contained = do.call(
    cloud_example_helpers$render_perlin_cloud_scene,
    c(args, list(crop_to_extent = TRUE))
  )$scene_elements
  expect_equal(contained, original)
  partial = do.call(
    cloud_example_helpers$render_perlin_cloud_scene,
    c(args, list(crop_to_extent = TRUE, center = c(500900, 700)))
  )$scene_elements
  size = partial$shape_info[[1]]$shape_properties$boxinfo
  expect_equal(
    size[1],
    original$shape_info[[1]]$shape_properties$boxinfo[1] / 2
  )
  expect_equal(partial$x + size[1] / 2, info$scene_bounds["x", "max"])
  expect_equal(size[3], original$shape_info[[1]]$shape_properties$boxinfo[3])
  expect_error(
    cloud_example_helpers$render_perlin_cloud_scene(
      crop_to_extent = NA,
      resolution = 24
    ),
    "crop_to_extent must be TRUE or FALSE"
  )
  expect_error(
    cloud_example_helpers$render_perlin_cloud_scene(
      crop_to_extent = TRUE,
      center = c(600000, 700),
      resolution = 24
    ),
    "does not overlap the terrain extent"
  )
})
