test_that("medium descriptions validate coefficients and fields", {
  expect_equal(homogeneous_medium(sigma_s = 2)$sigma_s, rep(2, 3))
  expect_error(homogeneous_medium(sigma_a = -1), "nonnegative")
  expect_error(homogeneous_medium(sigma_s = c(1, NA, 0)), "finite")
  expect_error(homogeneous_medium(g = 1), "strictly")
  expect_error(homogeneous_medium(g = -1), "strictly")
  expect_error(homogeneous_medium(temperature = 3000, emission = 1), "not both")
  expect_error(
    homogeneous_medium(medium_transform = matrix(0, 4, 4)),
    "transform"
  )
  expect_error(grid_medium(array(-1, c(2, 2, 2))), "density")
  expect_error(grid_medium(matrix(1, 2, 2)), "3D")
  expect_error(
    grid_medium(array(1, c(2, 2, 2)), temperature = array(1, c(3, 2, 2))),
    "dimensions"
  )
  expect_error(
    grid_medium(array(1, c(2, 2, 2)), emission = array(1, c(2, 2, 2, 4))),
    "dimensions"
  )
  expect_error(nanovdb_medium(tempfile()), "existing")
  expect_error(
    set_medium(xy_rect(), homogeneous_medium()),
    "boundary|closed|supported"
  )
})

test_that("medium metadata survives scene operations and serialization", {
  medium = grid_medium(array(1, c(2, 3, 4)), sigma_a = c(0, 0.2, 0.5))
  scene = set_medium(sphere(), medium) |> add_object(cube())
  restored = unserialize(serialize(scene, NULL))
  expect_identical(restored$shape_info[[1]]$medium, medium)
  expect_null(restored$shape_info[[2]]$medium)
  replaced = set_scene_material(scene, diffuse("red"))
  expect_identical(replaced$shape_info[[1]]$medium, medium)
  instanced = create_instances(restored, x = c(-2, 2))
  child = instanced$shape_info[[1]]$shape_properties$original_scene[[1]]
  expect_identical(child$shape_info[[1]]$medium, medium)
  expect_null(set_medium(sphere(), NULL)$shape_info[[1]]$medium)
  expect_error(
    render_scene(
      scene,
      integrator_type = "basic",
      width = 4,
      height = 4,
      samples = 1,
      preview = FALSE,
      plot_scene = FALSE
    ),
    'require.*nee'
  )
})

test_that("volume foreground converts to straight RGB once", {
  raw = list(
    r = matrix(c(0, 0.1, 0.6), 1),
    g = matrix(c(0, 0.2, 0.3), 1),
    b = matrix(c(0, 0.3, 0.15), 1),
    a = matrix(c(0, 0.25, 1), 1),
    premultiplied = TRUE
  )
  straight = rayrender:::straight_volume_rgb(raw)
  expect_equal(straight$r, matrix(c(0, 0.4, 0.6), 1))
  expect_equal(straight$a, raw$a)
  expect_false(straight$premultiplied)
})

# Measurements here diagnose transport; tools/codex/volume-validation.R provides
# the corresponding images that must be opened and reviewed before release.
medium_test_render = function(scene, samples = 64, ...) {
  render_scene(
    scene,
    integrator_type = "nee",
    width = 12,
    height = 12,
    samples = samples,
    lookfrom = c(0, 0, 3),
    lookat = c(0, 0, 0),
    fov = 0,
    ortho_dimensions = c(0.7, 0.7),
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black",
    min_variance = 0,
    clamp_value = Inf,
    tonemap = "raw",
    denoise = FALSE,
    bloom = FALSE,
    preview = FALSE,
    plot_scene = FALSE,
    parallel = FALSE,
    progress = FALSE,
    ...
  )
}

test_that("primary alpha measures world distance, nesting, and surface coverage", {
  a = c(0, 1, 2)
  medium = homogeneous_medium(sigma_a = a, sigma_s = 0)
  slab = set_medium(cube(width = 1), medium)
  expected = 1 - mean(exp(-a))
  image = medium_test_render(slab, samples = 4, transparent_background = TRUE)
  expect_equal(as.numeric(image[,, 4]), rep(expected, 144), tolerance = 2e-5)
  cavity = slab |>
    add_object(set_medium(cube(width = 0.5), homogeneous_medium(sigma_s = 0)))
  image = medium_test_render(cavity, samples = 4, transparent_background = TRUE)
  expect_equal(
    mean(image[5:8, 5:8, 4]),
    1 - mean(exp(-a * 0.5)),
    tolerance = 2e-5
  )
  glass = set_medium(cube(material = dielectric()), medium, keep_surface = TRUE)
  image = medium_test_render(glass, samples = 4, transparent_background = TRUE)
  expect_equal(as.numeric(image[,, 4]), rep(1, 144))
})

test_that("emitting slab and equivalent grid agree with radiative transfer", {
  set.seed(12)
  m = homogeneous_medium(sigma_a = 1, sigma_s = 0, emission = c(1, 0.5, 0.2))
  homogeneous = medium_test_render(set_medium(cube(), m), samples = 128)
  expected = (1 - exp(-1)) * c(1, 0.5, 0.2)
  expect_equal(apply(homogeneous[,, 1:3], 3, mean), expected, tolerance = 0.025)
  grid = grid_medium(
    array(1, c(4, 4, 4)),
    sigma_a = 1,
    sigma_s = 0,
    emission = c(1, 0.5, 0.2)
  )
  set.seed(12)
  dense = medium_test_render(set_medium(cube(), grid), samples = 128)
  expect_equal(apply(dense[,, 1:3], 3, mean), expected, tolerance = 0.025)
})

test_that("invalid NanoVDB files and majorant overflow fail before rendering", {
  malformed = tempfile(fileext = ".nvdb")
  on.exit(unlink(malformed))
  writeBin(as.raw(c(0, 1, 2)), malformed)
  expect_error(
    medium_test_render(
      set_medium(cube(), nanovdb_medium(malformed)),
      samples = 1
    ),
    "NanoVDB"
  )
  huge = grid_medium(array(1e30, c(2, 2, 2)), sigma_s = 1e30)
  expect_error(
    medium_test_render(set_medium(cube(), huge), samples = 1),
    "float range"
  )
})

test_that("sparse NanoVDB tiles load without expanding their distant bounding box", {
  fixture = test_path("fixtures", "volumes", "tiles.nvdb")
  medium = nanovdb_medium(fixture, sigma_a = 0.2, sigma_s = 0)
  transform = diag(4)
  transform[1:3, 4] = -3.5
  medium$medium_transform = transform
  image = medium_test_render(
    set_medium(cube(), medium),
    samples = 128,
    transparent_background = TRUE
  )
  expect_lt(abs(mean(image[,, 4]) - (1 - exp(-0.2))), 0.015)
  wrong = nanovdb_medium(test_path("fixtures", "volumes", "double.nvdb"))
  expect_error(
    medium_test_render(set_medium(cube(), wrong), samples = 1),
    "float"
  )
  missing = nanovdb_medium(fixture, density_grid = "missing")
  expect_error(
    medium_test_render(set_medium(cube(), missing), samples = 1),
    "not found|missing"
  )
})

test_that("raw NanoVDB buffers render like their file container", {
  original = readBin(
    test_path("fixtures", "volumes", "tiles.nvdb"),
    "raw",
    n = 1e7
  )
  # NanoVDB's fixed 16-byte file header, 176-byte metadata, and "density\0" name.
  raw_grid = original[-seq_len(200)]
  temporary = tempfile(fileext = ".nvdb")
  on.exit(unlink(temporary))
  writeBin(raw_grid, temporary)
  transform = diag(4)
  transform[1:3, 4] = -3.5
  raw_medium = nanovdb_medium(
    temporary,
    sigma_a = 0.2,
    sigma_s = 0,
    medium_transform = transform
  )
  set.seed(39)
  raw_image = medium_test_render(set_medium(cube(), raw_medium), samples = 16)
  file_medium = raw_medium
  file_medium$filename = normalizePath(test_path(
    "fixtures",
    "volumes",
    "tiles.nvdb"
  ))
  set.seed(39)
  file_image = medium_test_render(set_medium(cube(), file_medium), samples = 16)
  expect_identical(as.numeric(raw_image), as.numeric(file_image))
  raw_medium$density_grid = "missing"
  expect_error(
    medium_test_render(set_medium(cube(), raw_medium), samples = 1),
    "missing"
  )
  # The raw grid header claims a second grid that is absent from this allocation.
  raw_grid[29:32] = as.raw(c(2, 0, 0, 0))
  writeBin(raw_grid, temporary)
  expect_error(
    medium_test_render(set_medium(cube(), raw_medium), samples = 1),
    "Truncated"
  )
  raw_medium$density_grid = "density"
  raw_grid = original[-seq_len(200)]
  raw_grid[9:16] = as.raw(255) # An optional checksum must not be required for validation.
  raw_grid[457:464] = writeBin(Inf, raw(), size = 8, endian = "little")
  writeBin(raw_grid, temporary)
  expect_error(
    medium_test_render(set_medium(cube(), raw_medium), samples = 1),
    "transform"
  )
  raw_grid = original[-seq_len(200)]
  raw_grid[9:16] = as.raw(255)
  raw_grid[697:704] = as.raw(c(rep(0, 7), 127)) # Invalid tree-relative root offset.
  writeBin(raw_grid, temporary)
  expect_error(
    medium_test_render(set_medium(cube(), raw_medium), samples = 1),
    "root offset"
  )
})

test_that("alpha is independent of exposure and denoising", {
  scene = set_medium(
    sphere(),
    homogeneous_medium(sigma_s = 0.6, emission = 1, sigma_a = 0.2)
  )
  set.seed(17)
  a = medium_test_render(scene, samples = 8, transparent_background = TRUE)
  set.seed(17)
  b = medium_test_render(
    scene,
    samples = 8,
    transparent_background = TRUE,
    iso = 200
  )
  expect_identical(as.numeric(a[,, 4]), as.numeric(b[,, 4]))
})

test_that("compressed and truncated NanoVDB files report actionable errors", {
  original = readBin(
    test_path("fixtures", "volumes", "tiles.nvdb"),
    "raw",
    n = 1e7
  )
  temporary = tempfile(fileext = ".nvdb")
  on.exit(unlink(temporary))
  compressed = original
  compressed[15] = as.raw(1)
  writeBin(compressed, temporary)
  expect_error(
    medium_test_render(
      set_medium(cube(), nanovdb_medium(temporary)),
      samples = 1
    ),
    "Compressed"
  )
  writeBin(original[seq_len(length(original) %/% 2)], temporary)
  expect_error(
    medium_test_render(
      set_medium(cube(), nanovdb_medium(temporary)),
      samples = 1
    ),
    "Truncated|sizes"
  )
})

test_that("observed non-nested boundaries are diagnosed", {
  vacuum = homogeneous_medium(sigma_s = 0)
  overlap = set_medium(cube(x = -0.2, z = 0.3), vacuum) |>
    add_object(set_medium(cube(x = 0.2, z = -0.3), vacuum))
  expect_error(medium_test_render(overlap, samples = 1), "Non-nested")
})

test_that("glowing media inside glass include radiance transport across the interface", {
  medium = homogeneous_medium(sigma_a = 100, sigma_s = 0, emission = 1)
  scene = set_medium(
    cube(material = dielectric(refraction = 1.5)),
    medium,
    keep_surface = TRUE
  )
  set.seed(1402)
  image = medium_test_render(scene, samples = 128)
  fresnel = ((1.5 - 1) / (1.5 + 1))^2
  expect_lt(abs(mean(image[,, 1]) - (1 - fresnel) / 1.5^2), 0.01)
})
