load_rayrender <- function() {
  if (
    !identical(Sys.getenv("RAYRENDER_BENCH_USE_INSTALLED"), "true") &&
      requireNamespace("devtools", quietly = TRUE) &&
      file.exists(file.path(getwd(), "DESCRIPTION"))
  ) {
    suppressWarnings(devtools::load_all(".", quiet = TRUE, debug = FALSE))
  } else {
    library(rayrender)
  }
}

csg_case_specs <- function() {
  list(
    sphere = function() csg_sphere(radius = 0.9),
    plane_rotated = function() {
      csg_plane(normal = c(1, 2, 1), width_x = 2.4, width_z = 2.4)
    },
    box = function() csg_box(width = c(1.5, 1.2, 1.0)),
    rounded_box = function()
      csg_box(width = c(1.4, 1.2, 1.0), corner_radius = 0.18),
    torus = function() csg_torus(radius = 0.72, minor_radius = 0.22),
    capsule = function() {
      csg_capsule(
        start = c(-0.9, -0.35, 0),
        end = c(0.9, 0.35, 0),
        radius = 0.32
      )
    },
    cylinder = function() {
      csg_cylinder(
        start = c(-0.85, -0.45, 0),
        end = c(0.85, 0.45, 0),
        radius = 0.33
      )
    },
    rounded_cylinder = function() {
      csg_cylinder(
        start = c(-0.85, -0.45, 0),
        end = c(0.85, 0.45, 0),
        radius = 0.28,
        corner_radius = 0.12
      )
    },
    ellipsoid = function() csg_ellipsoid(axes = c(0.85, 0.55, 1.05)),
    rounded_cone = function() {
      csg_rounded_cone(
        start = c(-0.75, -0.45, 0),
        end = c(0.75, 0.45, 0),
        radius = 0.45,
        upper_radius = 0.2
      )
    },
    cone = function() {
      csg_cone(
        start = c(-0.75, -0.45, 0),
        end = c(0.75, 0.45, 0),
        radius = 0.55
      )
    },
    pyramid = function() csg_pyramid(y = -0.65, height = 1.3, base = 1.3),
    triangle = function() {
      csg_triangle(
        v1 = c(0, 0.9, 0),
        v2 = c(0.9, -0.55, 0),
        v3 = c(-0.9, -0.55, 0)
      )
    },
    elongate_fast = function() {
      csg_elongate(
        csg_sphere(radius = 0.42),
        elongate = c(0.55, 0.15, 0),
        robust = FALSE
      )
    },
    elongate_robust = function() {
      csg_elongate(csg_sphere(radius = 0.42), elongate = c(0.55, 0.15, 0.35))
    },
    round = function()
      csg_round(csg_box(width = c(1.1, 1.1, 1.1)), radius = 0.18),
    onion_subtract = function() {
      csg_combine(
        csg_onion(csg_sphere(radius = 0.9), thickness = 0.18),
        csg_box(y = 0.55, width = c(3, 1.1, 3)),
        operation = "subtract"
      )
    },
    scale = function()
      csg_scale(csg_pyramid(y = -0.45, height = 0.9, base = 0.9), scale = 1.35),
    rotate_angles = function() {
      csg_rotate(csg_box(width = c(1.4, 0.8, 1.0)), angles = c(18, 35, 12))
    },
    rotate_up = function() {
      csg_rotate(
        csg_torus(radius = 0.72, minor_radius = 0.2),
        up = c(1, 1, 0.35)
      )
    },
    translate = function()
      csg_translate(csg_sphere(radius = 0.7), x = 0.35, y = -0.2, z = 0.15),
    group = function() {
      csg_group(list(
        csg_sphere(x = -0.65, z = -0.65, radius = 0.38),
        csg_sphere(x = 0.65, z = -0.65, radius = 0.38),
        csg_sphere(x = -0.65, z = 0.65, radius = 0.38),
        csg_sphere(x = 0.65, z = 0.65, radius = 0.38)
      ))
    },
    union = function() {
      csg_combine(
        csg_sphere(x = -0.35, radius = 0.7),
        csg_sphere(x = 0.35, radius = 0.7),
        operation = "union"
      )
    },
    subtract = function() {
      csg_combine(
        csg_sphere(radius = 0.9),
        csg_sphere(x = 0.45, y = 0.1, z = 0.1, radius = 0.55),
        operation = "subtract"
      )
    },
    intersection = function() {
      csg_combine(
        csg_box(width = c(1.35, 1.35, 1.35)),
        csg_sphere(radius = 0.92),
        operation = "intersection"
      )
    },
    blend = function() {
      csg_combine(
        csg_sphere(x = -0.45, radius = 0.62),
        csg_sphere(x = 0.45, radius = 0.62),
        operation = "blend",
        radius = 0.28
      )
    },
    subtractblend = function() {
      csg_combine(
        csg_sphere(radius = 0.9),
        csg_sphere(x = 0.48, y = 0.1, radius = 0.55),
        operation = "subtractblend",
        radius = 0.22
      )
    },
    mix = function() {
      csg_combine(
        csg_box(width = c(1.25, 1.25, 1.25)),
        csg_torus(radius = 0.72, minor_radius = 0.24),
        operation = "mix"
      )
    },
    nested_classic = function() {
      csg_combine(
        csg_combine(
          csg_box(width = c(1.6, 1.6, 1.6)),
          csg_sphere(radius = 1.05),
          operation = "intersection"
        ),
        csg_group(list(
          csg_cylinder(
            start = c(-1.2, 0, 0),
            end = c(1.2, 0, 0),
            radius = 0.28
          ),
          csg_cylinder(
            start = c(0, -1.2, 0),
            end = c(0, 1.2, 0),
            radius = 0.28
          ),
          csg_cylinder(start = c(0, 0, -1.2), end = c(0, 0, 1.2), radius = 0.28)
        )),
        operation = "subtract"
      )
    }
  )
}

render_csg_case <- function(object, width, height) {
  scene <- csg_object(object, material = diffuse(color = "white"))
  set.seed(1)
  render_scene(
    scene,
    width = width,
    height = height,
    samples = 1,
    sample_method = "random",
    parallel = FALSE,
    preview = FALSE,
    progress = FALSE,
    verbose = FALSE,
    plot_scene = FALSE,
    debug_channel = "depth",
    lookfrom = c(4, 3, 6),
    lookat = c(0, 0, 0),
    fov = 22,
    aperture = 0,
    clamp_value = 10
  )
}

checksum_depth <- function(image) {
  finite <- is.finite(image)
  finite_values <- image[finite]
  list(
    finite_count = sum(finite),
    finite_sum = if (length(finite_values)) sum(finite_values) else 0
  )
}

benchmark_one_case <- function(case_name, factory, iterations, width, height) {
  object <- factory()
  image <- render_csg_case(object, width, height)
  checksum <- checksum_depth(image)

  elapsed <- numeric(iterations)
  for (i in seq_len(iterations)) {
    elapsed[[i]] <- unname(system.time(render_csg_case(object, width, height))[[
      "elapsed"
    ]])
  }

  data.frame(
    benchmark = case_name,
    iteration = seq_len(iterations),
    width = width,
    height = height,
    elapsed = elapsed,
    finite_count = checksum$finite_count,
    finite_sum = checksum$finite_sum
  )
}

run_csg_suite <- function(iterations = 5, width = 192, height = 144) {
  cases <- csg_case_specs()
  do.call(
    rbind,
    lapply(names(cases), function(case_name) {
      benchmark_one_case(
        case_name,
        cases[[case_name]],
        iterations,
        width,
        height
      )
    })
  )
}

summarize_results <- function(results) {
  aggregate(
    elapsed ~ benchmark + width + height + finite_count + finite_sum,
    results,
    median
  )
}

args <- commandArgs(trailingOnly = TRUE)
iterations <- if (length(args) >= 1) as.integer(args[[1]]) else 5L
width <- if (length(args) >= 2) as.integer(args[[2]]) else 192L
height <- if (length(args) >= 3) as.integer(args[[3]]) else 144L
output_file <- if (length(args) >= 4) args[[4]] else NA_character_

invisible(suppressPackageStartupMessages(load_rayrender()))
results <- run_csg_suite(iterations, width, height)
summary <- summarize_results(results)

print(summary[order(summary$benchmark), ], row.names = FALSE)
cat("suite_median_elapsed=", median(summary$elapsed), "\n", sep = "")

if (!is.na(output_file) && nzchar(output_file)) {
  write.csv(results, output_file, row.names = FALSE)
}
