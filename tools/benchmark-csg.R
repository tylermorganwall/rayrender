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

make_csg_scene <- function() {
  csg_shape <- csg_combine(
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
      csg_cylinder(
        start = c(0, 0, -1.2),
        end = c(0, 0, 1.2),
        radius = 0.28
      )
    )),
    operation = "subtract"
  )

  csg_object(csg_shape, material = diffuse(color = "white"))
}

render_csg_depth <- function(scene, width, height) {
  set.seed(1)
  invisible(render_scene(
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
  ))
}

run_csg_benchmark <- function(iterations = 5, width = 160, height = 120) {
  scene <- make_csg_scene()
  render_csg_depth(scene, width, height)

  elapsed <- numeric(iterations)
  for (i in seq_len(iterations)) {
    elapsed[[i]] <- unname(system.time(render_csg_depth(scene, width, height))[[
      "elapsed"
    ]])
  }

  data.frame(
    benchmark = "csg_depth_render",
    iteration = seq_len(iterations),
    width = width,
    height = height,
    elapsed = elapsed
  )
}

args <- commandArgs(trailingOnly = TRUE)
iterations <- if (length(args) >= 1) as.integer(args[[1]]) else 5L
width <- if (length(args) >= 2) as.integer(args[[2]]) else 160L
height <- if (length(args) >= 3) as.integer(args[[3]]) else 120L

suppressPackageStartupMessages(load_rayrender())
results <- run_csg_benchmark(iterations, width, height)
print(results)
cat("median_elapsed=", median(results$elapsed), "\n", sep = "")
