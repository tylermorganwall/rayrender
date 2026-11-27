# Run with the authored walking GLB from the Blender example generator:
# Rscript frozen-model-blur.R /path/to/rayrender-cycles output-directory
# The output contains ordinary PNG frames and a side-by-side blur comparison.
library(rayrender)

#' @param asset_dir Directory containing walk/walk.glb from the Blender example.
#' @return A scene with two independently phased walkers and a fast moving sphere.
#' @keywords internal
frozen_model_blur_scene = function(asset_dir) {
  source(
    system.file("animation", "blender-models.R", package = "rayrender"),
    local = TRUE
  )
  walkers = blender_walk_scene(asset_dir)
  # The moving sphere has no internal pose; its entire motion is sampled over
  # the shutter. The walkers retain only their outer placement motion blur.
  flyer = animation_keyframes(
    sphere(radius = .16, material = glossy("#FF8C28", gloss = .5)),
    animation_transform_track(
      c(0, .5, 1, 1.5, 2),
      position = rbind(
        c(-2, 1.3, .5),
        c(2, 1.7, .5),
        c(-2, 1.3, .5),
        c(2, 1.7, .5),
        c(-2, 1.3, .5)
      ),
      position_interpolation = "cubic"
    )
  )
  add_object(walkers$scene, flyer)
}

#' @param asset_dir Directory containing walk/walk.glb.
#' @param output Directory for PNG frames, comparison, and backend reports.
#' @return Invisibly, the rendered frames.
#' @keywords internal
render_frozen_model_blur = function(asset_dir, output) {
  scene = frozen_model_blur_scene(asset_dir)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  settings = list(
    scene = scene,
    width = 320,
    height = 240,
    samples = 32,
    preview = FALSE,
    progress = FALSE,
    plot_scene = FALSE,
    bloom = FALSE,
    denoise = FALSE,
    sample_method = "sobol",
    min_variance = 0
  )
  # Each output frame has a new internal walking pose. Within each exposure,
  # both walking poses are held while the outer paths and sphere keep moving.
  frames = do.call(
    render_animation,
    c(settings, list(times = (0:11) / 6, shutter_offsets = c(-.1, .1)))
  )
  for (i in seq_along(frames)) {
    rayimage::ray_write_image(
      frames[[i]],
      filename = file.path(output, sprintf("frame-%03d.png", i))
    )
  }
  sharp = do.call(
    render_scene,
    c(settings, list(time = 1, shutter_offsets = c(0, 0)))
  )
  blurred = frames[[7]]
  comparison = array(0, dim = c(dim(sharp)[1], 2 * dim(sharp)[2], 3))
  comparison[, seq_len(dim(sharp)[2]), ] = sharp[,, 1:3]
  comparison[, dim(sharp)[2] + seq_len(dim(sharp)[2]), ] = blurred[,, 1:3]
  rayimage::ray_write_image(
    comparison,
    filename = file.path(output, "sharp-and-placement-blur.png")
  )
  saveRDS(
    lapply(frames, attr, "wavefront"),
    file.path(output, "backend-reports.rds")
  )
  invisible(frames)
}

if (sys.nframe() == 0L) {
  args = commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) {
    stop("Supply the rayrender-cycles directory and an output directory.")
  }
  render_frozen_model_blur(args[1], args[2])
}
