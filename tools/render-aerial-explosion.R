# Run from the rayrender repository after installing rayrender and openvdbr.
# Optional arguments: scene.pbrt, output directory.
args = commandArgs(trailingOnly = TRUE)
scene_file = if (length(args)) {
  args[1]
} else {
  'tests/aerial_explosion/aerial_explosion.pbrt'
}
output_dir = if (length(args) > 1L) {
  args[2]
} else {
  'tests/aerial_explosion/openvdbr-render'
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
scene_file = normalizePath(scene_file, mustWork = TRUE)
output_dir = normalizePath(output_dir, mustWork = TRUE)
options(cores = 6L)
library(rayrender)

import_time = system.time(imported <- read_pbrt(scene_file, strict = FALSE))
write.csv(
  imported$diagnostics,
  file.path(output_dir, 'diagnostics.csv'),
  row.names = FALSE
)
settings = imported$render_args
settings$filename = file.path(output_dir, 'aerial-explosion.png')
settings$sample_method = 'sobol'
settings$plot_scene = FALSE
settings$verbose = TRUE
settings$progress = FALSE
settings$mode = 'image'
set.seed(20261008)
render_time = system.time(
  image <- do.call(render_scene, c(list(scene = imported$scene), settings))
)
saveRDS(
  list(
    import_seconds = unname(import_time['elapsed']),
    render_seconds = unname(render_time['elapsed']),
    cores = 6L,
    render_args = settings,
    path_warnings = attr(image, 'path_warnings'),
    openvdbr_messages = openvdbr::vdb_messages(),
    package_versions = c(
      rayrender = as.character(packageVersion('rayrender')),
      openvdbr = as.character(packageVersion('openvdbr'))
    )
  ),
  file.path(output_dir, 'results.rds')
)
cat(sprintf(
  'Import %.3f s; render %.3f s; output %s\n',
  import_time['elapsed'],
  render_time['elapsed'],
  settings$filename
))
