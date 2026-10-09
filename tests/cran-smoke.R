# Small end-to-end checks shipped with the package. The larger rendering and
# reference-image suite remains a developer check.
library(rayrender)
local({
  previous = options(cores = 2L)
  on.exit(options(previous), add = TRUE)
  file = tempfile(fileext = '.vdb')
  on.exit(unlink(file), add = TRUE)
  fog = openvdbr::vdb_fog(openvdbr::vdb_sphere(
    radius = 0.5,
    voxel_size = 0.15,
    name = 'density'
  ))
  openvdbr::vdb_write(fog, file)

  # Exercise both compiled shader headers and the runtime OpenVDB interface.
  surfaces = sphere(
    material = openpbr(
      base_color = texture_checker('coral', 'steelblue'),
      specular_roughness = 0.3
    )
  )
  volume = set_medium(
    cube(width = 2),
    vdb_medium(
      file,
      sigma_a = 0.2,
      sigma_s = 0.5,
      emission = 2
    )
  )
  for (scene in list(surfaces, volume)) {
    scene = add_infinite_light(
      scene,
      disk_light(
        direction = c(-1, 1, 2),
        angular_diameter = 30,
        intensity = 5
      )
    )
    set.seed(129)
    image = render_scene(
      scene,
      width = 12,
      height = 12,
      samples = 4,
      lookfrom = c(0, 0, 4),
      lookat = c(0, 0, 0),
      denoise = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      mode = 'image'
    )
    stopifnot(
      identical(dim(image)[1:2], c(12L, 12L)),
      all(is.finite(image)),
      max(image[,, 1:3]) > 0
    )
  }
})
