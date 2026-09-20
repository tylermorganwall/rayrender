# Run with rayrender's imgui branch and rayimgui >= 0.0.9 installed.
# Shift-click to select; Shift-click again to drill down. Tree labels select directly.
library(rayrender)

# Resolve the demo asset beside this script when source() is called from another
# working directory. Direct execution from the checkout uses its tools directory.
demo_sources = lapply(sys.frames(), function(frame) frame$ofile)
demo_sources = Filter(Negate(is.null), demo_sources)
demo_directory = if (length(demo_sources)) {
  dirname(normalizePath(tail(demo_sources, 1)[[1]], mustWork = TRUE))
} else {
  "tools"
}
roughness_map = normalizePath(
  file.path(demo_directory, "editor-assets", "watermask.00000.png"),
  mustWork = TRUE
)

# A nested display group exercises selection of the whole assembly, its plinth,
# and the imported mesh. Keep the logo's original material slots for inspection.
plinth = group_objects(rbind(
  cube(
    y = 0.25,
    xwidth = 4.8,
    ywidth = 0.5,
    zwidth = 2.4,
    material = diffuse(color = "#505c6a")
  ),
  cube(
    y = 0.65,
    xwidth = 4.1,
    ywidth = 0.3,
    zwidth = 1.8,
    material = diffuse(color = "#e1d8c8")
  )
))
assembly = group_objects(add_object(
  plinth,
  obj_model(r_obj(), y = 2.2652, scale_obj = 3.8)
))

# Each placement shares an R mesh and a base, with its own rotation and scale.
# The small offset centers the bundled letter and rests its feet on the base.
letter_display = group_objects(rbind(
  cube(
    y = 0.175,
    xwidth = 2.5,
    ywidth = 0.35,
    zwidth = 1.6,
    material = diffuse(color = "#c7cbd0")
  ),
  obj_model(
    r_obj(simple_r = TRUE),
    x = 0.1883,
    y = 0.4904,
    scale_obj = 1.3,
    load_material = FALSE,
    material = glossy(color = "#168f89", gloss = 0.8)
  )
))
copies = create_instances(
  letter_display,
  x = c(-5, 0, 5),
  z = c(3.5, 5.5, 4),
  angle_y = c(-15, 10, 30),
  scale_x = c(0.9, 1.15, 1),
  scale_y = c(0.9, 1.15, 1),
  scale_z = c(0.9, 1.15, 1)
)

# Different surface types make the material inspector and sky reflections useful.
material_displays = group_objects(rbind(
  group_objects(rbind(
    cube(
      x = -3.7,
      y = 0.2,
      z = -1.7,
      xwidth = 2.2,
      ywidth = 0.4,
      zwidth = 2.2,
      material = diffuse(color = "#e1d8c8")
    ),
    sphere(
      x = -3.7,
      y = 1.25,
      z = -1.7,
      radius = 0.85,
      # Dark mask pixels are polished; light pixels scatter the reflections.
      material = microfacet(
        color = "#d9a441",
        roughness = 0.18,
        roughness_texture = roughness_map,
        roughness_range = c(0.03, 0.45)
      )
    )
  )),
  group_objects(rbind(
    cube(
      x = 3.8,
      y = 0.2,
      z = -2,
      xwidth = 2.4,
      ywidth = 0.4,
      zwidth = 2.4,
      material = diffuse(color = "#505c6a")
    ),
    sphere(
      x = 3.8,
      y = 1.35,
      z = -2,
      radius = 0.95,
      material = dielectric(color = "#d9f1ff")
    )
  )),
  cube(
    x = 6,
    y = 0.7,
    z = 0.4,
    width = 1.4,
    angle = c(0, 25, 0),
    material = glossy(color = "#e36b4f", gloss = 0.7)
  )
))

# A broad, flat plaza meets the sky cleanly in the low camera view.
scene = xz_rect(
  y = 0,
  xwidth = 20000,
  zwidth = 20000,
  material = diffuse(
    color = "#c2c5c5",
    checkercolor = "#aab0b3",
    checkerperiod = 3
  )
) |>
  add_object(assembly) |>
  add_object(copies) |>
  add_object(material_displays) |>
  add_infinite_light(sky_light_image(
    lat = 40.7,
    long = -74,
    datetime = as.POSIXct("2026-06-21 16:00:00", tz = "UTC"),
    hosek = TRUE,
    resolution = 512,
    moon = FALSE,
    number_cores = 2
  ))

# The renderer reads its worker count from options(), not a function argument.
local({
  old_options = options(cores = 10L)
  on.exit(options(old_options), add = TRUE)
  render_scene(
    scene,
    gui = "imgui",
    interactive = TRUE,
    deferred_render = TRUE,
    width = 640,
    height = 480,
    samples = 128,
    auto_exposure = TRUE,
    # A wide, near-eye-level view keeps the horizon and sky above the sculptures.
    fov = 55,
    lookfrom = c(7, 3.4, -15),
    lookat = c(0, 1.8, 1),
    ambient_light = FALSE,
    denoise = has_denoiser(),
    parallel = TRUE,
    progress = FALSE
  )
})
