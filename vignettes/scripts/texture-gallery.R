## ---- gallery-tools ----
library(rayrender)

# Run from either the repository root or the vignette directory. When knitted,
# this script supplies named chunks; Rscript can also regenerate the whole set.
if (dir.exists("vignettes/scripts")) {
  setwd("vignettes")
}
gallery_dir = file.path("figures", "composable-textures")
rebuild_gallery = identical(Sys.getenv("RAYRENDER_REBUILD_TEXTURES"), "true")
gallery_groups = strsplit(
  Sys.getenv("RAYRENDER_TEXTURE_GROUPS"),
  ",",
  fixed = TRUE
)[[1]]
dir.create(gallery_dir, recursive = TRUE, showWarnings = FALSE)

#' @param name Stable figure name.
#' @param materials Named list of material descriptors, one per image.
#' @param layout Default `"sphere"`. Geometry/camera arrangement.
#' @param samples Default `64L`. Samples per pixel.
#' @param width Default `420L`. Width of each rendered image.
#' @param height Default `420L`. Height of each rendered image.
#' @return A knitr image reference to the assembled comparison.
texture_panel = function(
  name,
  materials,
  layout = "sphere",
  samples = 64L,
  width = 420L,
  height = 420L
) {
  filename = file.path(gallery_dir, paste0(name, ".png"))
  selected = !length(gallery_groups) || name %in% gallery_groups
  if (rebuild_gallery && selected) {
    old_options = options(cores = 6L)
    on.exit(options(old_options))
    images = character(length(materials))
    timings = vector("list", length(materials))
    for (i in seq_along(materials)) {
      images[i] = file.path(gallery_dir, sprintf("%s-%02d.png", name, i))
      scene = texture_gallery_scene(materials[[i]], layout)
      settings = if (layout == "pair") {
        list(
          lookfrom = c(0, 1.1, 6),
          lookat = c(0, -.1, 0),
          fov = 0,
          ortho_dimensions = c(4.9, 3.43)
        )
      } else if (layout == "plane") {
        list(lookfrom = c(2.6, 2.2, 5), lookat = c(0, -1, 0), fov = 46)
      } else {
        list(
          lookfrom = c(2.8, 1.4, 6),
          lookat = c(0, -.05, 0),
          fov = 0,
          ortho_dimensions = c(2.9, 2.9)
        )
      }
      set.seed(20261002)
      message("Rendering ", name, ": ", names(materials)[i])
      elapsed = system.time(do.call(
        render_scene,
        c(
          list(
            scene = scene,
            width = width,
            height = height,
            samples = samples,
            sample_method = "sobol",
            ambient_light = TRUE,
            backgroundhigh = "#d5e1ee",
            backgroundlow = "#d9d2c6",
            min_variance = 0,
            max_depth = 16,
            roulette_active_depth = 5,
            parallel = TRUE,
            preview = FALSE,
            progress = FALSE,
            plot_scene = FALSE,
            denoise = FALSE,
            bloom = FALSE,
            tonemap = "hbd",
            filename = images[i]
          ),
          settings
        )
      ))
      timings[[i]] = data.frame(
        group = name,
        example = names(materials)[i],
        samples = samples,
        cores = 6L,
        elapsed = unname(elapsed["elapsed"])
      )
    }
    utils::write.csv(
      do.call(rbind, timings),
      file.path(gallery_dir, paste0(name, "-timings.csv")),
      row.names = FALSE
    )
    # Captions are outside the rendered pixels, so every comparison uses the
    # exact same framing and exposure. Individual renders are retained too.
    grDevices::png(
      filename,
      width = width * length(images),
      height = height + 38
    )
    par(mar = c(0, 0, 0, 0), bg = "#f3f0e9", xaxs = "i", yaxs = "i")
    plot.new()
    plot.window(c(0, width * length(images)), c(0, height + 38))
    for (i in seq_along(images)) {
      rasterImage(
        png::readPNG(images[i]),
        (i - 1) * width,
        0,
        i * width,
        height
      )
      text(
        (i - .5) * width,
        height + 19,
        names(materials)[i],
        cex = 1.15,
        col = "#242d36"
      )
    }
    dev.off()
  }
  # Source archives use compact display copies; local regeneration keeps PNGs.
  compact = sub("[.]png$", ".jpg", filename)
  if (
    file.exists(compact) &&
      (!file.exists(filename) || file.mtime(compact) >= file.mtime(filename))
  ) {
    filename = compact
  }
  if (!file.exists(filename)) {
    stop(
      "Missing vignette render: ",
      filename,
      ". Run with RAYRENDER_REBUILD_TEXTURES=true."
    )
  }
  knitr::include_graphics(filename)
}

#' @param material Material for the demonstration geometry.
#' @param layout Arrangement selected by texture_panel().
#' @return Scene with consistent studio lighting.
texture_gallery_scene = function(material, layout) {
  floor_material = if (layout == "glass") {
    diffuse(texture_checker(
      "#e2d9c4",
      "#25444f",
      texture_coordinates("world", scale = 1.5)
    ))
  } else {
    diffuse("#b9b5ac")
  }
  scene = if (layout == "plane") {
    xz_rect(y = -1, xwidth = 9, zwidth = 9, material = material)
  } else {
    object = if (layout == "pair") {
      create_instances(
        sphere(material = material),
        x = c(-1.1, 1.1),
        scale_x = .85,
        scale_y = .85,
        scale_z = .85
      )
    } else {
      sphere(material = material)
    }
    object |>
      add_object(generate_ground(depth = -1.02, material = floor_material))
  }
  scene |>
    add_object(sphere(
      x = -3.5,
      y = 5,
      z = 4,
      radius = 1.5,
      material = light(color = "#fff4e6", intensity = 9)
    )) |>
    add_object(sphere(
      x = 4,
      y = 2.8,
      z = -2,
      radius = 1.1,
      material = light(color = "#d6e7ff", intensity = 5)
    ))
}

## ---- shared-field ----
grain = texture_noise(scale = 3, octaves = 5, seed = 14)
# The outer mix clamps its weight. Expanding the noise range increases contrast.
mask = texture_mix(-.65, 1.65, weight = grain)
pigment = texture_mix("#153f54", "#e1a34c", weight = mask)
roughness = texture_mix(.04, .7, weight = mask)

texture_panel(
  "shared-field",
  list(
    "Scalar mask" = diffuse(texture_mix(0, 1, mask)),
    "Mask into color" = openpbr(
      base_color = pigment,
      base_metalness = .7,
      specular_roughness = .25
    ),
    "Same mask into roughness" = openpbr(
      base_color = pigment,
      base_metalness = .7,
      specular_roughness = roughness
    )
  )
)

## ---- nested-patterns ----
fine_checker = texture_checker(
  "#17495a",
  "#f1d79c",
  coordinates = texture_coordinates("object", scale = 9)
)
nested_checker = texture_checker(
  fine_checker,
  "#cc674b",
  coordinates = texture_coordinates("object", scale = 2)
)
stripes = texture_checker(
  0,
  1,
  coordinates = texture_coordinates("object", scale = c(0, 7, 0))
)
striped_stone = texture_mix(
  texture_mix("#203d4a", "#a5c6bf", grain),
  "#e0b85d",
  weight = stripes
)

texture_panel(
  "nested-patterns",
  list(
    "Fine checker" = diffuse(fine_checker),
    "Checker inside checker" = diffuse(nested_checker),
    "Noise with horizontal bands" = openpbr(
      base_color = striped_stone,
      specular_roughness = texture_mix(.12, .65, stripes)
    )
  )
)

## ---- stretched-noise ----
stretched = texture_noise(
  octaves = 5,
  seed = 26,
  coordinates = texture_coordinates(
    "object",
    scale = c(1, 22, 2),
    rotation = 22
  )
)
wood_tones = texture_mix("#351c12", "#d99a49", texture_mix(-.6, 1.6, stretched))
stone = texture_mix(
  "#34464e",
  "#dce0d8",
  texture_mix(-.45, 1.45, texture_noise(scale = 5, seed = 42))
)
weathered = texture_mix(
  wood_tones,
  stone,
  texture_mix(-1, 2, texture_noise(scale = 1.6, seed = 2))
)

texture_panel(
  "stretched-noise",
  list(
    "Stretched, rotated noise" = openpbr(
      base_color = wood_tones,
      specular_roughness = texture_mix(.15, .6, stretched)
    ),
    "Stone-like color variation" = diffuse(stone),
    "Noise blends two recipes" = openpbr(
      base_color = weathered,
      specular_roughness = .45
    )
  )
)

## ---- directional ----
both_poles = texture_direction_mix("#eee2c2", "#26465a", direction = c(0, 1, 0))
upper_only = texture_direction_mix(
  "#eee2c2",
  "#26465a",
  direction = c(0, 1, 0),
  absolute = FALSE
)
up_mask = texture_direction_mix(1, 0, direction = c(0, 1, 0), absolute = FALSE)
broken_deposit = texture_scale(up_mask, texture_mix(-.6, 1.6, grain))

texture_panel(
  "directional",
  list(
    "Absolute: both poles" = diffuse(both_poles),
    "One-sided: upper surfaces" = diffuse(upper_only),
    "Orientation multiplied by noise" = openpbr(
      base_color = texture_mix("#26465a", "#eee2c2", broken_deposit),
      specular_roughness = texture_mix(.1, .8, broken_deposit)
    )
  )
)

## ---- coordinate-spaces ----
object_gradient = texture_gradient(
  "#1a6678",
  "#db8549",
  axis = "x",
  coordinates = texture_coordinates("object", scale = .5, offset = c(.5, 0, 0))
)
world_gradient = texture_gradient(
  "#1a6678",
  "#db8549",
  axis = "x",
  coordinates = texture_coordinates("world", scale = .5, offset = c(.5, 0, 0))
)

texture_panel(
  "coordinate-spaces",
  list(
    "Object: gradient repeats in each instance" = diffuse(object_gradient),
    "World: one gradient spans both instances" = diffuse(world_gradient)
  ),
  layout = "pair",
  width = 560,
  height = 392
)

## ---- uv-versus-object ----
object_tiles = texture_checker(
  "#15475c",
  "#edcc82",
  coordinates = texture_coordinates("object", scale = 4)
)
uv_tiles = texture_checker(
  "#15475c",
  "#edcc82",
  coordinates = texture_coordinates("uv", scale = c(16, 8))
)

texture_panel(
  "uv-versus-object",
  list(
    "3D checker cuts through the surface" = diffuse(object_tiles),
    "UV checker follows the sphere's UVs" = diffuse(uv_tiles)
  )
)

## ---- roughness-finishes ----
bands = texture_checker(
  0,
  1,
  coordinates = texture_coordinates("object", scale = c(0, 5, 0))
)
band_roughness = texture_mix(.03, .6, bands)
gold = texture_mix("#655035", "#ecc46d", grain)

texture_panel(
  "roughness-finishes",
  list(
    "Constant roughness" = openpbr(
      base_color = gold,
      base_metalness = 1,
      specular_roughness = .23
    ),
    "Polished and rough bands" = openpbr(
      base_color = gold,
      base_metalness = 1,
      specular_roughness = band_roughness
    ),
    "Legacy microfacet roughness" = microfacet(
      color = gold,
      eta = c(.18, .42, 1.37),
      kappa = c(3.42, 2.35, 1.77),
      roughness = band_roughness
    )
  ),
  samples = 128L
)

## ---- glass ----
view_axis = c(2.8, 1.4, 6)
glass_roughness = texture_direction_mix(.015, .28, direction = view_axis)
texture_panel(
  "glass",
  list(
    "Clear glass" = openpbr(
      transmission_weight = 1,
      specular_ior = 1.45,
      specular_roughness = .015
    ),
    "Roughness increases toward silhouette" = openpbr(
      transmission_weight = 1,
      specular_ior = 1.45,
      specular_roughness = glass_roughness
    )
  ),
  layout = "glass",
  samples = 128L
)

## ---- make-image ----
tile_file = file.path(gallery_dir, "source-tile.png")
if (rebuild_gallery || !file.exists(tile_file)) {
  # An original RGB test tile: rings in red, a ramp in green, diagonal bands in blue.
  uv = (seq_len(128) - .5) / 128
  xx = matrix(rep(uv, each = 128), 128, 128)
  yy = matrix(rep(rev(uv), 128), 128, 128)
  tile = array(0, c(128, 128, 3))
  tile[,, 1] = .1 + .8 * (sin(32 * sqrt((xx - .5)^2 + (yy - .5)^2)) > 0)
  tile[,, 2] = .15 + .7 * yy
  tile[,, 3] = .15 + .7 * (sin(15 * (xx + yy)) > 0)
  png::writePNG(tile, tile_file)
}

## ---- image-composition ----
image_color = texture_image(tile_file, encoding = "srgb")
red_mask = texture_channel(image_color, "r")
recolored = texture_mix("#123d51", "#dca554", red_mask)
image_roughness = texture_scale(
  texture_image(
    tile_file,
    type = "scalar",
    coordinates = texture_coordinates(
      "uv",
      scale = c(3, 2),
      rotation = 25,
      offset = c(.1, -.15)
    )
  ),
  .7
)

texture_panel(
  "image-composition",
  list(
    "RGB image" = openpbr(
      base_color = image_color,
      base_metalness = .65,
      specular_roughness = .3
    ),
    "Red channel into a new palette" = openpbr(
      base_color = recolored,
      base_metalness = .65,
      specular_roughness = .3
    ),
    "Independent roughness image mapping" = openpbr(
      base_color = recolored,
      base_metalness = .65,
      specular_roughness = image_roughness
    )
  ),
  samples = 128L
)

## ---- fine-patterns ----
coarse_floor = texture_checker(
  "#173b49",
  "#e3d9ba",
  texture_coordinates("world", scale = 3)
)
fine_floor = texture_checker(
  "#173b49",
  "#e3d9ba",
  texture_coordinates("world", scale = 90)
)

texture_panel(
  "fine-patterns",
  list(
    "Resolved checker cells" = diffuse(coarse_floor),
    "Subpixel cells: sampling artifacts" = diffuse(fine_floor)
  ),
  layout = "plane",
  samples = 64L
)
