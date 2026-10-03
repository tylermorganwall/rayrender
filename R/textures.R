#' Composable Surface Textures
#'
#' Build scalar and color texture graphs for `diffuse(color)`,
#' `microfacet(color, roughness)`, and `openpbr(base_color, specular_roughness)`.
#' Constants are accepted wherever a child texture is expected. Scalar values
#' broadcast to colors; use `texture_channel()` to convert a color to a scalar.
#' Graphs are evaluated in C++ at surface hits, without callbacks to R.
#'
#' @param value A finite scalar, linear RGB triple, R color, or texture.
#' @param space Default `"object"` for `texture_coordinates()` and `"world"` for
#' `texture_direction_mix()`. Coordinate space: `"object"`, `"world"`, or `"uv"`;
#' direction mixing supports only object and world space.
#' @param scale Default `1`. Scalar or three-component coordinate multiplier;
#' UV mappings also accept two components. For `texture_noise()`, a positive
#' scalar frequency multiplier applied in addition to its coordinate mapping.
#' @param offset Default `c(0, 0, 0)`. Coordinate translation with three components,
#' or two for UV mappings.
#' @param rotation Default `0`. Rotation in degrees about the texture z axis,
#' applied after scaling and before translation (UV rotation is in the UV plane).
#' @param coordinates Default `texture_coordinates()`, except for
#' `texture_image()`, which defaults to `texture_coordinates("uv")`.
#' Independent texture mapping.
#' @param a First value or texture.
#' @param b Second value or texture.
#' @param weight Default `0.5`. Scalar blend weight, clamped to `[0, 1]`.
#' @param direction Default `c(0, 1, 0)`. Nonzero direction in the selected space.
#' @param absolute Default `TRUE`. Use the absolute geometric-normal dot product.
#' If false, use the positive part, so opposite-facing surfaces receive zero weight.
#' @param octaves Default `4L`. Number of noise octaves, between one and sixteen.
#' @param seed Default `0L`. Integer seed; independent of the renderer's sampling seed.
#' @param filename Image filename.
#' @param type Default `"color"`. Output type: `"color"` or `"scalar"`.
#' @param encoding Default `"auto"`. `"auto"`, `"linear"`, or `"srgb"`.
#' Auto treats scalar images and HDR/EXR images as linear, other color images as sRGB.
#' @param wrap Default `"repeat"`. Image addressing: `"repeat"` or `"clamp"`.
#' @param channel Default `"luminance"`. Color reduction: `"r"`, `"g"`, `"b"`,
#' `"average"`, or linear RGB `"luminance"`.
#' @param factor Scalar multiplier or scalar texture. Multiplies the texture's
#' output values; does not change the size or spacing of its pattern.
#' @param axis Default `"y"`. Coordinate axis for a clamped zero-to-one gradient.
#' @details
#' Each operation produces a color or scalar value that can feed another
#' operation or a supported material input. For example, one spill mask can
#' blend a fabric's color toward brown and its roughness toward a glossy finish.
#' This composes individual inputs of a material, rather than blending two
#' complete materials or their scattering models.
#'
#' ## texture_constant(): reuse a uniform value
#'
#' Wrap a color or number that stays the same everywhere on a surface. Use it
#' for a shared paint color, a fixed roughness, or one endpoint of a blend.
#' Explicit wrapping is optional: other texture operations accept ordinary
#' colors and numbers and wrap them automatically. Passing an existing texture
#' returns that texture unchanged.
#'
#' ## texture_coordinates(): place and resize a pattern
#'
#' Control where a pattern is evaluated, independently of the object's geometry.
#' Object coordinates keep a pattern attached to an object as it moves or
#' rotates; world coordinates align patterns across separate objects; UV
#' coordinates follow the surface's UV layout, as needed for wrapping an image.
#' Scaling coordinates by two makes repeating features half as wide along those
#' axes. Unequal axis scales stretch noise into grain, while rotation and offset
#' orient and position the pattern. Transformations apply in the order scale,
#' rotation, then translation. This returns a mapping, not a color or scalar.
#'
#' ## texture_mix(): blend values using a mask
#'
#' Blend `a` and `b` with a constant amount or a spatially varying scalar mask.
#' The result is `(1-weight)*a + weight*b`: zero selects `a`, one selects `b`,
#' and intermediate weights blend them. Weights are clamped to `[0, 1]`. Use a
#' painted image mask for stains, or noise for mottled paint. Reuse the same
#' mask in color and roughness blends to make a spill both darker and shinier.
#'
#' ## texture_direction_mix(): vary a surface by its orientation
#'
#' Blend according to the angle between the geometric normal and a direction.
#' This is useful for snow-colored upward faces, moss on one side of an object,
#' or roughness that differs between horizontal and vertical surfaces.
#' Following PBRT, alignment weights **a**: parallel normals select `a`, and
#' perpendicular normals select `b`. With `absolute = TRUE`, opposite-facing
#' normals also select `a`; use `absolute = FALSE` for a one-sided effect such
#' as snow on top but not underneath. The default direction is world up.
#' This is an orientation mask, not a simulation of deposition or exposure;
#' occlusion and bump mapping do not change its weight.
#'
#' ## texture_noise(): introduce irregular variation
#'
#' Generate smooth scalar values in `[0, 1]` without an image file. Use noise as
#' a blend mask for mottling, as a roughness input for an uneven finish, or with
#' stretched coordinates for wood-like grain. Larger `scale` values make finer
#' features. More `octaves` add progressively finer detail with decreasing,
#' normalized weights; they do not increase the output range. `seed` chooses a
#' repeatable pattern independently of the render's sampling seed.
#'
#' ## texture_checker(): alternate between two patterns
#'
#' Select `a` or `b` in alternating unit cells of the mapped coordinates. Either
#' input can itself be a texture, so checks can alternate colors, finishes, or
#' more elaborate patterns. UV coordinates make a two-dimensional checker;
#' object and world coordinates make a three-dimensional checker through which
#' the surface passes. Setting two coordinate scales to zero produces stripes
#' along the remaining axis. Change cell size with `texture_coordinates()`.
#'
#' ## texture_gradient(): make a gradual spatial transition
#'
#' Blend from `a` to `b` along one mapped coordinate axis. Coordinates at or
#' below zero select `a`, coordinates at or above one select `b`, and values
#' between them blend linearly. Use this for a color fade, a gradual change in
#' polish, or a mask that changes with height. Coordinate scale and offset set
#' the transition's extent and starting position; it does not repeat.
#'
#' ## texture_image(): bring painted or measured data into a graph
#'
#' Read an image as surface color, a blend mask, or a roughness map, then combine
#' it with other operations. UV mapping is the default. `wrap = "repeat"` tiles
#' the image; `wrap = "clamp"` extends its edge values outside the image bounds.
#' Use sRGB decoding for ordinary color artwork, and `encoding = "linear"` for
#' masks or other numeric data. `type = "scalar"` averages the decoded RGB
#' channels; to extract one channel from a packed image, load it as color and
#' use `texture_channel()`.
#'
#' ## texture_channel(): extract a mask or scalar control
#'
#' Convert a color texture to a scalar for inputs such as blend weight or
#' roughness. Select `"r"`, `"g"`, or `"b"` to unpack separate masks stored in
#' one RGB image. Use `"luminance"` for a brightness-weighted combination of
#' linear RGB, or `"average"` to give the three channels equal weight. A scalar
#' input passes through unchanged. For packed numeric masks, load the image
#' with `encoding = "linear"` so color decoding does not alter the mask values.
#'
#' ## texture_scale(): adjust the strength of texture values
#'
#' Multiply the evaluated texture by `factor`, applying the same scalar factor
#' to all RGB channels for a color texture. Use `texture_scale(paint, 0.5)` to
#' halve its linear RGB reflectance, or `texture_scale(roughness_map, 0.25)` to
#' reduce roughness values while preserving their spatial pattern. A texture
#' factor provides a varying multiplier, for example to restrict a noise mask
#' to selected regions. This changes value amplitude, not pattern size: use
#' `texture_coordinates(scale = ...)` to resize features. It is not exposure
#' control, and halving linear RGB does not halve displayed sRGB brightness.
#' The operation itself does not clamp values; the receiving blend or material
#' input may clamp them.
#'
#' ## Color handling and current limitations
#'
#' Numeric colors are linear RGB, following rayrender's material convention.
#' Image color decoding happens before filtering and mixing. Scalar images use
#' the average of RGB channels. Image filtering is bilinear; procedural textures
#' currently evaluate a point sample, without ray-footprint antialiasing.
#' Object coordinates refer to the primitive's original coordinates, including
#' inside instances; world coordinates include the instance placement.
#' Noise is seeded smooth value noise with normalized, decreasing octave weights.
#' Roughness textures are clamped to `[0, 1]` at the material input. The microfacet
#' material uses the same conversion as its numeric roughness argument.
#' Spatial microfacet roughness uses a minimum alpha of `1e-6`; a constant zero
#' roughness uses the existing perfectly smooth material. Transmitting microfacet
#' color graphs require nonzero roughness.
#' Graph inputs cannot be combined with the legacy image/procedural settings for
#' the same input. Bump, alpha, emission, and volume inputs do not yet accept graphs.
#' @return A serializable `ray_texture` or `ray_texture_coordinates` descriptor.
#' @name textures
#' @md
#' @export
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # texture_constant(): one reusable color and one scalar polish value.
#' # Stack three lacquered toy blocks; constants do not vary over a surface.
#' paint = texture_constant("#168c89")
#' polish = texture_constant(.18)
#' toy = openpbr(base_color = paint, specular_roughness = polish, coat_weight = .5)
#' scene = generate_ground(depth = 0, material = diffuse("#efe1bf"))
#' for (i in 1:3) {
#'   scene = add_object(
#'     scene,
#'     cube(
#'       x = c(-.2, .15, -.1)[i],
#'       y = .4 + (i - 1) * .8,
#'       width = .8,
#'       angle = c(0, c(-15, 12, -8)[i], 0),
#'       material = toy
#'     )
#'   )
#' }
#' scene = scene |>
#'   add_infinite_light(disk_light(
#'     direction = c(-1, 1.2, 1),
#'     angular_diameter = 35,
#'     intensity = 12,
#'     color = "#fff0d8",
#'     name = "key"
#'   )) |>
#'   add_infinite_light(disk_light(
#'     direction = c(1, .4, -.7),
#'     angular_diameter = 45,
#'     intensity = .5,
#'     color = "#adcfff",
#'     name = "rim"
#'   ))
#' render_scene(
#'   scene,
#'   lookfrom = c(4, 3, 7),
#'   lookat = c(0, 1.15, 0),
#'   fov = 0,
#'   ortho_dimensions = c(3.8, 3.6),
#'   width = 420,
#'   height = 340,
#'   samples = 16,
#'   denoise = TRUE
#' )
texture_constant = function(value) {
  if (inherits(value, "ray_texture")) {
    return(value)
  }
  if (is.character(value)) {
    value = convert_color(value)
  }
  if (
    !is.numeric(value) ||
      !length(value) %in% c(1L, 3L) ||
      !is.null(dim(value)) ||
      any(!is.finite(value))
  ) {
    stop(
      "A texture constant must be a finite scalar or RGB triple.",
      call. = FALSE
    )
  }
  texture_node(
    "constant",
    if (length(value) == 1L) "scalar" else "color",
    value = as.numeric(value)
  )
}

#' @rdname textures
#' @export
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # texture_coordinates(): barber poles with local rings versus world-space bands.
#' # The left pattern tilts with the pole; the right stays level in the room.
#' local = texture_coordinates("object", scale = c(0, 5, 0))
#' world = texture_coordinates("world", scale = c(0, 5, 0), offset = c(0, .3, 0))
#' scene = generate_ground(depth = 0, material = diffuse("#e5ded2"))
#' for (i in 1:2) {
#'   stripes = texture_checker(
#'     "#fff4d9",
#'     "#d94759",
#'     coordinates = list(local, world)[[i]]
#'   )
#'   scene = add_object(
#'     scene,
#'     cylinder(
#'       x = c(-.9, .9)[i],
#'       y = 1.1,
#'       radius = .24,
#'       length = 2,
#'       angle = c(0, 0, 22),
#'       material = diffuse(stripes)
#'     )
#'   )
#' }
#' scene = scene |>
#'   add_infinite_light(disk_light(
#'     direction = c(-1, 1.2, 1),
#'     angular_diameter = 35,
#'     intensity = 12,
#'     color = "#fff0d8",
#'     name = "key"
#'   )) |>
#'   add_infinite_light(disk_light(
#'     direction = c(1, .4, -.7),
#'     angular_diameter = 45,
#'     intensity = .5,
#'     color = "#adcfff",
#'     name = "rim"
#'   ))
#' render_scene(
#'   scene,
#'   lookfrom = c(0, 2, 7),
#'   lookat = c(0, 1, 0),
#'   fov = 0,
#'   ortho_dimensions = c(4.8, 3.8),
#'   width = 420,
#'   height = 340,
#'   samples = 16,
#'   denoise = TRUE
#' )
texture_coordinates = function(
  space = "object",
  scale = 1,
  offset = c(0, 0, 0),
  rotation = 0
) {
  space = match.arg(space, c("object", "world", "uv"))
  if (space == "uv" && is.numeric(scale) && length(scale) == 2L) {
    scale = c(scale, 1)
  }
  if (space == "uv" && is.numeric(offset) && length(offset) == 2L) {
    offset = c(offset, 0)
  }
  if (
    !is.numeric(scale) ||
      !length(scale) %in% c(1L, 3L) ||
      any(!is.finite(scale)) ||
      !is.numeric(offset) ||
      length(offset) != 3L ||
      any(!is.finite(offset)) ||
      !is.numeric(rotation) ||
      length(rotation) != 1L ||
      !is.finite(rotation)
  ) {
    stop(
      "Invalid texture coordinate scale, offset, or rotation.",
      call. = FALSE
    )
  }
  structure(
    list(
      space = space,
      scale = rep(scale, length.out = 3L),
      offset = offset,
      rotation = rotation * pi / 180
    ),
    class = "ray_texture_coordinates"
  )
}

#' @rdname textures
#' @export
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # texture_mix(): coffee on a tightly woven rug.
#' # A grayscale spill mask blends both color and roughness: dry yarn versus wet coffee.
#' uv = seq(0, 1, length.out = 384)
#' u = matrix(rep(uv, each = length(uv)), length(uv))
#' v = t(u)
#' spill = pmin(
#'   1,
#'   pmax(
#'     0,
#'     5 *
#'       (exp(-((u - .62)^2 / .02 + (v - .55)^2 / .035)) +
#'         .65 * exp(-((u - .43)^2 / .008 + (v - .48)^2 / .01)) +
#'         .7 * exp(-((u - .28)^2 + (v - .7)^2) / .0007) +
#'         .5 * exp(-((u - .37)^2 + (v - .64)^2) / .0003) -
#'         .24)
#'   )
#' )
#' file = tempfile(fileext = ".png")
#' dim(spill) = dim(u)
#' png::writePNG(spill, file)
#' mask = texture_image(file, type = "scalar", encoding = "linear", wrap = "clamp")
#' # Alternating over/under threads form a fine height map, independent of the spill.
#' warp = cos(2 * pi * 48 * u)
#' weft = cos(2 * pi * 48 * v)
#' weave = .5 + .22 * warp + .22 * weft + .06 * warp * weft
#' fabric = texture_checker(
#'   "#d5ac78",
#'   "#698c95",
#'   texture_coordinates("uv", scale = c(24, 2))
#' )
#' rug = openpbr(
#'   base_color = texture_mix(fabric, "#35190c", mask),
#'   specular_roughness = texture_mix(.85, .18, mask)
#' )
#' rug_mesh = rayvertex::xz_rect_mesh(scale = c(3, 1, 2)) |>
#'   rayvertex::subdivide_mesh(subdivision_levels = 8, simple = TRUE) |>
#'   rayvertex::displace_mesh(
#'     array(weave, c(dim(weave), 3)),
#'     displacement_scale = .012,
#'     verbose = FALSE
#'   )
#' scene = generate_ground(depth = -.07, material = diffuse("#303843")) |>
#'   add_object(cube(
#'     y = -.025,
#'     xwidth = 3,
#'     ywidth = .05,
#'     zwidth = 2,
#'     material = diffuse("#b89668")
#'   )) |>
#'   add_object(raymesh_model(rug_mesh, y = .002, material = rug))
#' for (x in seq(-1.45, 1.45, length.out = 35)) {
#'   for (side in c(-1, 1)) {
#'     scene = add_object(
#'       scene,
#'       segment(
#'         start = c(x, -.015, side),
#'         end = c(x + .02 * sin(12 * x), -.035, side * 1.16),
#'         radius = .008,
#'         material = diffuse("#d5ac78")
#'       )
#'     )
#'   }
#' }
#'
#' scene = scene |>
#'   add_infinite_light(disk_light(
#'     direction = c(-1, .7, -.7),
#'     angular_diameter = 25,
#'     intensity = 16,
#'     color = "#fff0d8",
#'     name = "key"
#'   )) |>
#'   add_infinite_light(disk_light(
#'     direction = c(1, .4, -.7),
#'     angular_diameter = 45,
#'     intensity = 2,
#'     color = "#adcfff",
#'     name = "rim"
#'   ))
#' render_scene(
#'   scene,
#'   lookfrom = c(2.5, 3.7, 4.5),
#'   lookat = c(0, 0, 0),
#'   fov = 0,
#'   ortho_dimensions = c(4.1, 3.3),
#'   width = 500,
#'   height = 400,
#'   samples = 16,
#'   denoise = TRUE
#' )
#' unlink(file)
texture_mix = function(a, b, weight = 0.5) {
  a = texture_constant(a)
  b = texture_constant(b)
  weight = texture_scalar(weight)
  type = if (a$type == "color" || b$type == "color") "color" else "scalar"
  if (weight$op == "constant") {
    w = pmin(1, pmax(0, weight$value))
    if (a$op == "constant" && b$op == "constant") {
      return(texture_constant((1 - w) * a$value + w * b$value))
    }
  }
  texture_node("mix", type, a = a, b = b, weight = weight)
}

#' @rdname textures
#' @export
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # texture_direction_mix(): moss chooses a side of a tiny mushroom forest.
#' # The first color faces the direction; absolute=FALSE leaves the opposite side bare.
#' moss = texture_direction_mix(
#'   "#25752e",
#'   "#763713",
#'   direction = c(-1, .5, 1),
#'   absolute = FALSE,
#'   space = "world"
#' )
#' scene = generate_ground(depth = 0, material = diffuse("#d8c8a5"))
#' for (i in 1:3) {
#'   x = c(-1.1, 0, 1.1)[i]
#'   height = c(.7, 1.2, .85)[i]
#'   scene = scene |>
#'     add_object(cylinder(
#'       x = x,
#'       y = height / 2,
#'       radius = .12,
#'       length = height,
#'       material = diffuse("#f3dfb7")
#'     )) |>
#'     add_object(ellipsoid(
#'       x = x,
#'       y = height,
#'       a = .55,
#'       b = .28,
#'       c = .55,
#'       material = diffuse(moss)
#'     ))
#' }
#' scene = scene |>
#'   add_infinite_light(disk_light(
#'     direction = c(-1, 1.2, 1),
#'     angular_diameter = 35,
#'     intensity = 12,
#'     color = "#fff0d8",
#'     name = "key"
#'   )) |>
#'   add_infinite_light(disk_light(
#'     direction = c(1, .4, -.7),
#'     angular_diameter = 45,
#'     intensity = .5,
#'     color = "#adcfff",
#'     name = "rim"
#'   ))
#' render_scene(
#'   scene,
#'   lookfrom = c(4, 3, 7),
#'   lookat = c(0, 0.65, 0),
#'   fov = 0,
#'   ortho_dimensions = c(4.8, 3.8),
#'   width = 420,
#'   height = 340,
#'   samples = 16,
#'   denoise = TRUE
#' )
texture_direction_mix = function(
  a,
  b,
  direction = c(0, 1, 0),
  absolute = TRUE,
  space = "world"
) {
  space = match.arg(space, c("world", "object"))
  if (
    !is.numeric(direction) ||
      length(direction) != 3L ||
      any(!is.finite(direction)) ||
      !any(direction != 0)
  ) {
    stop("`direction` must be a finite nonzero three-vector.", call. = FALSE)
  }
  if (!is.logical(absolute) || length(absolute) != 1L || is.na(absolute)) {
    stop("`absolute` must be TRUE or FALSE.", call. = FALSE)
  }
  direction = direction / max(abs(direction))
  direction = direction / sqrt(sum(direction^2))
  weight = texture_node(
    "direction",
    "scalar",
    direction = direction,
    absolute = absolute,
    space = space
  )
  texture_mix(b, a, weight)
}

#' @rdname textures
#' @export
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # texture_noise(): a little wooden sailboat.
#' # Stretch the noise coordinates across the grain, leaving long fibers along x.
#' # A fixed seed reproduces the same piece of timber on every render.
#' grain = texture_noise(
#'   coordinates = texture_coordinates("object", scale = c(1.5, 35, 35)),
#'   octaves = 2,
#'   seed = 17
#' )
#' wood = openpbr(
#'   base_color = texture_mix("#251005", "#e6a25e", grain),
#'   specular_roughness = .32,
#'   coat_weight = .25
#' )
#' scene = generate_ground(depth = 0, material = diffuse("#386373")) |>
#'   add_object(ellipsoid(y = .24, a = 1.25, b = .24, c = .44, material = wood)) |>
#'   add_object(cylinder(
#'     y = 1.05,
#'     radius = .035,
#'     length = 1.5,
#'     material = wood
#'   )) |>
#'   add_object(triangle(
#'     v1 = c(.06, .6, 0),
#'     v2 = c(1.04, .6, 0),
#'     v3 = c(.06, 1.8, 0),
#'     material = diffuse("#fff0cb")
#'   )) |>
#'   add_object(triangle(
#'     v1 = c(-.06, .6, 0),
#'     v2 = c(-.75, .6, 0),
#'     v3 = c(-.06, 1.4, 0),
#'     material = diffuse("#bb4b3a")
#'   ))
#'
#' scene = scene |>
#'   add_infinite_light(disk_light(
#'     direction = c(-1, 1.2, 1),
#'     angular_diameter = 35,
#'     intensity = 12,
#'     color = "#fff0d8",
#'     name = "key"
#'   )) |>
#'   add_infinite_light(disk_light(
#'     direction = c(1, .4, -.7),
#'     angular_diameter = 45,
#'     intensity = .5,
#'     color = "#adcfff",
#'     name = "rim"
#'   ))
#' render_scene(
#'   scene,
#'   lookfrom = c(2.1, 1.9, 5),
#'   lookat = c(0, 0.8, 0),
#'   fov = 0,
#'   ortho_dimensions = c(3.5, 2.8),
#'   width = 500,
#'   height = 400,
#'   samples = 16,
#'   denoise = TRUE
#' )
texture_noise = function(
  scale = 1,
  coordinates = texture_coordinates(),
  octaves = 4L,
  seed = 0L
) {
  coordinates = texture_mapping(coordinates)
  if (
    !is.numeric(scale) ||
      length(scale) != 1L ||
      !is.finite(scale) ||
      scale <= 0 ||
      !is.numeric(octaves) ||
      length(octaves) != 1L ||
      !is.finite(octaves) ||
      octaves < 1 ||
      octaves > 16 ||
      octaves != floor(octaves) ||
      !is.numeric(seed) ||
      length(seed) != 1L ||
      !is.finite(seed) ||
      abs(seed) > .Machine$integer.max ||
      seed != floor(seed)
  ) {
    stop("Invalid noise scale, octave count, or integer seed.", call. = FALSE)
  }
  coordinates$scale = coordinates$scale * scale
  texture_node(
    "noise",
    "scalar",
    coordinates = coordinates,
    octaves = as.integer(octaves),
    seed = as.integer(seed)
  )
}

#' @rdname textures
#' @export
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # texture_checker(): a single checkerboard marble.
#' # UV coordinates wrap rows and columns around a curved surface.
#' checks = texture_checker(
#'   "#f8eaca",
#'   "#17494d",
#'   texture_coordinates("uv", scale = c(12, 6))
#' )
#' scene = generate_ground(depth = 0, material = diffuse("#d9bc98")) |>
#'   add_object(sphere(
#'     y = .85,
#'     radius = .85,
#'     material = openpbr(
#'       base_color = checks,
#'       specular_roughness = .3,
#'       coat_weight = .25
#'     )
#'   ))
#'
#' scene = scene |>
#'   add_infinite_light(disk_light(
#'     direction = c(-1, 1.2, 1),
#'     angular_diameter = 35,
#'     intensity = 12,
#'     color = "#fff0d8",
#'     name = "key"
#'   )) |>
#'   add_infinite_light(disk_light(
#'     direction = c(1, .4, -.7),
#'     angular_diameter = 45,
#'     intensity = .5,
#'     color = "#adcfff",
#'     name = "rim"
#'   ))
#' render_scene(
#'   scene,
#'   lookfrom = c(2, 1.8, 5),
#'   lookat = c(0, 0.7, 0),
#'   fov = 0,
#'   ortho_dimensions = c(2.8, 2.3),
#'   width = 500,
#'   height = 400,
#'   samples = 16,
#'   denoise = TRUE
#' )
texture_checker = function(a, b, coordinates = texture_coordinates()) {
  a = texture_constant(a)
  b = texture_constant(b)
  texture_node(
    "checker",
    if (a$type == "color" || b$type == "color") "color" else "scalar",
    a = a,
    b = b,
    coordinates = texture_mapping(coordinates)
  )
}

#' @rdname textures
#' @export
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # texture_gradient(): a sunset frozen into a popsicle.
#' # A world-space y mapping keeps the rounded top and block on the same gradient.
#' colors = texture_gradient(
#'   "#cc4f75",
#'   "#ffdb8a",
#'   axis = "y",
#'   coordinates = texture_coordinates("world", scale = .65, offset = c(0, -.1, 0))
#' )
#' scene = generate_ground(depth = -.1, material = diffuse("#bfd8d5")) |>
#'   add_object(cube(
#'     y = 1.15,
#'     xwidth = .95,
#'     ywidth = 1.2,
#'     zwidth = .32,
#'     material = diffuse(colors)
#'   )) |>
#'   add_object(ellipsoid(
#'     y = 1.75,
#'     a = .475,
#'     b = .3,
#'     c = .16,
#'     material = diffuse(colors)
#'   )) |>
#'   add_object(cube(
#'     y = .35,
#'     xwidth = .16,
#'     ywidth = .6,
#'     zwidth = .1,
#'     material = diffuse("#c99b65")
#'   ))
#' scene = scene |>
#'   add_infinite_light(disk_light(
#'     direction = c(-1, 1.2, 1),
#'     angular_diameter = 35,
#'     intensity = 12,
#'     color = "#fff0d8",
#'     name = "key"
#'   )) |>
#'   add_infinite_light(disk_light(
#'     direction = c(1, .4, -.7),
#'     angular_diameter = 45,
#'     intensity = .5,
#'     color = "#adcfff",
#'     name = "rim"
#'   ))
#' render_scene(
#'   scene,
#'   lookfrom = c(2, 1.6, 6),
#'   lookat = c(0, 1, 0),
#'   fov = 0,
#'   ortho_dimensions = c(3.3, 3.1),
#'   width = 420,
#'   height = 340,
#'   samples = 16,
#'   denoise = TRUE
#' )
texture_gradient = function(
  a,
  b,
  coordinates = texture_coordinates(),
  axis = "y"
) {
  axis = match.arg(axis, c("x", "y", "z"))
  weight = texture_node(
    "gradient",
    "scalar",
    coordinates = texture_mapping(coordinates),
    axis = match(axis, c("x", "y", "z")) - 1L
  )
  texture_mix(a, b, weight)
}

#' @rdname textures
#' @export
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # texture_image(): a framed reminder to Live, Laugh, Plot.
#' # An ordinary R plot becomes an sRGB texture, including its axes and script lettering.
#' file = tempfile(fileext = ".png")
#' grDevices::png(file, width = 800, height = 1000, pointsize = 28)
#' par(
#'   mar = c(7, 4, 3, 1),
#'   bg = "#fff8e8",
#'   col.axis = "#375665",
#'   col.lab = "#375665"
#' )
#' plot(
#'   mtcars$wt,
#'   mtcars$mpg,
#'   pch = 21,
#'   cex = 1.3,
#'   bg = "#438b98",
#'   col = "#244853",
#'   xlab = "Weight (1000 lbs)",
#'   ylab = "Miles per gallon",
#'   main = "Less is more"
#' )
#' abline(lm(mpg ~ wt, data = mtcars), col = "#d08152", lwd = 3)
#' mtext(
#'   "Live, Laugh, Plot",
#'   side = 1,
#'   line = 5,
#'   cex = 1.7,
#'   family = "HersheyScript"
#' )
#' dev.off()
#' print_texture = texture_image(file, encoding = "srgb", wrap = "clamp")
#' scene = generate_ground(depth = 0, material = diffuse("#d8d0bd")) |>
#'   add_object(cube(
#'     y = 1.26,
#'     z = -.06,
#'     xwidth = 2.02,
#'     ywidth = 2.52,
#'     zwidth = .12,
#'     material = diffuse("#744a2c")
#'   )) |>
#'   add_object(xy_rect(
#'     y = 1.26,
#'     z = .005,
#'     xwidth = 1.8,
#'     ywidth = 2.25,
#'     material = diffuse(print_texture)
#'   ))
#'
#' scene = scene |>
#'   add_infinite_light(disk_light(
#'     direction = c(-1, 1.2, 1),
#'     angular_diameter = 35,
#'     intensity = 12,
#'     color = "#fff0d8",
#'     name = "key"
#'   )) |>
#'   add_infinite_light(disk_light(
#'     direction = c(1, .4, -.7),
#'     angular_diameter = 45,
#'     intensity = .5,
#'     color = "#adcfff",
#'     name = "rim"
#'   ))
#' render_scene(
#'   scene,
#'   lookfrom = c(1, 2, 6),
#'   lookat = c(0, 1.3, 0),
#'   fov = 0,
#'   ortho_dimensions = c(3.5, 2.9),
#'   width = 500,
#'   height = 400,
#'   samples = 16,
#'   denoise = TRUE
#' )
#' unlink(file)
texture_image = function(
  filename,
  type = "color",
  encoding = "auto",
  coordinates = texture_coordinates("uv"),
  wrap = "repeat"
) {
  type = match.arg(type, c("color", "scalar"))
  encoding = match.arg(encoding, c("auto", "linear", "srgb"))
  wrap = match.arg(wrap, c("repeat", "clamp"))
  if (
    !is.character(filename) ||
      length(filename) != 1L ||
      is.na(filename) ||
      !file.exists(path.expand(filename))
  ) {
    stop("`filename` must name an existing texture image.", call. = FALSE)
  }
  if (encoding == "auto") {
    encoding = if (
      type == "scalar" ||
        tolower(tools::file_ext(filename)) %in% c("exr", "hdr")
    ) {
      "linear"
    } else {
      "srgb"
    }
  }
  texture_node(
    "image",
    type,
    filename = normalizePath(path.expand(filename)),
    encoding = encoding,
    coordinates = texture_mapping(coordinates),
    wrap = wrap
  )
}

#' @rdname textures
#' @export
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # texture_channel(): two secret stencils hidden in one RGB image.
#' # Red stores a heart; green stores a star. Extracting a channel gives a scalar mask.
#' a = seq(-1.25, 1.25, length.out = 256)
#' x = matrix(rep(a, each = 256), 256)
#' y = -t(x)
#' heart = (x^2 + y^2 - 1)^3 - x^2 * y^3 < 0
#' star = sqrt(x^2 + y^2) < .7 + .24 * cos(5 * (atan2(y, x) - pi / 2))
#' pixels = array(0, c(256, 256, 3))
#' pixels[,, 1] = heart
#' pixels[,, 2] = star
#' file = tempfile(fileext = ".png")
#' png::writePNG(pixels, file)
#' packed = texture_image(file, encoding = "linear", wrap = "clamp")
#' paints = list(
#'   packed,
#'   texture_mix("#143b52", "#ffdbc0", texture_channel(packed, "r")),
#'   texture_mix("#143b52", "#ffdbc0", texture_channel(packed, "g"))
#' )
#' scene = generate_ground(depth = 0, material = diffuse("#bec7c2"))
#' for (i in 1:3) {
#'   scene = scene |>
#'     add_object(cube(
#'       x = (i - 2) * 1.4,
#'       y = .8,
#'       xwidth = 1.2,
#'       ywidth = 1.4,
#'       zwidth = .1,
#'       material = diffuse("#e8d6b8")
#'     )) |>
#'     add_object(xy_rect(
#'       x = (i - 2) * 1.4,
#'       y = .85,
#'       z = .051,
#'       xwidth = 1.08,
#'       ywidth = 1.08,
#'       material = diffuse(paints[[i]])
#'     ))
#' }
#' labels = screen_text(
#'   c("RGB image", "Red: heart", "Green: star"),
#'   x = c(-1.4, 0, 1.4),
#'   y = .08,
#'   z = .1,
#'   size = 16,
#'   hjust = .5,
#'   color = "#152b3a"
#' )
#'
#' scene = scene |>
#'   add_infinite_light(disk_light(
#'     direction = c(-1, 1.2, 1),
#'     angular_diameter = 35,
#'     intensity = 12,
#'     color = "#fff0d8",
#'     name = "key"
#'   )) |>
#'   add_infinite_light(disk_light(
#'     direction = c(1, .4, -.7),
#'     angular_diameter = 45,
#'     intensity = .5,
#'     color = "#adcfff",
#'     name = "rim"
#'   ))
#' render_scene(
#'   scene,
#'   lookfrom = c(0, 2, 7),
#'   lookat = c(0, 0.7, 0),
#'   fov = 0,
#'   ortho_dimensions = c(4.9, 3.1),
#'   width = 500,
#'   height = 400,
#'   samples = 16,
#'   denoise = TRUE,
#'   screen_text = labels
#' )
#' unlink(file)
texture_channel = function(value, channel = "luminance") {
  value = texture_constant(value)
  channel = match.arg(channel, c("r", "g", "b", "average", "luminance"))
  if (value$type == "scalar") {
    return(value)
  }
  texture_node(
    "channel",
    "scalar",
    child = value,
    channel = match(channel, c("r", "g", "b", "average", "luminance")) - 1L
  )
}

#' @rdname textures
#' @export
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # texture_scale(): three candy-striped spinning tops, from dim paint to bright paint.
#' # Multiplication changes RGB values, not stripe spacing. Keep light and geometry fixed.
#' stripes = texture_checker(
#'   "#e78653",
#'   "#4c9fba",
#'   texture_coordinates("world", scale = c(0, 8, 0))
#' )
#' factors = c(.2, .5, 1)
#' scene = generate_ground(depth = 0, material = diffuse("#d8cbb6"))
#' for (i in 1:3) {
#'   x = (i - 2) * 1.35
#'   paint = diffuse(texture_scale(stripes, factors[i]))
#'   scene = scene |>
#'     add_object(cone(
#'       start = c(x, .7, 0),
#'       end = c(x, .02, 0),
#'       radius = .55,
#'       material = paint
#'     )) |>
#'     add_object(cylinder(
#'       x = x,
#'       y = .72,
#'       radius = .55,
#'       length = .06,
#'       material = paint
#'     )) |>
#'     add_object(cylinder(
#'       x = x,
#'       y = .94,
#'       radius = .06,
#'       length = .4,
#'       material = diffuse("#5a3621")
#'     ))
#' }
#' labels = screen_text(
#'   c("RGB x 0.2", "RGB x 0.5", "RGB x 1"),
#'   x = c(-1.35, 0, 1.35),
#'   y = 0,
#'   z = .7,
#'   size = 16,
#'   hjust = .5,
#'   color = "#20313e"
#' )
#'
#' scene = scene |>
#'   add_infinite_light(disk_light(
#'     direction = c(-1, 1.2, 1),
#'     angular_diameter = 35,
#'     intensity = 12,
#'     color = "#fff0d8",
#'     name = "key"
#'   )) |>
#'   add_infinite_light(disk_light(
#'     direction = c(1, .4, -.7),
#'     angular_diameter = 45,
#'     intensity = .5,
#'     color = "#adcfff",
#'     name = "rim"
#'   ))
#' render_scene(
#'   scene,
#'   lookfrom = c(1.5, 3, 7),
#'   lookat = c(0, 0.45, 0),
#'   fov = 0,
#'   ortho_dimensions = c(4.8, 3.1),
#'   width = 500,
#'   height = 400,
#'   samples = 16,
#'   denoise = TRUE,
#'   screen_text = labels
#' )
texture_scale = function(value, factor) {
  value = texture_constant(value)
  factor = texture_scalar(factor)
  if (value$op == "constant" && factor$op == "constant") {
    return(texture_constant(value$value * factor$value))
  }
  texture_node("scale", value$type, child = value, factor = factor)
}

#' @param op Operation name.
#' @param type Output type.
#' @param ... Node fields.
#' @return Typed texture descriptor.
#' @keywords internal
#' @noRd
texture_node = function(op, type, ...) {
  structure(c(list(op = op, type = type), list(...)), class = "ray_texture")
}

#' @param value Scalar value or texture.
#' @return Scalar texture, rejecting implicit color reduction.
#' @keywords internal
#' @noRd
texture_scalar = function(value) {
  value = texture_constant(value)
  if (value$type != "scalar") {
    stop(
      "Expected a scalar texture; use texture_channel() to reduce colors.",
      call. = FALSE
    )
  }
  value
}

#' @param coordinates Coordinate descriptor.
#' @return Validated coordinate descriptor.
#' @keywords internal
#' @noRd
texture_mapping = function(coordinates) {
  if (!inherits(coordinates, "ray_texture_coordinates")) {
    stop("Expected texture_coordinates().", call. = FALSE)
  }
  coordinates
}

#' @param image Legacy image input.
#' @param checker Legacy checker color.
#' @param noise Legacy noise scale.
#' @param gradient Legacy gradient color.
#' @return Nothing; rejects ambiguous color input combinations.
#' @keywords internal
#' @noRd
texture_check_legacy = function(image, checker, noise, gradient) {
  if (
    !identical(image, "") ||
      any(!is.na(checker)) ||
      noise != 0 ||
      any(!is.na(gradient))
  ) {
    stop(
      "Choose a color texture graph or legacy image/checker/noise/gradient arguments, not both.",
      call. = FALSE
    )
  }
}
