#' Texture a Mesh Without a UV Atlas Using Ptex
#'
#' Use Ptex when different parts of a mesh need unique surface detail and
#' arranging them into a shared 2D texture image would be cumbersome. For
#' example, localized weathering on an irregular sculpture or
#' a different painted pattern on each panel of a patchwork object. Ptex gives
#' each mesh face its own image, avoiding the need to unwrap and pack the whole
#' mesh into a UV atlas. This function reads those images from a Ptex file made
#' for the mesh, using the \pkg{ptex} package.
#'
#' @param filename Existing Ptex file containing one or three channels.
#' @param type Default `"color"`. Output type, either `"color"` or `"scalar"`.
#' Scalar output averages the decoded channels; single-channel files broadcast
#' to RGB when used as a color.
#' @param encoding Default `"linear"`. `"linear"`, `"srgb"`, or `"gamma N"`
#' with a positive exponent, for example `"gamma 2.2"`.
#' @param filter Default `"bspline"`. Ptex reconstruction filter: `"point"`,
#' `"bilinear"`, `"box"`, `"gaussian"`, `"bicubic"`, `"bspline"`,
#' `"catmullrom"`, or `"mitchell"`.
#' @details
#' Face images can have different resolutions, letting an artist allocate more
#' pixels to detailed areas and fewer to plain surfaces. Adjacency stored in the
#' file lets filtering sample neighboring faces across their shared edges,
#' avoiding the need to pad separate islands in an image atlas for filtering.
#'
#' Choose [texture_image()] for a conventional image mapped with mesh UVs,
#' such as a photograph, label, or repeating fabric pattern. `texture_ptex()`
#' reads an existing Ptex file made for the mesh's face layout; it does not
#' convert a PNG or JPEG into Ptex or automatically paint an untextured model.
#' Ptex removes the need for a global UV atlas, but rayrender still needs the
#' correct face IDs and coordinates within each face, as described below.
#'
#' Use this node in supported color and roughness inputs, or compose it with
#' [texture_mix()], [texture_scale()], and [texture_channel()]. Descriptors
#' contain filenames, not open handles, and can be serialized before rendering.
#'
#' Mesh UVs must be face-local coordinates in `[0, 1]` (and `u + v <= 1` for
#' triangular Ptex files). A rayvertex mesh can store a zero-based integer
#' `ptex_face_indices` vector in each `shapes` entry, with one value per triangle.
#' Triangles made from the same source quad share its Ptex face ID and UV chart.
#' Without explicit IDs, face zero is used, matching PBRT. Analytic primitives
#' have no Ptex face. Subdivision of an addressed mesh is currently rejected
#' because it would require propagating its original face parameterization.
#'
#' [read_pbrt()] preserves `faceIndices` and PLY `face_indices`, including in
#' lazy PLY storage, and imports Ptex color, roughness, and material bump maps.
#' PBRT's default encoding is `"gamma 2.2"`, including scalar files. Imported
#' nonlinear Ptex values use PBRT's 8-bit quantization before decoding; this
#' constructor decodes without that extra quantization. Material displacement
#' changes shading as a bump map; it does not displace the geometry.
#'
#' Filtering uses ray UV differentials when available. Footprints larger than a
#' face are uniformly reduced to fit the runtime's one-face limit. Lookup
#' centers are clamped to the face, including finite-difference bump probes;
#' filter taps still use Ptex adjacency. There is no additional image v flip.
#'
#' Each render shares one cache across its workers. Set
#' `options(rayrender.ptex_cache_memory = 256 * 1024^2,
#' rayrender.ptex_cache_files = 100)` to change its soft eviction limits (bytes
#' and open files). Workers never call R. Failed lookups return zero and produce
#' an aggregated warning after rendering; missing files fail during setup.
#' Verbose rendering reports cache statistics. The provider stays loaded for
#' the render; do not explicitly unload its DLL while rendering.
#' @return A composable `ray_texture` descriptor.
#' @importFrom ptex ptex_api
#' @export
#' @md
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # Two tiny quilt panels, each painted in its own Ptex face.
#' file = tempfile(fileext = ".ptx")
#' warm = cool = array(0, c(16, 16, 3))
#' grid = outer(1:16, 1:16, function(u, v) (u %/% 4 + v %/% 4) %% 2)
#' warm[,, 1] = .8
#' warm[,, 2] = .15 + .5 * grid
#' cool[,, 2] = .3 + .5 * grid
#' cool[,, 3] = .7
#' ptex::ptex_write(file, list(warm, cool))
#' mesh = rayvertex::construct_mesh(
#'   vertices = rbind(c(-1, -1, 0), c(0, -1, 0), c(0, 1, 0), c(-1, 1, 0),
#'                    c(0, -1, 0), c(1, -1, 0), c(1, 1, 0), c(0, 1, 0)),
#'   indices = rbind(c(0, 1, 2), c(0, 2, 3), c(4, 5, 6), c(4, 6, 7)),
#'   tex_indices = rbind(c(0, 1, 2), c(0, 2, 3), c(4, 5, 6), c(4, 6, 7)),
#'   texcoords = rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1),
#'                    c(0, 0), c(1, 0), c(1, 1), c(0, 1))
#' )
#' mesh$shapes[[1]]$ptex_face_indices = c(0L, 0L, 1L, 1L)
#' scene = raymesh_model(mesh, material = diffuse(texture_ptex(file))) |>
#'   add_infinite_light(disk_light(direction = c(-1, 1, 2),
#'                                angular_diameter = 30, intensity = 20))
#' render_scene(scene, lookfrom = c(2, 1, 5), lookat = c(0, 0, 0),
#'              fov = 30, samples = 16, denoise = TRUE)
#' unlink(file)
texture_ptex = function(
  filename,
  type = "color",
  encoding = "linear",
  filter = "bspline"
) {
  type = match.arg(type, c("color", "scalar"))
  filters = c(
    "point",
    "bilinear",
    "box",
    "gaussian",
    "bicubic",
    "bspline",
    "catmullrom",
    "mitchell"
  )
  filter = match.arg(filter, filters)
  if (
    !is.character(filename) ||
      length(filename) != 1L ||
      is.na(filename) ||
      !file.exists(path.expand(filename))
  ) {
    stop("`filename` must name an existing Ptex file.", call. = FALSE)
  }
  if (!is.character(encoding) || length(encoding) != 1L || is.na(encoding)) {
    stop(
      "`encoding` must be linear, srgb, or gamma followed by a positive exponent.",
      call. = FALSE
    )
  }
  encoding = tolower(trimws(encoding))
  gamma = if (encoding == "gamma") NA_real_ else 1
  if (grepl("^gamma[[:space:]]+", encoding)) {
    gamma = suppressWarnings(as.numeric(sub(
      "^gamma[[:space:]]+",
      "",
      encoding
    )))
    encoding = "gamma"
  }
  if (
    !encoding %in% c("linear", "srgb", "gamma") ||
      !is.finite(gamma) ||
      gamma <= 0
  ) {
    stop(
      "`encoding` must be linear, srgb, or gamma followed by a positive exponent.",
      call. = FALSE
    )
  }
  texture_node(
    "ptex",
    type,
    filename = normalizePath(path.expand(filename)),
    encoding = encoding,
    gamma = gamma,
    filter = match(filter, filters) - 1L,
    pbrt_encoding = FALSE
  )
}
