#' Pig Object
#'
#' Builds an earless, low-poly mesh pig from compact data included in rayrender.
#' The skiing outfit includes goggles, a scarf and skis; the spider variant has
#' eight legs, four eyes and procedural copper-colored hair curves.
#'
#' @md
#' @param x Default `0`. x-coordinate offset of the pig.
#' @param y Default `0`. y-coordinate offset of the pig.
#' @param z Default `0`. z-coordinate offset of the pig.
#' @param emotion Default `"neutral"`. Facial expression: `neutral`, `skeptical`,
#'   `worried`, `angry`, `surprised`, or `excited`. `cheerful` is an alias for
#'   `neutral`. When omitted with `ski = TRUE`, the expression is `excited`.
#' @param spider Default `FALSE`. Generate an eight-legged, four-eyed spider pig.
#' @param angle Default `c(0, 0, 0)`. Rotation in degrees about the x, y and z axes.
#' @param order_rotation Default `c(1, 2, 3)`. Order of the rotation axes.
#' @param scale Default `c(1, 1, 1)`. Uniform scalar or x/y/z scale factors.
#' @param diffuse_sigma Default `0`. Oren-Nayar roughness for diffuse surfaces,
#'   including the skin. Does not affect glossy eyes, lenses or hair.
#' @param ski Default `FALSE`. Add the ski outfit and skiing pose. If both `ski`
#'   and `spider` are `TRUE`, warns `"spiders can't ski"` and generates the spider
#'   without ski gear.
#' @param hair_count Default `12000`. Number of procedural spider hair curves.
#'   Set to zero for a bare spider mesh. Ignored for a regular or skiing pig.
#' @param hair_length Default `1`. Multiplier for spider hair lengths.
#' @param hair_seed Default `8`. Reproducible spider hair seed; preserves the
#'   caller's random-number state.
#' @param hair_material Default `NULL`. Use the built-in rust/copper hair palette.
#'   Supply a single [hair()] material to override the entire coat.
#'
#' @details Geometry, expression heads and the ski pose are assembled in memory;
#'   no OBJ paths, downloads, Blender or rgl installation are required. The
#'   character faces +X, uses Y up, and retains the rotation/scale pivot at
#'   `c(0, 1, 0)`. Feet or skis rest approximately at `y - 0.6` before rotation
#'   and scaling. The mask and goggles can obscure eyebrows in skiing poses.
#'   Hair uses tapered [bezier_curve()] objects with [hair()] materials. A
#'   finite bounce limit, such as `max_depth = 12`, is useful for dense fur.
#' @return A scene tibble containing mesh objects and, for a hairy spider, curves.
#' @export
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' generate_ground(depth = -0.6) |>
#'   add_object(pig()) |>
#'   render_scene(lookfrom = c(8, 5, 8), lookat = c(0, 1, 0), samples = 64)
#'
#' generate_ground(depth = -0.6, material = glossy(color="white")) |>
#'   add_object(pig(ski = TRUE)) |>
#'   render_scene(lookfrom = c(8, 5, 8),
#'                lookat = c(0, 1, 0),
#'                fov=40,samples = 64)
#'
#' generate_ground(depth = -0.6) |>
#'   add_object(pig(spider = TRUE)) |>
#'   render_scene(lookfrom = c(10, 7, 8), lookat = c(0, 1, 0),
#'                samples = 64, max_depth = 12)
pig = function(
  x = 0,
  y = 0,
  z = 0,
  emotion = "neutral",
  spider = FALSE,
  angle = c(0, 0, 0),
  order_rotation = c(1, 2, 3),
  scale = c(1, 1, 1),
  diffuse_sigma = 0,
  ski = FALSE,
  hair_count = 12000,
  hair_length = 1,
  hair_seed = 8,
  hair_material = NULL
) {
  scalar = function(value) {
    is.numeric(value) && length(value) == 1L && is.finite(value)
  }
  flag = function(value) {
    is.logical(value) && length(value) == 1L && !is.na(value)
  }
  if (!flag(ski) || !flag(spider)) {
    stop("ski and spider must each be TRUE or FALSE.", call. = FALSE)
  }
  if (ski && spider) {
    warning("spiders can't ski", call. = FALSE)
    ski = FALSE
  }
  if (missing(emotion) && ski) {
    emotion = "excited"
  }
  emotion = match.arg(
    emotion,
    c(
      "neutral",
      "skeptical",
      "worried",
      "angry",
      "surprised",
      "excited",
      "cheerful"
    )
  )
  if (emotion == "cheerful") {
    emotion = "neutral"
  }
  if (!all(vapply(list(x, y, z), scalar, logical(1)))) {
    stop("x, y and z must be finite numeric scalars.", call. = FALSE)
  }
  if (!is.numeric(angle) || length(angle) != 3L || any(!is.finite(angle))) {
    stop("angle must contain three finite numbers.", call. = FALSE)
  }
  if (
    !is.numeric(scale) ||
      !length(scale) %in% c(1L, 3L) ||
      any(!is.finite(scale)) ||
      any(scale == 0)
  ) {
    stop(
      "scale must contain one or three finite nonzero numbers.",
      call. = FALSE
    )
  }
  if (
    !is.numeric(order_rotation) ||
      length(order_rotation) != 3L ||
      !identical(sort(as.numeric(order_rotation)), c(1, 2, 3))
  ) {
    stop("order_rotation must be a permutation of c(1, 2, 3).", call. = FALSE)
  }
  if (!scalar(diffuse_sigma) || diffuse_sigma < 0) {
    stop("diffuse_sigma must be a finite nonnegative number.", call. = FALSE)
  }
  if (spider) {
    if (
      !scalar(hair_count) ||
        hair_count < 0 ||
        hair_count != floor(hair_count) ||
        hair_count > .Machine$integer.max
    ) {
      stop("hair_count must be a nonnegative integer.", call. = FALSE)
    }
    if (!scalar(hair_length) || hair_length <= 0) {
      stop("hair_length must be a finite positive number.", call. = FALSE)
    }
    if (
      !scalar(hair_seed) ||
        hair_seed < 0 ||
        hair_seed > .Machine$integer.max ||
        hair_seed != floor(hair_seed)
    ) {
      stop("hair_seed must be a nonnegative integer seed.", call. = FALSE)
    }
    if (
      !is.null(hair_material) &&
        (!inherits(hair_material, "ray_material") ||
          length(hair_material) != 1L ||
          !identical(hair_material[[1]]$type, hair()[[1]]$type))
    ) {
      stop("hair_material must be a single hair() material.", call. = FALSE)
    }
  }
  assets = pig_mesh_data()
  parts = pig_mesh_parts(assets, emotion, ski, spider)
  objects = list()
  for (part in parts) {
    for (id in unique(part$materials)) {
      mesh = structure(
        list(
          vb = t(part$vertices),
          normals = t(part$normals),
          it = t(part$faces[part$materials == id, , drop = FALSE])
        ),
        class = "mesh3d"
      )
      objects[[length(objects) + 1L]] = mesh3d_model(
        mesh,
        override_material = TRUE,
        material = pig_surface_material(
          part$name,
          assets$materials[[id]],
          diffuse_sigma
        )
      )
    }
  }
  if (spider && hair_count > 0) {
    guides = pig_fur_guides(
      parts,
      hair_count,
      hair_seed,
      hair_length,
      assets$floor
    )
    palette = if (is.null(hair_material)) {
      lapply(c("#B80301", "#D80702", "#E51804"), function(color) {
        hair(color = color, eta = 1.03, beta_m = .65, beta_n = .65)
      })
    } else {
      rep(list(hair_material), 3)
    }
    objects = c(
      objects,
      lapply(seq_len(nrow(guides)), function(i) {
        g = guides[i, ]
        bezier_curve(
          p1 = g[1:3],
          p2 = g[4:6],
          p3 = g[7:9],
          p4 = g[10:12],
          width = g[13],
          width_end = .00015,
          type = "flat",
          material = palette[[g[14]]]
        )
      })
    )
  }
  scene = do.call(vctrs::vec_rbind, objects)
  # All geometry is assembled in one coordinate system. Recycle one common
  # transform rather than repeatedly copying a tibble containing 12,000 curves.
  transform = group_objects(
    scene[1, ],
    translate = c(x, y, z),
    angle = angle,
    order_rotation = order_rotation,
    scale = scale,
    pivot_point = c(0, 1, 0)
  )$transforms
  scene$transforms = rep(transform, nrow(scene))
  scene
}

# Private, lazily loaded package data. Geometry pools share repeated topology,
# normals and vertices; ski clothes are stored separately from the posed body.
#' Read the shared mesh pig data
#' @keywords internal
#' @noRd
pig_mesh_data = local({
  cached = NULL
  function() {
    if (is.null(cached)) {
      path = system.file("extdata", "pig.rds", package = "rayrender")
      if (!nzchar(path)) {
        stop(
          "The installed rayrender package is missing its pig mesh data.",
          call. = FALSE
        )
      }
      cached <<- readRDS(path)
    }
    cached
  }
})

#' Assemble the selected pig meshes and expression
#' @param assets Decoded package mesh data.
#' @param emotion Facial expression.
#' @param ski Whether to use the skiing pose and outfit.
#' @param spider Whether to use the spider body and eyes.
#' @keywords internal
#' @noRd
pig_mesh_parts = function(assets, emotion, ski, spider) {
  decode = function(part) {
    part$vertices = assets$v[[part$v]] / 10000
    part$normals = assets$n[[part$n]] / 1000
    part$normals = part$normals / sqrt(rowSums(part$normals^2))
    part$faces = assets$f[[part$f]]
    part
  }
  head_name = switch(
    emotion,
    worried = "surprised",
    angry = "skeptical",
    emotion
  )
  body = if (spider) {
    assets$spider_body
  } else if (ski) {
    assets$ski_body
  } else {
    assets$body
  }
  if (spider && emotion == "neutral") {
    head = lapply(assets$spider_head, decode)
  } else {
    head = lapply(assets$heads[[head_name]], decode)
    if (spider) {
      secondary = Filter(
        function(p) grepl("secondary", p$name),
        assets$spider_head
      )
      head = c(head, lapply(secondary, decode))
      for (i in seq_along(head)) {
        p = head[[i]]
        # Reuse the exact spider material assignments where compatible.
        source_name = sub("^Spider_", "", p$name)
        reference = Filter(
          function(q) {
            sub("\\.[0-9]+$", "", sub("^Spider_", "", q$name)) ==
              sub("\\.[0-9]+$", "", source_name)
          },
          assets$spider_head
        )
        if (length(reference)) {
          if (length(unique(reference[[1]]$materials)) == 1L) {
            p$materials[] = reference[[1]]$materials[1]
          } else if (length(p$materials) == length(reference[[1]]$materials)) {
            p$materials = reference[[1]]$materials
          }
        }
        if (grepl("^Brow_", p$name)) {
          p$vertices[, 1] = p$vertices[, 1] - .13
          p$vertices[, 2] = p$vertices[, 2] - .035
          p$vertices[, 3] = p$vertices[, 3] + sign(mean(p$vertices[, 3])) * .13
        }
        head[[i]] = p
      }
    }
  }
  for (i in seq_along(head)) {
    p = head[[i]]
    if (grepl("Brow_", p$name) && emotion %in% c("worried", "angry")) {
      center = colMeans(p$vertices)
      slope = if (emotion == "worried") -.7 else .9
      slope = slope * sign(center[3])
      p$vertices[, 2] = p$vertices[, 2] +
        slope * (p$vertices[, 3] - center[3]) +
        if (emotion == "worried") .02 else -.08
      p$normals[, 3] = p$normals[, 3] - slope * p$normals[, 2]
      p$normals = p$normals / sqrt(rowSums(p$normals^2))
    }
    if (ski) {
      transform = assets$ski_head_transform
      p$vertices = t(transform[1:3, ] %*% rbind(t(p$vertices), 1))
      p$normals = p$normals %*% t(transform[1:3, 1:3])
    } else {
      p$vertices = sweep(p$vertices, 2, assets$pivot, "+")
    }
    head[[i]] = p
  }
  if (ski) {
    head = Filter(function(p) !grepl("Brow_", p$name), head)
  }
  c(lapply(body, decode), head, if (ski) lapply(assets$ski_gear, decode))
}

#' Assign native rayrender materials to pig surfaces
#' @param name Mesh component name.
#' @param definition Stored surface color and material properties.
#' @param sigma Diffuse roughness in degrees.
#' @keywords internal
#' @noRd
pig_surface_material = function(name, definition, sigma) {
  color = definition$Kd
  if (grepl("Goggles_amber_lens", name)) {
    dielectric(
      color = "white",
      refraction = 1.08,
      attenuation = c(.04, .11, .3),
      attenuation_intensity = .15
    )
  } else if (
    grepl("(Eye_.*(dome|iris|pupil)|secondary_(ruby_eye|iris|pupil))", name)
  ) {
    glossy(color = color, gloss = .75, reflectance = .045)
  } else if (grepl("(golden_frame|gold_release|strap_adjuster)", name)) {
    glossy(color = color, gloss = .7, reflectance = .2)
  } else {
    diffuse(color = color, sigma = sigma)
  }
}

#' Sample cubic hair guides on the selected spider surface
#' @param parts Assembled mesh components.
#' @param n Number of guides.
#' @param seed Local random seed.
#' @param length_scale Hair length multiplier.
#' @param floor Foot contact height.
#' @keywords internal
#' @noRd
pig_fur_guides = function(parts, n, seed, length_scale, floor) {
  surfaces = Filter(
    function(p) grepl("(^Spider_skin|(^|_)Head_earless_sculpt)", p$name),
    parts
  )
  triangles = lapply(surfaces, function(p) {
    f = p$faces
    cbind(
      p$vertices[f[, 1], ],
      p$vertices[f[, 2], ],
      p$vertices[f[, 3], ],
      p$normals[f[, 1], ],
      p$normals[f[, 2], ],
      p$normals[f[, 3], ],
      grepl("Head_", p$name)
    )
  })
  s = do.call(rbind, triangles)
  ab = s[, 4:6] - s[, 1:3]
  ac = s[, 7:9] - s[, 1:3]
  cross = cbind(
    ab[, 2] * ac[, 3] - ab[, 3] * ac[, 2],
    ab[, 3] * ac[, 1] - ab[, 1] * ac[, 3],
    ab[, 1] * ac[, 2] - ab[, 2] * ac[, 1]
  )
  area = sqrt(rowSums(cross^2)) / 2
  withr::with_seed(seed, {
    batches = list()
    total = 0L
    while (total < n) {
      count = max(100L, ceiling((n - total) * 1.35))
      index = sample.int(nrow(s), count, replace = TRUE, prob = area)
      u = sqrt(stats::runif(count))
      v = stats::runif(count)
      w1 = 1 - u
      w2 = u * (1 - v)
      w3 = u * v
      root = s[index, 1:3] * w1 + s[index, 4:6] * w2 + s[index, 7:9] * w3
      normal = s[index, 10:12] *
        w1 +
        s[index, 13:15] * w2 +
        s[index, 16:18] * w3
      normal = normal / sqrt(rowSums(normal^2))
      head = s[index, 19] == 1
      keep = root[, 2] > floor + .15 &
        (!head | root[, 1] < 1.6 | (abs(root[, 3]) > .64 & root[, 1] < 1.86))
      root = root[keep, , drop = FALSE]
      normal = normal[keep, , drop = FALSE]
      head = head[keep]
      count = nrow(root)
      if (!count) {
        next
      }
      strand_length = ifelse(
        head,
        .09,
        ifelse(abs(root[, 3]) > .98, .105, .17)
      ) *
        stats::runif(count, .65, 1.35) *
        length_scale
      guard = stats::runif(count) < .13
      strand_length[guard] = strand_length[guard] * 1.5
      tangent = matrix(stats::rnorm(count * 3), ncol = 3)
      tangent = tangent - normal * rowSums(tangent * normal)
      tangent = tangent / pmax(sqrt(rowSums(tangent^2)), 1e-8)
      sweep = tangent *
        .23 +
        matrix(rep(c(-.17, -.15, 0), each = count), ncol = 3)
      p1 = root + normal * .0008
      p2 = p1 + normal * strand_length * .34
      p3 = p1 + (normal * .72 + sweep * .45) * strand_length
      p4 = p1 + (normal + sweep) * strand_length
      width = stats::runif(count, .0032, .0055) * ifelse(guard, .85, 1)
      shade = sample.int(3, count, replace = TRUE, prob = c(.62, .33, .05))
      batches[[length(batches) + 1L]] = cbind(p1, p2, p3, p4, width, shade)
      total = total + count
    }
    do.call(rbind, batches)[seq_len(n), , drop = FALSE]
  })
}
