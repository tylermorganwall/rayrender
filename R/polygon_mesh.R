#' Signed polygon area in local coordinates
#' @keywords internal
#' @noRd
polygon_ring_area = function(p) {
  q = sweep(p, 2, p[1, ])
  size = max(abs(q))
  if (!is.finite(size) || size == 0) {
    return(0)
  }
  q = q / size
  j = c(seq_len(nrow(q))[-1L], 1L)
  sum(q[, 1] * q[j, 2] - q[j, 1] * q[, 2]) / 2
}

#' Normalize an extrusion ring without its closing vertex
#' @keywords internal
#' @noRd
polygon_ring = function(p) {
  if (
    !is.matrix(p) ||
      !is.numeric(p) ||
      ncol(p) < 2L ||
      nrow(p) < 3L ||
      any(!is.finite(p[, 1:2, drop = FALSE]))
  ) {
    stop(
      "Each polygon ring must contain at least three finite x/y vertices.",
      call. = FALSE
    )
  }
  p = unname(p[, 1:2, drop = FALSE])
  p = p[
    c(TRUE, rowSums(p[-1L, , drop = FALSE] != p[-nrow(p), , drop = FALSE]) > 0),
    ,
    drop = FALSE
  ]
  if (nrow(p) > 1L && all(p[1, ] == p[nrow(p), ])) {
    p = p[-nrow(p), , drop = FALSE]
  }
  if (nrow(p) < 3L) {
    stop(
      "Each polygon ring needs at least three distinct vertices.",
      call. = FALSE
    )
  }
  # Remove straight, forward-going redundant points from caps and walls alike.
  # A reversal along an edge is invalid, not a removable collinear vertex.
  repeat {
    n = nrow(p)
    q = sweep(p, 2, p[1, ])
    size = max(abs(q))
    if (!is.finite(size) || size == 0) {
      break
    }
    q = q / size
    a = q - q[c(n, seq_len(n - 1L)), , drop = FALSE]
    b = q[c(seq_len(n)[-1L], 1L), , drop = FALSE] - q
    straight = a[, 1] * b[, 2] == a[, 2] * b[, 1]
    if (any(straight & rowSums(a * b) < 0)) {
      stop(
        "Polygon rings must not self-intersect or retrace an edge.",
        call. = FALSE
      )
    }
    remove = straight & rowSums(a * b) > 0
    if (!any(remove) || sum(!remove) < 3L) {
      break
    }
    p = p[!remove, , drop = FALSE]
  }
  if (polygon_ring_area(p) == 0) {
    stop(
      "Polygon rings must have nonzero area and must not self-intersect.",
      call. = FALSE
    )
  }
  p
}

#' Point containment for disjoint simple polygon rings
#' @keywords internal
#' @noRd
polygon_contains = function(ring, point) {
  p = sweep(ring, 2, point)
  j = c(seq_len(nrow(p))[-1L], 1L)
  i = which((p[, 2] > 0) != (p[j, 2] > 0))
  if (!length(i)) {
    return(FALSE)
  }
  x = p[i, 1] + (p[j[i], 1] - p[i, 1]) * (-p[i, 2] / (p[j[i], 2] - p[i, 2]))
  sum(x > 0) %% 2L == 1L
}

#' Reject crossing or touching extrusion ring edges
#' @keywords internal
#' @noRd
polygon_check_edges = function(rings) {
  p = do.call(rbind, rings)
  p = sweep(p, 2, p[1, ])
  size = max(abs(p))
  if (!is.finite(size) || size == 0) {
    stop("Polygon coordinate extent must be finite and nonzero.", call. = FALSE)
  }
  p = p / size
  lengths = vapply(rings, nrow, 0L)
  starts = c(0L, utils::head(cumsum(lengths), -1L))
  next_id = unlist(Map(
    function(n, start) c(seq_len(n)[-1L], 1L) + start,
    lengths,
    starts
  ))
  q = p[next_id, , drop = FALSE]
  xmin = pmin(p[, 1], q[, 1])
  xmax = pmax(p[, 1], q[, 1])
  ymin = pmin(p[, 2], q[, 2])
  ymax = pmax(p[, 2], q[, 2])
  cross_sign = function(a, b) {
    left = a[, 1] * b[, 2]
    right = a[, 2] * b[, 1]
    value = left - right
    value[abs(value) <= 8 * .Machine$double.eps * (abs(left) + abs(right))] = 0
    sign(value)
  }
  active = integer()
  for (i in order(xmin)) {
    active = active[xmax[active] >= xmin[i]]
    candidates = active[
      ymin[active] <= ymax[i] &
        ymax[active] >= ymin[i] &
        next_id[active] != i &
        active != next_id[i]
    ]
    if (length(candidates)) {
      a = p[candidates, , drop = FALSE]
      b = q[candidates, , drop = FALSE]
      c = matrix(p[i, ], nrow(a), 2L, byrow = TRUE)
      d = matrix(q[i, ], nrow(a), 2L, byrow = TRUE)
      if (
        any(
          cross_sign(b - a, c - a) * cross_sign(b - a, d - a) <= 0 &
            cross_sign(d - c, a - c) * cross_sign(d - c, b - c) <= 0
        )
      ) {
        stop(
          "Polygon rings must not self-intersect, cross, or touch other rings.",
          call. = FALSE
        )
      }
    }
    active = c(active, i)
  }
  invisible(NULL)
}

#' Validate polygon holes and choose consistent ring winding
#' @keywords internal
#' @noRd
polygon_component = function(rings) {
  if (!length(rings)) {
    stop("Empty polygons cannot be extruded.", call. = FALSE)
  }
  rings = lapply(rings, polygon_ring)
  polygon_check_edges(rings)
  if (length(rings) > 1L) {
    for (i in 2:length(rings)) {
      if (!polygon_contains(rings[[1]], rings[[i]][1, ])) {
        stop(
          "Every hole must lie strictly inside its outer polygon.",
          call. = FALSE
        )
      }
      for (j in seq_len(i - 1L)[-1L]) {
        if (
          polygon_contains(rings[[j]], rings[[i]][1, ]) ||
            polygon_contains(rings[[i]], rings[[j]][1, ])
        ) {
          stop(
            "Polygon holes must not overlap or contain other holes.",
            call. = FALSE
          )
        }
      }
    }
  }
  for (i in seq_along(rings)) {
    if ((polygon_ring_area(rings[[i]]) > 0) != (i == 1L)) {
      rings[[i]] = rings[[i]][c(1L, nrow(rings[[i]]):2L), , drop = FALSE]
    }
  }
  rings
}

#' Resolve polygon input into components with per-feature heights
#' @keywords internal
#' @noRd
polygon_components = function(
  polygon,
  holes,
  top,
  bottom,
  data_column_top,
  data_column_bottom,
  scale_data
) {
  spatial = inherits(polygon, "SpatialPolygons")
  is_sf = inherits(polygon, "sf")
  if ((spatial || is_sf) && !is.null(holes)) {
    warning("`holes` is unused when input is sf or Spatial.", call. = FALSE)
  }
  components = list()
  append = function(rings, high, low) {
    components[[length(components) + 1L]] <<- list(
      rings = rings,
      top = high,
      bottom = low
    )
  }
  if (is_sf) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("The sf package is required for sf input.", call. = FALSE)
    }
    geometry = sf::st_geometry(polygon)
    if (!length(geometry)) {
      stop("Empty polygons cannot be extruded.", call. = FALSE)
    }
    height_column = function(name, fallback, argument) {
      if (is.null(name)) {
        return(rep(fallback, length(geometry)))
      }
      if (
        !is.character(name) ||
          length(name) != 1L ||
          is.na(name) ||
          !nzchar(name)
      ) {
        stop(
          sprintf("`%s` must be NULL or a single column name.", argument),
          call. = FALSE
        )
      }
      if (!name %in% names(polygon)) {
        warning(
          sprintf("Was not able to find %s `%s` in sf object.", argument, name),
          call. = FALSE
        )
        return(rep(fallback, length(geometry)))
      }
      values = polygon[[name]]
      if (
        !is.numeric(values) ||
          !is.null(dim(values)) ||
          length(values) != length(geometry) ||
          any(!is.finite(values)) ||
          any(!is.finite(values * scale_data))
      ) {
        stop(
          sprintf("`%s` must select a finite numeric height column.", argument),
          call. = FALSE
        )
      }
      values * scale_data
    }
    high = height_column(data_column_top, top, "data_column_top")
    low = height_column(data_column_bottom, bottom, "data_column_bottom")
    for (i in seq_along(geometry)) {
      g = geometry[[i]]
      if (inherits(g, "POLYGON")) {
        pieces = list(unclass(g))
      } else if (inherits(g, "MULTIPOLYGON")) {
        pieces = unclass(g)
      } else {
        stop(
          "sf geometry must contain only POLYGON or MULTIPOLYGON features.",
          call. = FALSE
        )
      }
      if (!length(pieces)) {
        stop("Empty polygons cannot be extruded.", call. = FALSE)
      }
      for (rings in pieces) {
        append(rings, high[i], low[i])
      }
    }
  } else if (spatial) {
    if (!length(polygon@polygons)) {
      stop("Empty polygons cannot be extruded.", call. = FALSE)
    }
    for (feature in polygon@polygons) {
      rings = lapply(feature@Polygons, function(r) polygon_ring(r@coords))
      hole = vapply(feature@Polygons, function(r) r@hole, TRUE)
      outer = which(!hole)
      if (!length(outer)) {
        stop("Spatial polygons require an outer ring.", call. = FALSE)
      }
      polygon_check_edges(rings)
      groups = lapply(outer, function(i) list(rings[[i]]))
      for (i in which(hole)) {
        parents = which(vapply(
          outer,
          function(j) polygon_contains(rings[[j]], rings[[i]][1, ]),
          TRUE
        ))
        # The innermost containing exterior owns the hole, independent of ring order.
        if (!length(parents)) {
          stop(
            "Every hole must lie strictly inside its outer polygon.",
            call. = FALSE
          )
        }
        depth = vapply(
          parents,
          function(j) {
            sum(vapply(
              outer[parents],
              function(k) {
                polygon_contains(rings[[k]], rings[[outer[j]]][1, ])
              },
              TRUE
            ))
          },
          0L
        )
        owner = parents[which.max(depth)]
        groups[[owner]][[length(groups[[owner]]) + 1L]] = rings[[i]]
      }
      for (rings in groups) {
        append(rings, top, bottom)
      }
    }
  } else {
    xy = tryCatch(grDevices::xy.coords(polygon), error = function(e) NULL)
    if (is.null(xy) || !length(xy$x)) {
      stop(
        "`polygon` must specify nonempty finite x/y coordinates.",
        call. = FALSE
      )
    }
    p = cbind(xy$x, xy$y)
    if (
      is.null(holes) ||
        (is.numeric(holes) && length(holes) == 1L && isTRUE(holes == 0))
    ) {
      starts = 1L
    } else {
      if (
        !is.numeric(holes) ||
          !length(holes) ||
          any(!is.finite(holes)) ||
          any(holes != floor(holes)) ||
          any(holes < 4L | holes > nrow(p)) ||
          any(diff(holes) <= 0)
      ) {
        stop(
          "`holes` must be NULL, zero, or strictly increasing integer start indices from 4 to the vertex count.",
          call. = FALSE
        )
      }
      starts = c(1L, as.integer(holes))
    }
    ends = c(starts[-1L] - 1L, nrow(p))
    append(
      Map(function(a, b) p[seq.int(a, b), , drop = FALSE], starts, ends),
      top,
      bottom
    )
  }
  components
}

#' Construct a shared-vertex polygon extrusion mesh
#' @keywords internal
#' @noRd
polygon_mesh_data = function(
  polygon,
  plane,
  top,
  bottom,
  holes,
  center,
  flip_horizontal,
  flip_vertical,
  data_column_top,
  data_column_bottom,
  scale_data,
  scale
) {
  for (name in c("top", "bottom", "scale_data")) {
    sweep_scalar(get(name), name)
  }
  for (name in c("center", "flip_horizontal", "flip_vertical")) {
    flag = get(name)
    if (!is.logical(flag) || length(flag) != 1L || is.na(flag)) {
      stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
    }
  }
  if (
    !is.character(plane) ||
      length(plane) != 1L ||
      is.na(plane) ||
      !tolower(plane) %in% c("xz", "zx", "xy", "yx", "zy", "yz")
  ) {
    stop("`plane` must be one of xz, zx, xy, yx, zy, or yz.", call. = FALSE)
  }
  plane = tolower(plane)
  if (
    !is.numeric(scale) ||
      !length(scale) %in% c(1L, 3L) ||
      any(!is.finite(scale)) ||
      any(scale == 0)
  ) {
    stop(
      "`scale` must contain one or three finite, nonzero numbers.",
      call. = FALSE
    )
  }
  scale = rep(scale, length.out = 3L)
  components = polygon_components(
    polygon,
    holes,
    top,
    bottom,
    data_column_top,
    data_column_bottom,
    scale_data
  )
  # Retain the established coordinate convention, including the first-axis mirror.
  for (i in seq_along(components)) {
    rings = lapply(components[[i]]$rings, function(p) {
      p = polygon_ring(p)
      p[, 1] = p[, 1] * if (flip_horizontal) 1 else -1
      if (flip_vertical) {
        p[, 2] = -p[, 2]
      }
      p
    })
    components[[i]]$rings = polygon_component(rings)
  }
  offset = c(0, 0)
  if (center) {
    all_points = do.call(
      rbind,
      lapply(components, function(c) do.call(rbind, c$rings))
    )
    offset = apply(all_points, 2, function(x) min(x) / 2 + max(x) / 2)
  }
  permutation = switch(
    plane,
    xz = 1:3,
    zx = c(3, 2, 1),
    xy = c(3, 1, 2),
    yx = c(1, 3, 2),
    zy = c(2, 3, 1),
    yz = c(2, 1, 3)
  )
  reflected = xor(plane %in% c("zx", "yx", "yz"), prod(sign(scale)) < 0)
  meshes = lapply(components, function(component) {
    rings = component$rings
    p = do.call(rbind, rings)
    n = nrow(p)
    lengths = vapply(rings, nrow, 0L)
    starts = c(0L, utils::head(cumsum(lengths), -1L))
    local = sweep(p, 2, p[1, ])
    local = local / max(abs(local))
    caps = matrix(
      decido::earcut(
        local,
        holes = if (length(starts) > 1L) starts[-1L] + 1L else 0L
      ),
      ncol = 3L,
      byrow = TRUE
    )
    if (!nrow(caps)) {
      stop("Polygon could not be triangulated.", call. = FALSE)
    }
    a = local[caps[, 2], , drop = FALSE] - local[caps[, 1], , drop = FALSE]
    b = local[caps[, 3], , drop = FALSE] - local[caps[, 1], , drop = FALSE]
    area = a[, 1] * b[, 2] - a[, 2] * b[, 1]
    if (any(area == 0)) {
      stop("Polygon triangulation produced a degenerate cap.", call. = FALSE)
    }
    caps[area < 0, ] = caps[area < 0, c(1, 3, 2), drop = FALSE]
    next_id = unlist(Map(
      function(n, start) c(seq_len(n)[-1L], 1L) + start,
      lengths,
      starts
    ))
    # Earcut can discard collinear vertices. Require the cap boundary to match
    # the wall rings exactly, so no T-junction can silently enter the mesh.
    edges = rbind(caps[, c(1, 2)], caps[, c(2, 3)], caps[, c(3, 1)])
    edge_key = function(a, b) paste(pmin(a, b), pmax(a, b))
    counts = table(edge_key(edges[, 1], edges[, 2]))
    boundary = edge_key(seq_len(n), next_id)
    if (
      anyNA(counts[boundary]) ||
        any(counts[boundary] != 1L) ||
        any(counts[!names(counts) %in% boundary] != 2L)
    ) {
      stop(
        "Polygon triangulation did not preserve a closed cap boundary.",
        call. = FALSE
      )
    }
    p = sweep(p, 2, offset)
    low = min(component$top, component$bottom)
    high = max(component$top, component$bottom)
    vertices = cbind(p[, 1], low, p[, 2])
    indices = caps
    if (low != high) {
      vertices = rbind(vertices, cbind(p[, 1], high, p[, 2]))
      i = seq_len(n)
      indices = rbind(
        caps,
        caps[, c(1, 3, 2)] + n,
        cbind(i, i + n, next_id + n),
        cbind(i, next_id + n, next_id)
      )
    }
    vertices = sweep(vertices[, permutation, drop = FALSE], 2, scale, "*")
    if (any(!is.finite(vertices))) {
      stop("Scaled polygon vertices must be finite.", call. = FALSE)
    }
    if (anyDuplicated(vertices)) {
      stop(
        "Scaling or centering collapsed distinct polygon vertices.",
        call. = FALSE
      )
    }
    if (reflected) {
      indices = indices[, c(1, 3, 2), drop = FALSE]
    }
    list(vertices = vertices, indices = indices)
  })
  counts = vapply(meshes, function(m) nrow(m$vertices), 0L)
  starts = c(0L, utils::head(cumsum(counts), -1L))
  list(
    vertices = do.call(rbind, lapply(meshes, `[[`, "vertices")),
    indices = do.call(
      rbind,
      Map(function(m, start) m$indices + start - 1L, meshes, starts)
    )
  )
}
