parse_args = function(args) {
  values = list()
  extra = list()
  i = 1
  while (i <= length(args)) {
    key = args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected positional argument: ", key)
    }
    name = sub("^--", "", key)
    name = gsub("-", "_", name)
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      value = "true"
      i = i + 1
    } else {
      value = args[[i + 1]]
      i = i + 2
    }
    if (
      name %in%
        c(
          "benchmark_name",
          "benchmark_source",
          "lib_path",
          "iteration",
          "seed",
          "width",
          "height",
          "samples",
          "threads",
          "settings_json",
          "output_json",
          "artifact_dir"
        )
    ) {
      values[[name]] = value
    } else {
      extra[[name]] = parse_scalar(value)
    }
  }
  values$extra = extra
  values
}

parse_scalar = function(value) {
  value_trim = trimws(value)
  value_lower = tolower(value_trim)
  if (value_lower %in% c("true", "false")) {
    return(value_lower == "true")
  }
  if (value_lower %in% c("null", "na")) {
    return(NA)
  }
  numeric_value = suppressWarnings(as.numeric(value_trim))
  if (
    !is.na(numeric_value) && grepl("^-?[0-9.]+([eE][-+]?[0-9]+)?$", value_trim)
  ) {
    return(numeric_value)
  }
  value_trim
}

parse_simple_json = function(text) {
  text = trimws(text)
  if (!nzchar(text) || text == "{}") {
    return(list())
  }
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    parsed = jsonlite::fromJSON(text, simplifyVector = FALSE)
    if (is.null(parsed)) {
      return(list())
    }
    return(parsed)
  }
  if (!startsWith(text, "{") || !endsWith(text, "}")) {
    stop("settings-json requires jsonlite for non-object JSON")
  }
  inner = substring(text, 2, nchar(text) - 1)
  if (!nzchar(trimws(inner))) {
    return(list())
  }
  parts = strsplit(inner, ",", fixed = TRUE)[[1]]
  output = list()
  for (part in parts) {
    key_value = strsplit(part, ":", fixed = TRUE)[[1]]
    if (length(key_value) < 2) {
      next
    }
    key = trimws(key_value[[1]])
    key = gsub('^"|"$', "", key)
    value = trimws(paste(key_value[-1], collapse = ":"))
    value = gsub('^"|"$', "", value)
    output[[key]] = parse_scalar(value)
  }
  output
}

json_escape = function(value) {
  value = as.character(value)
  value = gsub("\\\\", "\\\\\\\\", value)
  value = gsub("\"", "\\\\\"", value)
  value = gsub("\n", "\\\\n", value)
  value = gsub("\r", "\\\\r", value)
  value
}

json_value = function(value) {
  if (length(value) == 0 || is.null(value) || all(is.na(value))) {
    return("null")
  }
  if (is.logical(value)) {
    return(ifelse(isTRUE(value[[1]]), "true", "false"))
  }
  if (is.numeric(value)) {
    if (!is.finite(value[[1]])) {
      return("null")
    }
    return(format(value[[1]], scientific = FALSE, digits = 15))
  }
  paste0("\"", json_escape(value[[1]]), "\"")
}

write_json_object = function(path, values) {
  parent = dirname(path)
  if (!dir.exists(parent)) {
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  }
  fields = names(values)
  pieces = vapply(
    fields,
    function(name) {
      paste0("\"", json_escape(name), "\":", json_value(values[[name]]))
    },
    character(1)
  )
  writeLines(paste0("{", paste(pieces, collapse = ","), "}"), path)
}

hash_object = function(value, artifact_dir, write_artifact) {
  if (!dir.exists(artifact_dir)) {
    dir.create(artifact_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (is.array(value) || is.matrix(value) || is.atomic(value)) {
    payload = list(
      type = typeof(value),
      dim = dim(value),
      values = as.vector(value)
    )
  } else {
    payload = utils::capture.output(str(value))
  }
  temp_path = tempfile(
    pattern = "render-",
    fileext = ".rds",
    tmpdir = artifact_dir
  )
  saveRDS(payload, temp_path, version = 2)
  hash = unname(tools::md5sum(temp_path))
  artifact_path = NA_character_
  if (isTRUE(write_artifact)) {
    artifact_path = temp_path
  } else {
    unlink(temp_path)
  }
  list(hash = hash, artifact_path = artifact_path)
}

validate_render_output = function(value) {
  if (is.null(value)) {
    stop("Benchmark returned NULL")
  }
  if (length(value) == 0) {
    stop("Benchmark returned an empty object")
  }
  if (is.numeric(value) && all(is.na(value))) {
    stop("Benchmark returned only NA numeric values")
  }
  TRUE
}

metric_value = function(value, name, default = NA_real_) {
  if (is.list(value) && name %in% names(value)) {
    return(value[[name]])
  }
  if (is.data.frame(value) && name %in% names(value) && nrow(value) > 0) {
    return(value[[name]][[1]])
  }
  attr_value = attr(value, name, exact = TRUE)
  if (!is.null(attr_value)) {
    return(attr_value)
  }
  default
}

numeric_metric = function(value, name, default = NA_real_) {
  metric = suppressWarnings(as.numeric(metric_value(value, name, default)))
  if (length(metric) == 0 || is.na(metric[[1]])) {
    return(default)
  }
  metric[[1]]
}

benchmark_output_object = function(value) {
  if (is.list(value)) {
    for (name in c("image", "rendered", "output", "artifact")) {
      if (name %in% names(value) && !is.null(value[[name]])) {
        return(value[[name]])
      }
    }
  }
  value
}

args = parse_args(commandArgs(trailingOnly = TRUE))

required = c(
  "benchmark_name",
  "benchmark_source",
  "lib_path",
  "iteration",
  "seed",
  "width",
  "height",
  "samples",
  "threads",
  "settings_json",
  "output_json",
  "artifact_dir"
)
missing = required[!required %in% names(args)]
if (length(missing) > 0) {
  stop("Missing required arguments: ", paste(missing, collapse = ", "))
}

result = tryCatch(
  {
    lib_path = normalizePath(args$lib_path, mustWork = TRUE)
    .libPaths(c(lib_path, .libPaths()))
    library(rayrender)
    requireNamespace("rayimage", quietly = TRUE)

    source(args$benchmark_source, local = FALSE)
    if (!exists("run_benchmark", mode = "function")) {
      stop("Benchmark source must define run_benchmark(settings)")
    }

    settings = parse_simple_json(args$settings_json)
    settings$benchmark_name = args$benchmark_name
    settings$benchmark_source = args$benchmark_source
    settings$width = as.integer(args$width)
    settings$height = as.integer(args$height)
    settings$samples = as.integer(args$samples)
    settings$threads = as.integer(args$threads)
    settings$seed = as.integer(args$seed)
    settings$iteration = as.integer(args$iteration)
    settings$artifact_dir = args$artifact_dir
    for (name in names(args$extra)) {
      settings[[name]] = args$extra[[name]]
    }

    options(cores = settings$threads, Ncpus = settings$threads)
    set.seed(settings$seed + settings$iteration)

    elapsed = system.time({
      benchmark_result = run_benchmark(settings)
    })[["elapsed"]]
    rendered = benchmark_output_object(benchmark_result)
    render_seconds = numeric_metric(
      benchmark_result,
      "render_seconds",
      NA_real_
    )
    if (is.null(render_seconds) || is.na(render_seconds)) {
      render_seconds = elapsed
    }
    scene_build_seconds = numeric_metric(
      benchmark_result,
      "scene_build_seconds",
      NA_real_
    )
    bvh_build_seconds = numeric_metric(
      benchmark_result,
      "bvh_build_seconds",
      NA_real_
    )
    total_seconds = numeric_metric(benchmark_result, "total_seconds", elapsed)

    validate_render_output(rendered)
    hash_info = hash_object(
      rendered,
      settings$artifact_dir,
      isTRUE(settings$write_artifact)
    )

    list(
      render_seconds = as.numeric(render_seconds),
      scene_build_seconds = scene_build_seconds,
      bvh_build_seconds = bvh_build_seconds,
      total_seconds = total_seconds,
      status = "ok",
      error = NA_character_,
      output_hash = hash_info$hash,
      image_path = hash_info$artifact_path
    )
  },
  error = function(e) {
    list(
      render_seconds = NA_real_,
      scene_build_seconds = NA_real_,
      bvh_build_seconds = NA_real_,
      total_seconds = NA_real_,
      status = "render_failed",
      error = conditionMessage(e),
      output_hash = NA_character_,
      image_path = NA_character_
    )
  }
)

write_json_object(args$output_json, result)
if (!identical(result$status, "ok")) {
  quit(status = 1)
}
