#!/usr/bin/env Rscript

script_version = "pr4-spectral-assets-v1"
cie_y_integral = 106.856895
lambda_min = 360
lambda_max = 830

`%||%` = function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) {
    y
  } else {
    x
  }
}

parse_args = function(args) {
  values = list()
  i = 1
  while (i <= length(args)) {
    key = args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected positional argument: ", key, call. = FALSE)
    }
    name = gsub("-", "_", sub("^--", "", key))
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      values[[name]] = TRUE
      i = i + 1
    } else {
      values[[name]] = args[[i + 1]]
      i = i + 2
    }
  }
  values
}

script_path = function() {
  file_arg = grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
  }
  normalizePath(sys.frames()[[1]]$ofile, mustWork = TRUE)
}

repo_root = function() {
  normalizePath(file.path(dirname(script_path()), "..", ".."), mustWork = TRUE)
}

read_pbrt_spectrum_cpp = function(root) {
  path = file.path(root, "src", "pbrt", "util", "spectrum.cpp")
  if (!file.exists(path)) {
    stop("Unable to find pbrt spectrum.cpp at ", path, call. = FALSE)
  }
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

extract_array = function(source, name) {
  start_pattern = paste0(
    "const[[:space:]]+Float[[:space:]]+",
    name,
    "[[:space:]]*(?:\\[[^]]*\\])?[[:space:]]*="
  )
  start = regexpr(start_pattern, source, perl = TRUE)
  if (start < 0) {
    stop("Unable to find pbrt array: ", name, call. = FALSE)
  }
  tail = substring(source, start)
  open = regexpr("\\{", tail, perl = TRUE)
  close = regexpr("\\};", tail, perl = TRUE)
  if (open < 0 || close < 0 || close <= open) {
    stop("Malformed pbrt array: ", name, call. = FALSE)
  }
  body = substr(tail, open + 1, close - 1)
  body = gsub("//[^\n]*", "", body, perl = TRUE)
  tokens = regmatches(
    body,
    gregexpr("[-+]?(?:[0-9]*[.])?[0-9]+(?:[eE][-+]?[0-9]+)?", body, perl = TRUE)
  )[[1]]
  as.numeric(tokens)
}

as_pairs = function(samples) {
  if (length(samples) %% 2 != 0) {
    stop(
      "Interleaved samples must contain wavelength/value pairs",
      call. = FALSE
    )
  }
  data.frame(
    lambda_nm = samples[seq(1, length(samples), by = 2)],
    value = samples[seq(2, length(samples), by = 2)]
  )
}

extend_visible_interval = function(pairs) {
  if (pairs$lambda_nm[[1]] > lambda_min) {
    pairs = rbind(
      data.frame(lambda_nm = lambda_min - 1, value = pairs$value[[1]]),
      pairs
    )
  }
  if (pairs$lambda_nm[[nrow(pairs)]] < lambda_max) {
    pairs = rbind(
      pairs,
      data.frame(lambda_nm = lambda_max + 1, value = pairs$value[[nrow(pairs)]])
    )
  }
  row.names(pairs) = NULL
  pairs
}

eval_piecewise = function(pairs, lambda_nm) {
  stats::approx(
    x = pairs$lambda_nm,
    y = pairs$value,
    xout = lambda_nm,
    method = "linear",
    rule = 1,
    ties = "ordered"
  )$y
}

inner_product = function(lhs, rhs) {
  lambdas = seq(lambda_min, lambda_max)
  sum(eval_piecewise(lhs, lambdas) * eval_piecewise(rhs, lambdas))
}

normalize_to_y = function(pairs, cie_y) {
  y_integral = inner_product(pairs, cie_y)
  if (!is.finite(y_integral) || y_integral <= 0) {
    stop(
      "Cannot normalize spectrum with non-positive Y integral",
      call. = FALSE
    )
  }
  pairs$value = pairs$value * (cie_y_integral / y_integral)
  pairs
}

fmt = function(value) {
  formatC(value, digits = 15, format = "fg", flag = "#")
}

sha256_lines = function(lines) {
  temp = tempfile("spectral-asset-sha-")
  on.exit(unlink(temp), add = TRUE)
  writeLines(lines, temp, useBytes = TRUE)
  unname(tools::sha256sum(temp))
}

asset = function(
  name,
  semantic_type,
  pairs,
  value_unit,
  normalization,
  source_citation,
  source_license,
  storage
) {
  list(
    name = name,
    semantic_type = semantic_type,
    pairs = pairs,
    value_unit = value_unit,
    normalization = normalization,
    source_citation = source_citation,
    source_license = source_license,
    storage = storage
  )
}

source_name = "pbrt-v4 8c19f304558fd7681e2fef2c395a689d0106fb05 src/pbrt/util/spectrum.cpp"
pbrt_license = "Apache-2.0"

illuminant_specs = list(
  c("stdillum-A", "CIE_Illum_A"),
  c("stdillum-D50", "CIE_Illum_D5000"),
  c("illum-acesD60", "ACES_Illum_D60"),
  c("stdillum-D65", "CIE_Illum_D6500"),
  c("stdillum-F1", "CIE_Illum_F1"),
  c("stdillum-F2", "CIE_Illum_F2"),
  c("stdillum-F3", "CIE_Illum_F3"),
  c("stdillum-F4", "CIE_Illum_F4"),
  c("stdillum-F5", "CIE_Illum_F5"),
  c("stdillum-F6", "CIE_Illum_F6"),
  c("stdillum-F7", "CIE_Illum_F7"),
  c("stdillum-F8", "CIE_Illum_F8"),
  c("stdillum-F9", "CIE_Illum_F9"),
  c("stdillum-F10", "CIE_Illum_F10"),
  c("stdillum-F11", "CIE_Illum_F11"),
  c("stdillum-F12", "CIE_Illum_F12")
)

metal_specs = list(
  c("metal-Ag-eta", "Ag_eta", "optical_eta"),
  c("metal-Ag-k", "Ag_k", "optical_k"),
  c("metal-Al-eta", "Al_eta", "optical_eta"),
  c("metal-Al-k", "Al_k", "optical_k"),
  c("metal-Au-eta", "Au_eta", "optical_eta"),
  c("metal-Au-k", "Au_k", "optical_k"),
  c("metal-Cu-eta", "Cu_eta", "optical_eta"),
  c("metal-Cu-k", "Cu_k", "optical_k"),
  c("metal-CuZn-eta", "CuZn_eta", "optical_eta"),
  c("metal-CuZn-k", "CuZn_k", "optical_k"),
  c("metal-MgO-eta", "MgO_eta", "optical_eta"),
  c("metal-MgO-k", "MgO_k", "optical_k"),
  c("metal-TiO2-eta", "TiO2_eta", "optical_eta"),
  c("metal-TiO2-k", "TiO2_k", "optical_k")
)

glass_specs = list(
  c("glass-BK7", "GlassBK7_eta"),
  c("glass-BAF10", "GlassBAF10_eta"),
  c("glass-FK51A", "GlassFK51A_eta"),
  c("glass-LASF9", "GlassLASF9_eta"),
  c("glass-F5", "GlassSF5_eta"),
  c("glass-F10", "GlassSF10_eta"),
  c("glass-F11", "GlassSF11_eta")
)

sensor_specs = c(
  "canon_eos_100d",
  "canon_eos_1dx_mkii",
  "canon_eos_200d",
  "canon_eos_200d_mkii",
  "canon_eos_5d",
  "canon_eos_5d_mkii",
  "canon_eos_5d_mkiii",
  "canon_eos_5d_mkiv",
  "canon_eos_5ds",
  "canon_eos_m",
  "hasselblad_l1d_20c",
  "nikon_d810",
  "nikon_d850",
  "sony_ilce_6400",
  "sony_ilce_7m3",
  "sony_ilce_7rm3",
  "sony_ilce_9"
)

build_assets = function(source) {
  cie_lambda = extract_array(source, "CIE_lambda")
  cie_x = data.frame(
    lambda_nm = cie_lambda,
    value = extract_array(source, "CIE_X")
  )
  cie_y = data.frame(
    lambda_nm = cie_lambda,
    value = extract_array(source, "CIE_Y")
  )
  cie_z = data.frame(
    lambda_nm = cie_lambda,
    value = extract_array(source, "CIE_Z")
  )

  assets = list(
    asset(
      "cie-x",
      "color_matching_function",
      cie_x,
      "relative_response",
      "CIE 1931 2-degree observer, 1 nm samples",
      source_name,
      pbrt_license,
      "dense_1nm"
    ),
    asset(
      "cie-y",
      "color_matching_function",
      cie_y,
      "relative_response",
      paste0("CIE 1931 2-degree observer, integral ", cie_y_integral),
      source_name,
      pbrt_license,
      "dense_1nm"
    ),
    asset(
      "cie-z",
      "color_matching_function",
      cie_z,
      "relative_response",
      "CIE 1931 2-degree observer, 1 nm samples",
      source_name,
      pbrt_license,
      "dense_1nm"
    )
  )

  for (spec in illuminant_specs) {
    pairs = extend_visible_interval(as_pairs(extract_array(source, spec[[2]])))
    pairs = normalize_to_y(pairs, cie_y)
    assets[[length(assets) + 1]] = asset(
      spec[[1]],
      "illuminant",
      pairs,
      "relative_spd",
      paste0("pbrt luminance normalization to CIE_Y_integral=", cie_y_integral),
      source_name,
      pbrt_license,
      "piecewise_linear"
    )
  }

  for (spec in metal_specs) {
    pairs = extend_visible_interval(as_pairs(extract_array(source, spec[[2]])))
    assets[[length(assets) + 1]] = asset(
      spec[[1]],
      spec[[3]],
      pairs,
      "dimensionless",
      "none",
      source_name,
      pbrt_license,
      "piecewise_linear"
    )
  }

  for (spec in glass_specs) {
    pairs = extend_visible_interval(as_pairs(extract_array(source, spec[[2]])))
    assets[[length(assets) + 1]] = asset(
      spec[[1]],
      "glass_eta",
      pairs,
      "dimensionless",
      "none",
      "refractiveindex.info data embedded in pbrt-v4 spectrum.cpp",
      "CC0-1.0",
      "piecewise_linear"
    )
  }

  for (sensor in sensor_specs) {
    for (channel in c("r", "g", "b")) {
      array_name = paste(sensor, channel, sep = "_")
      pairs = extend_visible_interval(as_pairs(extract_array(
        source,
        array_name
      )))
      assets[[length(assets) + 1]] = asset(
        paste(sensor, channel, sep = "_"),
        "camera_sensor_response",
        pairs,
        "relative_response",
        "none",
        source_name,
        pbrt_license,
        "piecewise_linear"
      )
    }
  }

  assets
}

write_assets = function(assets, output_dir, manifest_path) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  data_lines = "name\tlambda_nm\tvalue"
  metadata_rows = data.frame(
    name = character(),
    semantic_type = character(),
    wavelength_unit = character(),
    value_unit = character(),
    normalization = character(),
    source_citation = character(),
    source_license = character(),
    asset_generation_script_version = character(),
    sha256 = character(),
    storage = character(),
    stringsAsFactors = FALSE
  )

  for (item in assets) {
    item_lines = sprintf(
      "%s\t%s\t%s",
      item$name,
      fmt(item$pairs$lambda_nm),
      fmt(item$pairs$value)
    )
    data_lines = c(data_lines, item_lines)
    metadata_rows = rbind(
      metadata_rows,
      data.frame(
        name = item$name,
        semantic_type = item$semantic_type,
        wavelength_unit = "nm",
        value_unit = item$value_unit,
        normalization = item$normalization,
        source_citation = item$source_citation,
        source_license = item$source_license,
        asset_generation_script_version = script_version,
        sha256 = sha256_lines(item_lines),
        storage = item$storage,
        stringsAsFactors = FALSE
      )
    )
  }

  data_path = file.path(output_dir, "named-spectra-v1.tsv")
  metadata_path = file.path(output_dir, "named-spectra-v1-metadata.tsv")

  writeLines(data_lines, data_path, useBytes = TRUE)
  utils::write.table(
    metadata_rows,
    metadata_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )

  manifest = data.frame(
    path = c(
      sub(
        paste0(normalizePath(repo_root()), "/"),
        "",
        normalizePath(data_path)
      ),
      sub(
        paste0(normalizePath(repo_root()), "/"),
        "",
        normalizePath(metadata_path)
      )
    ),
    kind = c("spectral_named_data", "spectral_named_metadata"),
    semantic_type = c("mixed", "metadata"),
    source = c(source_name, source_name),
    license = c("Apache-2.0 and CC0-1.0", "Apache-2.0 and CC0-1.0"),
    sha256 = c(
      unname(tools::sha256sum(data_path)),
      unname(tools::sha256sum(metadata_path))
    ),
    package_required = c("true", "true"),
    notes = c(
      "Generated by tools/spectral-tests/generate-pr4-spectral-assets.R",
      "Generated by tools/spectral-tests/generate-pr4-spectral-assets.R"
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(manifest, manifest_path, row.names = FALSE, quote = TRUE)
}

args = parse_args(commandArgs(trailingOnly = TRUE))
root = repo_root()
pbrt_root = normalizePath(
  args$pbrt_root %||% file.path(root, "mmp", "pbrt-v4"),
  mustWork = TRUE
)
output_dir = normalizePath(
  args$output_dir %||% file.path(root, "inst", "extdata", "spectral"),
  mustWork = FALSE
)
manifest_path = args$manifest %||%
  file.path(root, "docs", "spectral", "assets-manifest.csv")

source = read_pbrt_spectrum_cpp(pbrt_root)
assets = build_assets(source)
write_assets(assets, output_dir, manifest_path)
message("Generated ", length(assets), " spectral asset(s) in ", output_dir)
