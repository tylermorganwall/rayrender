#!/usr/bin/env Rscript

script_version = "pr5-rgb-spectrum-table-v1"
table_resolution = 64L
coefficient_count = 3L *
  table_resolution *
  table_resolution *
  table_resolution *
  3L
pbrt_commit = "8c19f304558fd7681e2fef2c395a689d0106fb05"

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

read_pbrt_rgb_table_cpp = function(path) {
  if (!file.exists(path)) {
    stop("Unable to find pbrt RGB spectrum table at ", path, call. = FALSE)
  }
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

extract_initializer = function(source, name) {
  start_pattern = paste0(name, "[[:space:]]*(?:\\[[^]]*\\])*[[:space:]]*=")
  start = regexpr(start_pattern, source, perl = TRUE)
  if (start < 0) {
    stop("Unable to find pbrt initializer: ", name, call. = FALSE)
  }
  tail = substring(source, start)
  open = regexpr("\\{", tail, perl = TRUE)
  close = regexpr("\\};", tail, perl = TRUE)
  if (open < 0 || close < 0 || close <= open) {
    stop("Malformed pbrt initializer: ", name, call. = FALSE)
  }
  body = substr(tail, open + 1, close - 1)
  gsub("//[^\n]*", "", body, perl = TRUE)
}

parse_numeric_initializer = function(source, name, expected_count) {
  body = extract_initializer(source, name)
  tokens = regmatches(
    body,
    gregexpr("[-+]?(?:[0-9]*[.])?[0-9]+(?:[eE][-+]?[0-9]+)?", body, perl = TRUE)
  )[[1]]
  values = as.numeric(tokens)
  if (length(values) != expected_count) {
    stop(
      "Unexpected value count for ",
      name,
      ": expected ",
      expected_count,
      ", got ",
      length(values),
      call. = FALSE
    )
  }
  values
}

write_uint32 = function(con, value) {
  if (!is.finite(value) || value < 0 || value > 4294967295) {
    stop("uint32 value out of range: ", value, call. = FALSE)
  }
  signed = if (value > 2147483647) value - 4294967296 else value
  writeBin(as.integer(signed), con, size = 4, endian = "little")
}

fixed_ascii = function(value, width) {
  bytes = charToRaw(value)
  if (length(bytes) > width) {
    stop("Fixed ASCII field is too long: ", value, call. = FALSE)
  }
  c(bytes, raw(width - length(bytes)))
}

adler32 = function(bytes) {
  modulus = 65521
  a = 1
  b = 0
  n = length(bytes)
  if (n == 0) {
    return(1)
  }
  chunk_size = 50000
  starts = seq(1, n, by = chunk_size)
  for (start in starts) {
    end = min(start + chunk_size - 1, n)
    values = as.integer(bytes[start:end])
    cumulative = cumsum(values)
    a_values = (a + cumulative) %% modulus
    b = (b + sum(a_values)) %% modulus
    a = (a + sum(values)) %% modulus
  }
  b * 65536 + a
}

payload_bytes = function(scale, coefficients) {
  con = rawConnection(raw(0), "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.numeric(scale), con, size = 4, endian = "little")
  writeBin(as.numeric(coefficients), con, size = 4, endian = "little")
  rawConnectionValue(con)
}

write_table = function(output, scale, coefficients) {
  payload = payload_bytes(scale, coefficients)
  checksum = adler32(payload)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)

  con = file(output, "wb")
  on.exit(close(con), add = TRUE)
  writeBin(fixed_ascii("RAYRGBSPECv1", 16), con, useBytes = TRUE)
  write_uint32(con, 1)
  write_uint32(con, 0x01020304)
  write_uint32(con, 4)
  write_uint32(con, table_resolution)
  write_uint32(con, 3)
  write_uint32(con, 3)
  write_uint32(con, table_resolution)
  write_uint32(con, coefficient_count)
  write_uint32(con, length(payload))
  write_uint32(con, checksum)
  writeBin(fixed_ascii("sRGB", 32), con, useBytes = TRUE)
  writeBin(payload, con, useBytes = TRUE)

  checksum
}

update_manifest = function(root, output) {
  manifest = file.path(root, "docs", "spectral", "assets-manifest.csv")
  assets = read.csv(
    manifest,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA"),
    colClasses = "character"
  )
  root_path = normalizePath(root)
  output_path = normalizePath(output)
  rel_path = substring(output_path, nchar(root_path) + 2)
  row = data.frame(
    path = rel_path,
    kind = "rgb_to_spectrum_table",
    semantic_type = "sRGB_reconstruction",
    source = paste0("pbrt-v4 ", pbrt_commit, " build/rgbspectrum_srgb.cpp"),
    license = "Apache-2.0",
    sha256 = unname(tools::sha256sum(output)),
    package_required = "true",
    notes = paste0(
      "Generated by tools/spectral-tests/generate-pr5-rgb-table.R (",
      script_version,
      ")"
    ),
    stringsAsFactors = FALSE
  )

  assets = assets[assets$path != rel_path, , drop = FALSE]
  assets = rbind(assets, row[names(assets)])
  write.csv(assets, manifest, row.names = FALSE, na = "")
}

args = parse_args(commandArgs(trailingOnly = TRUE))
root = repo_root()
source = args$source %||%
  file.path(root, "mmp", "pbrt-v4", "build", "rgbspectrum_srgb.cpp")
output = args$output %||%
  file.path(root, "inst", "extdata", "spectral", "rgb-to-spectrum-srgb-v1.bin")

source_text = read_pbrt_rgb_table_cpp(source)
scale = parse_numeric_initializer(
  source_text,
  "sRGBToSpectrumTable_Scale",
  table_resolution
)
coefficients = parse_numeric_initializer(
  source_text,
  "sRGBToSpectrumTable_Data",
  coefficient_count
)

if (!identical(scale[[1]], 0) || !identical(scale[[length(scale)]], 1)) {
  stop(
    "Unexpected pbrt scale endpoints for sRGB RGB-to-spectrum table",
    call. = FALSE
  )
}

checksum = write_table(output, scale, coefficients)
update_manifest(root, output)

message("Wrote ", output)
message("Payload Adler-32: ", sprintf("0x%08X", checksum))
