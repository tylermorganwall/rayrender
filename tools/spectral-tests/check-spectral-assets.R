#!/usr/bin/env Rscript

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
      stop("Unexpected positional argument: ", key)
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

is_true = function(value) {
  tolower(as.character(value)) %in% c("true", "t", "1", "yes", "y")
}

sha256_file = function(path) {
  unname(tools::sha256sum(path))
}

build_source_package = function(root) {
  oldwd = getwd()
  tmp = tempfile("spectral-assets-build-")
  dir.create(tmp, recursive = TRUE)
  on.exit(
    {
      setwd(oldwd)
      unlink(tmp, recursive = TRUE)
    },
    add = TRUE
  )

  setwd(tmp)
  output = system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "build", "--no-build-vignettes", "--no-manual", root),
    stdout = TRUE,
    stderr = TRUE
  )
  tarball = list.files(tmp, pattern = "[.]tar[.]gz$", full.names = TRUE)
  if (length(tarball) != 1) {
    stop(
      "R CMD build did not report a source tarball:\n",
      paste(output, collapse = "\n")
    )
  }
  persistent_tarball = tempfile("rayrender-source-", fileext = ".tar.gz")
  if (!file.copy(tarball[[1]], persistent_tarball, overwrite = TRUE)) {
    stop(
      "Unable to preserve source tarball for asset inspection",
      call. = FALSE
    )
  }
  normalizePath(persistent_tarball, mustWork = TRUE)
}

tar_entries = function(tarball) {
  entries = utils::untar(tarball, list = TRUE)
  sub("^[^/]+/", "", entries)
}

args = parse_args(commandArgs(trailingOnly = TRUE))
root = repo_root()
manifest = args$manifest %||%
  file.path(root, "docs", "spectral", "assets-manifest.csv")
manifest = normalizePath(manifest, mustWork = TRUE)

assets = read.csv(
  manifest,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA"),
  colClasses = "character"
)

required_columns = c(
  "path",
  "kind",
  "semantic_type",
  "source",
  "license",
  "sha256",
  "package_required",
  "notes"
)
missing_columns = setdiff(required_columns, names(assets))
if (length(missing_columns) > 0) {
  stop(
    "Asset manifest is missing columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (nrow(assets) == 0) {
  message("No spectral assets listed; manifest schema check passed.")
  quit(status = 0)
}

errors = character()
for (i in seq_len(nrow(assets))) {
  rel_path = assets$path[[i]]
  abs_path = file.path(root, rel_path)
  if (!file.exists(abs_path)) {
    errors = c(errors, sprintf("missing file: %s", rel_path))
    next
  }
  expected_hash = assets$sha256[[i]]
  if (!is.na(expected_hash) && nzchar(expected_hash)) {
    actual_hash = sha256_file(abs_path)
    if (!identical(tolower(actual_hash), tolower(expected_hash))) {
      errors = c(
        errors,
        sprintf(
          "checksum mismatch for %s: expected %s, got %s",
          rel_path,
          expected_hash,
          actual_hash
        )
      )
    }
  }
}

build_source = isTRUE(args$build_source)
if (build_source) {
  tarball = build_source_package(root)
  entries = tar_entries(tarball)
  package_assets = assets[is_true(assets$package_required), , drop = FALSE]
  for (i in seq_len(nrow(package_assets))) {
    rel_path = package_assets$path[[i]]
    if (!rel_path %in% entries) {
      errors = c(
        errors,
        sprintf("source package does not contain %s", rel_path)
      )
    }
  }
}

if (length(errors) > 0) {
  stop(paste(errors, collapse = "\n"))
}

message(sprintf("Checked %d spectral asset manifest row(s).", nrow(assets)))
