# rwasm runs configure in host R. Its LinkingTo packages provide native
# libraries, even though em++ is compiling the consumer for wasm32.
wasm_exr_paths = function(package, archive) {
  host = find.package(package)
  version = read.dcf(file.path(host, "DESCRIPTION"), fields = "Version")[[1]]
  r_series = paste(
    R.version$major,
    sub("\\..*", "", R.version$minor),
    sep = "."
  )
  target_library = Sys.getenv("RAYRENDER_WASM_RLIB")
  cache = file.path(".cache", "wasm", r_series, paste0(package, "_", version))

  is_target = function(root) {
    desc = file.path(root, "DESCRIPTION")
    if (!file.exists(desc)) {
      return(FALSE)
    }
    pkg = read.dcf(desc, fields = c("Package", "Version", "Built"))
    identical(pkg[[1]], package) &&
      identical(pkg[[2]], version) &&
      isTRUE(grepl("; wasm32-unknown-emscripten;", pkg[[3]], fixed = TRUE))
  }

  if (nzchar(target_library)) {
    root = file.path(target_library, package)
  } else if (is_target(host)) {
    root = host
  } else {
    root = file.path(cache, package)
    if (!is_target(root)) {
      url = sprintf(
        "https://tylermorganwall.r-universe.dev/bin/emscripten/contrib/%s/%s_%s.tgz",
        r_series,
        package,
        version
      )
      tgz = tempfile(fileext = ".tgz")
      on.exit(unlink(tgz))
      dir.create(cache, recursive = TRUE, showWarnings = FALSE)
      stage = tempfile(tmpdir = cache)
      dir.create(stage)
      on.exit(unlink(stage, recursive = TRUE), add = TRUE)
      if (utils::download.file(url, tgz, mode = "wb") != 0L) {
        stop("Could not download ", url)
      }
      # Unpack only the target headers, archives and package identity.
      entries = utils::untar(tgz, list = TRUE)
      entries = entries[!endsWith(entries, "/")]
      files = entries[
        startsWith(entries, paste0(package, "/include/")) |
          endsWith(entries, ".a") &
            startsWith(entries, paste0(package, "/lib/"))
      ]
      status = utils::untar(
        tgz,
        files = c(files, paste0(package, "/DESCRIPTION")),
        exdir = stage
      )
      unpacked = file.path(stage, package)
      if (isTRUE(status != 0L) || !is_target(unpacked)) {
        stop("Could not unpack the wasm32 build of ", package, " ", version)
      }
      if (!file.rename(unpacked, root) && !is_target(root)) {
        stop("Could not cache the wasm32 build of ", package, " at ", root)
      }
    }
  }
  if (!is_target(root)) {
    stop("Expected a wasm32 build of ", package, " ", version, " at ", root)
  }
  # Published packages can name this subdirectory after the build machine.
  libs = Sys.glob(file.path(root, "lib", "*", archive))
  if (length(libs) != 1L) {
    stop("Expected one ", archive, " under ", root)
  }
  list(
    include = normalizePath(file.path(root, "include"), mustWork = TRUE),
    lib = dirname(normalizePath(libs))
  )
}
