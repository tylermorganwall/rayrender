#' Preserve transport diagnostics, optionally writing a debug log
#' @keywords internal
save_path_warnings = function(diagnostics, filename) {
  if (is.null(diagnostics) || diagnostics$terminated_paths == 0) {
    return(NULL)
  }
  report = tolower(Sys.getenv("RAYRENDER_DEBUG_PATHS", "false")) %in%
    c("true", "1")
  if (!report) {
    return(list(diagnostics = diagnostics, log = NULL, report = FALSE))
  }
  path = if (is.na(filename)) {
    tempfile("rayrender-", fileext = ".warnings.log")
  } else {
    paste0(filename, ".warnings.log")
  }
  contents = c(
    "rayrender: terminated individual paths after transport failures",
    paste("UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("Package:", utils::packageVersion("rayrender")),
    paste("Terminated paths:", diagnostics$terminated_paths),
    "Accumulated light was retained; missing contributions can bias the image.",
    "Counts include all failures. Examples are bounded per failure category.",
    capture.output(dput(diagnostics))
  )
  write_log = function(path) {
    tryCatch(
      {
        suppressWarnings(writeLines(contents, path))
        normalizePath(path, winslash = "/", mustWork = TRUE)
      },
      error = function(e) NULL
    )
  }
  saved = write_log(path)
  if (is.null(saved)) {
    saved = write_log(tempfile("rayrender-", fileext = ".warnings.log"))
  }
  list(diagnostics = diagnostics, log = saved, report = TRUE)
}

#' Report a completed render's terminated paths on the R thread
#' @keywords internal
warn_path_failures = function(record) {
  if (is.null(record) || !isTRUE(record$report)) {
    return(invisible(NULL))
  }
  location = if (is.null(record$log)) {
    "Could not write a warning log; diagnostics are attached to the returned image."
  } else {
    paste0("Diagnostics: ", record$log)
  }
  warning(structure(
    list(
      message = paste0(
        "Terminated ",
        record$diagnostics$terminated_paths,
        " path(s) after transport failures; rendering continued with accumulated light. ",
        location
      ),
      call = NULL,
      diagnostics = record$diagnostics,
      log = record$log
    ),
    class = c("rayrender_path_warning", "warning", "condition")
  ))
  invisible(NULL)
}

#' Attach transport diagnostics to a returned image
#' @keywords internal
attach_path_warnings = function(image, record) {
  if (!is.null(record)) {
    attr(image, "path_warnings") = record$diagnostics
    attr(image, "path_warning_log") = record$log
  }
  image
}
