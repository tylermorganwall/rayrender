# Copyright 2017-2021  Kevin Ushey
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

# configure-database.R -------------------------------------------------------

#' Retrieve the Global Configuration Database
#'
#' Retrieve the global configuration database.
#' `db` is a helper alias for the database
#' returned by `configure_database()`.
#'
#' @export
configure_database <- local({
	database <- new.env(parent = emptyenv())
	class(database) <- "configure_database"
	function() database
})

#' @export
print.configure_database <- function(x, ...) {
	str.configure_database(x, ...)
}

#' @export
str.configure_database <- function(object, ...) {
	writeLines("<configure database>")
	objects <- mget(ls(envir = object, all.names = TRUE), object)
	output <- utils::capture.output(utils::str(objects, ...))
	writeLines(output[-1])
	invisible(output)
}


#' Define Variables for the Configuration Database
#'
#' Define variables to be used as part of the default configuration database.
#' These will be used by [configure_file()] when no configuration database
#' is explicitly supplied. [define()] is provided as a shorter alias for the
#' same function.
#'
#' @param ... A set of named arguments, mapping configuration names to values.
#'
#' @export
configure_define <- function(...) {
	envir <- configure_database()
	list2env(list(...), envir = envir)
}

#' @rdname configure_define
#' @export
define <- configure_define

#' @rdname configure_database
#' @export
db <- configure_database()


# utils.R --------------------------------------------------------------------

#' Configure a File
#'
#' Configure a file, replacing (by default) any instances of `@`-delimited
#' variables, e.g. `@VAR@`, with the value of the variable called `VAR` in the
#' associated `config` environment.
#'
#' @param source The file to be configured.
#' @param target The file to be generated.
#' @param config The configuration database.
#' @param lhs The left-hand side marker; defaults to `@`.
#' @param rhs The right-hand side marker; defaults to `@`.
#' @param verbose Boolean; report files as they are configured?
#'
#' @family configure
#'
#' @export
configure_file <- function(
	source,
	target = sub("[.]in$", "", source),
	config = configure_database(),
	lhs = "@",
	rhs = "@",
	verbose = configure_verbose()
) {
	# read source file
	contents <- readLines(source, warn = FALSE)

	# replace defined variables
	enumerate(config, function(key, val) {
		needle <- paste(lhs, key, rhs, sep = "")
		replacement <- val
		contents <<- gsub(needle, replacement, contents, fixed = TRUE)
	})

	ensure_directory(dirname(target))

	# write configured file to target location
	# prefer unix newlines for Makevars
	mode <- if (target %in% "Makevars") "wb" else "w"
	conn <- file(target, open = mode)
	on.exit(close(conn), add = TRUE)
	writeLines(contents, con = conn)

	# copy over source permissions
	info <- file.info(source)
	Sys.chmod(target, mode = info$mode)

	if (isTRUE(verbose)) {
		fmt <- "*** configured file: '%s' => '%s'"
		message(sprintf(fmt, source, target))
	}
}

#' Configure Files in a Directory
#'
#' This companion function to [configure_file()] can be used to
#' configure all `.in` files within a directory.
#'
#' @param path The path to a directory in which files should be configured.
#' @param config The configuration database to be used.
#' @param verbose Boolean; report files as they are configured?
#'
#' @family configure
#'
#' @export
configure_directory <- function(
	path = ".",
	config = configure_database(),
	verbose = configure_verbose()
) {
	files <- list.files(
		path = path,
		pattern = "[.]in$",
		full.names = TRUE
	)

	lapply(files, configure_file, config = config, verbose = verbose)
}

configure_auto <- function(type) {
	if (!isTRUE(getOption("configure.auto", default = TRUE)))
		return(invisible(FALSE))

	if (isTRUE(getOption("configure.common", default = TRUE)))
		configure_common(type = type)

	if (isTRUE(getOption("configure.platform", default = TRUE)))
		configure_platform(type = type)
}

configure_common <- function(type) {
	sources <- list.files(
		path = c("R", "src"),
		pattern = "[.]in$",
		full.names = TRUE
	)

	sources <- sub("[.]/", "", sources)

	if (type == "configure") {
		lapply(sources, configure_file)
	} else if (type == "cleanup") {
		targets <- sub("[.]in$", "", sources)
		lapply(targets, remove_file)
	}

	invisible(TRUE)
}

configure_platform <- function(type) {
	sysname <- tolower(Sys.info()[["sysname"]])

	subdirs <- sysname
	if (sysname != "windows") subdirs <- c("unix", subdirs)

	dirs <- c("R", "src")
	for (dir in dirs) {
		# list files (take care to remove directories)
		sources <- Filter(
			function(file) identical(file.info(file)$isdir, FALSE),
			list.files(file.path(dir, subdirs), full.names = TRUE)
		)

		# configure all discovered sources
		for (source in sources) {
			target <- file.path(dir, basename(source))
			switch(
				type,
				configure = configure_file(source, target),
				cleanup = remove_file(target)
			)
		}
	}
}

#' Execute R CMD config
#'
#' Read information about how \R is configured as through `R CMD config`.
#'
#' @param ... The names of potential configuration values.
#' @param simplify Boolean; simplify in the case where a single value was
#'   requested?
#'
#' @export
r_cmd_config <- function(..., simplify = TRUE) {
	R <- file.path(R.home("bin"), "R")

	# suppress cygwin path warnings for windows
	if (Sys.info()[["sysname"]] == "Windows") {
		CYGWIN <- Sys.getenv("CYGWIN")
		Sys.setenv(CYGWIN = "nodosfilewarning")
		on.exit(Sys.setenv(CYGWIN = CYGWIN), add = TRUE)
	}

	# loop through requested values and call R CMD config
	values <- unlist(list(...), recursive = TRUE)
	config <- lapply(values, function(value) {
		# execute it
		stdout <- tempfile("r-cmd-config-", fileext = ".txt")
		on.exit(unlink(stdout), add = TRUE)
		status <- system2(R, c("CMD", "config", value), stdout = stdout)

		# report failures as NULL (distinct from empty string)
		if (status) return(NULL)

		readLines(stdout)
	})

	names(config) <- values

	if (simplify && length(config) == 1) return(config[[1]])

	config
}

#' Read R Configuration for a Package
#'
#' Read the \R configuration, as through `R CMD config`.
#'
#' @param ... The \R configuration values to read (as a character vector).
#'   If empty, all values are read as through `R CMD config --all`).
#' @param package The path to the \R package's sources.
#' @param envir The environment in which the configuration information should
#'   be assigned. By default, the [configure_database()] is populated with the
#'   requested values.
#' @param verbose Boolean; notify the user as \R configuration is read?
#'
#' @export
read_r_config <- function(
	...,
	package = Sys.getenv("R_PACKAGE_DIR", unset = "."),
	envir = configure_database(),
	verbose = configure_verbose()
) {
	# move to requested directory
	owd <- setwd(package)
	on.exit(setwd(owd), add = TRUE)
	R <- file.path(R.home("bin"), "R")

	# suppress cygwin path warnings for windows
	if (Sys.info()[["sysname"]] == "Windows") {
		CYGWIN <- Sys.getenv("CYGWIN")
		Sys.setenv(CYGWIN = "nodosfilewarning")
		on.exit(Sys.setenv(CYGWIN = CYGWIN), add = TRUE)
	}

	values <- unlist(list(...), recursive = TRUE)
	if (length(values) == 0) {
		# R CMD config --all only available since R 3.4.0
		if (getRversion() < "3.4.0") {
			fmt <- "'R CMD config --all' not available in R version '%s'"
			stop(sprintf(fmt, getRversion()))
		}

		# execute action
		stdout <- tempfile("r-cmd-config-", fileext = ".txt")
		on.exit(unlink(stdout), add = TRUE)
		status <- system2(R, c("CMD", "config", "--all"), stdout = stdout)
		if (status) stop("failed to execute 'R CMD config --all'")

		# read and parse output
		output <- readLines(stdout, warn = FALSE)
		config <- parse_key_value(output)
	} else {
		# loop through requested values and call R CMD config
		config <- lapply(values, function(value) {
			# execute it
			stdout <- tempfile("r-cmd-config-", fileext = ".txt")
			on.exit(unlink(stdout), add = TRUE)
			status <- system2(R, c("CMD", "config", value), stdout = stdout)

			# report failures as NULL (distinct from empty string)
			if (status) return(NULL)

			readLines(stdout)
		})
		names(config) <- values
	}

	if (is.null(envir)) return(config)

	list2env(config, envir = envir)
}

#' Concatenate the Contents of a Set of Files
#'
#' Given a set of files, concatenate their contents into
#' a single file.
#'
#' @param sources An \R list of files
#' @param target The file to use for generation.
#' @param headers Headers to be used for each file copied.
#' @param preamble Text to be included at the beginning of the document.
#' @param postamble Text to be included at the end of the document.
#' @param verbose Boolean; inform the user when the requested file is created?
#'
#' @export
concatenate_files <- function(
	sources,
	target,
	headers = section_header(basename(sources)),
	preamble = NULL,
	postamble = NULL,
	verbose = configure_verbose()
) {
	pieces <- vapply(
		seq_along(sources),
		function(i) {
			source <- sources[[i]]
			header <- headers[[i]]
			contents <- trim_whitespace(read_file(source))
			paste(header, contents, "", sep = "\n\n")
		},
		character(1)
	)

	all <- c(preamble, pieces, postamble)

	ensure_directory(dirname(target))
	writeLines(all, con = target)

	if (verbose) {
		fmt <- "*** created file '%s'"
		message(sprintf(fmt, target))
	}

	TRUE
}

#' Add Configure Infrastructure to an R Package
#'
#' Add the infrastructure needed to configure an R package.
#'
#' @param package The path to the top-level directory of an \R package.
#' @export
use_configure <- function(package = ".") {
	# preserve working directory
	owd <- getwd()
	on.exit(setwd(owd), add = TRUE)

	# find resources
	package <- normalizePath(package, winslash = "/")
	resources <- system.file("resources", package = "configure")

	# copy into temporary directory
	dir <- tempfile("configure-")
	on.exit(unlink(dir, recursive = TRUE), add = TRUE)

	dir.create(dir)
	file.copy(resources, dir, recursive = TRUE)

	# rename resources directory
	setwd(dir)
	file.rename(basename(resources), basename(package))

	# now, copy these files back into the target directory
	file.copy(basename(package), dirname(package), recursive = TRUE)

	# ensure DESCRIPTION contains 'Biarch: TRUE' for Windows
	setwd(package)
	DESCRIPTION <- read_file("DESCRIPTION")
	if (!grepl("(?:^|\n)Biarch:", DESCRIPTION)) {
		DESCRIPTION <- paste(DESCRIPTION, "Biarch: TRUE", sep = "\n")
		DESCRIPTION <- gsub("\n{2,}", "\n", DESCRIPTION)
		cat(DESCRIPTION, file = "DESCRIPTION", sep = "\n")
	}

	# write placeholders for 'configure.R', 'cleanup.R' if none exist
	ensure_directory("tools/config")
	configure <- "tools/config/configure.R"
	if (!file.exists("tools/config/configure.R")) {
		text <- c(
			"# Prepare your package for installation here.",
			"# Use 'define()' to define configuration variables.",
			"# Use 'configure_file()' to substitute configuration values.",
			"",
			""
		)
		writeLines(text, con = configure)
	}

	cleanup <- "tools/config/cleanup.R"
	if (!file.exists("tools/config/cleanup.R")) {
		text <- c(
			"# Clean up files generated during configuration here.",
			"# Use 'remove_file()' to remove files generated during configuration.",
			"",
			""
		)
		writeLines(text, con = cleanup)
	}

	# notify the user what we did
	message("* Copied 'configure{.win}' and 'cleanup{.win}'.")
	message("* Updated 'tools/config.R'.")

	# open 'configure.R', 'cleanup.R' for editing if in RStudio
	rstudio <-
		!is.na(Sys.getenv("RSTUDIO", unset = NA)) &&
		requireNamespace("rstudioapi", quietly = TRUE)

	if (rstudio) {
		rstudioapi::navigateToFile("tools/config/configure.R", 5, 1)
		rstudioapi::navigateToFile("tools/config/cleanup.R", 4, 1)
	} else {
		message("* Use 'tools/config/configure.R' for package configuration.")
		message("* Use 'tools/config/cleanup.R' for package cleanup.")
	}
}

ensure_directory <- function(dir) {
	info <- file.info(dir)

	# no file exists at this location; try to make it
	if (is.na(info$isdir)) {
		dir.create(dir, recursive = TRUE, showWarnings = FALSE)
		if (!file.exists(dir)) stop("failed to create directory '", dir, "'")
		return(TRUE)
	}

	# a directory already exists
	if (isTRUE(info$isdir)) return(TRUE)

	# a file exists, but it's not a directory
	stop("file already exists at path '", dir, "'")
}

enumerate <- function(x, f, ...) {
	nms <- if (is.environment(x)) ls(envir = x) else names(x)
	lapply(nms, function(nm) {
		f(nm, x[[nm]], ...)
	})
}

read_file <- function(path) {
	paste(readLines(path, warn = FALSE), collapse = "\n")
}

remove_file <- function(
	path,
	verbose = configure_verbose()
) {
	info <- file.info(path)
	if (is.na(info$isdir)) return(TRUE)

	name <- if (info$isdir) "directory" else "file"

	unlink(path, recursive = isTRUE(info$isdir))
	if (file.exists(path)) {
		fmt <- "failed to remove %s '%s' (insufficient permissions?)"
		stop(sprintf(fmt, name, path))
	}

	if (verbose) {
		fmt <- "*** removed %s '%s'"
		message(sprintf(fmt, name, path))
	}

	TRUE
}

source_file <- function(
	path,
	envir = parent.frame()
) {
	contents <- read_file(path)
	invisible(eval(parse(text = contents), envir = envir))
}

trim_whitespace <- function(x) {
	gsub("^[[:space:]]*|[[:space:]]*$", "", x)
}

configure_verbose <- function() {
	getOption("configure.verbose", !interactive())
}

named <- function(object, nm) {
	names(object) <- nm
	object
}

parse_key_value <- function(
	text,
	separator = "=",
	trim = TRUE
) {
	# find the separator
	index <- regexpr(separator, text, fixed = TRUE)

	# split into parts
	keys <- substring(text, 1, index - 1)
	vals <- substring(text, index + 1)

	# trim if requested
	if (trim) {
		keys <- trim_whitespace(keys)
		vals <- trim_whitespace(vals)
	}

	# put together into R list
	named(as.list(vals), keys)
}

move_directory <- function(source, target) {
	# ensure we're trying to move a directory
	info <- file.info(source)
	if (is.na(info$isdir)) {
		fmt <- "no directory exists at path '%s'"
		stop(sprintf(fmt, source), call. = FALSE)
	}

	if (!info$isdir) {
		fmt <- "'%s' exists but is not a directory"
		stop(sprintf(fmt, source), call. = FALSE)
	}

	# good to go -- let's move it
	unlink(target, recursive = TRUE)
	file.rename(source, target)
	unlink(source, recursive = TRUE)
}

section_header <- function(
	label,
	prefix = "#",
	suffix = "-",
	length = 78L
) {
	# figure out length of full header
	n <- length - nchar(label) - nchar(prefix) - 2L
	n[n < 0] <- 0

	# generate '-' suffixes
	tail <- vapply(
		n,
		function(i) {
			paste(rep(suffix, i), collapse = "")
		},
		character(1)
	)

	# join it all together
	paste(prefix, label, tail)
}


# run.R ----------------------------------------------------------------------

if (!interactive()) {
	# extract path to install script
	args <- commandArgs(TRUE)
	type <- args[[1]]

	# preserve working directory
	owd <- getwd()

	on.exit(setwd(owd), add = TRUE)

	# switch working directory to the calling scripts's directory as set
	# by the shell, in case the R working directory was set to something else
	basedir <- Sys.getenv("PWD", unset = NA)
	if (!is.na(basedir)) setwd(basedir)
	# report start of execution
	package <- Sys.getenv("R_PACKAGE_NAME", unset = "<unknown>")
	fmt <- "** preparing to %s package '%s' ..."
	message(sprintf(fmt, type, package))

	# execute the requested script
	path <- sprintf("tools/config/%s.R", type)
	if (file.exists(path)) source_file(path)

	# perform automatic configuration
	configure_auto(type = type)

	# report end of execution
	fmt <- "** finished %s for package '%s'"
	message(sprintf(fmt, type, package))
}
