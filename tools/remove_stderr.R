#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# remove_stderr_base.R
#
# Delete every call of the form
#     fprintf ( stderr , ... );
# from C / C++ sources, even when the argument list spans
# several lines.  A “.bak” copy of each modified file is kept.
# ------------------------------------------------------------------

#' Strip \code{fprintf(stderr, …)} calls from a source tree
#'
#' @param root_dir Default `'.'`. Directory that contains the C / C++
#'   sources (searched recursively).
#' @param exts Default `c('c', 'cc', 'cpp', 'cxx', 'C')`. File‑name
#'   extensions that should be scanned.
#' @param backup_suffix Default `'.bak'`. Suffix for the backup copy
#'   written before each file is overwritten.
strip_stderr_calls = function(
	root_dir = ".",
	exts = c("c", "cc", "cpp", "cxx", "C"),
	backup_suffix = ".bak"
) {
	## PCRE pattern: (?s) turns on DOTALL so '.' matches newlines;
	## non‑greedy '.*?' stops at the first ');'
	pat = "(?s)fprintf\\s*\\(\\s*stderr\\s*,.*?\\);"

	files = list.files(
		root_dir,
		pattern = paste0("\\.(", paste(exts, collapse = "|"), ")$"),
		recursive = TRUE,
		full.names = TRUE
	)

	for (f in files) {
		size = file.info(f)$size
		txt = if (is.na(size) || size == 0) "" else
			readChar(f, size, useBytes = TRUE)

		new_txt = gsub(pat, ";", txt, perl = TRUE)

		if (!identical(txt, new_txt)) {
			file.copy(f, paste0(f, backup_suffix), overwrite = TRUE)
			writeChar(new_txt, f, eos = NULL, useBytes = TRUE)
			message("cleaned: ", f)
		}
	}
}

## When executed directly (not sourced) run on supplied directory,
## defaulting to '.'.
args = commandArgs(trailingOnly = TRUE)
strip_stderr_calls(root_dir = args[1])
unlink(list.files(pattern = "\\.bak", recursive = TRUE))
