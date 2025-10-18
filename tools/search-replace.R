#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
fn <- args[[1]]
stringval = args[[2]]
replace = args[[3]]
txt <- readLines(fn)
txt_new <- gsub(
	stringval,
	replace,
	txt
)
newlines = txt_new != txt
if (sum(newlines) > 0) {
	message(
		sprintf("Replaced the following lines in '%s':\n", fn),
		paste0(
			sprintf(" Old: '%s'\n New: '%s'", txt[newlines], txt_new[newlines]),
			collapse = "\n"
		)
	)
} else {
	message(
		sprintf("Did not find any changes to make in '%s'", fn)
	)
}
newline = if (.Platform$OS.type == "windows") "\r\n" else "\n"
writeLines(txt_new, fn, sep = newline)
