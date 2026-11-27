# Prepare the installed reference headers for Metal's source-string compiler.
# The C preprocessor selects the MSL backend and expands includes/guards once;
# no compiler is invoked at render time, and no OpenPBR source is vendored here.
write_openpbr_metal = function(compiler, output) {
  header = system.file('include/openpbr/openpbr.h', package = 'openpbr')
  if (!nzchar(header)) {
    stop('The openpbr headers are required for Metal shading.')
  }
  expanded = tempfile(fileext = '.metal')
  on.exit(unlink(expanded))
  status = system2(
    compiler[[1]],
    c(
      compiler[-1],
      '-E',
      '-P',
      '-CC',
      '-x',
      'c++',
      '-D__METAL_VERSION__=230',
      '-DOPENPBR_LANGUAGE_TARGET_MSL=1',
      shQuote(header)
    ),
    stdout = expanded
  )
  if (status != 0L) {
    stop('Could not preprocess the OpenPBR Metal headers.')
  }
  source = readLines(expanded, warn = FALSE)
  if (any(grepl(')OPENPBR"', source, fixed = TRUE))) {
    stop('Unexpected raw-string delimiter in the OpenPBR headers.')
  }
  writeLines(c('R"OPENPBR(', source, ')OPENPBR"'), output)
}
