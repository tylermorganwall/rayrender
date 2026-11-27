test_that('packed transform operands preserve double precision without text round trips', {
  values = c(
    '1.0000000000000002',
    '-0',
    '1e-200',
    '0',
    '0.10000000000000001',
    '1',
    '0',
    '0',
    '0',
    '0',
    '1',
    '0',
    '123456789.12345679',
    '-2.2250738585072014e-308',
    '1.7976931348623157e308',
    '1'
  )
  file = tempfile(fileext = '.pbrt')
  on.exit(unlink(file))
  for (directive in c('Transform', 'ConcatTransform')) {
    writeLines(
      sprintf('%s [%s]', directive, paste(values, collapse = ' ')),
      file
    )
    parsed = pbrt_parse_file(file)[[1]]
    expect_type(parsed$args, 'double')
    expect_identical(parsed$args, as.numeric(values))
    expect_identical(1 / parsed$args, 1 / as.numeric(values))
    writeLines(sprintf('%s %s', directive, paste(values, collapse = ' ')), file)
    unbracketed = pbrt_parse_file(file)[[1]]
    expect_identical(as.numeric(unbracketed$args), parsed$args)
  }
})

test_that('numeric matrix packing preserves arity and nonfinite errors', {
  file = tempfile(fileext = '.pbrt')
  on.exit(unlink(file))
  for (count in c(0L, 15L, 17L)) {
    writeLines(
      sprintf('Transform [%s]', paste(rep(1, count), collapse = ' ')),
      file
    )
    expect_error(pbrt_parse_file(file), '16 bracketed numbers')
  }
  values = as.character(as.vector(diag(4)))
  writeLines(sprintf('Transform [%s] 2', paste(values, collapse = ' ')), file)
  expect_error(pbrt_parse_file(file), 'Unexpected extra operands')
  values[13] = '1e999'
  writeLines(
    c('WorldBegin', sprintf('Transform [%s]', paste(values, collapse = ' '))),
    file
  )
  expect_error(read_pbrt(file), 'finite numbers')
})
