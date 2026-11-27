test_that("PBRT buffers preserve order and references across chunk boundaries", {
  for (size in c(0L, 1L, 255L, 256L, 257L, 513L)) {
    buffer = pbrt_row_buffer()
    rows = lapply(seq_len(size), function(i) list(index = i))
    for (row in rows) {
      pbrt_buffer_append(buffer, row)
    }
    pbrt_buffer_append(buffer, NULL)
    expect_identical(pbrt_buffer_rows(buffer), rows)
    expect_equal(buffer$chunks, size %/% 256L)
    expect_length(buffer$pending, size %% 256L)
  }
  reference = structure(
    list(command = list(args = 'part')),
    class = 'pbrt_instance_reference'
  )
  buffer = pbrt_row_buffer()
  pbrt_buffer_append(buffer, reference)
  expect_identical(buffer$new_instances, 1L)
  expect_identical(pbrt_buffer_rows(buffer), list(reference))
})

test_that("buffer assembly retains shared nested prototype identity", {
  skip_if_not(capabilities('profmem'))
  row = create_instances(create_instances(sphere()))
  prototype = row$shape_info[[1]]$shape_properties$original_scene[[1]]
  address = tracemem(prototype)
  buffer = pbrt_row_buffer()
  for (i in 1:257) {
    pbrt_buffer_append(buffer, row)
  }
  rows = pbrt_buffer_rows(buffer)
  expect_length(rows, 257L)
  for (i in c(1L, 256L, 257L)) {
    child = rows[[i]]$shape_info[[1]]$shape_properties$original_scene[[1]]
    expect_identical(tracemem(child), address)
    untracemem(child)
  }
  untracemem(prototype)
})
