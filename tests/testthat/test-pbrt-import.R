test_that("PBRT Import isolates graphics state and shares object definitions", {
  directory = tempfile('pbrt-import-')
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  root = file.path(directory, 'scene.pbrt')
  writeLines(
    c(
      'WorldBegin',
      'MakeNamedMaterial "shared" "string type" "diffuse" "rgb reflectance" [.2 .3 .4]',
      'NamedMaterial "shared" Translate 1 0 0',
      'Import "child.pbrt"',
      'Shape "sphere"'
    ),
    root
  )
  writeLines(
    c(
      'Translate 2 0 0 Shape "sphere"',
      'MakeNamedMaterial "private" "string type" "diffuse"',
      'ObjectBegin "local" Shape "sphere" ObjectEnd',
      'ObjectInstance "local"'
    ),
    file.path(directory, 'child.pbrt')
  )
  imported = read_pbrt(root)
  expect_equal(nrow(imported$scene), 3)
  matrices = lapply(imported$scene$transforms, function(x) {
    x$group_transform[[1]]
  })
  expect_equal(vapply(matrices, function(x) x[1, 4], numeric(1)), c(3, 3, 1))
  expect_length(imported$source_files, 2)
  write('NamedMaterial "private"', root, append = TRUE)
  expect_error(read_pbrt(root), 'Undefined material')
  lines = readLines(root)
  writeLines(c(head(lines, -1), 'ObjectInstance "local"'), root)
  expect_equal(nrow(read_pbrt(root)$scene), 4)
})

test_that("PBRT instances resolve forward and sibling-import definitions", {
  directory = tempfile('pbrt-forward-')
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  root = file.path(directory, 'scene.pbrt')
  writeLines('WorldBegin Import "uses.pbrt" Import "definitions.pbrt"', root)
  writeLines(
    'Translate 3 0 0 ObjectInstance "tree"',
    file.path(directory, 'uses.pbrt')
  )
  writeLines(
    'ObjectBegin "tree" Shape "sphere" ObjectEnd',
    file.path(directory, 'definitions.pbrt')
  )
  imported = read_pbrt(root)
  expect_equal(nrow(imported$scene), 1)
  expect_equal(imported$scene$transforms[[1]]$group_transform[[1]][1, 4], 3)
  writeLines(
    'ObjectBegin "tree" ObjectInstance "tree" ObjectEnd',
    file.path(directory, 'definitions.pbrt')
  )
  expect_error(read_pbrt(root), 'Cyclic object instance')
  writeLines(
    'ObjectBegin "other" Shape "sphere" ObjectEnd',
    file.path(directory, 'definitions.pbrt')
  )
  expect_error(read_pbrt(root), 'Undefined object: tree')
})

test_that("PBRT Import scopes and nested assets are checked", {
  directory = tempfile('pbrt-import-scope-')
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  root = file.path(directory, 'scene.pbrt')
  child = file.path(directory, 'child.pbrt')
  writeLines(
    'WorldBegin ObjectBegin "parts" Import "child.pbrt" ObjectEnd ObjectInstance "parts"',
    root
  )
  writeLines('Shape "sphere"', child)
  expect_equal(nrow(read_pbrt(root)$scene), 1)
  writeLines('AttributeBegin Shape "sphere"', child)
  expect_error(read_pbrt(root), 'Unclosed PBRT scope')
  writeLines('ObjectEnd', child)
  expect_error(read_pbrt(root), 'Unmatched ObjectEnd')
  writeLines('Import "child.pbrt"', child)
  expect_error(read_pbrt(root), 'Cyclic PBRT Include')
  writeLines('Import "child.pbrt" WorldBegin', root)
  expect_error(read_pbrt(root), 'Import must follow WorldBegin')
})

test_that("file-boundary compaction preserves unresolved forward references", {
  directory = tempfile('pbrt-incremental-')
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  root = file.path(directory, 'scene.pbrt')
  uses = file.path(directory, 'uses.pbrt')
  writeLines(
    c(rep('ObjectInstance "known"', 65536L), 'ObjectInstance "later"'),
    uses
  )
  writeLines(
    c(
      'WorldBegin ObjectBegin "known" Shape "sphere" ObjectEnd',
      'Include "uses.pbrt"',
      'ObjectBegin "later" Shape "sphere" "float radius" 2 ObjectEnd'
    ),
    root
  )
  scene = read_pbrt(root)$scene
  expect_equal(nrow(scene), 2)
  expect_equal(
    ncol(scene$shape_info[[1]]$shape_properties$instance_transforms),
    65536L
  )
  later = scene$shape_info[[2]]$shape_properties$original_scene[[1]]
  expect_equal(later$shape_info[[1]]$shape_properties$radius, 2)
  writeLines(head(readLines(root), -1), root)
  expect_error(
    read_pbrt(root),
    'uses[.]pbrt:65537: ObjectInstance: Undefined object: later'
  )
})
