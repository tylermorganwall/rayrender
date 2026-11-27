test_that("descriptor casts preserve records, names and metadata", {
  classes = c(
    "ray_material",
    "ray_shape_info",
    "ray_transform",
    "ray_animated_transform"
  )
  for (cl in classes) {
    x = vctrs::new_vctr(list(list(a = 1), list(b = "two")), class = cl)
    named = x
    names(named) = c("one", "two")
    extra = x
    attr(extra, "metadata") = list(source = "custom", units = "metres")
    subclass = x
    class(subclass) = c("custom_descriptor", class(x))
    for (value in list(x, named, extra, subclass, x[FALSE])) {
      prototype = vctrs::vec_ptype(value)
      expect_identical(vctrs::vec_cast(value, prototype), value)
      expect_identical(
        vctrs::vec_cast(value, prototype),
        vctrs::vec_default_cast(value, prototype)
      )
      expect_identical(
        vctrs::vec_ptype2(value, prototype),
        vctrs::vec_default_ptype2(value, prototype)
      )
    }
    expect_error(
      vctrs::vec_cast(extra, x),
      class = "vctrs_error_incompatible_type"
    )
    expect_error(
      vctrs::vec_ptype2(extra, x),
      class = "vctrs_error_incompatible_type"
    )
    other = vctrs::new_vctr(list(list(a = 1)), class = setdiff(classes, cl)[1])
    expect_error(
      vctrs::vec_cast(x, other),
      class = "vctrs_error_incompatible_type"
    )
    expect_error(
      vctrs::vec_ptype2(x, other),
      class = "vctrs_error_incompatible_type"
    )
  }
})

test_that("typed and inferred scene assembly retain mixed rows and empty rows", {
  rows = list(
    sphere(material = diffuse("red")),
    cube(material = metal("gold")),
    cylinder(material = dielectric())
  )
  # Base data-frame binding supplies an independent assembly reference.
  expected = do.call(rbind, rows)
  rows = append(rows, list(rows[[1]][FALSE, ]), after = 1L)
  expect_identical(
    vctrs::list_unchop(rows, ptype = expected[FALSE, ]),
    expected
  )
  expect_identical(vctrs::list_unchop(rows), expected)
  expect_identical(
    vctrs::list_unchop(list(expected[FALSE, ])),
    expected[FALSE, ]
  )
})
