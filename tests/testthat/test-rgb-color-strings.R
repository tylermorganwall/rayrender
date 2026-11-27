test_that("RGB strings use byte-scale channels while numeric vectors stay normalized", {
  expected = c(129, 123, 201) / 255
  for (color in c(
    "rgb(129, 123, 201)",
    "rgb(129,123,201)",
    " RGB ( 129 , 123 , 201 ) ",
    "rgb(129.0,\t123.0,201)"
  )) {
    expect_equal(convert_color(color), expected)
    expect_identical(toupper(convert_color(color, as_hex = TRUE)), "#817BC9")
  }
  expect_equal(convert_color(expected), expected)
  expect_equal(convert_color("rgb(0, 0, 0)"), c(0, 0, 0))
  expect_equal(convert_color("rgb(255, 255, 255)"), c(1, 1, 1))
  expect_equal(convert_color("rgb(1, 0, 0)"), c(1 / 255, 0, 0))
  expect_equal(
    convert_color("rgb(12.5, 0, 255)"),
    convert_color(grDevices::rgb(12.5, 0, 255, maxColorValue = 255))
  )
})

test_that("RGB normalization preserves mixed palettes and non-string inputs", {
  colors = c(
    css = "rgb(129, 123, 201)",
    named = "red",
    hex = "#33d",
    alpha = "#11223380",
    missing = NA_character_,
    clear = "transparent"
  )
  expected = colors
  expected[1] = "#817BC9"
  expect_identical(normalize_rgb_color(colors), expected)
  expect_identical(normalize_rgb_color(character()), character())
  for (value in list(NULL, NA, TRUE, c(0.1, 0.2, 0.3), matrix(0.5, 2, 3))) {
    expect_identical(normalize_rgb_color(value), value)
  }
  expect_equal(convert_color(c("rgb(255, 0, 0)", "blue")), c(1, 0, 0, 0, 0, 1))
  expect_error(convert_color("not-a-color"), "invalid color")
})

test_that("invalid RGB strings fail without evaluating their contents", {
  for (color in c(
    "rgb(256, 0, 0)",
    "rgb(-1, 0, 0)",
    "rgb(1, 2)",
    "rgb(1, 2, 3, 4)",
    "rgb(1, 2,)",
    "rgb(NA, 2, 3)",
    "rgb(Inf, 2, 3)",
    "rgb(1%, 2%, 3%)",
    "rgb(1, 2, 3) trailing",
    "rgb(stop('evaluated'), 2, 3)"
  )) {
    expect_error(convert_color(color), "RGB color")
  }
})

test_that("material constructors and texture graphs match equivalent hexadecimal colors", {
  color = "rgb(129, 123, 201)"
  hex = "#817BC9"
  for (material in list(
    diffuse,
    metal,
    dielectric,
    microfacet,
    glossy,
    light,
    subsurface,
    subsurface_diffusion
  )) {
    expect_equal(material(color = color), material(color = hex))
  }
  expect_equal(
    diffuse(
      color,
      checkercolor = color,
      noisecolor = color,
      gradient_color = color
    ),
    diffuse(hex, checkercolor = hex, noisecolor = hex, gradient_color = hex)
  )
  expect_equal(openpbr(base_color = color), openpbr(base_color = hex))
  expect_equal(texture_constant(color), texture_constant(hex))
  expect_equal(point_light(color = color), point_light(color = hex))
  expect_equal(disk_light(color = color), disk_light(color = hex))
})

test_that("screen overlay constructors and data frames accept RGB strings", {
  color = "rgb(129, 123, 201)"
  hex = "#817BC9"
  expect_equal(screen_line(color = color), screen_line(color = hex))
  expect_equal(
    screen_text(
      "text",
      color = color,
      halo_color = color,
      background_color = color
    ),
    screen_text("text", color = hex, halo_color = hex, background_color = hex)
  )
  line = screen_line(color = hex)
  line$color = color
  expect_identical(normalize_screen_line(line)$color, hex)
  expect_equal(
    prepare_screen_line_preview(line),
    prepare_screen_line_preview(screen_line(color = hex))
  )
  text = screen_text("text", color = hex)
  text$color = color
  expect_identical(normalize_screen_text(text)$color, hex)
  expect_equal(
    rasterize_screen_line(
      0,
      0,
      2,
      2,
      1,
      normalize_screen_line(line)$color,
      1,
      "round"
    ),
    rasterize_screen_line(0, 0, 2, 2, 1, hex, 1, "round")
  )
})
