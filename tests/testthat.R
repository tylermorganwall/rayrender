# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)
library(rayrender)

if (!isTRUE(as.logical(Sys.getenv("RAY_COLOR_DEBUG", "false")))) {
  if (isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    test_check("rayrender", filter = c("integrator"), invert = TRUE)
  } else {
    test_check("rayrender", filter = c("integrator|cpp"), invert = TRUE)
  }
} else {
  test_check("rayrender", filter = "integrator")
}
