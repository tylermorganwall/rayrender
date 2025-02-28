if(length(find.package("xml2",quiet=TRUE)) > 0 && isTRUE(as.logical(Sys.getenv("NOT_CRAN", "true")))) {
  run_cpp_tests("rayrender")
}