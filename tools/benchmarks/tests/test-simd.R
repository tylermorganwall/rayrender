testthat::test_that("scalar and supported SIMD backends compile and return safe integer signs", {
  compiler = Sys.which("clang++")
  if (!nzchar(compiler)) {
    compiler = Sys.which("g++")
  }
  testthat::skip_if(!nzchar(compiler), "No standalone C++ compiler")
  root = tempfile()
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  arch = Sys.info()[["machine"]]
  x86 = grepl("x86|amd64", arch)
  rosetta = Sys.info()[["sysname"]] == "Darwin" &&
    grepl("arm", arch) &&
    system2(
      "/usr/bin/arch",
      c("-x86_64", "/usr/bin/true"),
      stdout = FALSE,
      stderr = FALSE
    ) ==
      0
  modes = list(scalar = character())
  if (x86 || rosetta) {
    target = if (rosetta) c("-arch", "x86_64") else character()
    modes = c(
      modes,
      list(
        scalar_x86 = c(target, "-msse4.1"),
        sse2 = c(target, "-msse2", "-mno-sse4.1", "-DRAYSIMD", "-DHAS_SSE"),
        sse41 = c(target, "-msse4.1", "-DRAYSIMD", "-DHAS_SSE")
      )
    )
  }
  if (grepl("arm|aarch64", arch)) {
    modes$neon = c("-DRAYSIMD", "-DHAS_NEON")
  }
  for (mode in names(modes)) {
    exe = file.path(root, mode)
    args = c(
      "-std=c++20",
      "-O3",
      modes[[mode]],
      paste0("-I", benchmark_repo),
      file.path(benchmark_repo, "tools/benchmarks/tests/simd_sign.cpp"),
      "-o",
      exe
    )
    output = suppressWarnings(system2(
      compiler,
      shQuote(args),
      stdout = TRUE,
      stderr = TRUE
    ))
    testthat::expect_equal(
      attr(output, "status") %||% 0L,
      0L,
      info = paste(output, collapse = "\n")
    )
    if (file.exists(exe)) testthat::expect_equal(system2(exe), 0L, info = mode)
  }
})
