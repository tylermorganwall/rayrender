test_that("Test rendering with and without OIDN", {
	set.seed(1)
	rgb_output = generate_cornell() |>
		render_scene(
			samples = 16,
			denoise = FALSE,
			preview = FALSE,
			print_debug_info = TRUE,
			return_raw_array = TRUE
		)
	if (abs(sum(rgb_output) - 389644.66) > 300) {
		warning(sprintf("Sum of rgb_output is %f", sum(rgb_output)))
	}
	expect_true(abs(sum(rgb_output) - 389644.66) < 300)
	if (rayrender:::cppdef_HAS_OIDN()) {
		set.seed(1)
		rgb_output_oidn = generate_cornell() |>
			render_scene(samples = 16, preview = FALSE, return_raw_array = TRUE)
		if (abs(sum(rgb_output_oidn) - 389644.66) > 300) {
			warning(sprintf("Sum of rgb_output_oidn is %f", sum(rgb_output_oidn)))
		}
		expect_true(abs(sum(rgb_output_oidn) - 389644.66) < 300)
	}
})

test_that("Test rendering with and without OIDN, snapshot", {
	set.seed(1)
	path = tempfile(fileext = ".png")
	generate_cornell() |>
		render_scene(
			samples = 16,
			denoise = FALSE,
			preview = FALSE,
			filename = path
		)
	expect_snapshot_file(path, variant = Sys.info()[["sysname"]])
	if (rayrender:::cppdef_HAS_OIDN()) {
		set.seed(1)
		path_oidn = tempfile(pattern = "oidn", fileext = ".png")

		generate_cornell() |>
			render_scene(samples = 16, preview = FALSE, filename = path_oidn)
		expect_snapshot_file(path_oidn, variant = Sys.info()[["sysname"]])
	}
})

test_that("Test different rendering integrator types, snapshot", {
	set.seed(1)
	types = c("rtiow", "nee", "basic")
	paths = tempfile(pattern = types, fileext = ".png")
	for (i in seq_along(types)) {
		set.seed(1)
		generate_cornell() |>
			render_scene(
				samples = 16,
				denoise = FALSE,
				preview = FALSE,
				filename = paths[i],
				integrator_type = types[i]
			)
		expect_snapshot_file(paths[i], variant = Sys.info()[["sysname"]])
	}
	if (rayrender:::cppdef_HAS_OIDN()) {
		types_oidn = c("rtiow", "nee", "basic")
		paths_oidn = tempfile(
			pattern = paste0("denoise_", types_oidn),
			fileext = ".png"
		)

		for (i in seq_along(types)) {
			set.seed(1)
			generate_cornell() |>
				render_scene(
					samples = 16,
					denoise = TRUE,
					preview = FALSE,
					filename = paths_oidn[i],
					integrator_type = types_oidn[i]
				)
			expect_snapshot_file(paths_oidn[i], variant = Sys.info()[["sysname"]])
		}
	}
})
