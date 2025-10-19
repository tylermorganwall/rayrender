# Configure rayrender using an R-based workflow.

is_windows = identical(.Platform$OS.type, "windows")
sysname = Sys.info()[["sysname"]]
is_macos = identical(sysname, "Darwin")
target_arch = Sys.info()[["machine"]]
r_home = R.home()

split_flags = function(flags) {
	if (!nzchar(flags)) {
		return(character())
	}
	strsplit(flags, "[[:space:]]+")[[1]]
}

flag_with_path = function(flag, path) {
	if (!nzchar(path)) {
		return(character())
	}
	path = normalizePath(path, winslash = "/", mustWork = FALSE)
	if (grepl("\\s", path)) {
		c(flag, path)
	} else {
		paste0(flag, path)
	}
}

append_flags = function(existing, ...) {
	new_flags = unlist(list(...))
	new_flags = new_flags[nzchar(new_flags)]
	if (length(new_flags) == 0) {
		return(existing)
	}
	c(existing, new_flags)
}

append_unique_flags = function(existing, ...) {
	new_flags = unlist(list(...))
	new_flags = new_flags[nzchar(new_flags)]
	if (length(new_flags) == 0) {
		return(existing)
	}
	for (flag in new_flags) {
		if (!(flag %in% existing)) {
			existing = c(existing, flag)
		}
	}
	existing
}

collapse_flags = function(flags) {
	if (length(flags) == 0) {
		return("")
	}
	paste(flags, collapse = " ")
}

CXX = r_cmd_config("CXX20")
CXXFLAGS_BASE = split_flags(r_cmd_config("CXX20FLAGS"))
if (!nzchar(CXX)) {
	CXX = r_cmd_config("CXX17")
	CXXFLAGS_BASE = split_flags(r_cmd_config("CXX17FLAGS"))
}
if (!nzchar(CXX)) {
	stop("Failed to determine the C++ compiler via R CMD config.")
}
CPPFLAGS_BASE = split_flags(r_cmd_config("CPPFLAGS"))

split_command = function(command) {
	if (!nzchar(command)) {
		return(character())
	}
	tokens = strsplit(command, "[[:space:]]+")[[1]]
	tokens = tokens[nzchar(tokens)]
	if (length(tokens) == 0) {
		return(tokens)
	}
	# remove wrapping quotes if present
	gsub("(^['\"]|['\"]$)", "", tokens)
}

unwrap_compiler = function(tokens) {
	if (length(tokens) < 2) {
		return(tokens)
	}
	wrappers = c("ccache", "sccache", "icecc", "distcc")
	while (length(tokens) > 1 && tokens[[1]] %in% wrappers) {
		tokens = tokens[-1]
	}
	tokens
}

CXX_COMMAND = unwrap_compiler(split_command(CXX))
if (length(CXX_COMMAND) == 0) {
	stop("Failed to parse the C++ compiler command from R CMD config.")
}

build_command = function(base_tokens, extra_tokens = character()) {
	tokens = c(base_tokens, extra_tokens)
	tokens = tokens[nzchar(tokens)]
	if (!length(tokens)) {
		return("")
	}
	quote_type = if (is_windows) "cmd" else "sh"
	quoted = vapply(
		tokens,
		function(token) {
			if (!nzchar(token)) {
				return("")
			}
			shQuote(token, type = quote_type)
		},
		character(1),
		USE.NAMES = FALSE
	)
	paste(quoted[nzchar(quoted)], collapse = " ")
}

PKG_CPPFLAGS = character()
PKG_CXXFLAGS = append_unique_flags(
	character(),
	"-ffp-contract=off",
	"-fvisibility=hidden",
	"-fvisibility-inlines-hidden"
)
DEFINES = append_unique_flags(
	character(),
	"-DRAY_REPRODUCE_PERLIN",
	"-DSTRICT_R_HEADERS"
)
LINK_LIBS = character()
PKG_LIBS_ACC = character()
OIDN_CPPFLAGS = character()
OIDN_LDFLAGS = character()
OIDN_LIBS = character()

compile_test = function(
	code,
	extra_cppflags = character(),
	extra_cxxflags = character(),
	extra_ldflags = character(),
	link = FALSE,
	include_pkg_cppflags = TRUE,
	include_pkg_cxxflags = TRUE
) {
	src = tempfile(fileext = ".cpp")
	writeLines(code, src)
	on.exit(unlink(src), add = TRUE)

	output_ext = if (link) {
		if (is_windows) ".exe" else ""
	} else {
		".o"
	}
	output = tempfile(fileext = output_ext)
	on.exit(unlink(output), add = TRUE)

	args = c(
		CPPFLAGS_BASE,
		if (include_pkg_cppflags) PKG_CPPFLAGS else character(),
		extra_cppflags,
		CXXFLAGS_BASE,
		if (include_pkg_cxxflags) PKG_CXXFLAGS else character(),
		extra_cxxflags
	)

	if (link) {
		args = c(args, src, extra_ldflags, "-o", output)
	} else {
		args = c(args, "-c", src, "-o", output)
	}

	cmd = build_command(CXX_COMMAND, args)
	status = tryCatch(
		{
			system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
		},
		warning = function(w) 1L,
		error = function(e) 1L
	)

	isTRUE(status == 0)
}

if (is_windows) {
	DEFINES = append_unique_flags(
		DEFINES,
		"-DRAY_WINDOWS"
	)
	PKG_LIBS_ACC = append_unique_flags(
		PKG_LIBS_ACC,
		"-lgdi32"
	)
}

if (!is_windows) {
	pkgconfig = Sys.which("pkg-config")
	x11_cppflags = character()
	x11_ldflags = character()
	x11_found = FALSE

	if (nzchar(pkgconfig)) {
		status = tryCatch(
			system2(pkgconfig, c("--exists", "x11"), stdout = FALSE, stderr = FALSE),
			warning = function(w) 1L,
			error = function(e) 1L
		)
		if (isTRUE(status == 0)) {
			x11_cppflags = split_flags(paste(
				system2(pkgconfig, c("--cflags", "x11"), stdout = TRUE),
				collapse = " "
			))
			x11_ldflags = split_flags(paste(
				system2(pkgconfig, c("--libs", "x11"), stdout = TRUE),
				collapse = " "
			))
			x11_found = compile_test(
				code = "#include <X11/Xlib.h>\nint main() { XOpenDisplay(NULL); return 0; }\n",
				extra_cppflags = x11_cppflags,
				extra_ldflags = x11_ldflags,
				link = TRUE
			)
		}
	}

	if (!x11_found) {
		fallback_cpp = character()
		fallback_ld = character()
		if (is_macos) {
			fallback_cpp = append_flags(
				fallback_cpp,
				flag_with_path("-I", "/opt/X11/include")
			)
			fallback_ld = append_flags(
				fallback_ld,
				flag_with_path("-L", "/opt/X11/lib")
			)
		}
		fallback_ld = append_flags(fallback_ld, "-lX11")
		x11_found = compile_test(
			code = "#include <X11/Xlib.h>\nint main() { XOpenDisplay(NULL); return 0; }\n",
			extra_cppflags = fallback_cpp,
			extra_ldflags = fallback_ld,
			link = TRUE
		)
		if (x11_found) {
			x11_cppflags = fallback_cpp
			x11_ldflags = fallback_ld
		}
	}

	if (x11_found) {
		message("*** configure: found X11 headers and libraries")
		PKG_CPPFLAGS = append_flags(PKG_CPPFLAGS, x11_cppflags)
		LINK_LIBS = append_flags(LINK_LIBS, x11_ldflags)
		DEFINES = append_unique_flags(DEFINES, "-DRAY_HAS_X11")
		if (
			is_macos &&
				!any(grepl("/opt/X11/include", x11_cppflags, fixed = TRUE))
		) {
			PKG_CPPFLAGS = append_flags(
				PKG_CPPFLAGS,
				flag_with_path("-I", "/opt/X11/include")
			)
		}
	} else {
		message("*** configure: X11 not found; render preview will be disabled")
	}
}

oidn_path = Sys.getenv("OIDN_PATH", unset = "")
if (nzchar(oidn_path) && dir.exists(oidn_path)) {
	message(sprintf("*** configure: using Open Image Denoise at %s", oidn_path))
	OIDN_CPPFLAGS = append_flags(
		OIDN_CPPFLAGS,
		flag_with_path("-I", file.path(oidn_path, "include"))
	)
	lib_dir = file.path(oidn_path, "lib")
	if (dir.exists(lib_dir)) {
		OIDN_LDFLAGS = append_flags(OIDN_LDFLAGS, flag_with_path("-L", lib_dir))
		OIDN_LDFLAGS = append_flags(
			OIDN_LDFLAGS,
			sprintf(
				"-Wl,-rpath,%s",
				normalizePath(lib_dir, winslash = "/", mustWork = FALSE)
			)
		)
		OIDN_LIBS = append_flags(OIDN_LIBS, "-lOpenImageDenoise", "-ltbb")
	}
	DEFINES = append_unique_flags(DEFINES, "-DHAS_OIDN")
} else {
	message(
		"*** configure: Open Image Denoise (OIDN) not found; skipping denoiser support"
	)
}

openexr_lib_dir = tryCatch(
	normalizePath(
		system.file("lib", package = "libopenexr", mustWork = TRUE),
		winslash = "/",
		mustWork = TRUE
	),
	error = function(e) {
		stop(
			"OpenEXR not found; please install libopenexr or ensure it is available on the library path.",
			call. = FALSE
		)
	}
)
openexr_inc_dir = tryCatch(
	normalizePath(
		system.file("include", "OpenEXR", package = "libopenexr", mustWork = TRUE),
		winslash = "/",
		mustWork = TRUE
	),
	error = function(e) {
		stop("OpenEXR headers not found; please install libopenexr.", call. = FALSE)
	}
)
openexr_lib_arch = file.path(openexr_lib_dir, target_arch)
if (!dir.exists(openexr_lib_arch)) {
	stop(
		sprintf(
			"OpenEXR library directory for architecture '%s' not found at %s",
			target_arch,
			openexr_lib_arch
		),
		call. = FALSE
	)
}
PKG_CPPFLAGS = append_flags(
	PKG_CPPFLAGS,
	flag_with_path("-I", openexr_inc_dir)
)
PKG_LIBS_ACC = append_flags(
	PKG_LIBS_ACC,
	flag_with_path("-L", openexr_lib_arch),
	"-lIlmThread-3_4",
	"-lOpenEXR-3_4",
	"-lIex-3_4",
	"-lOpenEXRCore-3_4",
	"-lOpenEXRUtil-3_4"
)

imath_lib_dir = tryCatch(
	normalizePath(
		system.file("lib", package = "libimath", mustWork = TRUE),
		winslash = "/",
		mustWork = TRUE
	),
	error = function(e) {
		stop("Imath not found; please install libimath.", call. = FALSE)
	}
)
imath_inc_dir = tryCatch(
	normalizePath(
		system.file("include", "Imath", package = "libimath", mustWork = TRUE),
		winslash = "/",
		mustWork = TRUE
	),
	error = function(e) {
		stop("Imath headers not found; please install libimath.", call. = FALSE)
	}
)
imath_lib_arch = file.path(imath_lib_dir, target_arch)
if (!dir.exists(imath_lib_arch)) {
	stop(
		sprintf(
			"Imath library directory for architecture '%s' not found at %s",
			target_arch,
			imath_lib_arch
		),
		call. = FALSE
	)
}
PKG_CPPFLAGS = append_flags(
	PKG_CPPFLAGS,
	flag_with_path("-I", imath_inc_dir)
)
PKG_LIBS_ACC = append_flags(
	PKG_LIBS_ACC,
	flag_with_path("-L", imath_lib_arch),
	"-lImath-3_2"
)

libdeflate_lib_dir = tryCatch(
	normalizePath(
		system.file("lib", package = "libdeflate", mustWork = TRUE),
		winslash = "/",
		mustWork = TRUE
	),
	error = function(e) {
		stop("libdeflate not found; please install libdeflate.", call. = FALSE)
	}
)
libdeflate_inc_dir = tryCatch(
	normalizePath(
		system.file("include", package = "libdeflate", mustWork = TRUE),
		winslash = "/",
		mustWork = TRUE
	),
	error = function(e) {
		stop(
			"libdeflate headers not found; please install libdeflate.",
			call. = FALSE
		)
	}
)
libdeflate_lib_arch = file.path(libdeflate_lib_dir, target_arch)
if (!dir.exists(libdeflate_lib_arch)) {
	stop(
		sprintf(
			"libdeflate library directory for architecture '%s' not found at %s",
			target_arch,
			libdeflate_lib_arch
		),
		call. = FALSE
	)
}
PKG_CPPFLAGS = append_flags(
	PKG_CPPFLAGS,
	flag_with_path("-I", libdeflate_inc_dir)
)
PKG_LIBS_ACC = append_flags(
	PKG_LIBS_ACC,
	flag_with_path("-L", libdeflate_lib_arch),
	"-ldeflate"
)

sse_checked = FALSE
if (
	compile_test(
		"#include <smmintrin.h>\nint main() { __m128 v = _mm_dp_ps(_mm_set1_ps(1.0f), _mm_set1_ps(2.0f), 0xFF); (void)v; return 0; }\n"
	)
) {
	PKG_CXXFLAGS = append_unique_flags(PKG_CXXFLAGS, "-msse4.1")
	DEFINES = append_unique_flags(
		DEFINES,
		"-DHAS_SSE",
		"-DHAS_SSE41",
		"-DRAYSIMD",
		"-DRAYSIMDVECOFF"
	)
	sse_checked = TRUE
	message("*** configure: enabling SSE4.1 support")
}

if (
	!sse_checked &&
		compile_test(
			"#include <pmmintrin.h>\nint main() { __m128 v = _mm_hadd_ps(_mm_set1_ps(1.0f), _mm_set1_ps(1.0f)); (void)v; return 0; }\n"
		)
) {
	PKG_CXXFLAGS = append_unique_flags(PKG_CXXFLAGS, "-msse3")
	DEFINES = append_unique_flags(
		DEFINES,
		"-DHAS_SSE",
		"-DHAS_SSE3",
		"-DRAYSIMD",
		"-DRAYSIMDVECOFF"
	)
	sse_checked = TRUE
	message("*** configure: enabling SSE3 support")
}

if (
	!sse_checked &&
		compile_test(
			"#include <emmintrin.h>\nint main() { __m128d v = _mm_setzero_pd(); (void)v; return 0; }\n"
		)
) {
	PKG_CXXFLAGS = append_unique_flags(PKG_CXXFLAGS, "-msse2")
	DEFINES = append_unique_flags(
		DEFINES,
		"-DHAS_SSE",
		"-DHAS_SSE2",
		"-DRAYSIMD",
		"-DRAYSIMDVECOFF"
	)
	sse_checked = TRUE
	message("*** configure: enabling SSE2 support")
}

if (
	!sse_checked &&
		compile_test(
			"#include <xmmintrin.h>\nint main() { __m128 v = _mm_setzero_ps(); (void)v; return 0; }\n"
		)
) {
	PKG_CXXFLAGS = append_unique_flags(PKG_CXXFLAGS, "-msse")
	DEFINES = append_unique_flags(
		DEFINES,
		"-DHAS_SSE",
		"-DRAYSIMD",
		"-DRAYSIMDVECOFF"
	)
	sse_checked = TRUE
	message("*** configure: enabling SSE support")
}

if (!sse_checked) {
	message("*** configure: SSE intrinsics not available")
}

if (
	compile_test(
		"#if defined(__ARM_NEON) || defined(__ARM_NEON__)\n#include <arm_neon.h>\n#endif\nint main() { float32x4_t v = vdupq_n_f32(0.0f); (void)v; return 0; }\n"
	)
) {
	DEFINES = append_unique_flags(
		DEFINES,
		"-DHAS_NEON",
		"-DRAYSIMD",
		"-DRAYSIMDVECOFF"
	)
	if (is_macos) {
		PKG_CXXFLAGS = append_unique_flags(PKG_CXXFLAGS, "-march=armv8-a+simd")
	} else {
		PKG_CXXFLAGS = append_unique_flags(PKG_CXXFLAGS, "-mfpu=neon")
	}
	message("*** configure: enabling NEON support")
}

ray_color_debug = tolower(Sys.getenv("RAY_COLOR_DEBUG", unset = "false"))
if (identical(ray_color_debug, "true")) {
	message("*** configure: enabling RAY_COLOR_DEBUG")
	DEFINES = append_unique_flags(DEFINES, "-DRAY_COLOR_DEBUG")
} else {
	message("*** configure: RAY_COLOR_DEBUG disabled")
}

collect_sources = function(subdir, pattern) {
	if (!dir.exists(file.path("src", subdir))) {
		return(character())
	}
	files = list.files(
		file.path("src", subdir),
		pattern = pattern,
		full.names = FALSE,
		no.. = TRUE
	)
	if (length(files) == 0) {
		return(character())
	}
	file.path(subdir, files)
}

DIR_SOURCES = sort(list.files("src", pattern = "\\.cpp$", full.names = FALSE))
SUBDIR_SOURCES = sort(unlist(lapply(
	c("core", "hitables", "materials", "math", "utils"),
	collect_sources,
	pattern = "\\.cpp$"
)))
EXT_CPP_SOURCES = sort(unlist(lapply(
	c("ext/miniply"),
	collect_sources,
	pattern = "\\.cpp$"
)))
EXT_C_SOURCES = sort(list.files(
	"src",
	pattern = "\\.c$",
	recursive = TRUE,
	full.names = FALSE
))

PKG_CPPFLAGS_STR = collapse_flags(append_flags(PKG_CPPFLAGS, OIDN_CPPFLAGS))
PKG_LIBS_STR = collapse_flags(append_flags(
	LINK_LIBS,
	PKG_LIBS_ACC,
	OIDN_LDFLAGS,
	OIDN_LIBS
))
PKG_CXXFLAGS_STR = collapse_flags(PKG_CXXFLAGS)
DEFINES_STR = collapse_flags(DEFINES)

define(
	PKG_CPPFLAGS = PKG_CPPFLAGS_STR,
	PKG_LIBS = PKG_LIBS_STR,
	PKG_CXXFLAGS = PKG_CXXFLAGS_STR,
	DEFINES = DEFINES_STR,
	VISIBILITY_LDFLAGS = "",
	DIR_SOURCES = collapse_flags(DIR_SOURCES),
	SUBDIR_SOURCES = collapse_flags(SUBDIR_SOURCES),
	EXT_CPP_SOURCES = collapse_flags(EXT_CPP_SOURCES),
	EXT_C_SOURCES = collapse_flags(EXT_C_SOURCES)
)

configure_file("src/Makevars.in")
configure_file("src/Makevars.win.in")

message("--------------------------------------------------")
message("Configuration for rayrender")
message(sprintf(
	"  cppflags: %s %s",
	collapse_flags(CPPFLAGS_BASE),
	DEFINES_STR
))
message(sprintf("  cxxflags: %s", PKG_CXXFLAGS_STR))
message(sprintf("  includes: %s", PKG_CPPFLAGS_STR))
message(sprintf("  libs:     %s", PKG_LIBS_STR))
message("--------------------------------------------------")
