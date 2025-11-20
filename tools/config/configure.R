# Configure rayrender using an R-based workflow.

is_windows = identical(.Platform$OS.type, "windows")
sysname = Sys.info()[["sysname"]]
is_macos = identical(sysname, "Darwin")
is_linux = !is_windows && !is_macos
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

add_rpath = function(flags, dir) {
	flags = append_flags(flags, flag_with_path("-L", dir))
	flags = append_flags(
		flags,
		sprintf(
			"-Wl,-rpath,%s",
			normalizePath(dir, winslash = "/", mustWork = FALSE)
		)
	)
	if (is_linux) {
		flags = append_flags(
			flags,
			sprintf(
				"-Wl,-rpath-link,%s",
				normalizePath(dir, winslash = "/", mustWork = FALSE)
			)
		)
	}
	flags
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

# ---- TBB detection ----

detect_tbb = function(extra_lib_dirs = character()) {
	tbb_root = Sys.getenv("TBB_ROOT", unset = "")
	tbb_lib_env = Sys.getenv("TBB_LIB", unset = "")
	tbb_inc_env = Sys.getenv("TBB_INC", unset = "")

	candidate_lib_dirs = character()

	# 1) Explicit env hints first
	if (nzchar(tbb_lib_env)) {
		candidate_lib_dirs = c(
			candidate_lib_dirs,
			normalizePath(tbb_lib_env, winslash = "/", mustWork = FALSE)
		)
	}
	if (nzchar(tbb_root)) {
		candidate_lib_dirs = c(
			candidate_lib_dirs,
			normalizePath(
				file.path(tbb_root, "lib"),
				winslash = "/",
				mustWork = FALSE
			)
		)
	}

	# 2) Any library dirs provided by the caller (e.g. alongside OIDN)
	if (length(extra_lib_dirs)) {
		candidate_lib_dirs = c(
			candidate_lib_dirs,
			normalizePath(extra_lib_dirs, winslash = "/", mustWork = FALSE)
		)
	}

	# 3) Windows: always add Rtools candidates, regardless of what we already have
	if (is_windows) {
		gcc_path = normalizePath(
			Sys.which("gcc"),
			winslash = "/",
			mustWork = FALSE
		)
		if (nzchar(gcc_path)) {
			gcc_bin_dir = dirname(gcc_path) # .../x86_64-w64-mingw32.static.posix/bin
			triple_dir = dirname(gcc_bin_dir) # .../x86_64-w64-mingw32.static.posix
			rtools_root = dirname(triple_dir) # .../rtools45 (or similar)
			triple_name = basename(triple_dir) # x86_64-w64-mingw32.static.posix

			# Legacy ../../lib from gcc
			legacy_lib = normalizePath(
				file.path(gcc_bin_dir, "..", "..", "lib"),
				winslash = "/",
				mustWork = FALSE
			)

			# Rtools layout: <rtools_root>/<triple>/lib
			rtools_tbb_lib = normalizePath(
				file.path(rtools_root, triple_name, "lib"),
				winslash = "/",
				mustWork = FALSE
			)

			rtools_candidates = c(legacy_lib, rtools_tbb_lib)
			rtools_candidates = rtools_candidates[dir.exists(rtools_candidates)]

			candidate_lib_dirs = c(candidate_lib_dirs, rtools_candidates)

			# If we have no include yet, prefer <rtools_root>/<triple>/include
			if (!nzchar(tbb_inc_env)) {
				rtools_tbb_inc = normalizePath(
					file.path(rtools_root, triple_name, "include"),
					winslash = "/",
					mustWork = FALSE
				)
				if (dir.exists(rtools_tbb_inc)) {
					tbb_inc_env = rtools_tbb_inc
				}
			}
		}
	}

	# 4) Unix autodetect (optional); add, do NOT replace, candidates
	if (
		.Platform$OS.type == "unix" &&
			identical(Sys.getenv("TBB_AUTODETECT", unset = "FALSE"), "TRUE")
	) {
		sysname = Sys.info()[["sysname"]]
		homebrew_prefix = if (identical(sysname, "Darwin")) {
			"/opt/homebrew"
		} else {
			"/usr/local"
		}

		if (identical(sysname, "Darwin")) {
			hb_lib = normalizePath(
				file.path(homebrew_prefix, "opt", "tbb", "lib"),
				winslash = "/",
				mustWork = FALSE
			)
			if (dir.exists(hb_lib)) {
				candidate_lib_dirs = c(candidate_lib_dirs, hb_lib)
			}
			if (!nzchar(tbb_inc_env)) {
				hb_inc = normalizePath(
					file.path(homebrew_prefix, "opt", "tbb", "include"),
					winslash = "/",
					mustWork = FALSE
				)
				if (dir.exists(hb_inc)) {
					tbb_inc_env = hb_inc
				}
			}
		} else {
			tbb_lib_search = Sys.glob(c(
				"/usr/*/libtbb.so",
				"/usr/*/*/libtbb.so",
				"/usr/*/*/*/libtbb.so"
			))
			if (length(tbb_lib_search)) {
				candidate_lib_dirs = c(
					candidate_lib_dirs,
					dirname(tbb_lib_search[[1L]])
				)
			}
			if (!nzchar(tbb_inc_env)) {
				tbb_inc_search = Sys.glob(c(
					"/usr/include/tbb.h",
					"/usr/include/*/tbb.h"
				))
				if (length(tbb_inc_search)) {
					tbb_inc_env = dirname(tbb_inc_search[[1L]])
				}
			}
		}
	}

	candidate_lib_dirs = unique(candidate_lib_dirs)

	# Nothing to scan
	if (!length(candidate_lib_dirs)) {
		message(
			"*** configure: TBB not found; TBB-dependent features will be disabled"
		)
		return(list(found = FALSE, lib_dir = "", inc_dir = "", lib_name = ""))
	}

	# 5) Scan all candidate lib dirs in order; pick the first that actually has a libtbb*
	for (lib_dir in candidate_lib_dirs) {
		if (!dir.exists(lib_dir)) {
			next
		}

		files = list.files(lib_dir, pattern = "^libtbb", full.names = FALSE)
		if (!length(files)) {
			next
		}

		# Match libtbb.a, libtbb12.a, libtbb12_static.a, libtbb12.so, ...
		pattern = "^lib(tbb[0-9]*(?:_static)?)\\.(a|so|dylib|lib|dll)$"
		matches = grep(pattern, files, perl = TRUE, value = TRUE)
		if (!length(matches)) {
			next
		}

		chosen = matches[[1L]]
		lib_name = sub(pattern, "\\1", chosen, perl = TRUE)

		lib_dir_norm = normalizePath(
			lib_dir,
			winslash = "/",
			mustWork = FALSE
		)

		# Derive an include dir:
		# 1) explicit TBB_INC if valid
		inc_dir = ""
		if (nzchar(tbb_inc_env) && dir.exists(tbb_inc_env)) {
			inc_dir = normalizePath(tbb_inc_env, winslash = "/", mustWork = FALSE)
		} else if (nzchar(tbb_root) && dir.exists(file.path(tbb_root, "include"))) {
			# 2) TBB_ROOT/include
			inc_dir = normalizePath(
				file.path(tbb_root, "include"),
				winslash = "/",
				mustWork = FALSE
			)
		} else if (is_windows) {
			# 3) For Rtools-style layout, <triple>/include next to <triple>/lib
			win_inc_candidate = normalizePath(
				file.path(dirname(lib_dir_norm), "include"),
				winslash = "/",
				mustWork = FALSE
			)
			if (dir.exists(win_inc_candidate)) {
				inc_dir = win_inc_candidate
			}
		}

		# 4) Generic ../include fallback relative to lib dir
		if (!nzchar(inc_dir)) {
			generic_inc = normalizePath(
				file.path(lib_dir_norm, "..", "include"),
				winslash = "/",
				mustWork = FALSE
			)
			if (dir.exists(generic_inc)) {
				inc_dir = generic_inc
			}
		}

		message(sprintf(
			"*** configure: found TBB (%s) in %s",
			lib_name,
			lib_dir_norm
		))

		return(list(
			found = TRUE,
			lib_dir = lib_dir_norm,
			inc_dir = inc_dir,
			lib_name = lib_name
		))
	}

	message(
		"*** configure: TBB not found in any candidate locations; TBB-dependent features will be disabled"
	)
	list(found = FALSE, lib_dir = "", inc_dir = "", lib_name = "")
}


oidn_path = Sys.getenv("OIDN_PATH", unset = "")
if (nzchar(oidn_path) && dir.exists(oidn_path)) {
	lib_dir = file.path(oidn_path, "lib")

	# Prefer TBB bundled alongside OIDN if present; otherwise fall back
	tbb_info = detect_tbb(
		extra_lib_dirs = if (dir.exists(lib_dir)) lib_dir else character()
	)

	if (!tbb_info$found) {
		message(
			"*** configure: Open Image Denoise (OIDN) found but TBB is not available; skipping denoiser support"
		)
	} else {
		message(sprintf(
			"*** configure: using Open Image Denoise at %s (TBB %s in %s)",
			oidn_path,
			tbb_info$lib_name,
			tbb_info$lib_dir
		))

		# OIDN headers
		OIDN_CPPFLAGS = append_flags(
			OIDN_CPPFLAGS,
			flag_with_path("-I", file.path(oidn_path, "include"))
		)

		# TBB headers (if we know an include dir)
		if (nzchar(tbb_info$inc_dir)) {
			OIDN_CPPFLAGS = append_flags(
				OIDN_CPPFLAGS,
				flag_with_path("-I", tbb_info$inc_dir)
			)
		}

		# OIDN lib dir + rpath
		if (dir.exists(lib_dir)) {
			OIDN_LDFLAGS = append_flags(OIDN_LDFLAGS, flag_with_path("-L", lib_dir))
			OIDN_LDFLAGS = append_flags(
				OIDN_LDFLAGS,
				sprintf(
					"-Wl,-rpath,%s",
					normalizePath(lib_dir, winslash = "/", mustWork = FALSE)
				)
			)
		}

		# TBB lib dir + rpath (may be distinct from OIDN lib dir)
		if (nzchar(tbb_info$lib_dir) && dir.exists(tbb_info$lib_dir)) {
			OIDN_LDFLAGS = append_flags(
				OIDN_LDFLAGS,
				flag_with_path("-L", tbb_info$lib_dir)
			)
			OIDN_LDFLAGS = append_flags(
				OIDN_LDFLAGS,
				sprintf(
					"-Wl,-rpath,%s",
					normalizePath(tbb_info$lib_dir, winslash = "/", mustWork = FALSE)
				)
			)
		}

		# Link against OIDN + detected TBB name (tbb, tbb12, tbb12_static, ...)
		tbb_lib_flag = paste0(
			"-l",
			if (nzchar(tbb_info$lib_name)) tbb_info$lib_name else "tbb"
		)

		OIDN_LIBS = append_flags(
			OIDN_LIBS,
			"-lOpenImageDenoise",
			tbb_lib_flag
		)

		DEFINES = append_unique_flags(DEFINES, "-DHAS_OIDN")
	}
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
	flag_with_path("-I", (openexr_inc_dir))
)

PKG_LIBS_ACC = add_rpath(PKG_LIBS_ACC, openexr_lib_arch)
PKG_LIBS_ACC = append_flags(
	PKG_LIBS_ACC,
	"-lOpenEXR-3_4",
	"-lOpenEXRUtil-3_4",
	"-lOpenEXRCore-3_4",
	"-lIex-3_4",
	"-lIlmThread-3_4"
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
	flag_with_path("-I", (imath_inc_dir))
)
PKG_LIBS_ACC = add_rpath(PKG_LIBS_ACC, imath_lib_arch)
PKG_LIBS_ACC = append_flags(PKG_LIBS_ACC, "-lImath-3_2")


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
	flag_with_path("-I", (libdeflate_inc_dir))
)
PKG_LIBS_ACC = add_rpath(PKG_LIBS_ACC, libdeflate_lib_arch)
PKG_LIBS_ACC = append_flags(PKG_LIBS_ACC, "-ldeflate")

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

	# Use target_arch to decide which (if any) extra flags are safe:
	# - aarch64 / arm64: NEON is mandatory; use -march=armv8-a+simd
	# - 32-bit ARM: use -mfpu=neon
	# - others: no extra flags (the compile_test wouldn't have succeeded anyway)
	if (grepl("aarch64|arm64", target_arch)) {
		PKG_CXXFLAGS = append_unique_flags(
			PKG_CXXFLAGS,
			"-march=armv8-a+simd"
		)
	} else if (grepl("^arm", target_arch)) {
		PKG_CXXFLAGS = append_unique_flags(
			PKG_CXXFLAGS,
			"-mfpu=neon"
		)
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
