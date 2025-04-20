AC_INIT([rayrender])

# Find the compiler and compiler flags used by R.
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`

if test `uname` = "Darwin" ; then
  darwin="yes"
  cmd=`echo $CC $CFLAGS | grep -E 'x86_64|ppc64|-m64'`
  if test -n "$cmd"; then
    have_64bit="yes"
  else
    have_64bit="no"
  fi 
else
  darwin="no"
fi

CXX=`"${R_HOME}/bin/R" CMD config CXX20`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXX20FLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
DEFINES="-D RAY_REPRODUCE_PERLIN -DSTRICT_R_HEADERS"

AC_LANG(C++)
AC_PROG_CPP

PKG_CXXFLAGS=""

AC_MSG_CHECKING([for hidden‑visibility support])
PKG_CXXFLAGS="${PKG_CXXFLAGS} -fvisibility=hidden -fvisibility-inlines-hidden"
CFLAGS="${CFLAGS}  -fvisibility=hidden"
AC_MSG_RESULT([enabled])

dnl Detect the target architecture
AC_MSG_CHECKING([target architecture])
TARGET_ARCH=`${R_HOME}/bin/Rscript -e 'cat(Sys.info()[["machine"]])'`
AC_MSG_RESULT([${TARGET_ARCH}])


if test "$OS" = "Windows_NT"; then
  DEFINES="${DEFINES} -DPBRT_IS_WINDOWS -DNOMINMAX -DPBRT_HAVE__ALIGNED_MALLOC"
elif test "$(uname)" = "Darwin"; then
  DEFINES="${DEFINES} -DPBRT_IS_OSX -DPBRT_HAVE_POSIX_MEMALIGN"
elif test "$(uname)" = "Linux"; then
  DEFINES="${DEFINES} -DPBRT_IS_LINUX -DPBRT_HAVE_POSIX_MEMALIGN"
fi

AC_PATH_X
if test x$no_x = xyes ; then
  AC_MSG_WARN([X11 not found, compiling without render preview capability.])
else
  if test -n "${x_includes}"; then
    PKG_CPPFLAGS="${PKG_CPPFLAGS} -I${x_includes}"
  fi
  if test -n "${x_libraries}"; then
    LIBS="${LIBS} -L${x_libraries} -lX11"
  else
    LIBS="${LIBS} -lX11"
  fi
  DEFINES="${DEFINES} -DRAY_HAS_X11"
  if test $darwin = yes; then
    PKG_CPPFLAGS="${PKG_CPPFLAGS} -I/opt/X11/include"
  fi
fi

AC_MSG_CHECKING([for Open Image Denoise (OIDN)])

OIDN_CPPFLAGS=""
OIDN_LDFLAGS=""
OIDN_LIBS=""

if test -d "$OIDN_PATH"; then
  AC_MSG_RESULT([found at $OIDN_PATH])
  OIDN_CPPFLAGS="-I'$OIDN_PATH/include'"
  OIDN_LDFLAGS="-L$OIDN_PATH/lib"
  OIDN_LIBS="-lOpenImageDenoise -ltbb"
  DEFINES="${DEFINES} -DHAS_OIDN"
  OIDN_LDFLAGS="${OIDN_LDFLAGS} -Wl,-rpath,$OIDN_PATH/lib"
else
  AC_MSG_WARN([OIDN not found. Skipping OIDN support.])
  OIDN_CPPFLAGS=""
  OIDN_LDFLAGS=""
  OIDN_LIBS=""
fi

AC_SUBST([OIDN_CPPFLAGS])
AC_SUBST([OIDN_LDFLAGS])
AC_SUBST([OIDN_LIBS])

########################################################
# Locate + link to OpenEXR (static)
########################################################
AC_MSG_CHECKING([OpenEXR package])

# Directories from the R package 'libopenexr'
OPENEXR_LIB_DIR=`${R_HOME}/bin/Rscript -e 'cat(system.file("lib", package="libopenexr", mustWork=TRUE))'`
OPENEXR_INC_DIR=`${R_HOME}/bin/Rscript -e 'cat(system.file("include", "OpenEXR", package="libopenexr", mustWork=TRUE))'`

# Detect architecture subfolder (e.g. "x86_64", "arm64", etc.)
OPENEXR_LIB_ARCH="${OPENEXR_LIB_DIR}/${TARGET_ARCH}"

if test -d "${OPENEXR_LIB_ARCH}"; then
  AC_MSG_RESULT([found at ${OPENEXR_LIB_ARCH}])
  # Add -I for headers
  PKG_CPPFLAGS="${PKG_CPPFLAGS} -I${OPENEXR_INC_DIR}"
  # Add -L (the library folder) plus -l (the library name)
  PKG_LIBS="${PKG_LIBS} -L${OPENEXR_LIB_ARCH}  -lIlmThread-3_4 -lOpenEXR-3_4 -lIex-3_4 -lOpenEXRCore-3_4 -lOpenEXRUtil-3_4"
else
  AC_MSG_RESULT([not found])
  echo ${OPENEXR_LIB_ARCH}
  AC_MSG_ERROR([OpenEXR not found; please install libopenexr or make sure it is on the library path])
fi


########################################################
# Locate + link to Imath (static)
########################################################
AC_MSG_CHECKING([Imath package])

IMATH_LIB_DIR=`${R_HOME}/bin/Rscript -e 'cat(system.file("lib", package="libimath", mustWork=TRUE))'`
IMATH_INC_DIR=`${R_HOME}/bin/Rscript -e 'cat(system.file("include", "Imath", package="libimath", mustWork=TRUE))'`
IMATH_LIB_ARCH="${IMATH_LIB_DIR}/${TARGET_ARCH}"

if test -d "${IMATH_LIB_ARCH}"; then
  AC_MSG_RESULT([found at ${IMATH_LIB_ARCH}])
  PKG_CPPFLAGS="${PKG_CPPFLAGS} -I${IMATH_INC_DIR}"
  PKG_LIBS="${PKG_LIBS} -L${IMATH_LIB_ARCH} -lImath-3_2"
else
  AC_MSG_RESULT([not found])
  AC_MSG_ERROR([Imath not found; please install libimath])
fi

########################################################
# Locate + link to Ptex (static)
########################################################
AC_MSG_CHECKING([Ptex package])

PTEX_LIB_DIR=`${R_HOME}/bin/Rscript -e 'cat(system.file("lib", package="libptex", mustWork=TRUE))'`
PTEX_INC_DIR=`${R_HOME}/bin/Rscript -e 'cat(system.file("include", package="libptex", mustWork=TRUE))'`
PTEX_LIB_ARCH="${PTEX_LIB_DIR}/${TARGET_ARCH}/ptex"

if test -d "${PTEX_LIB_ARCH}"; then
  AC_MSG_RESULT([found at ${PTEX_LIB_ARCH}])
  PKG_CPPFLAGS="${PKG_CPPFLAGS} -I${PTEX_INC_DIR}"
  PKG_LIBS="${PKG_LIBS} -L${PTEX_LIB_ARCH} -lPtex"
else
  AC_MSG_RESULT([not found])
  AC_MSG_ERROR([Ptex not found; please install libptex])
fi

########################################################
# Locate + link to libdeflate (static)
########################################################
AC_MSG_CHECKING([libdeflate package])

LIBDEFLATE_LIB_DIR=`${R_HOME}/bin/Rscript -e 'cat(system.file("lib", package="libdeflate", mustWork=TRUE))'`
LIBDEFLATE_INC_DIR=`${R_HOME}/bin/Rscript -e 'cat(system.file("include", package="libdeflate", mustWork=TRUE))'`
LIBDEFLATE_LIB_ARCH="${LIBDEFLATE_LIB_DIR}/${TARGET_ARCH}"

if test -d "${LIBDEFLATE_LIB_ARCH}"; then
  AC_MSG_RESULT([found at ${LIBDEFLATE_LIB_ARCH}])
  PKG_CPPFLAGS="${PKG_CPPFLAGS} -I${LIBDEFLATE_INC_DIR}"
  PKG_LIBS="${PKG_LIBS} -L${LIBDEFLATE_LIB_ARCH} -ldeflate"
else
  AC_MSG_RESULT([not found])
  AC_MSG_ERROR([libdeflate not found; please install libdeflate])
fi

########################################################
# Locate + link (header-only) NanoVDB
########################################################
AC_MSG_CHECKING([NanoVDB headers])

NANOVDB_INC_DIR=`${R_HOME}/bin/Rscript -e 'cat(system.file("include", package="libnanovdb", mustWork=TRUE))'`
NANOVDB_LIB_DIR=`${R_HOME}/bin/Rscript -e 'cat(system.file("include", package="libnanovdb", mustWork=TRUE))'`

if test -d "${NANOVDB_INC_DIR}"; then
  AC_MSG_RESULT([found at ${NANOVDB_INC_DIR}])
  # Headers only => just need -I
  PKG_CPPFLAGS="${PKG_CPPFLAGS} -I${NANOVDB_INC_DIR}"
  PKG_LIBS="${PKG_LIBS} -L${NANOVDB_LIB_DIR}"
else
  AC_MSG_RESULT([not found])
  AC_MSG_WARN([NanoVDB headers not found; please install libnanovdb or remove LinkingTo: libnanovdb if not needed.])
fi

# Check for SSE/NEON support
AC_MSG_CHECKING([for SSE4.1 support])
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([
    #include <smmintrin.h>
  ], [
    __m128 test = _mm_dp_ps(_mm_set1_ps(1.0f), _mm_set1_ps(2.0f), 0xFF);
  ])],
  [AC_MSG_RESULT([yes])
   DEFINES="${DEFINES} -DHAS_SSE -DHAS_SSE41 -DRAYSIMD -DRAYSIMDVECOFF"
   PKG_CXXFLAGS="${PKG_CXXFLAGS} -msse4.1"
  ],
  [AC_MSG_RESULT([no])
   
   # If SSE4.1 not found, check SSE3
   AC_MSG_CHECKING([for SSE3 support])
   AC_COMPILE_IFELSE(
     [AC_LANG_PROGRAM([
       #include <pmmintrin.h>
     ], [
       __m128 test = _mm_hadd_ps(_mm_set1_ps(1.0f), _mm_set1_ps(1.0f));
     ])],
     [AC_MSG_RESULT([yes])
      DEFINES="${DEFINES} -DHAS_SSE -DHAS_SSE3 -DRAYSIMD -DRAYSIMDVECOFF"
      PKG_CXXFLAGS="${PKG_CXXFLAGS} -msse3"
     ],
     [AC_MSG_RESULT([no])
      
      # If SSE3 not found, fall back to SSE2
      AC_MSG_CHECKING([for SSE2 support])
      AC_COMPILE_IFELSE(
        [AC_LANG_PROGRAM([
          #include <emmintrin.h>
        ], [
          __m128d test = _mm_setzero_pd();
        ])],
        [AC_MSG_RESULT([yes])
         DEFINES="${DEFINES} -DHAS_SSE -DHAS_SSE2 -DRAYSIMD -DRAYSIMDVECOFF"
         PKG_CXXFLAGS="${PKG_CXXFLAGS} -msse2"
        ],
        [AC_MSG_RESULT([no])
         # If SSE2 not found, try SSE or nothing
         AC_MSG_CHECKING([for SSE support])
         AC_COMPILE_IFELSE(
           [AC_LANG_PROGRAM([
             #include <xmmintrin.h>
           ], [
             __m128 test = _mm_setzero_ps();
           ])],
           [AC_MSG_RESULT([yes])
            DEFINES="${DEFINES} -DHAS_SSE -DRAYSIMD -DRAYSIMDVECOFF"
            PKG_CXXFLAGS="${PKG_CXXFLAGS} -msse"
           ],
           [AC_MSG_RESULT([no])
            # SSE not found at all; no vector intrinsics
           ])
        ])
     ])
  ])

AC_MSG_CHECKING([for NEON support])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#endif
]], [[float32x4_t test = vdupq_n_f32(0.0);]])],
  [AC_MSG_RESULT([yes])
   DEFINES="${DEFINES} -DHAS_NEON -DRAYSIMD -DRAYSIMDVECOFF"
   if test "$darwin" = "yes"; then
     PKG_CXXFLAGS="${PKG_CXXFLAGS} -march=armv8-a+simd"
   else
     PKG_CXXFLAGS="${PKG_CXXFLAGS} -mfpu=neon"
   fi],
  [AC_MSG_RESULT([no])])

AC_MSG_CHECKING([If test debugging should be turned on:])
AC_MSG_NOTICE([RAY_COLOR_DEBUG is "$RAY_COLOR_DEBUG"])

if test "$RAY_COLOR_DEBUG" = "true"; then
  AC_MSG_RESULT([Yes, enabling debug])
  DEFINES="${DEFINES} -DRAY_COLOR_DEBUG"
else
  AC_MSG_RESULT([No, not enabling debug])
fi

SUBDIR_SOURCES="$(cd src/ && ls {core,hitables,materials,math,utils}/*.cpp | tr '\n' ' ')"
EXT_CPP_SOURCES="$(cd src/ && ls {ext/rply,ext/lodepng,ext/sky,ext/miniply}/*.cpp | tr '\n' ' ')"

# EXT_CPP_SOURCES="$(cd src/ && ls {ext/rply,ext/lodepng,ext/OpenEXR,ext/ptex,ext/sky}/*.cpp | tr '\n' ' ')"
# PBRT_SOURCES="pbrt/util/colorspace.cpp pbrt/util/spectrum.cpp pbrt/util/file.cpp pbrt/util/error.cpp pbrt/util/parallel.cpp pbrt/util/color.cpp pbrt/util/stats.cpp pbrt/util/log.cpp pbrt/util/string.cpp pbrt/util/pstd.cpp pbrt/util/mesh.cpp pbrt/pbrt.cpp pbrt/shapes.cpp pbrt/scene.cpp pbrt/cmd/pbrt.cpp"
PBRT_SOURCES="$(cd src/ && ls {ext/pbrt,ext/pbrt/base,ext/pbrt/cpu,ext/pbrt/util,ext/pbrt/spectrums}/*.cpp | tr '\n' ' ')"

DIR_SOURCES="$(cd src/ && ls *.cpp | tr '\n' ' ')"

# Supply it as a variable
AC_SUBST(SUBDIR_SOURCES)
AC_SUBST(EXT_C_SOURCES)
AC_SUBST(EXT_CPP_SOURCES)
AC_SUBST(PBRT_SOURCES)
AC_SUBST(DIR_SOURCES)
AC_SUBST(PKG_CXXFLAGS)


AC_SUBST([PKG_CPPFLAGS], ["${PKG_CPPFLAGS} ${OIDN_CPPFLAGS}"])
AC_SUBST([PKG_LIBS], ["${LIBS} ${PKG_LIBS} ${OIDN_LDFLAGS} ${OIDN_LIBS} "])
AC_SUBST([DEFINES], ["${DEFINES}"])

AC_CONFIG_FILES([
  src/Makevars
  src/Makevars.win
])
AC_OUTPUT

echo "
  --------------------------------------------------
  Configuration for ${PACKAGE_NAME} ${PACKAGE_VERSION}

    cppflags: ${CPPFLAGS} ${DEFINES}
    libs:     ${PKG_LIBS}
    includes: ${PKG_CPPFLAGS}
    cxxflags: ${PKG_CXXFLAGS}

  --------------------------------------------------
"
