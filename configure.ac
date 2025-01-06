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
   CXXFLAGS="${CXXFLAGS} -msse4.1"
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
      CXXFLAGS="${CXXFLAGS} -msse3"
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
         CXXFLAGS="${CXXFLAGS} -msse2"
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
            CXXFLAGS="${CXXFLAGS} -msse"
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
     CXXFLAGS="${CXXFLAGS} -march=armv8-a+simd"
   else
     CXXFLAGS="${CXXFLAGS} -mfpu=neon"
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

AC_SUBST([PKG_CPPFLAGS], ["${PKG_CPPFLAGS} ${OIDN_CPPFLAGS}"])
AC_SUBST([PKG_LIBS], ["${LIBS} ${PKG_LIBS} ${OIDN_LDFLAGS} ${OIDN_LIBS}"])
AC_SUBST([DEFINES], ["${DEFINES}"])

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

echo "
  --------------------------------------------------
  Configuration for ${PACKAGE_NAME} ${PACKAGE_VERSION}

    cppflags: ${CPPFLAGS} ${DEFINES}
    libs:     ${PKG_LIBS}
    includes: ${PKG_CPPFLAGS}
    cxxflags: ${CXXFLAGS}

  --------------------------------------------------
"
