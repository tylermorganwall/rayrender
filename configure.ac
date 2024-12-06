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

# Check for SSE2 support first
AC_MSG_CHECKING([for SSE2 support])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([#include <emmintrin.h>], [__m128d test = _mm_setzero_pd();])],
  [AC_MSG_RESULT([yes])
   DEFINES="${DEFINES} -DHAS_SSE -DRAYSIMD -DRAYSIMDVEC"
   CXXFLAGS="${CXXFLAGS} -msse2"],
  [AC_MSG_RESULT([no])
   # If SSE2 not found, try SSE
   AC_MSG_CHECKING([for SSE support])
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([#include <xmmintrin.h>], [__m128 test = _mm_setzero_ps();])],
     [AC_MSG_RESULT([yes])
      DEFINES="${DEFINES} -DHAS_SSE -DRAYSIMD -DRAYSIMDVEC"
      CXXFLAGS="${CXXFLAGS} -msse"],
     [AC_MSG_RESULT([no])])])

AC_MSG_CHECKING([for NEON support])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#endif
]], [[float32x4_t test = vdupq_n_f32(0.0);]])],
  [AC_MSG_RESULT([yes])
   DEFINES="${DEFINES} -DHAS_NEON -DRAYSIMD -DRAYSIMDVEC"
   CXXFLAGS="${CXXFLAGS} -mfpu=neon"],
  [AC_MSG_RESULT([no])])

AC_SUBST([PKG_CPPFLAGS], ["${PKG_CPPFLAGS}"])
AC_SUBST([PKG_LIBS], ["${LIBS} ${PKG_LIBS}"])
AC_SUBST([DEFINES], ["${DEFINES}"])

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

echo "
  --------------------------------------------------
  Configuration for ${PACKAGE_NAME} ${PACKAGE_VERSION}

    cppflags: ${CPPFLAGS} ${DEFINES}
    libs:     ${PKG_LIBS}
    includes: ${PKG_CPPFLAGS}

  --------------------------------------------------
"
