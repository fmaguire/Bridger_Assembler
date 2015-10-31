#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
m4_include([boost.m4])
m4_include([ax_boost_base.m4])

#AC_PREREQ([2.65])
AC_INIT([Bridger], [r2013-03-21], [changzmaths@gmail.com])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_AWK
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LIBTOOL
AC_SUBST(LIBTOOL_DEPS)
AX_BOOST_BASE([1.47.0])

# Checks for libraries.
BOOST_REQUIRE([1.47.0])
BOOST_FILESYSTEM
BOOST_SERIALIZATION

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_MODE_T
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T

# Checks for library functions.
AC_FUNC_FORK
AC_CHECK_FUNCS([floor mkdir sqrt strdup])

# set CFLAGS and CXXFLAGS
user_CFLAGS=${CFLAGS}
generic_CFLAGS="-Wall -Wno-deprecated -g -gdwarf-2 -Wunused -Wuninitialized"
ext_CFLAGS=""
debug_CFLAGS=""
#echo "${host_cpu}-${host_os}"
case "${host_cpu}-${host_os}" in
  i*86-*linux*)
            ext_CFLAGS="-march=i686";;
  i*86-darwin*)
        CFLAGS="-m64"
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM], [ext_CFLAGS="-arch x86_64"], []);;
  *)
    AC_MSG_CHECKING([if gcc accepts -m64])
    CFLAGS="-m64"
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM], [ext_CFLAGS="-m64"; AC_MSG_RESULT([yes])],
                                          [ext_CFLAGS="-D_FILE_OFFSET_BITS=64"; AC_MSG_RESULT([no])]);;
esac

AC_ARG_ENABLE(vectorize,    [  --enable-vectorize        Enable GCC auto-vectorization],
                          [ext_CFLAGS="${ext_CFLAGS} -ftree-vectorize -msse3 -ffast-math -ftree-vectorizer-verbose=99"], [])

AC_ARG_ENABLE(intel64,      [  --enable-intel64        optimize for Intel64 CPU such as Xeon and Core2],
                          [ext_CFLAGS="${ext_CFLAGS} -march=nocona"], [])

AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
                              [enable debugging info (default is no)])],
              [], [enable_debug=no])


AC_ARG_ENABLE([optim],
              [AS_HELP_STRING([--enable-optim@<:@=0|1|2|3@:>@],
                              [set optimization level (default is 3)])],
              [if test "x$enable_optim" = xyes; then enable_optim=3; fi],
              [enable_optim=3])

AS_IF([test "x$enable_optim" != xno], [ext_CFLAGS="$ext_CFLAGS -O$enable_optim"])

AS_IF([test "x$enable_debug" = xyes],
      [debug_CFLAGS="-DDEBUG"],
      [debug_CFLAGS="-DNDEBUG"])

AC_ARG_ENABLE(profiling,      [  --enable-profiling        enable profiling with google-perftools],
          [ext_LDFLAGS="-lprofiler -ltcmalloc"], [])

CFLAGS="${generic_CFLAGS} ${ext_CFLAGS} ${user_CFLAGS}"
CXXFLAGS="$CFLAGS"
CXXFLAGS="${CXXFLAGS} ${BOOST_CPPFLAGS}"
LDFLAGS="$ext_LDFLAGS"


AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT

# dump some configuration confirmations
echo \
"
-- ${PACKAGE_STRING} Configuration Results --
  C++ compiler:        ${CXX} ${CXXFLAGS} ${LDFLAGS}"

if test x"${GCC}" = x"yes" ; then
   gcc_version=`${CC} --version | head -n 1`
   echo "  GCC version:         ${gcc_version}"
else
   gcc_version=''
fi

echo \
"  Host System type:    ${host}
  Install prefix:      ${prefix}
  Install eprefix:     ${exec_prefix}

  See config.h for further configuration information.
  Email <${PACKAGE_BUGREPORT}> with questions and bug reports.
"

