#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([MrBayes],[pre-3.2.7])

AC_CONFIG_AUX_DIR([build])

AM_INIT_AUTOMAKE([foreign])
AM_MAINTAINER_MODE

AC_CONFIG_SRCDIR([src/bayes.c])
AC_CONFIG_HEADERS([src/config.h])

AC_CANONICAL_HOST
AC_DEFINE_UNQUOTED([HOST_TYPE], ["$host"], [Host triplet])
AC_DEFINE_UNQUOTED([HOST_CPU], ["$host_cpu"], [Host CPU type])

# NB: DEFAULT CFLAGS:
#     Using --enable-debug sets CFLAGS to "-g"  ($debug_cflags)
#     Not using --enable-debug sets CFLAGS to "-O3 -DNDEBUG" ($release_cflags)
#     "-std=c99 -pedantic" will always be added ($default_cflags)
#
# (warning flags are added by the AX_CFLAGS_WARN_ALL macro further down)

debug_cflags="-g"
release_cflags="-O3 -DNDEBUG"
default_cflags="-std=c99 -pedantic"

# Optional feature: debug

AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug],
        [Enable debugging @<:@default=no@:>@])],
    [enable_debug=$enableval],
    [enable_debug=no])

# Handle CFLAGS at once, before AX_PROG_CC_MPI further down.

AS_IF([test "x$enable_debug" = "xyes"],
    [our_cflags="$debug_cflags"],
    [our_cflags="$release_cflags"])

# Use our_cflags, the default_cflags,
# plus any CFLAGS given to us by the user.
CFLAGS="$our_cflags $default_cflags $CFLAGS"

# Optional features: SSE, AVX, FMA

AC_ARG_ENABLE([sse],
    [AS_HELP_STRING([--enable-sse],
        [Enable use of code written using SSE (Streaming SIMD Extensions)
         @<:@default=yes@:>@ (if supported)])],
    [enable_sse=$enableval],
    [enable_sse=yes])

AC_ARG_ENABLE([avx],
    [AS_HELP_STRING([--enable-avx],
        [Enable use of code written using AVX (Advanced Vector Extensions)
         @<:@default=yes@:>@ (if supported)])],
    [enable_avx=$enableval],
    [enable_avx=yes])

AC_ARG_ENABLE([fma],
    [AS_HELP_STRING([--enable-fma],
        [Enable use of code written using FMA (Fused multiply-add SIMD Extensions)
         @<:@default=yes@:>@ (if supported)])],
    [enable_fma=$enableval],
    [enable_fma=yes])

# Optional external libraries:  MPI, Readline, Beagle

AC_ARG_WITH([mpi],
    [AS_HELP_STRING([--with-mpi],
        [Build with MPI parallelization support @<:@default=check@:>@])],
    [],
    [with_mpi=no])

AC_ARG_WITH([readline],
    [AS_HELP_STRING([--with-readline],
        [Build with the Readline library @<:@default=check@:>@])],
    [],
    [with_readline=check])

AC_ARG_WITH([beagle],
    [AS_HELP_STRING([--with-beagle[[=PATH]]],
        [Build with the Beagle library @<:@default=check@:>@.
        The PATH argument should be the prefix path where Beagle
        was installed, e.g. @S|@HOME/local, /opt/local etc.
        The path /usr/local is checked if no path is specified.])],
    [],
    [with_beagle=check])

# Checks for programs.

AC_PROG_INSTALL

# Determine if we can (and want) to use an MPI C compiler.

AX_PROG_CC_MPI([test "x$with_mpi" != "xno"],
    [use_mpi=yes
     AC_DEFINE([MPI_ENABLED], [1], [Define if you want to enable MPI])],
    [use_mpi=no
     AS_IF([test "x$with_mpi" = "xyes"],
        [AC_MSG_ERROR([--with-mpi was given, but test for MPI failed])
    ]) dnl AS_IF
]) dnl AX_PROG_CC_MPI

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])

# Check for the Beagle library.

# Strategy:  Use AC_CHECK_LIB + AC_CHECK_HEADER to look for the library
# and headers with the given prefix, or in /usr/local if no prefix path
# was given, and if that fails use PKG_CHECK_MODULES (and trust its
# result).

AS_IF([test "x$with_beagle" != "xno"], [
    old_cppflags="$CPPFLAGS"
    old_ldflags="$LDFLAGS"

    beagle_subinc="include/libhmsbeagle-1"
    beagle_sublib="lib"

    AS_IF([test -d "$with_beagle"],
          [beagle_prefix="$with_beagle"],
          [beagle_prefix="/usr/local"])

    AC_MSG_NOTICE([Trying to locate Beagle in '$beagle_prefix'])

    CPPFLAGS="-I$beagle_prefix/$beagle_subinc"
    LDFLAGS="-L$beagle_prefix/$beagle_sublib -lhmsbeagle"

    AC_CHECK_HEADER([libhmsbeagle/beagle.h],
        [AC_SUBST([BEAGLE_CFLAGS], [$CPPFLAGS])], [])
    AC_CHECK_LIB(hmsbeagle, beagleCreateInstance,
        [AC_SUBST([BEAGLE_LIBS], [$LDFLAGS])], [])

    CPPFLAGS="$old_cppflags"
    LDFLAGS="$old_ldflags"

    AS_IF([test -n "$BEAGLE_CFLAGS" -a -n "$BEAGLE_LIBS"],
        [AC_MSG_NOTICE([Found Beagle in '$beagle_prefix'])
         beagle_found="yes"],
        [AC_MSG_WARN([Failed to locate Beagle, falling back on pkg-config])
         PKG_CHECK_MODULES([BEAGLE], [hmsbeagle-1],
            [AC_MSG_NOTICE([Found Beagle with pkg-config])
             beagle_found="yes"],
            [beagle_found="no"])
    ]) dnl AS_IF

    AS_IF([test "x$beagle_found" != "xyes"],
        [AS_IF([test "x$with_beagle" != "xcheck"],
            [AC_MSG_ERROR(
                [--with-beagle was given, but test for Beagle failed])]
        )], dnl AS_IF
        [AC_MSG_NOTICE([Using BEAGLE_CFLAGS="$BEAGLE_CFLAGS"])
         AC_MSG_NOTICE([Using BEAGLE_LIBS="$BEAGLE_LIBS"])
         AC_DEFINE([BEAGLE_ENABLED], [1],
            [Define if you want to enable Beagle])])
]) dnl AS_IF

# Check for readline, but only if we don't have MPI, because aparantly
# it won't work.

AS_IF([test "x$use_mpi" = "xno"],
    [AS_IF([test "x$with_readline" != "xno"], [AX_LIB_READLINE])],
    [AS_IF([test "x$with_readline" = "xyes"],
        [AC_MSG_WARN(
            [--with-readline was given, but MPI support requires
            readline to be disabled.  Retry with --without-mpi if you
            *really* need readline support])
    ]) dnl AS_IF
]) dnl AS_IF

# Checks for header files.
AC_CHECK_HEADERS([float.h limits.h memory.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_CHECK_HEADER_STDBOOL

# Get compiler's vendor and version
AX_COMPILER_VENDOR
AX_COMPILER_VERSION

AC_DEFINE_UNQUOTED([COMPILER_VENDOR],
    ["$ax_cv_c_compiler_vendor"], [Compiler vendor])
AC_DEFINE_UNQUOTED([COMPILER_VERSION],
    ["$ax_cv_c_compiler_version"], [Compiler version])

# Add warning flags to compile command line.
AX_CFLAGS_WARN_ALL

# Check for SIMD (etc.) extensions supported by hardware and compiler
# (-msse2 etc.)
AX_EXT

AS_IF([test "x$enable_sse" = "xno"],
    [AC_DEFINE([DISABLE_SSE], [1],
        [SSE code path should not be compiled.])])
AS_IF([test "x$enable_avx" = "xno"],
    [AC_DEFINE([DISABLE_AVX], [1],
        [AVX code path should not be compiled.])])
AS_IF([test "x$enable_fma" = "xno"],
    [AC_DEFINE([DISABLE_FMA], [1],
        [FMA code path should not be compiled.])])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor memset pow sqrt strchr])

AC_CONFIG_FILES([Makefile doc/Makefile doc/tutorial/Makefile examples/Makefile src/Makefile])
AC_OUTPUT
