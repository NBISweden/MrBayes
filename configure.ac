#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([MrBayes],[3.2.7a])

AC_CONFIG_AUX_DIR([autotools])
AC_CONFIG_MACRO_DIRS([autotools])

AC_REQUIRE_AUX_FILE([tap-driver.sh])

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

# Optional feature: Enable/disable building of documentation under doc/
# (enabled by default)

AC_ARG_ENABLE([doc],
    [AS_HELP_STRING([--enable-doc],
        [Enable building supplied documentation from sources in doc/
         @<:@default=yes@:>@])],
    [enable_doc=$enableval],
    [enable_doc=yes])

# Optional external libraries:  MPI, Readline, Beagle

AC_ARG_WITH([mpi],
    [AS_HELP_STRING([--with-mpi],
        [Build with MPI parallelization support @<:@default=no@:>@])],
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
        The paths /usr/local and /usr are checked if no path is specified.])],
    [],
    [with_beagle=check])

# Checks for programs.

AC_PROG_INSTALL
PKG_PROG_PKG_CONFIG

AM_CONDITIONAL([BUILD_DOC], [test "$enable_doc" = "yes"])
AM_COND_IF([BUILD_DOC], [

    # Test whether we can rebuild the documentation.
    # This does not affect the pre-built PDF documents unless these are deleted.
    AC_PATH_PROG(PDFLATEX, [pdflatex])
    AC_PATH_PROG(BIBTEX, [bibtex])
    AC_PATH_PROG(LATEXMK, [latexmk])

]) dnl AM_COND_IF([BUILD_DOC])

AM_CONDITIONAL([HAVE_LATEX], [test -n "$PDFLATEX" && test -n "$BIBTEX"])
AM_COND_IF([BUILD_DOC], [
    AM_COND_IF([HAVE_LATEX], [],
               [AC_MSG_WARN([pdflatex and/or bibtex is missing, documentation will not be rebuilt])])
])
AM_CONDITIONAL([HAVE_LATEXMK], [test -n "$LATEXMK"])

# Determine if we can (and want) to use an MPI C compiler.
AX_PROG_CC_MPI([test "x$with_mpi" != "xno"],
    [use_mpi=yes
     AC_DEFINE([MPI_ENABLED], [1], [Define if you want to enable MPI])],
    [use_mpi=no
     AS_IF([test "x$with_mpi" = "xyes"],
        [AC_MSG_ERROR([--with-mpi was given, but test for MPI failed])
    ]) dnl AS_IF
]) dnl AX_PROG_CC_MPI
TESTS_USE_MPI=$use_mpi
AC_SUBST([TESTS_USE_MPI])

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])

# Check for the Beagle library.

# Strategy:  Use AC_CHECK_LIB + AC_CHECK_HEADER to look for the library
# and headers with the given prefix, or in /usr/local followed by /usr
# if no prefix path was given, and if that fails use PKG_CHECK_MODULES
# (and trust its result).

# Note: We assume that Beagle's top-level include directory is located under
# include/libhmsbeagle-1 and nowhere else.

AS_IF([test "x$with_beagle" != "xno"], [
    old_cppflags="$CPPFLAGS"
    old_ldflags="$LDFLAGS"
    old_ldlibs="$LDLIBS"

    beagle_pkg_name="hmsbeagle-1"
    beagle_subinc="include/lib$beagle_pkg_name"
    beagle_sublib="lib"

    AS_IF([test -n "$BEAGLE_CFLAGS" || test -n "$BEAGLE_LIBS"],
        [AC_MSG_NOTICE([Using the supplied values of BEAGLE_CFLAGS and BEAGLE_LIBS])
         beagle_found="yes"],
        [AS_IF([test -d "$with_beagle"],
              [beagle_prefix_list="$with_beagle"],
              [beagle_prefix_list="/usr/local /usr"])])

    for beagle_prefix in $beagle_prefix_list; do

        AC_MSG_NOTICE([Trying to locate Beagle ($beagle_pkg_name) in '$beagle_prefix'])

        CPPFLAGS="-I$beagle_prefix/$beagle_subinc"
        LDFLAGS="-L$beagle_prefix/$beagle_sublib -Wl,-rpath,$beagle_prefix/$beagle_sublib"
        LDLIBS="-lhmsbeagle"

        AC_CHECK_HEADER([libhmsbeagle/beagle.h],
            [AC_SUBST([BEAGLE_CFLAGS], [$CPPFLAGS])], [])
        AC_CHECK_LIB(hmsbeagle, beagleCreateInstance,
            [AC_SUBST([BEAGLE_LIBS], ["$LDFLAGS $LDLIBS"])], [])

        CPPFLAGS="$old_cppflags"
        LDFLAGS="$old_ldflags"
        LDLIBS="$old_ldlibs"

        AS_IF([test -n "$BEAGLE_CFLAGS" && test -n "$BEAGLE_LIBS"],
            [AC_MSG_NOTICE([Found Beagle in '$beagle_prefix'])
             beagle_found="yes"])

        AS_IF([test "x$beagle_found" = "xyes"], [break],
              [AS_UNSET(BEAGLE_LIBS)
               AS_UNSET(BEAGLE_CFLAGS)
               AX_RESET_HEADERS_CACHE([libhmsbeagle_beagle_h])])
    done

    AS_IF([test "x$beagle_found" != "xyes"], [
        AC_MSG_WARN([Failed to locate Beagle, falling back on pkg-config])

             PKG_CHECK_MODULES([BEAGLE], [$beagle_pkg_name],
                [AC_MSG_NOTICE([Found Beagle with pkg-config])
                 beagle_found="yes"
                 # Augument the BEAGLE_LIBS value with -Wl,-rpath,...
                BEAGLE_LIBS="-Wl,-rpath,$(pkg-config --libs-only-L $beagle_pkg_name | sed 's/-L//') $BEAGLE_LIBS"],
                [beagle_found="no"])
    ]) dnl AS_IF

    AS_IF([test "x$beagle_found" != "xyes"],
        [AS_IF([test "x$with_beagle" != "xcheck"],
            [AC_MSG_ERROR(
                [--with-beagle was given, but test for Beagle failed])]
        )], dnl AS_IF - ELSE
        [
            AC_DEFINE([BEAGLE_ENABLED], [1],
                [Define if you want to enable Beagle])
            # Test whether we should enable Beagle3 support.
            old_ldflags=$LDFLAGS
            LDFLAGS="$BEAGLE_LIBS"
            AC_CHECK_LIB(hmsbeagle, beagleSetCPUThreadCount,
                [AC_DEFINE([BEAGLE_V3_ENABLED], [1],
                    [Define to enable Beagle v3.1 extensions])
                BEAGLE_CFLAGS="$BEAGLE_CFLAGS -pthread"])
            LDFLAGS=$old_ldflags
            AS_UNSET(old_ldflags)
            AC_MSG_NOTICE([Using BEAGLE_CFLAGS="$BEAGLE_CFLAGS"])
            AC_MSG_NOTICE([Using BEAGLE_LIBS="$BEAGLE_LIBS"])
        ])
]) dnl AS_IF

# Check for readline, but only if we don't have MPI, because apparently
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
AC_CHECK_TYPE([long long])
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

# Check for SIMD (etc.) extensions supported by hardware and compiler.
AX_EXT

# Act on the above tests (which sets the "ax_cv_support_*_ext"
# variables), and add appropriate flags to CFLAGS, or define
# DISABLE_SSE, DISABLE_AVX, DISABLE_FMA if the extensions are not
# supported or disabled by the user.
AS_IF([test "x$enable_sse" = "xno" || test "x$ax_cv_support_sse2_ext" != "xyes"],
    [AC_DEFINE([DISABLE_SSE], [1], [SSE code path should not be compiled.])],
    [CFLAGS="$CFLAGS -msse"])
AS_IF([test "x$enable_avx" = "xno" || test "x$ax_cv_support_avx_ext" != "xyes"],
    [AC_DEFINE([DISABLE_AVX], [1], [AVX code path should not be compiled.])],
    [CFLAGS="$CFLAGS -mavx"])
AS_IF([test "x$enable_fma" = "xno" || test "x$ax_cv_support_fma3_ext" != "xyes"],
    [AC_DEFINE([DISABLE_FMA], [1], [FMA code path should not be compiled.])],
    [CFLAGS="$CFLAGS -mfma"])

# Checks for library functions.
AC_CHECK_FUNCS([floor memset pow sqrt strchr])

AC_CONFIG_FILES([Makefile
                 doc/Makefile
                 doc/SMME_paper/Makefile
                 doc/manual/Makefile
                 examples/Makefile
                 src/Makefile
                 testing/Makefile
])
AC_CONFIG_FILES([testing/runtests.sh], [chmod +x testing/runtests.sh])
AC_OUTPUT
