#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([MrBayes],[3.2.7-dev])

AC_CONFIG_AUX_DIR([build])

AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([bayes.c])
AC_CONFIG_HEADERS([config.h])

AC_CANONICAL_HOST
AC_DEFINE_UNQUOTED([HOST_TYPE], ["$host"], [Host triplet])
AC_DEFINE_UNQUOTED([HOST_CPU], ["$host_cpu"], [Host CPU type])

# Checks for programs.

# Determine if we can (and want) to use an MPI C compiler.

AC_ARG_WITH([mpi],
    [AS_HELP_STRING([--with-mpi],
        [Build with MPI parallelization support @<:@default=check@:>@])],
    [],
    [with_mpi=check])

AX_PROG_CC_MPI([test "x$with_mpi" != "xno"],
    [use_mpi=yes
     AC_DEFINE([MPI_ENABLED], [1], [Define if you want to enable MPI])],
    [use_mpi=no
     AS_IF([test "x$with_mpi" = "xyes"],
        [AC_MSG_ERROR([--with-mpi was given, but test for MPI failed])
    ]) dnl AS_IF
]) dnl AX_PROG_CC_MPI

# Get compiler's vendor and version
AX_COMPILER_VENDOR
AX_COMPILER_VERSION

AC_DEFINE_UNQUOTED([COMPILER_VENDOR],
    ["$ax_cv_c_compiler_vendor"], [Compiler vendor])
AC_DEFINE_UNQUOTED([COMPILER_VERSION],
    ["$ax_cv_c_compiler_version"], [Compiler version])

AX_CFLAGS_WARN_ALL

# Check for SIMD extensions supported by hardware and compiler (-msse2 etc.)
AX_EXT

AC_PROG_INSTALL

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])

# Check for readline, but only if we don't have MPI, because aparantly
# it won't work.

AC_ARG_WITH([readline],
    [AS_HELP_STRING([--with-readline],
        [Build with the Readline library @<:@default=check@:>@])],
    [],
    [with_readline=check])

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

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor memset pow sqrt strchr strdup strncasecmp])

# Check for the Beagle library.

AC_ARG_WITH([beagle],
    [AS_HELP_STRING([--with-beagle],
        [Build with the Beagle library @<:@default=check@:>@])],
    [],
    [with_beagle=check])

AS_IF([test "x$with_beagle" != "xno"],
    [PKG_CHECK_MODULES([BEAGLE], [hmsbeagle-1],
        [AC_DEFINE([BEAGLE_ENABLED], [1],
            [Define if you want to enable Beagle])],
        [AS_IF([test "x$with_beagle" != "xcheck"],
            [AC_MSG_ERROR(
                [--with-beagle was given, but test for Beagle failed])
        ]) dnl AS_IF
    ]) dnl PKG_CHECK_MODULES
]) dnl AS_IF

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
