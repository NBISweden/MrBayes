#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([MrBayes],[3.2.7-dev])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([bayes.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.

AC_ARG_WITH([mpi],
            [AS_HELP_STRING([--with-mpi],
                [Build with MPI parallelization support @<:@default=check@:>@])
            ],
            [],
            [with_mpi=check])

AX_PROG_CC_MPI([test "x$with_mpi" != "xno"],
               [AC_DEFINE([MPI_ENABLED], [1], [Define if you want to enable MPI])],
       [if test "x$with_mpi" != "xno"; then
           AC_MSG_FAILURE(
               [--with-mpi was given, but test for MPI failed])
        fi
       ])

AC_PROG_INSTALL

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])

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
                [Build with the Beagle library @<:@default=check@:>@])
            ],
            [],
            [with_beagle=check])

if test "x$with_beagle" != "xno"; then
   PKG_CHECK_MODULES([BEAGLE], [hmsbeagle-1],
                     AC_DEFINE([BEAGLE_ENABLED], [1],
                               [Define if you want to enable Beagle]),
                     [if test "x$with_beagle" != "xcheck"; then
                        AC_MSG_FAILURE(
                           [--with-beagle was given, but test for Beagle failed])
                      fi
                     ])
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
