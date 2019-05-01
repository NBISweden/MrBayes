# SYNOPSIS
#
# AX_RESET_HEADERS_CACHE(headers ...)
#
# DESCRIPTION
#
# This macro invalidates the headers cache variables created by previous AC_CHECK_HEADER/AC_CHECK_HEADERS checks.
#
AC_DEFUN([AX_RESET_HEADERS_CACHE], [
    AS_FOR([AX_var], [ax_var], [$1], [
        dnl You can replace "ac_cv_header_" with any prefix from http://www.gnu.org/software/autoconf/manual/html_node/Cache-Variable-Index.html
        AS_VAR_PUSHDEF([ax_Var], [ac_cv_header_${ax_var}])
        AS_UNSET([ax_Var])
        AS_VAR_POPDEF([ax_Var])
    ])
]) # AX_RESET_HEADERS_CACHE
