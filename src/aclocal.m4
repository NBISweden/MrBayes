# generated automatically by aclocal 1.15 -*- Autoconf -*-

# Copyright (C) 1996-2014 Free Software Foundation, Inc.

# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

m4_ifndef([AC_CONFIG_MACRO_DIRS], [m4_defun([_AM_CONFIG_MACRO_DIRS], [])m4_defun([AC_CONFIG_MACRO_DIRS], [_AM_CONFIG_MACRO_DIRS($@)])])
m4_ifndef([AC_AUTOCONF_VERSION],
  [m4_copy([m4_PACKAGE_VERSION], [AC_AUTOCONF_VERSION])])dnl
m4_if(m4_defn([AC_AUTOCONF_VERSION]), [2.69],,
[m4_warning([this file was generated for autoconf 2.69.
You have another version of autoconf.  It may work, but is not guaranteed to.
If you have problems, you may need to regenerate the build system entirely.
To do so, use the procedure documented by the package, typically 'autoreconf'.])])

# ===========================================================================
#   http://www.gnu.org/software/autoconf-archive/ax_check_compile_flag.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CHECK_COMPILE_FLAG(FLAG, [ACTION-SUCCESS], [ACTION-FAILURE], [EXTRA-FLAGS], [INPUT])
#
# DESCRIPTION
#
#   Check whether the given FLAG works with the current language's compiler
#   or gives an error.  (Warnings, however, are ignored)
#
#   ACTION-SUCCESS/ACTION-FAILURE are shell commands to execute on
#   success/failure.
#
#   If EXTRA-FLAGS is defined, it is added to the current language's default
#   flags (e.g. CFLAGS) when the check is done.  The check is thus made with
#   the flags: "CFLAGS EXTRA-FLAGS FLAG".  This can for example be used to
#   force the compiler to issue an error when a bad flag is given.
#
#   INPUT gives an alternative input source to AC_COMPILE_IFELSE.
#
#   NOTE: Implementation based on AX_CFLAGS_GCC_OPTION. Please keep this
#   macro in sync with AX_CHECK_{PREPROC,LINK}_FLAG.
#
# LICENSE
#
#   Copyright (c) 2008 Guido U. Draheim <guidod@gmx.de>
#   Copyright (c) 2011 Maarten Bosmans <mkbosmans@gmail.com>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 4

AC_DEFUN([AX_CHECK_COMPILE_FLAG],
[AC_PREREQ(2.64)dnl for _AC_LANG_PREFIX and AS_VAR_IF
AS_VAR_PUSHDEF([CACHEVAR],[ax_cv_check_[]_AC_LANG_ABBREV[]flags_$4_$1])dnl
AC_CACHE_CHECK([whether _AC_LANG compiler accepts $1], CACHEVAR, [
  ax_check_save_flags=$[]_AC_LANG_PREFIX[]FLAGS
  _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $4 $1"
  AC_COMPILE_IFELSE([m4_default([$5],[AC_LANG_PROGRAM()])],
    [AS_VAR_SET(CACHEVAR,[yes])],
    [AS_VAR_SET(CACHEVAR,[no])])
  _AC_LANG_PREFIX[]FLAGS=$ax_check_save_flags])
AS_VAR_IF(CACHEVAR,yes,
  [m4_default([$2], :)],
  [m4_default([$3], :)])
AS_VAR_POPDEF([CACHEVAR])dnl
])dnl AX_CHECK_COMPILE_FLAGS

# ===========================================================================
#    http://www.gnu.org/software/autoconf-archive/ax_compiler_vendor.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_COMPILER_VENDOR
#
# DESCRIPTION
#
#   Determine the vendor of the C/C++ compiler, e.g., gnu, intel, ibm, sun,
#   hp, borland, comeau, dec, cray, kai, lcc, metrowerks, sgi, microsoft,
#   watcom, etc. The vendor is returned in the cache variable
#   $ax_cv_c_compiler_vendor for C and $ax_cv_cxx_compiler_vendor for C++.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2008 Matteo Frigo
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 15

AC_DEFUN([AX_COMPILER_VENDOR],
[AC_CACHE_CHECK([for _AC_LANG compiler vendor], ax_cv_[]_AC_LANG_ABBREV[]_compiler_vendor,
  dnl Please add if possible support to ax_compiler_version.m4
  [# note: don't check for gcc first since some other compilers define __GNUC__
  vendors="intel:     __ICC,__ECC,__INTEL_COMPILER
           ibm:       __xlc__,__xlC__,__IBMC__,__IBMCPP__
           pathscale: __PATHCC__,__PATHSCALE__
           clang:     __clang__
           cray:      _CRAYC
           fujitsu:   __FUJITSU
           gnu:       __GNUC__
           sun:       __SUNPRO_C,__SUNPRO_CC
           hp:        __HP_cc,__HP_aCC
           dec:       __DECC,__DECCXX,__DECC_VER,__DECCXX_VER
           borland:   __BORLANDC__,__CODEGEARC__,__TURBOC__
           comeau:    __COMO__
           kai:       __KCC
           lcc:       __LCC__
           sgi:       __sgi,sgi
           microsoft: _MSC_VER
           metrowerks: __MWERKS__
           watcom:    __WATCOMC__
           portland:  __PGI
	   tcc:       __TINYC__
           unknown:   UNKNOWN"
  for ventest in $vendors; do
    case $ventest in
      *:) vendor=$ventest; continue ;;
      *)  vencpp="defined("`echo $ventest | sed 's/,/) || defined(/g'`")" ;;
    esac
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM(,[
      #if !($vencpp)
        thisisanerror;
      #endif
    ])], [break])
  done
  ax_cv_[]_AC_LANG_ABBREV[]_compiler_vendor=`echo $vendor | cut -d: -f1`
 ])
])

# ===========================================================================
#    http://www.gnu.org/software/autoconf-archive/ax_compiler_version.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_COMPILER_VERSION
#
# DESCRIPTION
#
#   This macro retrieves the compiler version and returns it in the cache
#   variable $ax_cv_c_compiler_version for C and $ax_cv_cxx_compiler_version
#   for C++.
#
#   Version is returned as epoch:major.minor.patchversion
#
#   Epoch is used in order to have an increasing version number in case of
#   marketing change.
#
#   Epoch use: * borland compiler use chronologically 0turboc for turboc
#   era,
#
#     1borlanc BORLANC++ before 5, 2cppbuilder for cppbuilder era,
#     3borlancpp for return of BORLANC++ (after version 5.5),
#     4cppbuilder for cppbuilder with year version,
#     and 5xe for XE era.
#
#   An empty string is returned otherwise.
#
# LICENSE
#
#   Copyright (c) 2014 Bastien ROUCARIES <roucaries.bastien+autoconf@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 4

# for intel
AC_DEFUN([_AX_COMPILER_VERSION_INTEL],
  [ dnl
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_major,
    [__INTEL_COMPILER/100],,
    AC_MSG_FAILURE([[[$0]] unknown intel compiler version]))
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor,
    [(__INTEL_COMPILER%100)/10],,
    AC_MSG_FAILURE([[[$0]] unknown intel compiler version]))
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch,
    [(__INTEL_COMPILER%10)],,
    AC_MSG_FAILURE([[[$0]] unknown intel compiler version]))
  ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="$_ax_[]_AC_LANG_ABBREV[]_compiler_version_major.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch"
  ])

# for IBM
AC_DEFUN([_AX_COMPILER_VERSION_IBM],
  [ dnl
  dnl check between z/OS C/C++  and XL C/C++
  AC_COMPILE_IFELSE([
    AC_LANG_PROGRAM([],
      [
        #if defined(__COMPILER_VER__)
        choke me;
        #endif
      ])],
    [
      AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_major,
        [__xlC__/100],,
      	AC_MSG_FAILURE([[[$0]] unknown IBM compiler major version]))
      AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor,
        [__xlC__%100],,
      	AC_MSG_FAILURE([[[$0]] unknown IBM compiler minor version]))
      AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch,
        [__xlC_ver__/0x100],,
      	AC_MSG_FAILURE([[[$0]] unknown IBM compiler patch version]))
      AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_build,
        [__xlC_ver__%0x100],,
      	AC_MSG_FAILURE([[[$0]] unknown IBM compiler build version]))
      ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="$_ax_[]_AC_LANG_ABBREV[]_compiler_version_major.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_build"
    ],
    [
      AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch,
        [__xlC__%1000],,
      	AC_MSG_FAILURE([[[$0]] unknown IBM compiler patch version]))
      AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor,
        [(__xlC__/10000)%10],,
      	AC_MSG_FAILURE([[[$0]] unknown IBM compiler minor version]))
      AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_major,
        [(__xlC__/100000)%10],,
      	AC_MSG_FAILURE([[[$0]] unknown IBM compiler major version]))
      ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="$_ax_[]_AC_LANG_ABBREV[]_compiler_version_major.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch"
    ])
])

# for pathscale
AC_DEFUN([_AX_COMPILER_VERSION_PATHSCALE],[
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_major,
    __PATHCC__,,
    AC_MSG_FAILURE([[[$0]] unknown pathscale major]))
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor,
    __PATHCC_MINOR__,,
    AC_MSG_FAILURE([[[$0]] unknown pathscale minor]))
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch,
    [__PATHCC_PATCHLEVEL__],,
    AC_MSG_FAILURE([[[$0]] unknown pathscale patch level]))
  ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="$_ax_[]_AC_LANG_ABBREV[]_compiler_version_major.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch"
  ])

# for clang
AC_DEFUN([_AX_COMPILER_VERSION_CLANG],[
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_major,
    __clang_major__,,
    AC_MSG_FAILURE([[[$0]] unknown clang major]))
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor,
    __clang_minor__,,
    AC_MSG_FAILURE([[[$0]] unknown clang minor]))
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch,
    [__clang_patchlevel__],,0)
  ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="$_ax_[]_AC_LANG_ABBREV[]_compiler_version_major.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch"
  ])

# for crayc
AC_DEFUN([_AX_COMPILER_VERSION_CRAY],[
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_major,
    _RELEASE,,
    AC_MSG_FAILURE([[[$0]] unknown crayc release]))
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor,
    _RELEASE_MINOR,,
    AC_MSG_FAILURE([[[$0]] unknown crayc minor]))
  ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="$_ax_[]_AC_LANG_ABBREV[]_compiler_version_major.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor"
  ])

# for fujitsu
AC_DEFUN([_AX_COMPILER_VERSION_FUJITSU],[
  AC_COMPUTE_INT(ax_cv_[]_AC_LANG_ABBREV[]_compiler_version,
                 __FCC_VERSION,,
		 AC_MSG_FAILURE([[[$0]]unknown fujitsu release]))
  ])

# for GNU
AC_DEFUN([_AX_COMPILER_VERSION_GNU],[
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_major,
    __GNUC__,,
    AC_MSG_FAILURE([[[$0]] unknown gcc major]))
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor,
    __GNUC_MINOR__,,
    AC_MSG_FAILURE([[[$0]] unknown gcc minor]))
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch,
    [__GNUC_PATCHLEVEL__],,
    AC_MSG_FAILURE([[[$0]] unknown gcc patch level]))
  ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="$_ax_[]_AC_LANG_ABBREV[]_compiler_version_major.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch"
  ])

# For sun
AC_DEFUN([_AX_COMPILER_VERSION_SUN],[
  m4_define([_AX_COMPILER_VERSION_SUN_NUMBER],
            [
	     #if defined(__SUNPRO_CC)
	     __SUNPRO_CC
	     #else
	     __SUNPRO_C
	     #endif
	    ])
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_until59,
    !!(_AX_COMPILER_VERSION_SUN_NUMBER < 0x1000),,
    AC_MSG_FAILURE([[[$0]] unknown sun release version]))
  AS_IF([test "X$_ax_[]_AC_LANG_ABBREV[]_compiler_version_until59" = X1],
    [dnl
      AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch,
        _AX_COMPILER_VERSION_SUN_NUMBER % 0x10,,
	AC_MSG_FAILURE([[[$0]] unknown sun patch version]))
      AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor,
        (_AX_COMPILER_VERSION_SUN_NUMBER / 0x10) % 0x10,,
        AC_MSG_FAILURE([[[$0]] unknown sun minor version]))
      AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_major,
        (_AX_COMPILER_VERSION_SUN_NUMBER / 0x100),,
        AC_MSG_FAILURE([[[$0]] unknown sun major version]))
    ],
    [dnl
      AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch,
        _AX_COMPILER_VERSION_SUN_NUMBER % 0x10,,
        AC_MSG_FAILURE([[[$0]] unknown sun patch version]))
      AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor,
        (_AX_COMPILER_VERSION_SUN_NUMBER / 0x100) % 0x100,,
        AC_MSG_FAILURE([[[$0]] unknown sun minor version]))
      AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_major,
        (_AX_COMPILER_VERSION_SUN_NUMBER / 0x1000),,
        AC_MSG_FAILURE([[[$0]] unknown sun major version]))
    ])
  ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="$_ax_[]_AC_LANG_ABBREV[]_compiler_version_major.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch"
])

AC_DEFUN([_AX_COMPILER_VERSION_HP],[
  m4_define([_AX_COMPILER_VERSION_HP_NUMBER],
            [
	     #if defined(__HP_cc)
	     __HP_cc
	     #else
	     __HP_aCC
	     #endif
	    ])
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_untilA0121,
    !!(_AX_COMPILER_VERSION_HP_NUMBER <= 1),,
    AC_MSG_FAILURE([[[$0]] unknown hp release version]))
  AS_IF([test "X$_ax_[]_AC_LANG_ABBREV[]_compiler_version_untilA0121" = X1],
    [dnl By default output last version with this behavior.
     dnl it is so old
      ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="01.21.00"
    ],
    [dnl
      AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch,
        (_AX_COMPILER_VERSION_HP_NUMBER % 100),,
        AC_MSG_FAILURE([[[$0]] unknown hp release version]))
      AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor,
        ((_AX_COMPILER_VERSION_HP_NUMBER / 100)%100),,
        AC_MSG_FAILURE([[[$0]] unknown hp minor version]))
      AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_major,
        ((_AX_COMPILER_VERSION_HP_NUMBER / 10000)%100),,
        AC_MSG_FAILURE([[[$0]] unknown hp major version]))
      ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="$_ax_[]_AC_LANG_ABBREV[]_compiler_version_major.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch"
    ])
])

AC_DEFUN([_AX_COMPILER_VERSION_DEC],[dnl
  m4_define([_AX_COMPILER_VERSION_DEC_NUMBER],
            [
	     #if defined(__DECC_VER)
	     __DECC_VER
	     #else
	     __DECCXX_VER
	     #endif
	    ])
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch,
    (_AX_COMPILER_VERSION_DEC_NUMBER % 10000),,
    AC_MSG_FAILURE([[[$0]] unknown dec release version]))
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor,
    ((_AX_COMPILER_VERSION_DEC_NUMBER / 100000UL)%100),,
    AC_MSG_FAILURE([[[$0]] unknown dec minor version]))
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_major,
    ((_AX_COMPILER_VERSION_DEC_NUMBER / 10000000UL)%100),,
    AC_MSG_FAILURE([[[$0]] unknown dec major version]))
  ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="$_ax_[]_AC_LANG_ABBREV[]_compiler_version_major.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch"
  ])

# borland
AC_DEFUN([_AX_COMPILER_VERSION_BORLAND],[dnl
  m4_define([_AX_COMPILER_VERSION_TURBOC_NUMBER],
            [
	     #if defined(__TURBOC__)
	     __TURBOC__
	     #else
	     choke me
	     #endif
	    ])
  m4_define([_AX_COMPILER_VERSION_BORLANDC_NUMBER],
            [
	     #if defined(__BORLANDC__)
	     __BORLANDC__
	     #else
	     __CODEGEARC__
	     #endif
	    ])
 AC_COMPILE_IFELSE(
   [AC_LANG_PROGRAM(,
     _AX_COMPILER_VERSION_TURBOC_NUMBER)],
   [dnl TURBOC
     AC_COMPUTE_INT(
       _ax_[]_AC_LANG_ABBREV[]_compiler_version_turboc_raw,
       _AX_COMPILER_VERSION_TURBOC_NUMBER,,
       AC_MSG_FAILURE([[[$0]] unknown turboc version]))
     AS_IF(
       [test $_ax_[]_AC_LANG_ABBREV[]_compiler_version_turboc_raw -lt 661 || test $_ax_[]_AC_LANG_ABBREV[]_compiler_version_turboc_raw -gt 1023],
       [dnl compute normal version
        AC_COMPUTE_INT(
	  _ax_[]_AC_LANG_ABBREV[]_compiler_version_minor,
	  _AX_COMPILER_VERSION_TURBOC_NUMBER % 0x100,,
	  AC_MSG_FAILURE([[[$0]] unknown turboc minor version]))
	AC_COMPUTE_INT(
	  _ax_[]_AC_LANG_ABBREV[]_compiler_version_major,
	  (_AX_COMPILER_VERSION_TURBOC_NUMBER/0x100)%0x100,,
	  AC_MSG_FAILURE([[[$0]] unknown turboc major version]))
	ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="0turboc:$_ax_[]_AC_LANG_ABBREV[]_compiler_version_major.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor"],
      [dnl special version
       AS_CASE([$_ax_[]_AC_LANG_ABBREV[]_compiler_version_turboc_raw],
         [661],[ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="0turboc:1.00"],
	 [662],[ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="0turboc:1.01"],
         [663],[ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="0turboc:2.00"],
	 [
	 AC_MSG_WARN([[[$0]] unknown turboc version between 0x295 and 0x400 please report bug])
	 ax_cv_[]_AC_LANG_ABBREV[]_compiler_version=""
	 ])
      ])
    ],
    # borlandc
    [
    AC_COMPUTE_INT(
      _ax_[]_AC_LANG_ABBREV[]_compiler_version_borlandc_raw,
      _AX_COMPILER_VERSION_BORLANDC_NUMBER,,
      AC_MSG_FAILURE([[[$0]] unknown borlandc version]))
    AS_CASE([$_ax_[]_AC_LANG_ABBREV[]_compiler_version_borlandc_raw],
      dnl BORLANC++ before 5.5
      [512] ,[ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="1borlanc:2.00"],
      [1024],[ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="1borlanc:3.00"],
      [1024],[ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="1borlanc:3.00"],
      [1040],[ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="1borlanc:3.1"],
      [1106],[ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="1borlanc:4.0"],
      [1280],[ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="1borlanc:5.0"],
      [1312],[ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="1borlanc:5.02"],
      dnl C++ Builder era
      [1328],[ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="2cppbuilder:3.0"],
      [1344],[ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="2cppbuilder:4.0"],
      dnl BORLANC++ after 5.5
      [1360],[ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="3borlancpp:5.5"],
      [1361],[ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="3borlancpp:5.51"],
      [1378],[ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="3borlancpp:5.6.4"],
      dnl C++ Builder with year number
      [1392],[ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="4cppbuilder:2006"],
      [1424],[ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="4cppbuilder:2007"],
      [1555],[ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="4cppbuilder:2009"],
      [1569],[ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="4cppbuilder:2010"],
      dnl XE version
      [1584],[ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="5xe"],
      [1600],[ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="5xe:2"],
      [1616],[ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="5xe:3"],
      [1632],[ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="5xe:4"],
      [
      AC_MSG_WARN([[[$0]] Unknow borlanc compiler version $_ax_[]_AC_LANG_ABBREV[]_compiler_version_borlandc_raw please report bug])
      ])
    ])
  ])

# COMO
AC_DEFUN([_AX_COMPILER_VERSION_COMEAU],
  [ dnl
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor,
    [__COMO_VERSION__%100],,
    AC_MSG_FAILURE([[[$0]] unknown comeau compiler minor version]))
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_major,
    [(__COMO_VERSION__/100)%10],,
    AC_MSG_FAILURE([[[$0]] unknown comeau compiler major version]))
  ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="$_ax_[]_AC_LANG_ABBREV[]_compiler_version_major.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor"
  ])

# KAI
AC_DEFUN([_AX_COMPILER_VERSION_KAI],[
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch,
    [__KCC_VERSION%100],,
    AC_MSG_FAILURE([[[$0]] unknown kay compiler patch version]))
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor,
    [(__KCC_VERSION/100)%10],,
    AC_MSG_FAILURE([[[$0]] unknown kay compiler minor version]))
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_major,
    [(__KCC_VERSION/1000)%10],,
    AC_MSG_FAILURE([[[$0]] unknown kay compiler major version]))
  ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="$_ax_[]_AC_LANG_ABBREV[]_compiler_version_major.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch"
  ])

dnl LCC
dnl LCC does not output version...

# SGI
AC_DEFUN([_AX_COMPILER_VERSION_SGI],[
   m4_define([_AX_COMPILER_VERSION_SGI_NUMBER],
            [
	     #if defined(_COMPILER_VERSION)
	     _COMPILER_VERSION
	     #else
	     _SGI_COMPILER_VERSION
	     #endif
	    ])
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch,
    [_AX_COMPILER_VERSION_SGI_NUMBER%10],,
    AC_MSG_FAILURE([[[$0]] unknown SGI compiler patch version]))
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor,
    [(_AX_COMPILER_VERSION_SGI_NUMBER/10)%10],,
    AC_MSG_FAILURE([[[$0]] unknown SGI compiler minor version]))
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_major,
    [(_AX_COMPILER_VERSION_SGI_NUMBER/100)%10],,
    AC_MSG_FAILURE([[[$0]] unknown SGI compiler major version]))
  ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="$_ax_[]_AC_LANG_ABBREV[]_compiler_version_major.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch"
  ])

# microsoft
AC_DEFUN([_AX_COMPILER_VERSION_MICROSOFT],[
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor,
    _MSC_VER%100,,
    AC_MSG_FAILURE([[[$0]] unknown microsoft compiler minor version]))
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_major,
    (_MSC_VER/100)%100,,
    AC_MSG_FAILURE([[[$0]] unknown microsoft compiler major version]))
  dnl could be overriden
  _ax_[]_AC_LANG_ABBREV[]_compiler_version_patch=0
  _ax_[]_AC_LANG_ABBREV[]_compiler_version_build=0
  # special case for version 6
  AS_IF([test "X$_ax_[]_AC_LANG_ABBREV[]_compiler_version_major" = "X12"],
    [AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch,
       _MSC_FULL_VER%1000,,
       _ax_[]_AC_LANG_ABBREV[]_compiler_version_patch=0)])
  # for version 7
  AS_IF([test "X$_ax_[]_AC_LANG_ABBREV[]_compiler_version_major" = "X13"],
    [AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch,
       _MSC_FULL_VER%1000,,
       AC_MSG_FAILURE([[[$0]] unknown microsoft compiler patch version]))
    ])
  # for version > 8
 AS_IF([test $_ax_[]_AC_LANG_ABBREV[]_compiler_version_major -ge 14],
    [AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch,
       _MSC_FULL_VER%10000,,
       AC_MSG_FAILURE([[[$0]] unknown microsoft compiler patch version]))
    ])
 AS_IF([test $_ax_[]_AC_LANG_ABBREV[]_compiler_version_major -ge 15],
    [AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_build,
       _MSC_BUILD,,
       AC_MSG_FAILURE([[[$0]] unknown microsoft compiler build version]))
    ])
 ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="$_ax_[]_AC_LANG_ABBREV[]_compiler_version_major.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_build"
 ])

# for metrowerks
AC_DEFUN([_AX_COMPILER_VERSION_METROWERKS],[dnl
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch,
    __MWERKS__%0x100,,
    AC_MSG_FAILURE([[[$0]] unknown metrowerks compiler patch version]))
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor,
    (__MWERKS__/0x100)%0x10,,
    AC_MSG_FAILURE([[[$0]] unknown metrowerks compiler minor version]))
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_major,
    (__MWERKS__/0x1000)%0x10,,
    AC_MSG_FAILURE([[[$0]] unknown metrowerks compiler major version]))
  ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="$_ax_[]_AC_LANG_ABBREV[]_compiler_version_major.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch"
  ])

# for watcom
AC_DEFUN([_AX_COMPILER_VERSION_WATCOM],[dnl
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor,
    __WATCOMC__%100,,
    AC_MSG_FAILURE([[[$0]] unknown watcom compiler minor version]))
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_major,
    (__WATCOMC__/100)%100,,
    AC_MSG_FAILURE([[[$0]] unknown watcom compiler major version]))
  ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="$_ax_[]_AC_LANG_ABBREV[]_compiler_version_major.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor"
  ])

# for PGI
AC_DEFUN([_AX_COMPILER_VERSION_PORTLAND],[
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_major,
    __PGIC__,,
    AC_MSG_FAILURE([[[$0]] unknown pgi major]))
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor,
    __PGIC_MINOR__,,
    AC_MSG_FAILURE([[[$0]] unknown pgi minor]))
  AC_COMPUTE_INT(_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch,
    [__PGIC_PATCHLEVEL__],,
    AC_MSG_FAILURE([[[$0]] unknown pgi patch level]))
  ax_cv_[]_AC_LANG_ABBREV[]_compiler_version="$_ax_[]_AC_LANG_ABBREV[]_compiler_version_major.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_minor.$_ax_[]_AC_LANG_ABBREV[]_compiler_version_patch"
  ])

# tcc
AC_DEFUN([_AX_COMPILER_VERSION_TCC],[
  ax_cv_[]_AC_LANG_ABBREV[]_compiler_version=[`tcc -v | $SED 's/^[ ]*tcc[ ]\+version[ ]\+\([0-9.]\+\).*/\1/g'`]
  ])
# main entry point
AC_DEFUN([AX_COMPILER_VERSION],[dnl
  AC_REQUIRE([AX_COMPILER_VENDOR])
  AC_REQUIRE([AC_PROG_SED])
  AC_CACHE_CHECK([for _AC_LANG compiler version],
    ax_cv_[]_AC_LANG_ABBREV[]_compiler_version,
    [ dnl
      AS_CASE([$ax_cv_[]_AC_LANG_ABBREV[]_compiler_vendor],
        [intel],[_AX_COMPILER_VERSION_INTEL],
	[ibm],[_AX_COMPILER_VERSION_IBM],
	[pathscale],[_AX_COMPILER_VERSION_PATHSCALE],
	[clang],[_AX_COMPILER_VERSION_CLANG],
	[cray],[_AX_COMPILER_VERSION_CRAY],
	[fujitsu],[_AX_COMPILER_VERSION_FUJITSU],
        [gnu],[_AX_COMPILER_VERSION_GNU],
	[sun],[_AX_COMPILER_VERSION_SUN],
	[hp],[_AX_COMPILER_VERSION_HP],
	[dec],[_AX_COMPILER_VERSION_DEC],
	[borland],[_AX_COMPILER_VERSION_BORLAND],
	[comeau],[_AX_COMPILER_VERSION_COMEAU],
	[kai],[_AX_COMPILER_VERSION_KAI],
	[sgi],[_AX_COMPILER_VERSION_SGI],
	[microsoft],[_AX_COMPILER_VERSION_MICROSOFT],
	[metrowerks],[_AX_COMPILER_VERSION_METROWERKS],
	[watcom],[_AX_COMPILER_VERSION_WATCOM],
	[portland],[_AX_COMPILER_VERSION_PORTLAND],
	[tcc],[_AX_COMPILER_VERSION_TCC],
  	[ax_cv_[]_AC_LANG_ABBREV[]_compiler_version=""])
    ])
])

# ===========================================================================
#          http://www.gnu.org/software/autoconf-archive/ax_ext.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_EXT
#
# DESCRIPTION
#
#   Find supported SIMD extensions by requesting cpuid. When a SIMD
#   extension is found, the -m"simdextensionname" is added to SIMD_FLAGS if
#   compiler supports it. For example, if "sse2" is available then "-msse2"
#   is added to SIMD_FLAGS.
#
#   Find other supported CPU extensions by requesting cpuid. When a
#   processor extension is found, the -m"extensionname" is added to
#   CPUEXT_FLAGS if compiler supports it. For example, if "bmi2" is
#   available then "-mbmi2" is added to CPUEXT_FLAGS.
#
#   This macro calls:
#
#     AC_SUBST(SIMD_FLAGS)
#     AC_SUBST(CPUEXT_FLAGS)
#
#   And defines:
#
#     HAVE_RDRND / HAVE_BMI1 / HAVE_BMI2 / HAVE_ADX / HAVE_MPX
#     HAVE_PREFETCHWT1 / HAVE_ABM / HAVE_MMX / HAVE_SSE / HAVE_SSE2
#     HAVE_SSE3 / HAVE_SSSE3 / HAVE_SSE4_1 / HAVE_SSE4_2 / HAVE_SSE4a
#     HAVE_SHA / HAVE_AES / HAVE_AVX / HAVE_FMA3 / HAVE_FMA4 / HAVE_XOP
#     HAVE_AVX2 / HAVE_AVX512_F / HAVE_AVX512_CD / HAVE_AVX512_PF
#     HAVE_AVX512_ER / HAVE_AVX512_VL / HAVE_AVX512_BW / HAVE_AVX512_DQ
#     HAVE_AVX512_IFMA / HAVE_AVX512_VBMI
#
# LICENSE
#
#   Copyright (c) 2007 Christophe Tournayre <turn3r@users.sourceforge.net>
#   Copyright (c) 2013,2015 Michael Petch <mpetch@capp-sysware.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 15

AC_DEFUN([AX_EXT],
[
  AC_REQUIRE([AC_CANONICAL_HOST])
  AC_REQUIRE([AC_PROG_CC])

  CPUEXT_FLAGS=""
  SIMD_FLAGS=""

  case $host_cpu in
    powerpc*)
      AC_CACHE_CHECK([whether altivec is supported], [ax_cv_have_altivec_ext],
          [
            if test `/usr/sbin/sysctl -a 2>/dev/null| grep -c hw.optional.altivec` != 0; then
                if test `/usr/sbin/sysctl -n hw.optional.altivec` = 1; then
                  ax_cv_have_altivec_ext=yes
                fi
            fi
          ])

          if test "$ax_cv_have_altivec_ext" = yes; then
            AC_DEFINE(HAVE_ALTIVEC,,[Support Altivec instructions])
            AX_CHECK_COMPILE_FLAG(-faltivec, SIMD_FLAGS="$SIMD_FLAGS -faltivec", [])
          fi
    ;;

    i[[3456]]86*|x86_64*|amd64*)

      AC_REQUIRE([AX_GCC_X86_CPUID])
      AC_REQUIRE([AX_GCC_X86_CPUID_COUNT])
      AC_REQUIRE([AX_GCC_X86_AVX_XGETBV])

      eax_cpuid0=0
      AX_GCC_X86_CPUID(0x00000000)
      if test "$ax_cv_gcc_x86_cpuid_0x00000000" != "unknown";
      then
        eax_cpuid0=`echo $ax_cv_gcc_x86_cpuid_0x00000000 | cut -d ":" -f 1`
      fi

      eax_cpuid80000000=0
      AX_GCC_X86_CPUID(0x80000000)
      if test "$ax_cv_gcc_x86_cpuid_0x80000000" != "unknown";
      then
        eax_cpuid80000000=`echo $ax_cv_gcc_x86_cpuid_0x80000000 | cut -d ":" -f 1`
      fi

      ecx_cpuid1=0
      edx_cpuid1=0
      if test "$((0x$eax_cpuid0))" -ge 1 ; then
        AX_GCC_X86_CPUID(0x00000001)
        if test "$ax_cv_gcc_x86_cpuid_0x00000001" != "unknown";
        then
          ecx_cpuid1=`echo $ax_cv_gcc_x86_cpuid_0x00000001 | cut -d ":" -f 3`
          edx_cpuid1=`echo $ax_cv_gcc_x86_cpuid_0x00000001 | cut -d ":" -f 4`
        fi
      fi

      ebx_cpuid7=0
      ecx_cpuid7=0
      if test "$((0x$eax_cpuid0))" -ge 7 ; then
        AX_GCC_X86_CPUID_COUNT(0x00000007, 0x00)
        if test "$ax_cv_gcc_x86_cpuid_0x00000007" != "unknown";
        then
          ebx_cpuid7=`echo $ax_cv_gcc_x86_cpuid_0x00000007 | cut -d ":" -f 2`
          ecx_cpuid7=`echo $ax_cv_gcc_x86_cpuid_0x00000007 | cut -d ":" -f 3`
        fi
      fi

      ecx_cpuid80000001=0
      edx_cpuid80000001=0
      if test "$((0x$eax_cpuid80000000))" -ge "$((0x80000001))" ; then
        AX_GCC_X86_CPUID(0x80000001)
        if test "$ax_cv_gcc_x86_cpuid_0x80000001" != "unknown";
        then
          ecx_cpuid80000001=`echo $ax_cv_gcc_x86_cpuid_0x80000001 | cut -d ":" -f 3`
          edx_cpuid80000001=`echo $ax_cv_gcc_x86_cpuid_0x80000001 | cut -d ":" -f 4`
        fi
      fi

      AC_CACHE_VAL([ax_cv_have_mmx_os_support_ext],
      [
        ax_cv_have_mmx_os_support_ext=yes
      ])

      ax_cv_have_none_os_support_ext=yes

      AC_CACHE_VAL([ax_cv_have_sse_os_support_ext],
      [
        ax_cv_have_sse_os_support_ext=no,
        if test "$((0x$edx_cpuid1>>25&0x01))" = 1; then
          AC_LANG_PUSH([C])
          AC_TRY_RUN([
#include <signal.h>
#include <stdlib.h>
            /* No way at ring1 to ring3 in protected mode to check the CR0 and CR4
               control registers directly. Execute an SSE instruction.
               If it raises SIGILL then OS doesn't support SSE based instructions */
            void sig_handler(int signum){ exit(1); }
            int main(){
              signal(SIGILL, sig_handler);
              /* SSE instruction xorps  %xmm0,%xmm0 */
              __asm__ __volatile__ (".byte 0x0f, 0x57, 0xc0");
              return 0;
            }],
            ax_cv_have_sse_os_support_ext=yes,
            ax_cv_have_sse_os_support_ext=no,
            ax_cv_have_sse_os_support_ext=no)
          AC_LANG_POP([C])
        fi
      ])

      xgetbv_eax=0
      if test "$((0x$ecx_cpuid1>>28&0x01))" = 1; then
        AX_GCC_X86_AVX_XGETBV(0x00000000)

        if test x"$ax_cv_gcc_x86_avx_xgetbv_0x00000000" != x"unknown"; then
          xgetbv_eax=`echo $ax_cv_gcc_x86_avx_xgetbv_0x00000000 | cut -d ":" -f 1`
        fi

        AC_CACHE_VAL([ax_cv_have_avx_os_support_ext],
        [
          ax_cv_have_avx_os_support_ext=no
          if test "$((0x$ecx_cpuid1>>27&0x01))" = 1; then
            if test "$((0x$xgetbv_eax&0x6))" = 6; then
              ax_cv_have_avx_os_support_ext=yes
            fi
          fi
        ])
      fi

      AC_CACHE_VAL([ax_cv_have_avx512_os_support_ext],
      [
        ax_cv_have_avx512_os_support_ext=no
        if test "$ax_cv_have_avx_os_support_ext" = yes; then
          if test "$((0x$xgetbv_eax&0xe6))" = "$((0xe6))"; then
            ax_cv_have_avx512_os_support_ext=yes
          fi
        fi
      ])

      for ac_instr_info dnl
      in "none;rdrnd;RDRND;ecx_cpuid1,30;-mrdrnd;HAVE_RDRND;CPUEXT_FLAGS" dnl
         "none;bmi1;BMI1;ebx_cpuid7,3;-mbmi;HAVE_BMI1;CPUEXT_FLAGS" dnl
         "none;bmi2;BMI2;ebx_cpuid7,8;-mbmi2;HAVE_BMI2;CPUEXT_FLAGS" dnl
         "none;adx;ADX;ebx_cpuid7,19;-madx;HAVE_ADX;CPUEXT_FLAGS" dnl
         "none;mpx;MPX;ebx_cpuid7,14;-mmpx;HAVE_MPX;CPUEXT_FLAGS" dnl
         "none;prefetchwt1;PREFETCHWT1;ecx_cpuid7,0;-mprefetchwt1;HAVE_PREFETCHWT1;CPUEXT_FLAGS" dnl
         "none;abm;ABM;ecx_cpuid80000001,5;-mabm;HAVE_ABM;CPUEXT_FLAGS" dnl
         "mmx;mmx;MMX;edx_cpuid1,23;-mmmx;HAVE_MMX;SIMD_FLAGS" dnl
         "sse;sse;SSE;edx_cpuid1,25;-msse;HAVE_SSE;SIMD_FLAGS" dnl
         "sse;sse2;SSE2;edx_cpuid1,26;-msse2;HAVE_SSE2;SIMD_FLAGS" dnl
         "sse;sse3;SSE3;ecx_cpuid1,1;-msse3;HAVE_SSE3;SIMD_FLAGS" dnl
         "sse;ssse3;SSSE3;ecx_cpuid1,9;-mssse3;HAVE_SSSE3;SIMD_FLAGS" dnl
         "sse;sse41;SSE4.1;ecx_cpuid1,19;-msse4.1;HAVE_SSE4_1;SIMD_FLAGS" dnl
         "sse;sse42;SSE4.2;ecx_cpuid1,20;-msse4.2;HAVE_SSE4_2;SIMD_FLAGS" dnl
         "sse;sse4a;SSE4a;ecx_cpuid80000001,6;-msse4a;HAVE_SSE4a;SIMD_FLAGS" dnl
         "sse;sha;SHA;ebx_cpuid7,29;-msha;HAVE_SHA;SIMD_FLAGS" dnl
         "sse;aes;AES;ecx_cpuid1,25;-maes;HAVE_AES;SIMD_FLAGS" dnl
         "avx;avx;AVX;ecx_cpuid1,28;-mavx;HAVE_AVX;SIMD_FLAGS" dnl
         "avx;fma3;FMA3;ecx_cpuid1,12;-mfma;HAVE_FMA3;SIMD_FLAGS" dnl
         "avx;fma4;FMA4;ecx_cpuid80000001,16;-mfma4;HAVE_FMA4;SIMD_FLAGS" dnl
         "avx;xop;XOP;ecx_cpuid80000001,11;-mxop;HAVE_XOP;SIMD_FLAGS" dnl
         "avx;avx2;AVX2;ebx_cpuid7,5;-mavx2;HAVE_AVX2;SIMD_FLAGS" dnl
         "avx512;avx512f;AVX512-F;ebx_cpuid7,16;-mavx512f;HAVE_AVX512_F;SIMD_FLAGS" dnl
         "avx512;avx512cd;AVX512-CD;ebx_cpuid7,28;-mavx512cd;HAVE_AVX512_CD;SIMD_FLAGS" dnl
         "avx512;avx512pf;AVX512-PF;ebx_cpuid7,26;-mavx512pf;HAVE_AVX512_PF;SIMD_FLAGS" dnl
         "avx512;avx512er;AVX512-ER;ebx_cpuid7,27;-mavx512er;HAVE_AVX512_ER;SIMD_FLAGS" dnl
         "avx512;avx512vl;AVX512-VL;ebx_cpuid7,31;-mavx512vl;HAVE_AVX512_VL;SIMD_FLAGS" dnl
         "avx512;avx512bw;AVX512-BW;ebx_cpuid7,30;-mavx512bw;HAVE_AVX512_BW;SIMD_FLAGS" dnl
         "avx512;avx512dq;AVX512-DQ;ebx_cpuid7,17;-mavx512dq;HAVE_AVX512_DQ;SIMD_FLAGS" dnl
         "avx512;avx512ifma;AVX512-IFMA;ebx_cpuid7,21;-mavx512ifma;HAVE_AVX512_IFMA;SIMD_FLAGS" dnl
         "avx512;avx512vbmi;AVX512-VBMI;ecx_cpuid7,1;-mavx512vbmi;HAVE_AVX512_VBMI;SIMD_FLAGS" dnl
         #
      do ac_instr_os_support=$(eval echo \$ax_cv_have_$(echo $ac_instr_info | cut -d ";" -f 1)_os_support_ext)
         ac_instr_acvar=$(echo $ac_instr_info | cut -d ";" -f 2)
         ac_instr_shortname=$(echo $ac_instr_info | cut -d ";" -f 3)
         ac_instr_chk_loc=$(echo $ac_instr_info | cut -d ";" -f 4)
         ac_instr_chk_reg=0x$(eval echo \$$(echo $ac_instr_chk_loc | cut -d "," -f 1))
         ac_instr_chk_bit=$(echo $ac_instr_chk_loc | cut -d "," -f 2)
         ac_instr_compiler_flags=$(echo $ac_instr_info | cut -d ";" -f 5)
         ac_instr_have_define=$(echo $ac_instr_info | cut -d ";" -f 6)
         ac_instr_flag_type=$(echo $ac_instr_info | cut -d ";" -f 7)

         AC_CACHE_CHECK([whether ${ac_instr_shortname} is supported by the processor], [ax_cv_have_${ac_instr_acvar}_cpu_ext],
         [
           eval ax_cv_have_${ac_instr_acvar}_cpu_ext=no
           if test "$((${ac_instr_chk_reg}>>${ac_instr_chk_bit}&0x01))" = 1 ; then
             eval ax_cv_have_${ac_instr_acvar}_cpu_ext=yes
           fi
         ])

         if test x"$(eval echo \$ax_cv_have_${ac_instr_acvar}_cpu_ext)" = x"yes"; then
           AC_CACHE_CHECK([whether ${ac_instr_shortname} is supported by the processor and OS], [ax_cv_have_${ac_instr_acvar}_ext],
           [
             eval ax_cv_have_${ac_instr_acvar}_ext=no
             if test x"${ac_instr_os_support}" = x"yes"; then
               eval ax_cv_have_${ac_instr_acvar}_ext=yes
             fi
           ])

           if test "$(eval echo \$ax_cv_have_${ac_instr_acvar}_ext)" = yes; then
             AX_CHECK_COMPILE_FLAG(${ac_instr_compiler_flags}, eval ax_cv_support_${ac_instr_acvar}_ext=yes,
                                                               eval ax_cv_support_${ac_instr_acvar}_ext=no)
             if test x"$(eval echo \$ax_cv_support_${ac_instr_acvar}_ext)" = x"yes"; then
               eval ${ac_instr_flag_type}=\"\$${ac_instr_flag_type} ${ac_instr_compiler_flags}\"
               AC_DEFINE_UNQUOTED([${ac_instr_have_define}])
             else
               AC_MSG_WARN([Your processor and OS supports ${ac_instr_shortname} instructions but not your compiler, can you try another compiler?])
             fi
           else
             if test x"${ac_instr_os_support}" = x"no"; then
               AC_CACHE_VAL(ax_cv_support_${ac_instr_acvar}_ext, eval ax_cv_support_${ac_instr_acvar}_ext=no)
               AC_MSG_WARN([Your processor supports ${ac_instr_shortname}, but your OS doesn't])
             fi
           fi
         else
           AC_CACHE_VAL(ax_cv_have_${ac_instr_acvar}_ext, eval ax_cv_have_${ac_instr_acvar}_ext=no)
           AC_CACHE_VAL(ax_cv_support_${ac_instr_acvar}_ext, eval ax_cv_support_${ac_instr_acvar}_ext=no)
         fi
      done
  ;;
  esac

  AH_TEMPLATE([HAVE_RDRND],[Define to 1 to support Digital Random Number Generator])
  AH_TEMPLATE([HAVE_BMI1],[Define to 1 to support Bit Manipulation Instruction Set 1])
  AH_TEMPLATE([HAVE_BMI2],[Define to 1 to support Bit Manipulation Instruction Set 2])
  AH_TEMPLATE([HAVE_ADX],[Define to 1 to support Multi-Precision Add-Carry Instruction Extensions])
  AH_TEMPLATE([HAVE_MPX],[Define to 1 to support Memory Protection Extensions])
  AH_TEMPLATE([HAVE_PREFETCHWT1],[Define to 1 to support Prefetch Vector Data Into Caches WT1])
  AH_TEMPLATE([HAVE_ABM],[Define to 1 to support Advanced Bit Manipulation])
  AH_TEMPLATE([HAVE_MMX],[Define to 1 to support Multimedia Extensions])
  AH_TEMPLATE([HAVE_SSE],[Define to 1 to support Streaming SIMD Extensions])
  AH_TEMPLATE([HAVE_SSE2],[Define to 1 to support Streaming SIMD Extensions])
  AH_TEMPLATE([HAVE_SSE3],[Define to 1 to support Streaming SIMD Extensions 3])
  AH_TEMPLATE([HAVE_SSSE3],[Define to 1 to support Supplemental Streaming SIMD Extensions 3])
  AH_TEMPLATE([HAVE_SSE4_1],[Define to 1 to support Streaming SIMD Extensions 4.1])
  AH_TEMPLATE([HAVE_SSE4_2],[Define to 1 to support Streaming SIMD Extensions 4.2])
  AH_TEMPLATE([HAVE_SSE4a],[Define to 1 to support AMD Streaming SIMD Extensions 4a])
  AH_TEMPLATE([HAVE_SHA],[Define to 1 to support Secure Hash Algorithm Extension])
  AH_TEMPLATE([HAVE_AES],[Define to 1 to support Advanced Encryption Standard New Instruction Set (AES-NI)])
  AH_TEMPLATE([HAVE_AVX],[Define to 1 to support Advanced Vector Extensions])
  AH_TEMPLATE([HAVE_FMA3],[Define to 1 to support  Fused Multiply-Add Extensions 3])
  AH_TEMPLATE([HAVE_FMA4],[Define to 1 to support Fused Multiply-Add Extensions 4])
  AH_TEMPLATE([HAVE_XOP],[Define to 1 to support eXtended Operations Extensions])
  AH_TEMPLATE([HAVE_AVX2],[Define to 1 to support Advanced Vector Extensions 2])
  AH_TEMPLATE([HAVE_AVX512_F],[Define to 1 to support AVX-512 Foundation Extensions])
  AH_TEMPLATE([HAVE_AVX512_CD],[Define to 1 to support AVX-512 Conflict Detection Instructions])
  AH_TEMPLATE([HAVE_AVX512_PF],[Define to 1 to support AVX-512 Conflict Prefetch Instructions])
  AH_TEMPLATE([HAVE_AVX512_ER],[Define to 1 to support AVX-512 Exponential & Reciprocal Instructions])
  AH_TEMPLATE([HAVE_AVX512_VL],[Define to 1 to support AVX-512 Vector Length Extensions])
  AH_TEMPLATE([HAVE_AVX512_BW],[Define to 1 to support AVX-512 Byte and Word Instructions])
  AH_TEMPLATE([HAVE_AVX512_DQ],[Define to 1 to support AVX-512 Doubleword and Quadword Instructions])
  AH_TEMPLATE([HAVE_AVX512_IFMA],[Define to 1 to support AVX-512 Integer Fused Multiply Add Instructions])
  AH_TEMPLATE([HAVE_AVX512_VBMI],[Define to 1 to support AVX-512 Vector Byte Manipulation Instructions])
  AC_SUBST(SIMD_FLAGS)
  AC_SUBST(CPUEXT_FLAGS)
])

# ===========================================================================
#   http://www.gnu.org/software/autoconf-archive/ax_gcc_x86_avx_xgetbv.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_GCC_X86_AVX_XGETBV
#
# DESCRIPTION
#
#   On later x86 processors with AVX SIMD support, with gcc or a compiler
#   that has a compatible syntax for inline assembly instructions, run a
#   small program that executes the xgetbv instruction with input OP. This
#   can be used to detect if the OS supports AVX instruction usage.
#
#   On output, the values of the eax and edx registers are stored as
#   hexadecimal strings as "eax:edx" in the cache variable
#   ax_cv_gcc_x86_avx_xgetbv.
#
#   If the xgetbv instruction fails (because you are running a
#   cross-compiler, or because you are not using gcc, or because you are on
#   a processor that doesn't have this instruction),
#   ax_cv_gcc_x86_avx_xgetbv_OP is set to the string "unknown".
#
#   This macro mainly exists to be used in AX_EXT.
#
# LICENSE
#
#   Copyright (c) 2013 Michael Petch <mpetch@capp-sysware.com>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 2

AC_DEFUN([AX_GCC_X86_AVX_XGETBV],
[AC_REQUIRE([AC_PROG_CC])
AC_LANG_PUSH([C])
AC_CACHE_CHECK(for x86-AVX xgetbv $1 output, ax_cv_gcc_x86_avx_xgetbv_$1,
 [AC_RUN_IFELSE([AC_LANG_PROGRAM([#include <stdio.h>], [
     int op = $1, eax, edx;
     FILE *f;
      /* Opcodes for xgetbv */
      __asm__ __volatile__ (".byte 0x0f, 0x01, 0xd0"
        : "=a" (eax), "=d" (edx)
        : "c" (op));
     f = fopen("conftest_xgetbv", "w"); if (!f) return 1;
     fprintf(f, "%x:%x\n", eax, edx);
     fclose(f);
     return 0;
])],
     [ax_cv_gcc_x86_avx_xgetbv_$1=`cat conftest_xgetbv`; rm -f conftest_xgetbv],
     [ax_cv_gcc_x86_avx_xgetbv_$1=unknown; rm -f conftest_xgetbv],
     [ax_cv_gcc_x86_avx_xgetbv_$1=unknown])])
AC_LANG_POP([C])
])

# ===========================================================================
#     http://www.gnu.org/software/autoconf-archive/ax_gcc_x86_cpuid.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_GCC_X86_CPUID(OP)
#   AX_GCC_X86_CPUID_COUNT(OP, COUNT)
#
# DESCRIPTION
#
#   On Pentium and later x86 processors, with gcc or a compiler that has a
#   compatible syntax for inline assembly instructions, run a small program
#   that executes the cpuid instruction with input OP. This can be used to
#   detect the CPU type. AX_GCC_X86_CPUID_COUNT takes an additional COUNT
#   parameter that gets passed into register ECX before calling cpuid.
#
#   On output, the values of the eax, ebx, ecx, and edx registers are stored
#   as hexadecimal strings as "eax:ebx:ecx:edx" in the cache variable
#   ax_cv_gcc_x86_cpuid_OP.
#
#   If the cpuid instruction fails (because you are running a
#   cross-compiler, or because you are not using gcc, or because you are on
#   a processor that doesn't have this instruction), ax_cv_gcc_x86_cpuid_OP
#   is set to the string "unknown".
#
#   This macro mainly exists to be used in AX_GCC_ARCHFLAG.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2008 Matteo Frigo
#   Copyright (c) 2015 Michael Petch <mpetch@capp-sysware.com>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 9

AC_DEFUN([AX_GCC_X86_CPUID],
[AX_GCC_X86_CPUID_COUNT($1, 0)
])

AC_DEFUN([AX_GCC_X86_CPUID_COUNT],
[AC_REQUIRE([AC_PROG_CC])
AC_LANG_PUSH([C])
AC_CACHE_CHECK(for x86 cpuid $1 output, ax_cv_gcc_x86_cpuid_$1,
 [AC_RUN_IFELSE([AC_LANG_PROGRAM([#include <stdio.h>], [
     int op = $1, level = $2, eax, ebx, ecx, edx;
     FILE *f;
      __asm__ __volatile__ ("xchg %%ebx, %1\n"
        "cpuid\n"
        "xchg %%ebx, %1\n"
        : "=a" (eax), "=r" (ebx), "=c" (ecx), "=d" (edx)
        : "a" (op), "2" (level));

     f = fopen("conftest_cpuid", "w"); if (!f) return 1;
     fprintf(f, "%x:%x:%x:%x\n", eax, ebx, ecx, edx);
     fclose(f);
     return 0;
])],
     [ax_cv_gcc_x86_cpuid_$1=`cat conftest_cpuid`; rm -f conftest_cpuid],
     [ax_cv_gcc_x86_cpuid_$1=unknown; rm -f conftest_cpuid],
     [ax_cv_gcc_x86_cpuid_$1=unknown])])
AC_LANG_POP([C])
])

# ===========================================================================
#      http://www.gnu.org/software/autoconf-archive/ax_lib_readline.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_LIB_READLINE
#
# DESCRIPTION
#
#   Searches for a readline compatible library. If found, defines
#   `HAVE_LIBREADLINE'. If the found library has the `add_history' function,
#   sets also `HAVE_READLINE_HISTORY'. Also checks for the locations of the
#   necessary include files and sets `HAVE_READLINE_H' or
#   `HAVE_READLINE_READLINE_H' and `HAVE_READLINE_HISTORY_H' or
#   'HAVE_HISTORY_H' if the corresponding include files exists.
#
#   The libraries that may be readline compatible are `libedit',
#   `libeditline' and `libreadline'. Sometimes we need to link a termcap
#   library for readline to work, this macro tests these cases too by trying
#   to link with `libtermcap', `libcurses' or `libncurses' before giving up.
#
#   Here is an example of how to use the information provided by this macro
#   to perform the necessary includes or declarations in a C file:
#
#     #ifdef HAVE_LIBREADLINE
#     #  if defined(HAVE_READLINE_READLINE_H)
#     #    include <readline/readline.h>
#     #  elif defined(HAVE_READLINE_H)
#     #    include <readline.h>
#     #  else /* !defined(HAVE_READLINE_H) */
#     extern char *readline ();
#     #  endif /* !defined(HAVE_READLINE_H) */
#     char *cmdline = NULL;
#     #else /* !defined(HAVE_READLINE_READLINE_H) */
#       /* no readline */
#     #endif /* HAVE_LIBREADLINE */
#
#     #ifdef HAVE_READLINE_HISTORY
#     #  if defined(HAVE_READLINE_HISTORY_H)
#     #    include <readline/history.h>
#     #  elif defined(HAVE_HISTORY_H)
#     #    include <history.h>
#     #  else /* !defined(HAVE_HISTORY_H) */
#     extern void add_history ();
#     extern int write_history ();
#     extern int read_history ();
#     #  endif /* defined(HAVE_READLINE_HISTORY_H) */
#       /* no history */
#     #endif /* HAVE_READLINE_HISTORY */
#
# LICENSE
#
#   Copyright (c) 2008 Ville Laurikari <vl@iki.fi>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 6

AU_ALIAS([VL_LIB_READLINE], [AX_LIB_READLINE])
AC_DEFUN([AX_LIB_READLINE], [
  AC_CACHE_CHECK([for a readline compatible library],
                 ax_cv_lib_readline, [
    ORIG_LIBS="$LIBS"
    for readline_lib in readline edit editline; do
      for termcap_lib in "" termcap curses ncurses; do
        if test -z "$termcap_lib"; then
          TRY_LIB="-l$readline_lib"
        else
          TRY_LIB="-l$readline_lib -l$termcap_lib"
        fi
        LIBS="$ORIG_LIBS $TRY_LIB"
        AC_TRY_LINK_FUNC(readline, ax_cv_lib_readline="$TRY_LIB")
        if test -n "$ax_cv_lib_readline"; then
          break
        fi
      done
      if test -n "$ax_cv_lib_readline"; then
        break
      fi
    done
    if test -z "$ax_cv_lib_readline"; then
      ax_cv_lib_readline="no"
    fi
    LIBS="$ORIG_LIBS"
  ])

  if test "$ax_cv_lib_readline" != "no"; then
    LIBS="$LIBS $ax_cv_lib_readline"
    AC_DEFINE(HAVE_LIBREADLINE, 1,
              [Define if you have a readline compatible library])
    AC_CHECK_HEADERS(readline.h readline/readline.h)
    AC_CACHE_CHECK([whether readline supports history],
                   ax_cv_lib_readline_history, [
      ax_cv_lib_readline_history="no"
      AC_TRY_LINK_FUNC(add_history, ax_cv_lib_readline_history="yes")
    ])
    if test "$ax_cv_lib_readline_history" = "yes"; then
      AC_DEFINE(HAVE_READLINE_HISTORY, 1,
                [Define if your readline library has \`add_history'])
      AC_CHECK_HEADERS(history.h readline/history.h)
    fi
  fi
])dnl

# ===========================================================================
#      http://www.gnu.org/software/autoconf-archive/ax_prog_cc_mpi.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_PROG_CC_MPI([MPI-WANTED-TEST[, ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   This macro tries to find out how to compile C programs that use MPI
#   (Message Passing Interface), a standard API for parallel process
#   communication (see http://www-unix.mcs.anl.gov/mpi/). The macro has to
#   be used instead of the standard macro AC_PROG_CC and will replace the
#   standard variable CC with the found compiler.
#
#   MPI-WANTED-TEST is used to test whether MPI is actually wanted by the
#   user. If MPI-WANTED_TEST is omitted or if it succeeds, the macro will
#   try to find out how to use MPI, if it fails, the macro will call
#   AC_PROG_CC to find a standard C compiler instead.
#
#   When MPI is found, ACTION-IF-FOUND will be executed, if MPI is not found
#   (or MPI-WANTED-TEST fails) ACTION-IF-NOT-FOUND is executed. If
#   ACTION-IF-FOUND is not set, the macro will define HAVE_MPI.
#
#   The following example demonstrates usage of the macro:
#
#     # If --with-mpi=auto is used, try to find MPI, but use standard C compiler if it is not found.
#     # If --with-mpi=yes is used, try to find MPI and fail if it isn't found.
#     # If --with-mpi=no is used, use a standard C compiler instead.
#     AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],
#         [compile with MPI (parallelization) support. If none is found,
#         MPI is not used. Default: auto])
#     ],,[with_mpi=auto])
#     #
#     AX_PROG_CC_MPI([test x"$with_mpi" != xno],[use_mpi=yes],[
#       use_mpi=no
#       if test x"$with_mpi" = xyes; then
#         AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
#       else
#         AC_MSG_WARN([No MPI compiler found, won't use MPI.])
#       fi
#     ])
#
# LICENSE
#
#   Copyright (c) 2010,2011 Olaf Lenz <olenz@icp.uni-stuttgart.de>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 1

AC_DEFUN([AX_PROG_CC_MPI], [
AC_PREREQ(2.50)

# Check for compiler
# Needs to be split off into an extra macro to ensure right expansion
# order.
AC_REQUIRE([_AX_PROG_CC_MPI],[_AX_PROG_CC_MPI([$1])])

AS_IF([test x"$_ax_prog_cc_mpi_mpi_wanted" = xno],
  [ _ax_prog_cc_mpi_mpi_found=no ],
  [
    AC_LANG_PUSH([C])
    # test whether MPI_Init is available
    # We do not use AC_SEARCH_LIBS here, as it caches its outcome and
    # thus disallows corresponding calls in the other AX_PROG_*_MPI
    # macros.
    for lib in NONE mpi mpich; do
      save_LIBS=$LIBS
      if test x"$lib" = xNONE; then
        AC_MSG_CHECKING([for function MPI_Init])
      else
        AC_MSG_CHECKING([for function MPI_Init in -l$lib])
        LIBS="-l$lib $LIBS"
      fi
      AC_LINK_IFELSE([AC_LANG_CALL([],[MPI_Init])],
        [ _ax_prog_cc_mpi_mpi_found=yes ],
        [ _ax_prog_cc_mpi_mpi_found=no ])
      AC_MSG_RESULT($_ax_prog_cc_mpi_mpi_found)
      if test "x$_ax_prog_cc_mpi_mpi_found" = "xyes"; then
        break;
      fi
      LIBS=$save_LIBS
    done

    # Check for header
    AS_IF([test x"$_ax_prog_cc_mpi_mpi_found" = xyes], [
      AC_MSG_CHECKING([for mpi.h])
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([#include <mpi.h>])],
        [ AC_MSG_RESULT(yes)],
        [ AC_MSG_RESULT(no)
         _ax_prog_cc_mpi_mpi_found=no
      ])
    ])
    AC_LANG_POP([C])
])

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
AS_IF([test x"$_ax_prog_cc_mpi_mpi_found" = xyes], [
        ifelse([$2],,[AC_DEFINE(HAVE_MPI,1,[Define if you have the MPI library.])],[$2])
        :
],[
        $3
        :
])

])dnl AX_PROG_CC_MPI

dnl _AX_PROG_CC_MPI is an internal macro required by AX_PROG_CC_MPI.
dnl To ensure the right expansion order, the main function AX_PROG_CC_MPI
dnl has to be split into two parts.
dnl
dnl Known MPI C compilers:
dnl  mpicc
dnl  mpixlc_r
dnl  mpixlc
dnl  hcc
dnl  mpxlc_r
dnl  mpxlc
dnl  sxmpicc  NEC SX
dnl  mpifcc   Fujitsu
dnl  mpgcc
dnl  mpcc
dnl  cmpicc
dnl  cc
dnl
AC_DEFUN([_AX_PROG_CC_MPI], [
  AC_ARG_VAR(MPICC,[MPI C compiler command])
  ifelse([$1],,[_ax_prog_cc_mpi_mpi_wanted=yes],[
    AC_MSG_CHECKING([whether to compile using MPI])
    if $1; then
      _ax_prog_cc_mpi_mpi_wanted=yes
    else
      _ax_prog_cc_mpi_mpi_wanted=no
    fi
    AC_MSG_RESULT($_ax_prog_cc_mpi_mpi_wanted)
  ])
  if test x"$_ax_prog_cc_mpi_mpi_wanted" = xyes; then
    if test -z "$CC" && test -n "$MPICC"; then
      CC="$MPICC"
    else
      AC_CHECK_TOOLS([CC], [mpicc mpixlc_r mpixlc hcc mpxlc_r mpxlc sxmpicc mpifcc mpgcc mpcc cmpicc cc gcc])
    fi
  fi
  AC_PROG_CC
])dnl _AX_PROG_CC_MPI

# pkg.m4 - Macros to locate and utilise pkg-config.            -*- Autoconf -*-
# serial 1 (pkg-config-0.24)
# 
# Copyright © 2004 Scott James Remnant <scott@netsplit.com>.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
#
# As a special exception to the GNU General Public License, if you
# distribute this file as part of a program that contains a
# configuration script generated by Autoconf, you may include it under
# the same distribution terms that you use for the rest of that program.

# PKG_PROG_PKG_CONFIG([MIN-VERSION])
# ----------------------------------
AC_DEFUN([PKG_PROG_PKG_CONFIG],
[m4_pattern_forbid([^_?PKG_[A-Z_]+$])
m4_pattern_allow([^PKG_CONFIG(_(PATH|LIBDIR|SYSROOT_DIR|ALLOW_SYSTEM_(CFLAGS|LIBS)))?$])
m4_pattern_allow([^PKG_CONFIG_(DISABLE_UNINSTALLED|TOP_BUILD_DIR|DEBUG_SPEW)$])
AC_ARG_VAR([PKG_CONFIG], [path to pkg-config utility])
AC_ARG_VAR([PKG_CONFIG_PATH], [directories to add to pkg-config's search path])
AC_ARG_VAR([PKG_CONFIG_LIBDIR], [path overriding pkg-config's built-in search path])

if test "x$ac_cv_env_PKG_CONFIG_set" != "xset"; then
	AC_PATH_TOOL([PKG_CONFIG], [pkg-config])
fi
if test -n "$PKG_CONFIG"; then
	_pkg_min_version=m4_default([$1], [0.9.0])
	AC_MSG_CHECKING([pkg-config is at least version $_pkg_min_version])
	if $PKG_CONFIG --atleast-pkgconfig-version $_pkg_min_version; then
		AC_MSG_RESULT([yes])
	else
		AC_MSG_RESULT([no])
		PKG_CONFIG=""
	fi
fi[]dnl
])# PKG_PROG_PKG_CONFIG

# PKG_CHECK_EXISTS(MODULES, [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
#
# Check to see whether a particular set of modules exists.  Similar
# to PKG_CHECK_MODULES(), but does not set variables or print errors.
#
# Please remember that m4 expands AC_REQUIRE([PKG_PROG_PKG_CONFIG])
# only at the first occurence in configure.ac, so if the first place
# it's called might be skipped (such as if it is within an "if", you
# have to call PKG_CHECK_EXISTS manually
# --------------------------------------------------------------
AC_DEFUN([PKG_CHECK_EXISTS],
[AC_REQUIRE([PKG_PROG_PKG_CONFIG])dnl
if test -n "$PKG_CONFIG" && \
    AC_RUN_LOG([$PKG_CONFIG --exists --print-errors "$1"]); then
  m4_default([$2], [:])
m4_ifvaln([$3], [else
  $3])dnl
fi])

# _PKG_CONFIG([VARIABLE], [COMMAND], [MODULES])
# ---------------------------------------------
m4_define([_PKG_CONFIG],
[if test -n "$$1"; then
    pkg_cv_[]$1="$$1"
 elif test -n "$PKG_CONFIG"; then
    PKG_CHECK_EXISTS([$3],
                     [pkg_cv_[]$1=`$PKG_CONFIG --[]$2 "$3" 2>/dev/null`
		      test "x$?" != "x0" && pkg_failed=yes ],
		     [pkg_failed=yes])
 else
    pkg_failed=untried
fi[]dnl
])# _PKG_CONFIG

# _PKG_SHORT_ERRORS_SUPPORTED
# -----------------------------
AC_DEFUN([_PKG_SHORT_ERRORS_SUPPORTED],
[AC_REQUIRE([PKG_PROG_PKG_CONFIG])
if $PKG_CONFIG --atleast-pkgconfig-version 0.20; then
        _pkg_short_errors_supported=yes
else
        _pkg_short_errors_supported=no
fi[]dnl
])# _PKG_SHORT_ERRORS_SUPPORTED


# PKG_CHECK_MODULES(VARIABLE-PREFIX, MODULES, [ACTION-IF-FOUND],
# [ACTION-IF-NOT-FOUND])
#
#
# Note that if there is a possibility the first call to
# PKG_CHECK_MODULES might not happen, you should be sure to include an
# explicit call to PKG_PROG_PKG_CONFIG in your configure.ac
#
#
# --------------------------------------------------------------
AC_DEFUN([PKG_CHECK_MODULES],
[AC_REQUIRE([PKG_PROG_PKG_CONFIG])dnl
AC_ARG_VAR([$1][_CFLAGS], [C compiler flags for $1, overriding pkg-config])dnl
AC_ARG_VAR([$1][_LIBS], [linker flags for $1, overriding pkg-config])dnl

pkg_failed=no
AC_MSG_CHECKING([for $1])

_PKG_CONFIG([$1][_CFLAGS], [cflags], [$2])
_PKG_CONFIG([$1][_LIBS], [libs], [$2])

m4_define([_PKG_TEXT], [Alternatively, you may set the environment variables $1[]_CFLAGS
and $1[]_LIBS to avoid the need to call pkg-config.
See the pkg-config man page for more details.])

if test $pkg_failed = yes; then
   	AC_MSG_RESULT([no])
        _PKG_SHORT_ERRORS_SUPPORTED
        if test $_pkg_short_errors_supported = yes; then
	        $1[]_PKG_ERRORS=`$PKG_CONFIG --short-errors --print-errors --cflags --libs "$2" 2>&1`
        else 
	        $1[]_PKG_ERRORS=`$PKG_CONFIG --print-errors --cflags --libs "$2" 2>&1`
        fi
	# Put the nasty error message in config.log where it belongs
	echo "$$1[]_PKG_ERRORS" >&AS_MESSAGE_LOG_FD

	m4_default([$4], [AC_MSG_ERROR(
[Package requirements ($2) were not met:

$$1_PKG_ERRORS

Consider adjusting the PKG_CONFIG_PATH environment variable if you
installed software in a non-standard prefix.

_PKG_TEXT])[]dnl
        ])
elif test $pkg_failed = untried; then
     	AC_MSG_RESULT([no])
	m4_default([$4], [AC_MSG_FAILURE(
[The pkg-config script could not be found or is too old.  Make sure it
is in your PATH or set the PKG_CONFIG environment variable to the full
path to pkg-config.

_PKG_TEXT

To get pkg-config, see <http://pkg-config.freedesktop.org/>.])[]dnl
        ])
else
	$1[]_CFLAGS=$pkg_cv_[]$1[]_CFLAGS
	$1[]_LIBS=$pkg_cv_[]$1[]_LIBS
        AC_MSG_RESULT([yes])
	$3
fi[]dnl
])# PKG_CHECK_MODULES

# PKG_INSTALLDIR(DIRECTORY)
# -------------------------
# Substitutes the variable pkgconfigdir as the location where a module
# should install pkg-config .pc files. By default the directory is
# $libdir/pkgconfig, but the default can be changed by passing
# DIRECTORY. The user can override through the --with-pkgconfigdir
# parameter.
AC_DEFUN([PKG_INSTALLDIR],
[m4_pushdef([pkg_default], [m4_default([$1], ['${libdir}/pkgconfig'])])
m4_pushdef([pkg_description],
    [pkg-config installation directory @<:@]pkg_default[@:>@])
AC_ARG_WITH([pkgconfigdir],
    [AS_HELP_STRING([--with-pkgconfigdir], pkg_description)],,
    [with_pkgconfigdir=]pkg_default)
AC_SUBST([pkgconfigdir], [$with_pkgconfigdir])
m4_popdef([pkg_default])
m4_popdef([pkg_description])
]) dnl PKG_INSTALLDIR


# PKG_NOARCH_INSTALLDIR(DIRECTORY)
# -------------------------
# Substitutes the variable noarch_pkgconfigdir as the location where a
# module should install arch-independent pkg-config .pc files. By
# default the directory is $datadir/pkgconfig, but the default can be
# changed by passing DIRECTORY. The user can override through the
# --with-noarch-pkgconfigdir parameter.
AC_DEFUN([PKG_NOARCH_INSTALLDIR],
[m4_pushdef([pkg_default], [m4_default([$1], ['${datadir}/pkgconfig'])])
m4_pushdef([pkg_description],
    [pkg-config arch-independent installation directory @<:@]pkg_default[@:>@])
AC_ARG_WITH([noarch-pkgconfigdir],
    [AS_HELP_STRING([--with-noarch-pkgconfigdir], pkg_description)],,
    [with_noarch_pkgconfigdir=]pkg_default)
AC_SUBST([noarch_pkgconfigdir], [$with_noarch_pkgconfigdir])
m4_popdef([pkg_default])
m4_popdef([pkg_description])
]) dnl PKG_NOARCH_INSTALLDIR


# PKG_CHECK_VAR(VARIABLE, MODULE, CONFIG-VARIABLE,
# [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# -------------------------------------------
# Retrieves the value of the pkg-config variable for the given module.
AC_DEFUN([PKG_CHECK_VAR],
[AC_REQUIRE([PKG_PROG_PKG_CONFIG])dnl
AC_ARG_VAR([$1], [value of $3 for $2, overriding pkg-config])dnl

_PKG_CONFIG([$1], [variable="][$3]["], [$2])
AS_VAR_COPY([$1], [pkg_cv_][$1])

AS_VAR_IF([$1], [""], [$5], [$4])dnl
])# PKG_CHECK_VAR

# PKG_WITH_MODULES(VARIABLE-PREFIX, MODULES,
#                  [ACTION-IF-FOUND],[ACTION-IF-NOT-FOUND],
#                  [DESCRIPTION], [DEFAULT])
#
# Prepare a "--with-" configure option using the lowercase [VARIABLE-PREFIX]
# name, merging the behaviour of AC_ARG_WITH and PKG_CHECK_MODULES in a single
# macro
#
# --------------------------------------------------------------
AC_DEFUN([PKG_WITH_MODULES],
[
m4_pushdef([with_arg], m4_tolower([$1]))

m4_pushdef([description],
           [m4_default([$5], [build with ]with_arg[ support])])

m4_pushdef([def_arg], [m4_default([$6], [auto])])
m4_pushdef([def_action_if_found], [AS_TR_SH([with_]with_arg)=yes])
m4_pushdef([def_action_if_not_found], [AS_TR_SH([with_]with_arg)=no])

m4_case(def_arg,
            [yes],[m4_pushdef([with_without], [--without-]with_arg)],
            [m4_pushdef([with_without],[--with-]with_arg)])

AC_ARG_WITH(with_arg,
     AS_HELP_STRING(with_without, description[ @<:@default=]def_arg[@:>@]),,
    [AS_TR_SH([with_]with_arg)=def_arg])

AS_CASE([$AS_TR_SH([with_]with_arg)],
            [yes],[PKG_CHECK_MODULES([$1],[$2],$3,$4)],
            [auto],[PKG_CHECK_MODULES([$1],[$2],
                                        [m4_n([def_action_if_found]) $3],
                                        [m4_n([def_action_if_not_found]) $4])])

m4_popdef([with_arg])
m4_popdef([description])
m4_popdef([def_arg])

]) dnl PKG_WITH_MODULES

# PKG_HAVE_WITH_MODULES(VARIABLE-PREFIX, MODULES,
#                       [DESCRIPTION], [DEFAULT])
#
# Convenience macro to trigger AM_CONDITIONAL after
# PKG_WITH_MODULES check.
#
# HAVE_[VARIABLE-PREFIX] is exported as make variable.
#
# --------------------------------------------------------------
AC_DEFUN([PKG_HAVE_WITH_MODULES],
[
PKG_WITH_MODULES([$1],[$2],,,[$3],[$4])

AM_CONDITIONAL([HAVE_][$1],
               [test "$AS_TR_SH([with_]m4_tolower([$1]))" = "yes"])
])

# PKG_HAVE_DEFINE_WITH_MODULES(VARIABLE-PREFIX, MODULES,
#                              [DESCRIPTION], [DEFAULT])
#
# Convenience macro to run AM_CONDITIONAL and AC_DEFINE after
# PKG_WITH_MODULES check.
#
# HAVE_[VARIABLE-PREFIX] is exported as make and preprocessor variable.
#
# --------------------------------------------------------------
AC_DEFUN([PKG_HAVE_DEFINE_WITH_MODULES],
[
PKG_HAVE_WITH_MODULES([$1],[$2],[$3],[$4])

AS_IF([test "$AS_TR_SH([with_]m4_tolower([$1]))" = "yes"],
        [AC_DEFINE([HAVE_][$1], 1, [Enable ]m4_tolower([$1])[ support])])
])

# Copyright (C) 2002-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_AUTOMAKE_VERSION(VERSION)
# ----------------------------
# Automake X.Y traces this macro to ensure aclocal.m4 has been
# generated from the m4 files accompanying Automake X.Y.
# (This private macro should not be called outside this file.)
AC_DEFUN([AM_AUTOMAKE_VERSION],
[am__api_version='1.15'
dnl Some users find AM_AUTOMAKE_VERSION and mistake it for a way to
dnl require some minimum version.  Point them to the right macro.
m4_if([$1], [1.15], [],
      [AC_FATAL([Do not call $0, use AM_INIT_AUTOMAKE([$1]).])])dnl
])

# _AM_AUTOCONF_VERSION(VERSION)
# -----------------------------
# aclocal traces this macro to find the Autoconf version.
# This is a private macro too.  Using m4_define simplifies
# the logic in aclocal, which can simply ignore this definition.
m4_define([_AM_AUTOCONF_VERSION], [])

# AM_SET_CURRENT_AUTOMAKE_VERSION
# -------------------------------
# Call AM_AUTOMAKE_VERSION and AM_AUTOMAKE_VERSION so they can be traced.
# This function is AC_REQUIREd by AM_INIT_AUTOMAKE.
AC_DEFUN([AM_SET_CURRENT_AUTOMAKE_VERSION],
[AM_AUTOMAKE_VERSION([1.15])dnl
m4_ifndef([AC_AUTOCONF_VERSION],
  [m4_copy([m4_PACKAGE_VERSION], [AC_AUTOCONF_VERSION])])dnl
_AM_AUTOCONF_VERSION(m4_defn([AC_AUTOCONF_VERSION]))])

# AM_AUX_DIR_EXPAND                                         -*- Autoconf -*-

# Copyright (C) 2001-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# For projects using AC_CONFIG_AUX_DIR([foo]), Autoconf sets
# $ac_aux_dir to '$srcdir/foo'.  In other projects, it is set to
# '$srcdir', '$srcdir/..', or '$srcdir/../..'.
#
# Of course, Automake must honor this variable whenever it calls a
# tool from the auxiliary directory.  The problem is that $srcdir (and
# therefore $ac_aux_dir as well) can be either absolute or relative,
# depending on how configure is run.  This is pretty annoying, since
# it makes $ac_aux_dir quite unusable in subdirectories: in the top
# source directory, any form will work fine, but in subdirectories a
# relative path needs to be adjusted first.
#
# $ac_aux_dir/missing
#    fails when called from a subdirectory if $ac_aux_dir is relative
# $top_srcdir/$ac_aux_dir/missing
#    fails if $ac_aux_dir is absolute,
#    fails when called from a subdirectory in a VPATH build with
#          a relative $ac_aux_dir
#
# The reason of the latter failure is that $top_srcdir and $ac_aux_dir
# are both prefixed by $srcdir.  In an in-source build this is usually
# harmless because $srcdir is '.', but things will broke when you
# start a VPATH build or use an absolute $srcdir.
#
# So we could use something similar to $top_srcdir/$ac_aux_dir/missing,
# iff we strip the leading $srcdir from $ac_aux_dir.  That would be:
#   am_aux_dir='\$(top_srcdir)/'`expr "$ac_aux_dir" : "$srcdir//*\(.*\)"`
# and then we would define $MISSING as
#   MISSING="\${SHELL} $am_aux_dir/missing"
# This will work as long as MISSING is not called from configure, because
# unfortunately $(top_srcdir) has no meaning in configure.
# However there are other variables, like CC, which are often used in
# configure, and could therefore not use this "fixed" $ac_aux_dir.
#
# Another solution, used here, is to always expand $ac_aux_dir to an
# absolute PATH.  The drawback is that using absolute paths prevent a
# configured tree to be moved without reconfiguration.

AC_DEFUN([AM_AUX_DIR_EXPAND],
[AC_REQUIRE([AC_CONFIG_AUX_DIR_DEFAULT])dnl
# Expand $ac_aux_dir to an absolute path.
am_aux_dir=`cd "$ac_aux_dir" && pwd`
])

# AM_CONDITIONAL                                            -*- Autoconf -*-

# Copyright (C) 1997-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_CONDITIONAL(NAME, SHELL-CONDITION)
# -------------------------------------
# Define a conditional.
AC_DEFUN([AM_CONDITIONAL],
[AC_PREREQ([2.52])dnl
 m4_if([$1], [TRUE],  [AC_FATAL([$0: invalid condition: $1])],
       [$1], [FALSE], [AC_FATAL([$0: invalid condition: $1])])dnl
AC_SUBST([$1_TRUE])dnl
AC_SUBST([$1_FALSE])dnl
_AM_SUBST_NOTMAKE([$1_TRUE])dnl
_AM_SUBST_NOTMAKE([$1_FALSE])dnl
m4_define([_AM_COND_VALUE_$1], [$2])dnl
if $2; then
  $1_TRUE=
  $1_FALSE='#'
else
  $1_TRUE='#'
  $1_FALSE=
fi
AC_CONFIG_COMMANDS_PRE(
[if test -z "${$1_TRUE}" && test -z "${$1_FALSE}"; then
  AC_MSG_ERROR([[conditional "$1" was never defined.
Usually this means the macro was only invoked conditionally.]])
fi])])

# Copyright (C) 1999-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.


# There are a few dirty hacks below to avoid letting 'AC_PROG_CC' be
# written in clear, in which case automake, when reading aclocal.m4,
# will think it sees a *use*, and therefore will trigger all it's
# C support machinery.  Also note that it means that autoscan, seeing
# CC etc. in the Makefile, will ask for an AC_PROG_CC use...


# _AM_DEPENDENCIES(NAME)
# ----------------------
# See how the compiler implements dependency checking.
# NAME is "CC", "CXX", "OBJC", "OBJCXX", "UPC", or "GJC".
# We try a few techniques and use that to set a single cache variable.
#
# We don't AC_REQUIRE the corresponding AC_PROG_CC since the latter was
# modified to invoke _AM_DEPENDENCIES(CC); we would have a circular
# dependency, and given that the user is not expected to run this macro,
# just rely on AC_PROG_CC.
AC_DEFUN([_AM_DEPENDENCIES],
[AC_REQUIRE([AM_SET_DEPDIR])dnl
AC_REQUIRE([AM_OUTPUT_DEPENDENCY_COMMANDS])dnl
AC_REQUIRE([AM_MAKE_INCLUDE])dnl
AC_REQUIRE([AM_DEP_TRACK])dnl

m4_if([$1], [CC],   [depcc="$CC"   am_compiler_list=],
      [$1], [CXX],  [depcc="$CXX"  am_compiler_list=],
      [$1], [OBJC], [depcc="$OBJC" am_compiler_list='gcc3 gcc'],
      [$1], [OBJCXX], [depcc="$OBJCXX" am_compiler_list='gcc3 gcc'],
      [$1], [UPC],  [depcc="$UPC"  am_compiler_list=],
      [$1], [GCJ],  [depcc="$GCJ"  am_compiler_list='gcc3 gcc'],
                    [depcc="$$1"   am_compiler_list=])

AC_CACHE_CHECK([dependency style of $depcc],
               [am_cv_$1_dependencies_compiler_type],
[if test -z "$AMDEP_TRUE" && test -f "$am_depcomp"; then
  # We make a subdir and do the tests there.  Otherwise we can end up
  # making bogus files that we don't know about and never remove.  For
  # instance it was reported that on HP-UX the gcc test will end up
  # making a dummy file named 'D' -- because '-MD' means "put the output
  # in D".
  rm -rf conftest.dir
  mkdir conftest.dir
  # Copy depcomp to subdir because otherwise we won't find it if we're
  # using a relative directory.
  cp "$am_depcomp" conftest.dir
  cd conftest.dir
  # We will build objects and dependencies in a subdirectory because
  # it helps to detect inapplicable dependency modes.  For instance
  # both Tru64's cc and ICC support -MD to output dependencies as a
  # side effect of compilation, but ICC will put the dependencies in
  # the current directory while Tru64 will put them in the object
  # directory.
  mkdir sub

  am_cv_$1_dependencies_compiler_type=none
  if test "$am_compiler_list" = ""; then
     am_compiler_list=`sed -n ['s/^#*\([a-zA-Z0-9]*\))$/\1/p'] < ./depcomp`
  fi
  am__universal=false
  m4_case([$1], [CC],
    [case " $depcc " in #(
     *\ -arch\ *\ -arch\ *) am__universal=true ;;
     esac],
    [CXX],
    [case " $depcc " in #(
     *\ -arch\ *\ -arch\ *) am__universal=true ;;
     esac])

  for depmode in $am_compiler_list; do
    # Setup a source with many dependencies, because some compilers
    # like to wrap large dependency lists on column 80 (with \), and
    # we should not choose a depcomp mode which is confused by this.
    #
    # We need to recreate these files for each test, as the compiler may
    # overwrite some of them when testing with obscure command lines.
    # This happens at least with the AIX C compiler.
    : > sub/conftest.c
    for i in 1 2 3 4 5 6; do
      echo '#include "conftst'$i'.h"' >> sub/conftest.c
      # Using ": > sub/conftst$i.h" creates only sub/conftst1.h with
      # Solaris 10 /bin/sh.
      echo '/* dummy */' > sub/conftst$i.h
    done
    echo "${am__include} ${am__quote}sub/conftest.Po${am__quote}" > confmf

    # We check with '-c' and '-o' for the sake of the "dashmstdout"
    # mode.  It turns out that the SunPro C++ compiler does not properly
    # handle '-M -o', and we need to detect this.  Also, some Intel
    # versions had trouble with output in subdirs.
    am__obj=sub/conftest.${OBJEXT-o}
    am__minus_obj="-o $am__obj"
    case $depmode in
    gcc)
      # This depmode causes a compiler race in universal mode.
      test "$am__universal" = false || continue
      ;;
    nosideeffect)
      # After this tag, mechanisms are not by side-effect, so they'll
      # only be used when explicitly requested.
      if test "x$enable_dependency_tracking" = xyes; then
	continue
      else
	break
      fi
      ;;
    msvc7 | msvc7msys | msvisualcpp | msvcmsys)
      # This compiler won't grok '-c -o', but also, the minuso test has
      # not run yet.  These depmodes are late enough in the game, and
      # so weak that their functioning should not be impacted.
      am__obj=conftest.${OBJEXT-o}
      am__minus_obj=
      ;;
    none) break ;;
    esac
    if depmode=$depmode \
       source=sub/conftest.c object=$am__obj \
       depfile=sub/conftest.Po tmpdepfile=sub/conftest.TPo \
       $SHELL ./depcomp $depcc -c $am__minus_obj sub/conftest.c \
         >/dev/null 2>conftest.err &&
       grep sub/conftst1.h sub/conftest.Po > /dev/null 2>&1 &&
       grep sub/conftst6.h sub/conftest.Po > /dev/null 2>&1 &&
       grep $am__obj sub/conftest.Po > /dev/null 2>&1 &&
       ${MAKE-make} -s -f confmf > /dev/null 2>&1; then
      # icc doesn't choke on unknown options, it will just issue warnings
      # or remarks (even with -Werror).  So we grep stderr for any message
      # that says an option was ignored or not supported.
      # When given -MP, icc 7.0 and 7.1 complain thusly:
      #   icc: Command line warning: ignoring option '-M'; no argument required
      # The diagnosis changed in icc 8.0:
      #   icc: Command line remark: option '-MP' not supported
      if (grep 'ignoring option' conftest.err ||
          grep 'not supported' conftest.err) >/dev/null 2>&1; then :; else
        am_cv_$1_dependencies_compiler_type=$depmode
        break
      fi
    fi
  done

  cd ..
  rm -rf conftest.dir
else
  am_cv_$1_dependencies_compiler_type=none
fi
])
AC_SUBST([$1DEPMODE], [depmode=$am_cv_$1_dependencies_compiler_type])
AM_CONDITIONAL([am__fastdep$1], [
  test "x$enable_dependency_tracking" != xno \
  && test "$am_cv_$1_dependencies_compiler_type" = gcc3])
])


# AM_SET_DEPDIR
# -------------
# Choose a directory name for dependency files.
# This macro is AC_REQUIREd in _AM_DEPENDENCIES.
AC_DEFUN([AM_SET_DEPDIR],
[AC_REQUIRE([AM_SET_LEADING_DOT])dnl
AC_SUBST([DEPDIR], ["${am__leading_dot}deps"])dnl
])


# AM_DEP_TRACK
# ------------
AC_DEFUN([AM_DEP_TRACK],
[AC_ARG_ENABLE([dependency-tracking], [dnl
AS_HELP_STRING(
  [--enable-dependency-tracking],
  [do not reject slow dependency extractors])
AS_HELP_STRING(
  [--disable-dependency-tracking],
  [speeds up one-time build])])
if test "x$enable_dependency_tracking" != xno; then
  am_depcomp="$ac_aux_dir/depcomp"
  AMDEPBACKSLASH='\'
  am__nodep='_no'
fi
AM_CONDITIONAL([AMDEP], [test "x$enable_dependency_tracking" != xno])
AC_SUBST([AMDEPBACKSLASH])dnl
_AM_SUBST_NOTMAKE([AMDEPBACKSLASH])dnl
AC_SUBST([am__nodep])dnl
_AM_SUBST_NOTMAKE([am__nodep])dnl
])

# Generate code to set up dependency tracking.              -*- Autoconf -*-

# Copyright (C) 1999-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.


# _AM_OUTPUT_DEPENDENCY_COMMANDS
# ------------------------------
AC_DEFUN([_AM_OUTPUT_DEPENDENCY_COMMANDS],
[{
  # Older Autoconf quotes --file arguments for eval, but not when files
  # are listed without --file.  Let's play safe and only enable the eval
  # if we detect the quoting.
  case $CONFIG_FILES in
  *\'*) eval set x "$CONFIG_FILES" ;;
  *)   set x $CONFIG_FILES ;;
  esac
  shift
  for mf
  do
    # Strip MF so we end up with the name of the file.
    mf=`echo "$mf" | sed -e 's/:.*$//'`
    # Check whether this is an Automake generated Makefile or not.
    # We used to match only the files named 'Makefile.in', but
    # some people rename them; so instead we look at the file content.
    # Grep'ing the first line is not enough: some people post-process
    # each Makefile.in and add a new line on top of each file to say so.
    # Grep'ing the whole file is not good either: AIX grep has a line
    # limit of 2048, but all sed's we know have understand at least 4000.
    if sed -n 's,^#.*generated by automake.*,X,p' "$mf" | grep X >/dev/null 2>&1; then
      dirpart=`AS_DIRNAME("$mf")`
    else
      continue
    fi
    # Extract the definition of DEPDIR, am__include, and am__quote
    # from the Makefile without running 'make'.
    DEPDIR=`sed -n 's/^DEPDIR = //p' < "$mf"`
    test -z "$DEPDIR" && continue
    am__include=`sed -n 's/^am__include = //p' < "$mf"`
    test -z "$am__include" && continue
    am__quote=`sed -n 's/^am__quote = //p' < "$mf"`
    # Find all dependency output files, they are included files with
    # $(DEPDIR) in their names.  We invoke sed twice because it is the
    # simplest approach to changing $(DEPDIR) to its actual value in the
    # expansion.
    for file in `sed -n "
      s/^$am__include $am__quote\(.*(DEPDIR).*\)$am__quote"'$/\1/p' <"$mf" | \
	 sed -e 's/\$(DEPDIR)/'"$DEPDIR"'/g'`; do
      # Make sure the directory exists.
      test -f "$dirpart/$file" && continue
      fdir=`AS_DIRNAME(["$file"])`
      AS_MKDIR_P([$dirpart/$fdir])
      # echo "creating $dirpart/$file"
      echo '# dummy' > "$dirpart/$file"
    done
  done
}
])# _AM_OUTPUT_DEPENDENCY_COMMANDS


# AM_OUTPUT_DEPENDENCY_COMMANDS
# -----------------------------
# This macro should only be invoked once -- use via AC_REQUIRE.
#
# This code is only required when automatic dependency tracking
# is enabled.  FIXME.  This creates each '.P' file that we will
# need in order to bootstrap the dependency handling code.
AC_DEFUN([AM_OUTPUT_DEPENDENCY_COMMANDS],
[AC_CONFIG_COMMANDS([depfiles],
     [test x"$AMDEP_TRUE" != x"" || _AM_OUTPUT_DEPENDENCY_COMMANDS],
     [AMDEP_TRUE="$AMDEP_TRUE" ac_aux_dir="$ac_aux_dir"])
])

# Do all the work for Automake.                             -*- Autoconf -*-

# Copyright (C) 1996-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This macro actually does too much.  Some checks are only needed if
# your package does certain things.  But this isn't really a big deal.

dnl Redefine AC_PROG_CC to automatically invoke _AM_PROG_CC_C_O.
m4_define([AC_PROG_CC],
m4_defn([AC_PROG_CC])
[_AM_PROG_CC_C_O
])

# AM_INIT_AUTOMAKE(PACKAGE, VERSION, [NO-DEFINE])
# AM_INIT_AUTOMAKE([OPTIONS])
# -----------------------------------------------
# The call with PACKAGE and VERSION arguments is the old style
# call (pre autoconf-2.50), which is being phased out.  PACKAGE
# and VERSION should now be passed to AC_INIT and removed from
# the call to AM_INIT_AUTOMAKE.
# We support both call styles for the transition.  After
# the next Automake release, Autoconf can make the AC_INIT
# arguments mandatory, and then we can depend on a new Autoconf
# release and drop the old call support.
AC_DEFUN([AM_INIT_AUTOMAKE],
[AC_PREREQ([2.65])dnl
dnl Autoconf wants to disallow AM_ names.  We explicitly allow
dnl the ones we care about.
m4_pattern_allow([^AM_[A-Z]+FLAGS$])dnl
AC_REQUIRE([AM_SET_CURRENT_AUTOMAKE_VERSION])dnl
AC_REQUIRE([AC_PROG_INSTALL])dnl
if test "`cd $srcdir && pwd`" != "`pwd`"; then
  # Use -I$(srcdir) only when $(srcdir) != ., so that make's output
  # is not polluted with repeated "-I."
  AC_SUBST([am__isrc], [' -I$(srcdir)'])_AM_SUBST_NOTMAKE([am__isrc])dnl
  # test to see if srcdir already configured
  if test -f $srcdir/config.status; then
    AC_MSG_ERROR([source directory already configured; run "make distclean" there first])
  fi
fi

# test whether we have cygpath
if test -z "$CYGPATH_W"; then
  if (cygpath --version) >/dev/null 2>/dev/null; then
    CYGPATH_W='cygpath -w'
  else
    CYGPATH_W=echo
  fi
fi
AC_SUBST([CYGPATH_W])

# Define the identity of the package.
dnl Distinguish between old-style and new-style calls.
m4_ifval([$2],
[AC_DIAGNOSE([obsolete],
             [$0: two- and three-arguments forms are deprecated.])
m4_ifval([$3], [_AM_SET_OPTION([no-define])])dnl
 AC_SUBST([PACKAGE], [$1])dnl
 AC_SUBST([VERSION], [$2])],
[_AM_SET_OPTIONS([$1])dnl
dnl Diagnose old-style AC_INIT with new-style AM_AUTOMAKE_INIT.
m4_if(
  m4_ifdef([AC_PACKAGE_NAME], [ok]):m4_ifdef([AC_PACKAGE_VERSION], [ok]),
  [ok:ok],,
  [m4_fatal([AC_INIT should be called with package and version arguments])])dnl
 AC_SUBST([PACKAGE], ['AC_PACKAGE_TARNAME'])dnl
 AC_SUBST([VERSION], ['AC_PACKAGE_VERSION'])])dnl

_AM_IF_OPTION([no-define],,
[AC_DEFINE_UNQUOTED([PACKAGE], ["$PACKAGE"], [Name of package])
 AC_DEFINE_UNQUOTED([VERSION], ["$VERSION"], [Version number of package])])dnl

# Some tools Automake needs.
AC_REQUIRE([AM_SANITY_CHECK])dnl
AC_REQUIRE([AC_ARG_PROGRAM])dnl
AM_MISSING_PROG([ACLOCAL], [aclocal-${am__api_version}])
AM_MISSING_PROG([AUTOCONF], [autoconf])
AM_MISSING_PROG([AUTOMAKE], [automake-${am__api_version}])
AM_MISSING_PROG([AUTOHEADER], [autoheader])
AM_MISSING_PROG([MAKEINFO], [makeinfo])
AC_REQUIRE([AM_PROG_INSTALL_SH])dnl
AC_REQUIRE([AM_PROG_INSTALL_STRIP])dnl
AC_REQUIRE([AC_PROG_MKDIR_P])dnl
# For better backward compatibility.  To be removed once Automake 1.9.x
# dies out for good.  For more background, see:
# <http://lists.gnu.org/archive/html/automake/2012-07/msg00001.html>
# <http://lists.gnu.org/archive/html/automake/2012-07/msg00014.html>
AC_SUBST([mkdir_p], ['$(MKDIR_P)'])
# We need awk for the "check" target (and possibly the TAP driver).  The
# system "awk" is bad on some platforms.
AC_REQUIRE([AC_PROG_AWK])dnl
AC_REQUIRE([AC_PROG_MAKE_SET])dnl
AC_REQUIRE([AM_SET_LEADING_DOT])dnl
_AM_IF_OPTION([tar-ustar], [_AM_PROG_TAR([ustar])],
	      [_AM_IF_OPTION([tar-pax], [_AM_PROG_TAR([pax])],
			     [_AM_PROG_TAR([v7])])])
_AM_IF_OPTION([no-dependencies],,
[AC_PROVIDE_IFELSE([AC_PROG_CC],
		  [_AM_DEPENDENCIES([CC])],
		  [m4_define([AC_PROG_CC],
			     m4_defn([AC_PROG_CC])[_AM_DEPENDENCIES([CC])])])dnl
AC_PROVIDE_IFELSE([AC_PROG_CXX],
		  [_AM_DEPENDENCIES([CXX])],
		  [m4_define([AC_PROG_CXX],
			     m4_defn([AC_PROG_CXX])[_AM_DEPENDENCIES([CXX])])])dnl
AC_PROVIDE_IFELSE([AC_PROG_OBJC],
		  [_AM_DEPENDENCIES([OBJC])],
		  [m4_define([AC_PROG_OBJC],
			     m4_defn([AC_PROG_OBJC])[_AM_DEPENDENCIES([OBJC])])])dnl
AC_PROVIDE_IFELSE([AC_PROG_OBJCXX],
		  [_AM_DEPENDENCIES([OBJCXX])],
		  [m4_define([AC_PROG_OBJCXX],
			     m4_defn([AC_PROG_OBJCXX])[_AM_DEPENDENCIES([OBJCXX])])])dnl
])
AC_REQUIRE([AM_SILENT_RULES])dnl
dnl The testsuite driver may need to know about EXEEXT, so add the
dnl 'am__EXEEXT' conditional if _AM_COMPILER_EXEEXT was seen.  This
dnl macro is hooked onto _AC_COMPILER_EXEEXT early, see below.
AC_CONFIG_COMMANDS_PRE(dnl
[m4_provide_if([_AM_COMPILER_EXEEXT],
  [AM_CONDITIONAL([am__EXEEXT], [test -n "$EXEEXT"])])])dnl

# POSIX will say in a future version that running "rm -f" with no argument
# is OK; and we want to be able to make that assumption in our Makefile
# recipes.  So use an aggressive probe to check that the usage we want is
# actually supported "in the wild" to an acceptable degree.
# See automake bug#10828.
# To make any issue more visible, cause the running configure to be aborted
# by default if the 'rm' program in use doesn't match our expectations; the
# user can still override this though.
if rm -f && rm -fr && rm -rf; then : OK; else
  cat >&2 <<'END'
Oops!

Your 'rm' program seems unable to run without file operands specified
on the command line, even when the '-f' option is present.  This is contrary
to the behaviour of most rm programs out there, and not conforming with
the upcoming POSIX standard: <http://austingroupbugs.net/view.php?id=542>

Please tell bug-automake@gnu.org about your system, including the value
of your $PATH and any error possibly output before this message.  This
can help us improve future automake versions.

END
  if test x"$ACCEPT_INFERIOR_RM_PROGRAM" = x"yes"; then
    echo 'Configuration will proceed anyway, since you have set the' >&2
    echo 'ACCEPT_INFERIOR_RM_PROGRAM variable to "yes"' >&2
    echo >&2
  else
    cat >&2 <<'END'
Aborting the configuration process, to ensure you take notice of the issue.

You can download and install GNU coreutils to get an 'rm' implementation
that behaves properly: <http://www.gnu.org/software/coreutils/>.

If you want to complete the configuration process using your problematic
'rm' anyway, export the environment variable ACCEPT_INFERIOR_RM_PROGRAM
to "yes", and re-run configure.

END
    AC_MSG_ERROR([Your 'rm' program is bad, sorry.])
  fi
fi
dnl The trailing newline in this macro's definition is deliberate, for
dnl backward compatibility and to allow trailing 'dnl'-style comments
dnl after the AM_INIT_AUTOMAKE invocation. See automake bug#16841.
])

dnl Hook into '_AC_COMPILER_EXEEXT' early to learn its expansion.  Do not
dnl add the conditional right here, as _AC_COMPILER_EXEEXT may be further
dnl mangled by Autoconf and run in a shell conditional statement.
m4_define([_AC_COMPILER_EXEEXT],
m4_defn([_AC_COMPILER_EXEEXT])[m4_provide([_AM_COMPILER_EXEEXT])])

# When config.status generates a header, we must update the stamp-h file.
# This file resides in the same directory as the config header
# that is generated.  The stamp files are numbered to have different names.

# Autoconf calls _AC_AM_CONFIG_HEADER_HOOK (when defined) in the
# loop where config.status creates the headers, so we can generate
# our stamp files there.
AC_DEFUN([_AC_AM_CONFIG_HEADER_HOOK],
[# Compute $1's index in $config_headers.
_am_arg=$1
_am_stamp_count=1
for _am_header in $config_headers :; do
  case $_am_header in
    $_am_arg | $_am_arg:* )
      break ;;
    * )
      _am_stamp_count=`expr $_am_stamp_count + 1` ;;
  esac
done
echo "timestamp for $_am_arg" >`AS_DIRNAME(["$_am_arg"])`/stamp-h[]$_am_stamp_count])

# Copyright (C) 2001-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_PROG_INSTALL_SH
# ------------------
# Define $install_sh.
AC_DEFUN([AM_PROG_INSTALL_SH],
[AC_REQUIRE([AM_AUX_DIR_EXPAND])dnl
if test x"${install_sh+set}" != xset; then
  case $am_aux_dir in
  *\ * | *\	*)
    install_sh="\${SHELL} '$am_aux_dir/install-sh'" ;;
  *)
    install_sh="\${SHELL} $am_aux_dir/install-sh"
  esac
fi
AC_SUBST([install_sh])])

# Copyright (C) 2003-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# Check whether the underlying file-system supports filenames
# with a leading dot.  For instance MS-DOS doesn't.
AC_DEFUN([AM_SET_LEADING_DOT],
[rm -rf .tst 2>/dev/null
mkdir .tst 2>/dev/null
if test -d .tst; then
  am__leading_dot=.
else
  am__leading_dot=_
fi
rmdir .tst 2>/dev/null
AC_SUBST([am__leading_dot])])

# Check to see how 'make' treats includes.	            -*- Autoconf -*-

# Copyright (C) 2001-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_MAKE_INCLUDE()
# -----------------
# Check to see how make treats includes.
AC_DEFUN([AM_MAKE_INCLUDE],
[am_make=${MAKE-make}
cat > confinc << 'END'
am__doit:
	@echo this is the am__doit target
.PHONY: am__doit
END
# If we don't find an include directive, just comment out the code.
AC_MSG_CHECKING([for style of include used by $am_make])
am__include="#"
am__quote=
_am_result=none
# First try GNU make style include.
echo "include confinc" > confmf
# Ignore all kinds of additional output from 'make'.
case `$am_make -s -f confmf 2> /dev/null` in #(
*the\ am__doit\ target*)
  am__include=include
  am__quote=
  _am_result=GNU
  ;;
esac
# Now try BSD make style include.
if test "$am__include" = "#"; then
   echo '.include "confinc"' > confmf
   case `$am_make -s -f confmf 2> /dev/null` in #(
   *the\ am__doit\ target*)
     am__include=.include
     am__quote="\""
     _am_result=BSD
     ;;
   esac
fi
AC_SUBST([am__include])
AC_SUBST([am__quote])
AC_MSG_RESULT([$_am_result])
rm -f confinc confmf
])

# Fake the existence of programs that GNU maintainers use.  -*- Autoconf -*-

# Copyright (C) 1997-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_MISSING_PROG(NAME, PROGRAM)
# ------------------------------
AC_DEFUN([AM_MISSING_PROG],
[AC_REQUIRE([AM_MISSING_HAS_RUN])
$1=${$1-"${am_missing_run}$2"}
AC_SUBST($1)])

# AM_MISSING_HAS_RUN
# ------------------
# Define MISSING if not defined so far and test if it is modern enough.
# If it is, set am_missing_run to use it, otherwise, to nothing.
AC_DEFUN([AM_MISSING_HAS_RUN],
[AC_REQUIRE([AM_AUX_DIR_EXPAND])dnl
AC_REQUIRE_AUX_FILE([missing])dnl
if test x"${MISSING+set}" != xset; then
  case $am_aux_dir in
  *\ * | *\	*)
    MISSING="\${SHELL} \"$am_aux_dir/missing\"" ;;
  *)
    MISSING="\${SHELL} $am_aux_dir/missing" ;;
  esac
fi
# Use eval to expand $SHELL
if eval "$MISSING --is-lightweight"; then
  am_missing_run="$MISSING "
else
  am_missing_run=
  AC_MSG_WARN(['missing' script is too old or missing])
fi
])

# Helper functions for option handling.                     -*- Autoconf -*-

# Copyright (C) 2001-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# _AM_MANGLE_OPTION(NAME)
# -----------------------
AC_DEFUN([_AM_MANGLE_OPTION],
[[_AM_OPTION_]m4_bpatsubst($1, [[^a-zA-Z0-9_]], [_])])

# _AM_SET_OPTION(NAME)
# --------------------
# Set option NAME.  Presently that only means defining a flag for this option.
AC_DEFUN([_AM_SET_OPTION],
[m4_define(_AM_MANGLE_OPTION([$1]), [1])])

# _AM_SET_OPTIONS(OPTIONS)
# ------------------------
# OPTIONS is a space-separated list of Automake options.
AC_DEFUN([_AM_SET_OPTIONS],
[m4_foreach_w([_AM_Option], [$1], [_AM_SET_OPTION(_AM_Option)])])

# _AM_IF_OPTION(OPTION, IF-SET, [IF-NOT-SET])
# -------------------------------------------
# Execute IF-SET if OPTION is set, IF-NOT-SET otherwise.
AC_DEFUN([_AM_IF_OPTION],
[m4_ifset(_AM_MANGLE_OPTION([$1]), [$2], [$3])])

# Copyright (C) 1999-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# _AM_PROG_CC_C_O
# ---------------
# Like AC_PROG_CC_C_O, but changed for automake.  We rewrite AC_PROG_CC
# to automatically call this.
AC_DEFUN([_AM_PROG_CC_C_O],
[AC_REQUIRE([AM_AUX_DIR_EXPAND])dnl
AC_REQUIRE_AUX_FILE([compile])dnl
AC_LANG_PUSH([C])dnl
AC_CACHE_CHECK(
  [whether $CC understands -c and -o together],
  [am_cv_prog_cc_c_o],
  [AC_LANG_CONFTEST([AC_LANG_PROGRAM([])])
  # Make sure it works both with $CC and with simple cc.
  # Following AC_PROG_CC_C_O, we do the test twice because some
  # compilers refuse to overwrite an existing .o file with -o,
  # though they will create one.
  am_cv_prog_cc_c_o=yes
  for am_i in 1 2; do
    if AM_RUN_LOG([$CC -c conftest.$ac_ext -o conftest2.$ac_objext]) \
         && test -f conftest2.$ac_objext; then
      : OK
    else
      am_cv_prog_cc_c_o=no
      break
    fi
  done
  rm -f core conftest*
  unset am_i])
if test "$am_cv_prog_cc_c_o" != yes; then
   # Losing compiler, so override with the script.
   # FIXME: It is wrong to rewrite CC.
   # But if we don't then we get into trouble of one sort or another.
   # A longer-term fix would be to have automake use am__CC in this case,
   # and then we could set am__CC="\$(top_srcdir)/compile \$(CC)"
   CC="$am_aux_dir/compile $CC"
fi
AC_LANG_POP([C])])

# For backward compatibility.
AC_DEFUN_ONCE([AM_PROG_CC_C_O], [AC_REQUIRE([AC_PROG_CC])])

# Copyright (C) 2001-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_RUN_LOG(COMMAND)
# -------------------
# Run COMMAND, save the exit status in ac_status, and log it.
# (This has been adapted from Autoconf's _AC_RUN_LOG macro.)
AC_DEFUN([AM_RUN_LOG],
[{ echo "$as_me:$LINENO: $1" >&AS_MESSAGE_LOG_FD
   ($1) >&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD
   ac_status=$?
   echo "$as_me:$LINENO: \$? = $ac_status" >&AS_MESSAGE_LOG_FD
   (exit $ac_status); }])

# Check to make sure that the build environment is sane.    -*- Autoconf -*-

# Copyright (C) 1996-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_SANITY_CHECK
# ---------------
AC_DEFUN([AM_SANITY_CHECK],
[AC_MSG_CHECKING([whether build environment is sane])
# Reject unsafe characters in $srcdir or the absolute working directory
# name.  Accept space and tab only in the latter.
am_lf='
'
case `pwd` in
  *[[\\\"\#\$\&\'\`$am_lf]]*)
    AC_MSG_ERROR([unsafe absolute working directory name]);;
esac
case $srcdir in
  *[[\\\"\#\$\&\'\`$am_lf\ \	]]*)
    AC_MSG_ERROR([unsafe srcdir value: '$srcdir']);;
esac

# Do 'set' in a subshell so we don't clobber the current shell's
# arguments.  Must try -L first in case configure is actually a
# symlink; some systems play weird games with the mod time of symlinks
# (eg FreeBSD returns the mod time of the symlink's containing
# directory).
if (
   am_has_slept=no
   for am_try in 1 2; do
     echo "timestamp, slept: $am_has_slept" > conftest.file
     set X `ls -Lt "$srcdir/configure" conftest.file 2> /dev/null`
     if test "$[*]" = "X"; then
	# -L didn't work.
	set X `ls -t "$srcdir/configure" conftest.file`
     fi
     if test "$[*]" != "X $srcdir/configure conftest.file" \
	&& test "$[*]" != "X conftest.file $srcdir/configure"; then

	# If neither matched, then we have a broken ls.  This can happen
	# if, for instance, CONFIG_SHELL is bash and it inherits a
	# broken ls alias from the environment.  This has actually
	# happened.  Such a system could not be considered "sane".
	AC_MSG_ERROR([ls -t appears to fail.  Make sure there is not a broken
  alias in your environment])
     fi
     if test "$[2]" = conftest.file || test $am_try -eq 2; then
       break
     fi
     # Just in case.
     sleep 1
     am_has_slept=yes
   done
   test "$[2]" = conftest.file
   )
then
   # Ok.
   :
else
   AC_MSG_ERROR([newly created file is older than distributed files!
Check your system clock])
fi
AC_MSG_RESULT([yes])
# If we didn't sleep, we still need to ensure time stamps of config.status and
# generated files are strictly newer.
am_sleep_pid=
if grep 'slept: no' conftest.file >/dev/null 2>&1; then
  ( sleep 1 ) &
  am_sleep_pid=$!
fi
AC_CONFIG_COMMANDS_PRE(
  [AC_MSG_CHECKING([that generated files are newer than configure])
   if test -n "$am_sleep_pid"; then
     # Hide warnings about reused PIDs.
     wait $am_sleep_pid 2>/dev/null
   fi
   AC_MSG_RESULT([done])])
rm -f conftest.file
])

# Copyright (C) 2009-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_SILENT_RULES([DEFAULT])
# --------------------------
# Enable less verbose build rules; with the default set to DEFAULT
# ("yes" being less verbose, "no" or empty being verbose).
AC_DEFUN([AM_SILENT_RULES],
[AC_ARG_ENABLE([silent-rules], [dnl
AS_HELP_STRING(
  [--enable-silent-rules],
  [less verbose build output (undo: "make V=1")])
AS_HELP_STRING(
  [--disable-silent-rules],
  [verbose build output (undo: "make V=0")])dnl
])
case $enable_silent_rules in @%:@ (((
  yes) AM_DEFAULT_VERBOSITY=0;;
   no) AM_DEFAULT_VERBOSITY=1;;
    *) AM_DEFAULT_VERBOSITY=m4_if([$1], [yes], [0], [1]);;
esac
dnl
dnl A few 'make' implementations (e.g., NonStop OS and NextStep)
dnl do not support nested variable expansions.
dnl See automake bug#9928 and bug#10237.
am_make=${MAKE-make}
AC_CACHE_CHECK([whether $am_make supports nested variables],
   [am_cv_make_support_nested_variables],
   [if AS_ECHO([['TRUE=$(BAR$(V))
BAR0=false
BAR1=true
V=1
am__doit:
	@$(TRUE)
.PHONY: am__doit']]) | $am_make -f - >/dev/null 2>&1; then
  am_cv_make_support_nested_variables=yes
else
  am_cv_make_support_nested_variables=no
fi])
if test $am_cv_make_support_nested_variables = yes; then
  dnl Using '$V' instead of '$(V)' breaks IRIX make.
  AM_V='$(V)'
  AM_DEFAULT_V='$(AM_DEFAULT_VERBOSITY)'
else
  AM_V=$AM_DEFAULT_VERBOSITY
  AM_DEFAULT_V=$AM_DEFAULT_VERBOSITY
fi
AC_SUBST([AM_V])dnl
AM_SUBST_NOTMAKE([AM_V])dnl
AC_SUBST([AM_DEFAULT_V])dnl
AM_SUBST_NOTMAKE([AM_DEFAULT_V])dnl
AC_SUBST([AM_DEFAULT_VERBOSITY])dnl
AM_BACKSLASH='\'
AC_SUBST([AM_BACKSLASH])dnl
_AM_SUBST_NOTMAKE([AM_BACKSLASH])dnl
])

# Copyright (C) 2001-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_PROG_INSTALL_STRIP
# ---------------------
# One issue with vendor 'install' (even GNU) is that you can't
# specify the program used to strip binaries.  This is especially
# annoying in cross-compiling environments, where the build's strip
# is unlikely to handle the host's binaries.
# Fortunately install-sh will honor a STRIPPROG variable, so we
# always use install-sh in "make install-strip", and initialize
# STRIPPROG with the value of the STRIP variable (set by the user).
AC_DEFUN([AM_PROG_INSTALL_STRIP],
[AC_REQUIRE([AM_PROG_INSTALL_SH])dnl
# Installed binaries are usually stripped using 'strip' when the user
# run "make install-strip".  However 'strip' might not be the right
# tool to use in cross-compilation environments, therefore Automake
# will honor the 'STRIP' environment variable to overrule this program.
dnl Don't test for $cross_compiling = yes, because it might be 'maybe'.
if test "$cross_compiling" != no; then
  AC_CHECK_TOOL([STRIP], [strip], :)
fi
INSTALL_STRIP_PROGRAM="\$(install_sh) -c -s"
AC_SUBST([INSTALL_STRIP_PROGRAM])])

# Copyright (C) 2006-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# _AM_SUBST_NOTMAKE(VARIABLE)
# ---------------------------
# Prevent Automake from outputting VARIABLE = @VARIABLE@ in Makefile.in.
# This macro is traced by Automake.
AC_DEFUN([_AM_SUBST_NOTMAKE])

# AM_SUBST_NOTMAKE(VARIABLE)
# --------------------------
# Public sister of _AM_SUBST_NOTMAKE.
AC_DEFUN([AM_SUBST_NOTMAKE], [_AM_SUBST_NOTMAKE($@)])

# Check how to create a tarball.                            -*- Autoconf -*-

# Copyright (C) 2004-2014 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# _AM_PROG_TAR(FORMAT)
# --------------------
# Check how to create a tarball in format FORMAT.
# FORMAT should be one of 'v7', 'ustar', or 'pax'.
#
# Substitute a variable $(am__tar) that is a command
# writing to stdout a FORMAT-tarball containing the directory
# $tardir.
#     tardir=directory && $(am__tar) > result.tar
#
# Substitute a variable $(am__untar) that extract such
# a tarball read from stdin.
#     $(am__untar) < result.tar
#
AC_DEFUN([_AM_PROG_TAR],
[# Always define AMTAR for backward compatibility.  Yes, it's still used
# in the wild :-(  We should find a proper way to deprecate it ...
AC_SUBST([AMTAR], ['$${TAR-tar}'])

# We'll loop over all known methods to create a tar archive until one works.
_am_tools='gnutar m4_if([$1], [ustar], [plaintar]) pax cpio none'

m4_if([$1], [v7],
  [am__tar='$${TAR-tar} chof - "$$tardir"' am__untar='$${TAR-tar} xf -'],

  [m4_case([$1],
    [ustar],
     [# The POSIX 1988 'ustar' format is defined with fixed-size fields.
      # There is notably a 21 bits limit for the UID and the GID.  In fact,
      # the 'pax' utility can hang on bigger UID/GID (see automake bug#8343
      # and bug#13588).
      am_max_uid=2097151 # 2^21 - 1
      am_max_gid=$am_max_uid
      # The $UID and $GID variables are not portable, so we need to resort
      # to the POSIX-mandated id(1) utility.  Errors in the 'id' calls
      # below are definitely unexpected, so allow the users to see them
      # (that is, avoid stderr redirection).
      am_uid=`id -u || echo unknown`
      am_gid=`id -g || echo unknown`
      AC_MSG_CHECKING([whether UID '$am_uid' is supported by ustar format])
      if test $am_uid -le $am_max_uid; then
         AC_MSG_RESULT([yes])
      else
         AC_MSG_RESULT([no])
         _am_tools=none
      fi
      AC_MSG_CHECKING([whether GID '$am_gid' is supported by ustar format])
      if test $am_gid -le $am_max_gid; then
         AC_MSG_RESULT([yes])
      else
        AC_MSG_RESULT([no])
        _am_tools=none
      fi],

  [pax],
    [],

  [m4_fatal([Unknown tar format])])

  AC_MSG_CHECKING([how to create a $1 tar archive])

  # Go ahead even if we have the value already cached.  We do so because we
  # need to set the values for the 'am__tar' and 'am__untar' variables.
  _am_tools=${am_cv_prog_tar_$1-$_am_tools}

  for _am_tool in $_am_tools; do
    case $_am_tool in
    gnutar)
      for _am_tar in tar gnutar gtar; do
        AM_RUN_LOG([$_am_tar --version]) && break
      done
      am__tar="$_am_tar --format=m4_if([$1], [pax], [posix], [$1]) -chf - "'"$$tardir"'
      am__tar_="$_am_tar --format=m4_if([$1], [pax], [posix], [$1]) -chf - "'"$tardir"'
      am__untar="$_am_tar -xf -"
      ;;
    plaintar)
      # Must skip GNU tar: if it does not support --format= it doesn't create
      # ustar tarball either.
      (tar --version) >/dev/null 2>&1 && continue
      am__tar='tar chf - "$$tardir"'
      am__tar_='tar chf - "$tardir"'
      am__untar='tar xf -'
      ;;
    pax)
      am__tar='pax -L -x $1 -w "$$tardir"'
      am__tar_='pax -L -x $1 -w "$tardir"'
      am__untar='pax -r'
      ;;
    cpio)
      am__tar='find "$$tardir" -print | cpio -o -H $1 -L'
      am__tar_='find "$tardir" -print | cpio -o -H $1 -L'
      am__untar='cpio -i -H $1 -d'
      ;;
    none)
      am__tar=false
      am__tar_=false
      am__untar=false
      ;;
    esac

    # If the value was cached, stop now.  We just wanted to have am__tar
    # and am__untar set.
    test -n "${am_cv_prog_tar_$1}" && break

    # tar/untar a dummy directory, and stop if the command works.
    rm -rf conftest.dir
    mkdir conftest.dir
    echo GrepMe > conftest.dir/file
    AM_RUN_LOG([tardir=conftest.dir && eval $am__tar_ >conftest.tar])
    rm -rf conftest.dir
    if test -s conftest.tar; then
      AM_RUN_LOG([$am__untar <conftest.tar])
      AM_RUN_LOG([cat conftest.dir/file])
      grep GrepMe conftest.dir/file >/dev/null 2>&1 && break
    fi
  done
  rm -rf conftest.dir

  AC_CACHE_VAL([am_cv_prog_tar_$1], [am_cv_prog_tar_$1=$_am_tool])
  AC_MSG_RESULT([$am_cv_prog_tar_$1])])

AC_SUBST([am__tar])
AC_SUBST([am__untar])
]) # _AM_PROG_TAR

