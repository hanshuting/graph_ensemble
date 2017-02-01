#-------------------------------------------------------------------------------
#
#                       configure.ac for l1_logreg
#
#-------------------------------------------------------------------------------

AC_PREREQ(2.59)
AC_INIT([l1_logreg],[0.8.2],[deneb1@stanford.edu])
AC_COPYRIGHT([
Copyright (C) 2006   Kwangmoo Koh  <deneb1@stanford.edu>, SeungJean Kim <sjkim@stanford.edu>, Stephen Boyd  <boyd@stanford.edu>])

AC_CONFIG_HEADER([src_c/config.h])
AC_CONFIG_SRCDIR([src_c/l1_logreg.c])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([])

# Checks for static linking option.

AC_MSG_CHECKING(whether to link statically)
AC_ARG_WITH(static-linking,
    [AC_HELP_STRING([--with-static-linking], [link binaries statically (default no)])],
    [with_static_linking=$withval],
    [with_static_linking=no])
AC_MSG_RESULT($with_static_linking)

if test "x$with_static_linking" = xyes; then
    LDFLAGS="$LDFLAGS -static"
fi

# Checks for fortran compiler usage.
# Necessary for naming convention, but if skipped, use blind linking.

AC_MSG_CHECKING(whether to use fortran compiler)
AC_ARG_WITH(fortran,
    [AC_HELP_STRING([--with-fortran], [use fortran compiler (default no)])],
    [with_fortran=$withval],
    [with_fortran=no])
AC_MSG_RESULT($with_fortran)

AC_CANONICAL_HOST
case $host in
    *i*86*-linux*)
        try_cc="gcc icc pgcc cc"
        try_ff="g77 f77 gfortran ifort ifc pgf90 f90"
        ;;
    *ia*64*-linux*)
        try_cc="gcc icc cc"
        try_ff="g77 f77 gfortran ifort ifc pgf90 f90"
        ;;
    *x86*64*-linux*)
        try_cc="gcc icc pgcc cc"
        try_ff="g77 f77 gfortran pgf90 ifort ifc f90"
        ;;
    *alpha*-linux*)
        try_cc="ccc cc gcc"
        try_ff="g77 f77 gfortran fort f90"
        ;;
    *dec*-osf*)
        try_cc="ccc cc gcc"
        try_ff="f77 fort f90"
        ;;
    *ibm*-aix*)
        try_cc="cc gcc"
        try_ff="f77 xlf90 f90"
        ;;
    *sgi*-irix*)
        try_cc="cc gcc"
        try_ff="f77 f90"
        ;;
    *hppa*-hpux*)
        try_cc="cc gcc"
        try_ff="f77 f90"
        ;;
    *sparc*-solaris*)
        try_cc="cc gcc"
        try_ff="f77 f90"
        ;;
    *)
        AC_MSG_WARN([Unknown architecture! This might not work...])
        try_cc="cc gcc ccc cl"
        try_ff="f77 f90 xlf90 ifc pgf90 epcf90 f95 xlf95 lf95 fort"
        build_type="unknown"
        ;;
esac


# Checks for programs.
AC_LANG([C])
AC_PROG_CC([$try_cc])
AX_COMPILER_VENDOR
AX_CC_MAXOPT


if test "x$with_fortran" = xyes ; then
    AC_PROG_F77([$try_ff])
    AC_F77_WRAPPERS
    AC_F77_FUNC([sgemm])
    AC_F77_FUNC([dsyev])
else
    F77=""
    FLIBS=""
    sgemm="sgemm_"
    dsyev="dsyev_"
fi

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])
AC_CHECK_LIB([g2c], [main])
AC_CHECK_LIB([gfortran], [main])

# Checks for BLAS and LAPACK libraries:

sinclude(config/acx_blas.m4)
sinclude(config/acx_lapack.m4)

ACX_BLAS($sgemm, [], [AC_MSG_ERROR([BLAS was not found!])])
ACX_LAPACK($dsyev,[], [AC_MSG_ERROR([LAPACK was not found!])])

LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

# Checks for headers.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h malloc.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_HEADER_TIME

# Checks for functions.
AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([sqrt memset memcpy])
AC_CHECK_FUNCS([log1p gettimeofday])

# Checks for external programs.
# doxygen: a documentation tool
AC_CHECK_PROG([doxygen_ok],[doxygen],[yes],[no])
AM_CONDITIONAL([DOXYGEN_OK],[test "x$doxygen_ok" = xyes])


AC_CONFIG_FILES([Makefile               \
                 src_c/Makefile         \
                 doc/Makefile           \
                 util_m/Makefile        \
                 examples/Makefile      ])

AC_OUTPUT

AC_MSG_RESULT([

------------------------------------------------------------------------------
$PACKAGE-$VERSION configuration:
------------------------------------------------------------------------------

System:
  - CPU ............... : ${arch}
  - host type ......... : ${host_cpu}
  - vendor ............ : ${host_vendor}
  - OS ................ : ${host_os}

External libs:
  - BLAS .............. : ${acx_blas_ok}
  - LAPACK ............ : ${acx_lapack_ok}
                        : ${acx_blas_name}

Directories:
  - prefix ............ : ${prefix}
  - exec_prefix ....... : ${exec_prefix}
  - includedir ........ : ${includedir}


Compiler/linker flags libs:
  - CC ................ : ${CC}
  - F77 ............... : ${F77}
  - CFLAGS ............ : ${CFLAGS}
  - LDFLAGS ........... : ${LDFLAGS}
  - LIBS .............. : ${LIBS}
  - FLIBS ............. : ${FLIBS}

Build options:

  - with-static-linking : ${with_static_linking}
  - with-fortran        : ${with_fortran}

------------------------------------------------------------------------------


Type 'make' to build binaries.

])

