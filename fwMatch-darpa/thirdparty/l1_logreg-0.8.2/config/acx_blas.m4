dnl @synopsis ACX_BLAS(FUNC_NAME[, ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the BLAS
dnl linear-algebra interface (see http://www.netlib.org/blas/).
dnl On success, it sets the BLAS_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with BLAS, you should link with:
dnl
dnl 	$BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order.  FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS),
dnl and is sometimes necessary in order to link with F77 libraries.
dnl Users will also need to use AC_F77_DUMMY_MAIN (see the autoconf
dnl manual), for the same reason.
dnl
dnl Many libraries are searched for, from MKL, ACML, ATLAS, CXML to ESSL.
dnl The user may also use --with-blas=<lib> in order to use some
dnl specific BLAS library <lib>.  In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as
dnl was used to compile the BLAS library.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a BLAS
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_BLAS.
dnl
dnl This macro requires autoconf 2.50 or later.
dnl
dnl @version $Id: acx_f77blas.m4,v 1.2 2005/06/10 12:08:55 juselius Exp $
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl Modified by Jonas Juselius <jonas@iki.fi>
dnl
dnl 2006-10-17
dnl Modified by Kwangmoo Koh <deneb1@stanford.edu>

AC_DEFUN([ACX_BLAS], [
AC_PREREQ(2.59)

acx_blas_name=none
acx_blas_ok=no
acx_blas_static=no
acx_blas_save_LIBS="$LIBS"
acx_blas_save_LDFLAGS="$LDFLAGS"
acx_blas_libs=""
acx_blas_dir=""

#   Checks for static linking.
for tmpval in $LDFLAGS ; do
    if test "x$tmpval" = "x-static" ; then
        acx_blas_static=yes
    fi
done

AC_ARG_WITH(blas,
    [AC_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>])])
case $with_blas in
    yes | "") ;;
    no) acx_blas_ok=disable ;;
    -l* | */* | *.a | *.so | *.so.* | *.o) acx_blas_libs="$with_blas" ;;
    *) acx_blas_libs="-l$with_blas" ;;
esac

AC_ARG_WITH(blas_dir,
    [AC_HELP_STRING([--with-blas-dir=<dir>], [look for BLAS library in <dir>])])

case $with_blas_dir in
    yes | no | "") ;;
    -L*) LDFLAGS="$with_blas_dir $LDFLAGS" 
         acx_blas_dir="$with_blas_dir" ;;
    *)   LDFLAGS="-L$with_blas_dir $LDFLAGS" 
         acx_blas_dir="-L$with_blas_dir" ;;
esac

# If --with-blas is defined, then look for THIS AND ONLY THIS blas lib
if test $acx_blas_ok = no; then
case $with_blas in
    ""|yes) ;;
    *) save_LIBS="$LIBS"; LIBS="$acx_blas_libs $LIBS"
        AC_MSG_CHECKING([for $1 in $acx_blas_libs])
        AC_LINK_IFELSE([AC_LANG_CALL([], [$1])], 
            [acx_blas_ok=yes; acx_blas_name="user-defined"],
            [acx_blas_ok=none; acx_blas_name="no-user-defined"])
        AC_MSG_RESULT($acx_blas_ok)
        LIBS="$save_LIBS"
    ;;
esac
fi

# Checks for libraries.
AC_CHECK_LIB([m],[sqrt],[acx_blas_m=yes],[acx_blas_m=no])
AC_CHECK_LIB([g2c],[main],[acx_blas_g2c=yes],[acx_blas_g2c=no])
AC_CHECK_LIB([gfortran],[main],[acx_blas_gfortran=yes],[acx_blas_gfortran=no])

# Check BLAS_LIBS environment variable
if test $acx_blas_ok = no; then
if test "x$BLAS_LIBS" != x; then
    save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
    AC_MSG_CHECKING([for $1 in $BLAS_LIBS])
    AC_LINK_IFELSE([AC_LANG_CALL([], [$1])], 
        [acx_blas_ok=yes; acx_blas_name=environment; acx_blas_libs=$BLAS_LIBS])
    AC_MSG_RESULT($acx_blas_ok)
    LIBS="$save_LIBS"
fi
fi

# BLAS linked to by default?  (happens on some supercomputers)
if test $acx_blas_ok = no; then
    AC_MSG_CHECKING([for builtin $1])
    AC_LINK_IFELSE([AC_LANG_CALL([], [$1])], 
        [acx_blas_ok=yes; acx_blas_name="default"])
    AC_MSG_RESULT($acx_blas_ok)
fi

#
#   Intel Math Kernel Library (MKL)
#

if test "x$acx_blas_static" != xyes; then 
#   Intel MKL dynamic
AC_MSG_NOTICE([----- Intel MKL dynamic])

    if test $acx_blas_ok = no; then
    if test $acx_blas_m = yes; then
        AC_CHECK_LIB(mkl, $1, 
            [acx_blas_ok=yes; acx_blas_name=mkl; 
            acx_blas_libs="-lmkl -lguide -lpthread -lm"],
            [],[ -lguide -lpthread -lm])
    fi
    fi

else
#   Intel MKL static
AC_MSG_NOTICE([----- Intel MKL static])

    if test $acx_blas_ok = no; then
    if test $acx_blas_m = yes; then
        AC_CHECK_LIB(mkl_ia32, $1, 
            [acx_blas_ok=yes; acx_blas_name=mkl_ia32; 
                acx_blas_libs="-lmkl_ia32 -lguide -lpthread -lm"],
            [],[ -lguide -lpthread -lm])
    fi
    fi
    if test $acx_blas_ok = no; then
    if test $acx_blas_m = yes; then
        AC_CHECK_LIB(mkl_em64t, $1, 
            [acx_blas_ok=yes; acx_blas_name=mkl_em64t;
                 acx_blas_libs="-lmkl_em64t -lguide -lpthread -lm"],
            [],[ -lguide -lpthread -lm])
    fi
    fi
    if test $acx_blas_ok = no; then
    if test $acx_blas_m = yes; then
        AC_CHECK_LIB(mkl_ipf, $1, 
            [acx_blas_ok=yes; acx_blas_name=mkl_ipf;
                acx_blas_libs="-lmkl_ipf -lguide -lpthread -lm"],
            [],[ -lguide -lpthread -lm])
    fi
    fi
fi

#
#   AMD Core Math Library (ACML)
#
if test "x$acx_blas_static" != xyes; then 
#   ACML dynamic
    if test $acx_blas_ok = no; then
    AC_MSG_NOTICE([----- ACML dynamic])
    AC_CHECK_LIB(acml, $1,
        [acx_blas_ok=yes; acx_blas_name=acml;
            acx_blas_libs="-lacml"],
        [],[])
    fi
else
#   ACML static
    if test $acx_blas_ok = no; then
    if test $acx_blas_g2c = yes; then
    AC_MSG_NOTICE([----- ACML static])
    AC_CHECK_LIB(acml, $1,
        [acx_blas_ok=yes; acx_blas_name=acml;
            acx_blas_libs="-lacml -lg2c"],
        [],[-lg2c])
    fi
    fi
    if test $acx_blas_ok = no; then
    if test $acx_blas_gfortran = yes; then
    AC_MSG_NOTICE([----- ACML static])
    AC_CHECK_LIB(acml, $1,
        [acx_blas_ok=yes; acx_blas_name=acml;
            acx_blas_libs="-lacml -lgfortran"],
        [],[-lgfortran])
    fi
    fi
fi

#
#   Automatically Tuned Linear Algebra Software (ATLAS)
#
if test $acx_blas_ok = no; then
    AC_MSG_NOTICE([----- ATLAS])
    AC_CHECK_LIB(f77blas, $1,
        [acx_blas_ok=yes; acx_blas_name=atlas;
            acx_blas_libs="-lf77blas -latlas"], [], [-latlas])
fi

#
#   HP Mathematical software LIBrary (MLIB)
#
if test $acx_blas_ok = no; then
    AC_MSG_NOTICE([----- HP MLIB])
    AC_CHECK_LIB(veclib, $1, 
        [acx_blas_ok=yes; acx_blas_name=veclib; acx_blas_libs="-lveclib"])
fi

#
#   PhiPACK
#
if test $acx_blas_ok = no; then
    AC_MSG_NOTICE([----- PhiPACK])
    AC_CHECK_LIB(blas, $1,
        [AC_CHECK_LIB(dgemm, dgemm,
            [AC_CHECK_LIB(sgemm, $1,
                [acx_blas_ok=yes; acx_blas_name=phipack;
                    acx_blas_libs="-lsgemm -ldgemm -lblas"],
                [], [-lblas])],
            [], [-lblas])
        ])
fi

#
#   Compaq CXML
#
if test $acx_blas_ok = no; then
    AC_MSG_NOTICE([----- Compaq CXML])
    AC_CHECK_LIB(cxml, $1, 
        [acx_blas_ok=yes; acx_blas_name=cxml; acx_blas_libs="-lcxml"])
fi

#   Compaq DXML ( old version of CXML)
if test $acx_blas_ok = no; then
    AC_CHECK_LIB(dxml, $1, 
        [acx_blas_ok=yes; acx_blas_name=dxml; acx_blas_libs="-ldxml"])
fi

#
#   SUN Performance Library
#
if test $acx_blas_ok = no; then
    if test "x$GCC" != xyes; then # only works with Sun CC
        AC_MSG_NOTICE([----- SUNperf])
        AC_CHECK_LIB(sunmath, acosp,
            [AC_CHECK_LIB(sunperf, $1,
                [acx_blas_libs="-xlic_lib=sunperf -lsunmath"
                    acx_blas_ok=yes; acx_blas_name=sunperf],[],[-lsunmath])
            ])
    fi
fi

#
#   SGI/Cray Scientific Library
#
if test $acx_blas_ok = no; then
    AC_MSG_NOTICE([----- SGI/Cray])
    AC_CHECK_LIB(scs, $1,
        [acx_blas_ok=yes; acx_blas_name=scsl; acx_blas_libs="-lscs"])
fi


#
#   SGIMATH library
#
if test $acx_blas_ok = no; then
    AC_MSG_NOTICE([----- SGIMATH])
    AC_CHECK_LIB(complib.sgimath, $1,
        [acx_blas_ok=yes; acx_blas_name=sgimath;
            acx_blas_libs="-lcomplib.sgimath"])
fi

#
#   IBM ESSL (requires generic BLAS lib, too)
#
if test $acx_blas_ok = no; then
    AC_MSG_NOTICE([----- IBM ESSL])
    unset ac_cv_lib_blas_sgemm
    AC_CHECK_LIB(blas, $1,
        [AC_CHECK_LIB(essl, $1,
            [acx_blas_ok=yes; acx_blas_name=essl;
                acx_blas_libs="-lessl -lblas"],
            [], [-lblas])
	])
fi

#
# Generic BLAS library
#
if test $acx_blas_ok = no; then
    AC_MSG_NOTICE([----- Generic BLAS])
    unset ac_cv_lib_blas_sgemm
    AC_CHECK_LIB(blas, $1,
        [acx_blas_ok=yes; acx_blas_name=generic; acx_blas_libs="-lblas"])
fi

#
#   SGI/Cray other
#   (it must be located after Generic BLAS).

if test $acx_blas_ok = no; then
    unset ac_cv_lib_blas_sgemm
    AC_CHECK_LIB(blas, $1, 
	[acx_blas_ok=yes; acx_blas_name=craylib;
            acx_blas_libs="-lblas -lcraylibs"],[],[-lcraylibs])
fi

BLAS_LIBS="$acx_blas_libs"
AC_SUBST(BLAS_LIBS)

LIBS="$acx_blas_save_LIBS"
LDFLAGS="$acx_blas_save_LDFLAGS $acx_blas_dir"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test "x$acx_blas_ok" = xyes; then
    ifelse([$2],,AC_DEFINE(HAVE_BLAS,1,[Define if you have a BLAS library.]),[$2])
    :
else
    acx_blas_ok=no
    $3
fi
])dnl ACX_BLAS

