AC_INIT([gridsstools], [1.0])
AC_CONFIG_SRCDIR([gridsstools.c])
AC_CONFIG_HEADERS([config.h])

m4_include([m4/ax_with_htslib.m4])

AC_PROG_CC

AX_WITH_HTSLIB
if test "$ax_cv_htslib" != yes; then
  AC_MSG_ERROR([HTSlib development files not found])
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
