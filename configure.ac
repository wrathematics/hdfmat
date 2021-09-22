AC_PREREQ([2.69])
AC_INIT
AC_CONFIG_SRCDIR([DESCRIPTION])

m4_include([tools/ax_lib_hdf5.m4])
AX_LIB_HDF5()

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`

AC_PROG_CC

AC_OPENMP
if test -n "${OPENMP_CFLAGS}"; then
  OMPFLAGS_C="\$(SHLIB_OPENMP_CFLAGS)"
else
  OMPFLAGS_C=""
fi

AC_SUBST(OMPFLAGS_C)
AC_SUBST(HDF5_LDFLAGS)
AC_SUBST(HDF5_CPPFLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
