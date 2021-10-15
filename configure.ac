AC_PREREQ([2.69])
AC_INIT
AC_CONFIG_SRCDIR([DESCRIPTION])
AC_LANG_PUSH([C++])

m4_include([tools/ax_lib_hdf5.m4])
AX_LIB_HDF5()

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`

R_SCMD="${R_HOME}/bin/Rscript -e"
FLOAT_LIBS=`${R_SCMD} "float:::ldflags()"`

AC_OPENMP
if test -n "${OPENMP_CXXFLAGS}"; then
  USE_OMP="TRUE"
  OMPFLAGS_CXX="\$(SHLIB_OPENMP_CXXFLAGS)"
else
  USE_OMP="FALSE"
  OMPFLAGS_CXX=""
fi

echo " "
echo "****************** Results of hdfmat package configure ******************"
echo "* R Report:"
echo "    >> FLOAT_LIBS=${FLOAT_LIBS}"
echo "* OpenMP Report:"
if test "X${USE_OMP}" = "XTRUE"; then
  echo "    >> OMPFLAGS_CXX=${OMPFLAGS_CXX}"
else
  echo "    >> not used"
fi
echo "* HDF5 report:"
echo "    >> HDF5_CPPFLAGS=${HDF5_CPPFLAGS}"
echo "    >> HDF5_LDFLAGS=${HDF5_LDFLAGS}"
echo "*************************************************************************"
echo " "

AC_SUBST(OMPFLAGS_CXX)
AC_SUBST(FLOAT_LIBS)
AC_SUBST(HDF5_CPPFLAGS)
AC_SUBST(HDF5_LDFLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
