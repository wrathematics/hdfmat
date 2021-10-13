/* Automatically generated. Do not edit by hand. */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>


extern SEXP R_hdfmat_cp(SEXP x, SEXP ds);
extern SEXP R_hdfmat_eigen_sym(SEXP k_, SEXP n_, SEXP ds, SEXP type);
extern SEXP R_hdfmat_fill(SEXP ds, SEXP x, SEXP row_offset_, SEXP type);
extern SEXP R_hdfmat_fill_diag(SEXP m_, SEXP n_, SEXP ds, SEXP val_, SEXP type);
extern SEXP R_hdfmat_fill_linspace(SEXP m_, SEXP n_, SEXP ds, SEXP start_, SEXP stop_, SEXP type);
extern SEXP R_hdfmat_fill_runif(SEXP m_, SEXP n_, SEXP ds, SEXP min_, SEXP max_, SEXP type);
extern SEXP R_hdfmat_fill_val(SEXP m_, SEXP n_, SEXP ds, SEXP val_, SEXP type);
extern SEXP R_hdfmat_finalize(SEXP fp, SEXP ds);
extern SEXP R_hdfmat_inherit(SEXP fp, SEXP name);
extern SEXP R_hdfmat_init(SEXP fp, SEXP name, SEXP nrows, SEXP ncols, SEXP type);
extern SEXP R_hdfmat_open(SEXP filename, SEXP mode);
extern SEXP R_hdfmat_read(SEXP m_, SEXP n_, SEXP ds, SEXP type);
extern SEXP R_hdfmat_scale(SEXP m_, SEXP n_, SEXP ds, SEXP val_, SEXP type);
extern SEXP R_hdfmat_svd(SEXP k_, SEXP m_, SEXP n_, SEXP ds, SEXP type);
extern SEXP R_hdfmat_tcp(SEXP x, SEXP ds);

static const R_CallMethodDef CallEntries[] = {
  {"R_hdfmat_cp", (DL_FUNC) &R_hdfmat_cp, 2},
  {"R_hdfmat_eigen_sym", (DL_FUNC) &R_hdfmat_eigen_sym, 4},
  {"R_hdfmat_fill", (DL_FUNC) &R_hdfmat_fill, 4},
  {"R_hdfmat_fill_diag", (DL_FUNC) &R_hdfmat_fill_diag, 5},
  {"R_hdfmat_fill_linspace", (DL_FUNC) &R_hdfmat_fill_linspace, 6},
  {"R_hdfmat_fill_runif", (DL_FUNC) &R_hdfmat_fill_runif, 6},
  {"R_hdfmat_fill_val", (DL_FUNC) &R_hdfmat_fill_val, 5},
  {"R_hdfmat_finalize", (DL_FUNC) &R_hdfmat_finalize, 2},
  {"R_hdfmat_inherit", (DL_FUNC) &R_hdfmat_inherit, 2},
  {"R_hdfmat_init", (DL_FUNC) &R_hdfmat_init, 5},
  {"R_hdfmat_read", (DL_FUNC) &R_hdfmat_read, 4},
  {"R_hdfmat_open", (DL_FUNC) &R_hdfmat_open, 2},
  {"R_hdfmat_scale", (DL_FUNC) &R_hdfmat_scale, 5},
  {"R_hdfmat_svd", (DL_FUNC) &R_hdfmat_svd, 5},
  {"R_hdfmat_tcp", (DL_FUNC) &R_hdfmat_tcp, 2},
  {NULL, NULL, 0}
};

void R_init_coop(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
