/* Automatically generated. Do not edit by hand. */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>


extern SEXP R_hdfmat_cp(SEXP x, SEXP ds);
extern SEXP R_hdfmat_fill(SEXP ds, SEXP x);
extern SEXP R_hdfmat_init(SEXP fp, SEXP name, SEXP nrows, SEXP ncols, SEXP type);
extern SEXP R_hdfmat_open(SEXP filename);
extern SEXP R_hdfmat_read(SEXP m_, SEXP n_, SEXP ds);
extern SEXP R_hdfmat_scale(SEXP m_, SEXP n_, SEXP ds, SEXP val_);
extern SEXP R_hdfmat_set_diag(SEXP m_, SEXP n_, SEXP ds, SEXP val_);

static const R_CallMethodDef CallEntries[] = {
  {"R_hdfmat_cp", (DL_FUNC) &R_hdfmat_cp, 2},
  {"R_hdfmat_fill", (DL_FUNC) &R_hdfmat_fill, 2},
  {"R_hdfmat_init", (DL_FUNC) &R_hdfmat_init, 5},
  {"R_hdfmat_read", (DL_FUNC) &R_hdfmat_read, 3},
  {"R_hdfmat_open", (DL_FUNC) &R_hdfmat_open, 1},
  {"R_hdfmat_scale", (DL_FUNC) &R_hdfmat_scale, 4},
  {"R_hdfmat_set_diag", (DL_FUNC) &R_hdfmat_set_diag, 4},
  {NULL, NULL, 0}
};

void R_init_coop(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
