#include <cstdlib>

#include <float/float32.h>

#include "hdfmat.h"
#include "extptr.h"
#include "types.h"


template <typename T>
static inline void write(const hsize_t m, const hsize_t n,
  const hsize_t row_offset, T *x, H5::DataSet *dataset, H5::PredType h5type)
{
  hsize_t slice[2];
  slice[0] = m;
  slice[1] = n;
  
  H5::DataSpace mem_space(2, slice, NULL);
  H5::DataSpace data_space = dataset->getSpace();
  
  hsize_t offset[2];
  offset[0] = row_offset;
  offset[1] = 0;
  
  data_space.selectHyperslab(H5S_SELECT_SET, slice, offset);
  dataset->write(x, h5type, mem_space, data_space);
}

extern "C" SEXP R_hdfmat_fill(SEXP ds, SEXP x, SEXP row_offset_, SEXP type)
{
  H5::DataSet *dataset = (H5::DataSet*) getRptr(ds);
  
  const hsize_t m = (hsize_t) nrows(x);
  const hsize_t n = (hsize_t) ncols(x);
  
  const hsize_t row_offset = (hsize_t) DBL(row_offset_);
  
  if (INT(type) == TYPE_DOUBLE)
  {
    TRY_CATCH( write(m, n, row_offset, REAL(x), dataset, H5::PredType::IEEE_F64LE) );
  }
  else // if (INT(type) == TYPE_FLOAT)
  {
    TRY_CATCH( write(m, n, row_offset, FLOAT(x), dataset, H5::PredType::IEEE_F32LE) );
  }
  
  return R_NilValue;
}



template <typename T>
static inline void read(const hsize_t row_start, const hsize_t row_stop,
  const hsize_t col_start, const hsize_t col_stop,
  const hsize_t m, const hsize_t n, T *x, H5::DataSet *dataset,
  H5::PredType h5type)
{
  hsize_t slice[2];
  slice[0] = row_stop - row_start + 1;
  slice[1] = col_stop - col_start + 1;
  
  H5::DataSpace mem_space(2, slice, NULL);
  H5::DataSpace data_space = dataset->getSpace();
  
  hsize_t offset[2];
  offset[0] = row_start;
  offset[1] = col_start;
  
  data_space.selectHyperslab(H5S_SELECT_SET, slice, offset);
  dataset->read(x, h5type, mem_space, data_space);
}

extern "C" SEXP R_hdfmat_read(SEXP row_start_, SEXP row_stop_, SEXP col_start_, SEXP col_stop_, SEXP m_, SEXP n_, SEXP ds, SEXP type)
{
  SEXP ret;
  
  H5::DataSet *dataset = (H5::DataSet*) getRptr(ds);
  
  const hsize_t m = (hsize_t) DBL(m_);
  const hsize_t n = (hsize_t) DBL(n_);
  
  const hsize_t row_start = (hsize_t) DBL(row_start_);
  const hsize_t row_stop = (hsize_t) DBL(row_stop_);
  const hsize_t row_len = row_stop - row_start + 1;
  
  const hsize_t col_start = (hsize_t) DBL(col_start_);
  const hsize_t col_stop = (hsize_t) DBL(col_stop_);
  const hsize_t col_len = col_stop - col_start + 1;
  
  if (INT(type) == TYPE_DOUBLE)
  {
    PROTECT(ret = allocMatrix(REALSXP, col_len, row_len));
    TRY_CATCH( read(row_start, row_stop, col_start, col_stop, m, n, REAL(ret), dataset, H5::PredType::IEEE_F64LE) );
  }
  else // if (INT(type) == TYPE_FLOAT)
  {
    PROTECT(ret = allocMatrix(INTSXP, col_len, row_len));
    TRY_CATCH( read(row_start, row_stop, col_start, col_stop, m, n, FLOAT(ret), dataset, H5::PredType::IEEE_F32LE) );
  }
  
  UNPROTECT(1);
  return ret;
}
