#include <cstdlib>
#include <cstring>

#include <float/float32.h>

#include <fml/src/fml/cpu/linalg/blas.hh>

#include "hdfmat.h"
#include "extptr.h"
#include "types.h"


template <typename T>
static inline void cp(const int m, const int n, const T *x,
  H5::DataSet *dataset, H5::PredType h5type)
{
  T *y = (T*) std::malloc(m * sizeof(*y));
  T *cp_row = (T*) std::malloc(n * sizeof(*cp_row));
  
  hsize_t slice[2];
  slice[0] = (hsize_t) 1;
  slice[1] = (hsize_t) n;
  
  H5::DataSpace mem_space(2, slice, NULL);
  H5::DataSpace data_space = dataset->getSpace();
  
  hsize_t offset[2];
  offset[1] = 0;
  
  for (int i=0; i<n; i++)
  {
    std::memcpy(y, x+i*m, m*sizeof(T));
    fml::blas::gemm('T', 'N', n, 1, m, (T)1, x, m, y, m, (T)0, cp_row, n);
    
    offset[0] = (hsize_t) i;
    data_space.selectHyperslab(H5S_SELECT_SET, slice, offset);
    
    dataset->write(cp_row, h5type, mem_space, data_space);
  }
  
  std::free(y);
  std::free(cp_row);
}

extern "C" SEXP R_hdfmat_cp(SEXP x, SEXP ds, SEXP type)
{
  // H5::Exception::dontPrint();
  H5::DataSet *dataset = (H5::DataSet*) getRptr(ds);
  
  const int m = nrows(x);
  const int n = ncols(x);
  
  if (INT(type) == TYPE_DOUBLE)
    cp(m, n, REAL(x), dataset, H5::PredType::IEEE_F64LE);
  else // if (INT(type) == TYPE_FLOAT)
    cp(m, n, FLOAT(x), dataset, H5::PredType::IEEE_F32LE);
  
  return R_NilValue;
}
