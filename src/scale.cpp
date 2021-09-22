#include <cstdlib>

#include "hdfmat.h"
#include "extptr.h"


template <typename T>
static inline void scale(const T val, const hsize_t m, const hsize_t n,
  H5::DataSet *dataset, H5::PredType h5type)
{
  T *x = (T*) malloc(n * sizeof(*x));
  
  hsize_t slice[2];
  slice[0] = 1;
  slice[1] = n;
  
  H5::DataSpace mem_space(2, slice, NULL);
  H5::DataSpace data_space = dataset->getSpace();
  
  hsize_t offset[2];
  offset[1] = 0;
  
  for (hsize_t i=0; i<m; i++)
  {
    offset[0] = i;
    data_space.selectHyperslab(H5S_SELECT_SET, slice, offset);
    
    dataset->read(x, h5type, mem_space, data_space);
    for (hsize_t j=0; j<n; j++)
      x[j] *= val;
    
    dataset->write(x, h5type, mem_space, data_space);
  }
  
  free(x);
}



extern "C" SEXP R_hdfmat_scale(SEXP m_, SEXP n_, SEXP ds, SEXP val_)
{
  // H5::Exception::dontPrint();
  H5::DataSet *dataset = (H5::DataSet*) getRptr(ds);
  
  const hsize_t m = (hsize_t) REAL(m_)[0];
  const hsize_t n = (hsize_t) REAL(n_)[0];
  const double val = REAL(val_)[0];
  
  scale(val, m, n, dataset, H5::PredType::IEEE_F64LE);
  
  return R_NilValue;
}
