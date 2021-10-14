#include <cstdlib>

#include "hdfmat.h"
#include "extptr.h"
#include "omp.h"
#include "types.h"


template <typename T>
static inline void fill_val(const T v, const hsize_t m, const hsize_t n,
  H5::DataSet *dataset, H5::PredType h5type)
{
  T *x = (T*) std::malloc(n * sizeof(*x));
  
  #pragma omp for simd if(n > OMP_MIN_LEN)
  for (hsize_t j=0; j<n; j++)
    x[j] = v;
  
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
    dataset->write(x, h5type, mem_space, data_space);
  }
  
  std::free(x);
}

extern "C" SEXP R_hdfmat_fill_val(SEXP m_, SEXP n_, SEXP ds, SEXP val_, SEXP type)
{
  // H5::Exception::dontPrint();
  H5::DataSet *dataset = (H5::DataSet*) getRptr(ds);
  
  const hsize_t m = (hsize_t) REAL(m_)[0];
  const hsize_t n = (hsize_t) REAL(n_)[0];
  const double v = DBL(val_);
  
  if (INT(type) == TYPE_DOUBLE)
    fill_val(v, m, n, dataset, H5::PredType::IEEE_F64LE);
  else // if (INT(type) == TYPE_FLOAT)
    fill_val((float)v, m, n, dataset, H5::PredType::IEEE_F32LE);
  
  return R_NilValue;
}



template <typename T>
static inline void fill_linspace(const T start, const T stop, const hsize_t m,
  const hsize_t n, H5::DataSet *dataset, H5::PredType h5type)
{
  T *x = (T*) std::malloc(n * sizeof(*x));
  
  hsize_t slice[2];
  slice[0] = 1;
  slice[1] = n;
  
  H5::DataSpace mem_space(2, slice, NULL);
  H5::DataSpace data_space = dataset->getSpace();
  
  hsize_t offset[2];
  offset[1] = 0;
  
  const T v = (stop-start)/((T) (m*n - 1));
  
  for (hsize_t i=0; i<m; i++)
  {
    offset[0] = i;
    
    #pragma omp for simd if(n > OMP_MIN_LEN)
    for (hsize_t j=0; j<n; j++)
      x[j] = v * ((T) i + m*j) + start;
    
    data_space.selectHyperslab(H5S_SELECT_SET, slice, offset);
    dataset->write(x, h5type, mem_space, data_space);
  }
  
  std::free(x);
}

extern "C" SEXP R_hdfmat_fill_linspace(SEXP m_, SEXP n_, SEXP ds, SEXP start_, SEXP stop_, SEXP type)
{
  // H5::Exception::dontPrint();
  H5::DataSet *dataset = (H5::DataSet*) getRptr(ds);
  
  const hsize_t m = (hsize_t) REAL(m_)[0];
  const hsize_t n = (hsize_t) REAL(n_)[0];
  
  const double start = DBL(start_);
  const double stop = DBL(stop_);
  
  if (INT(type) == TYPE_DOUBLE)
    fill_linspace(start, stop, m, n, dataset, H5::PredType::IEEE_F64LE);
  else // if (INT(type) == TYPE_FLOAT)
    fill_linspace((float)start, (float)stop, m, n, dataset, H5::PredType::IEEE_F32LE);
  
  return R_NilValue;
}



template <typename T>
static inline void fill_runif(const T min, const T max, const hsize_t m,
  const hsize_t n, H5::DataSet *dataset, H5::PredType h5type)
{
  T *x = (T*) std::malloc(n * sizeof(*x));
  
  hsize_t slice[2];
  slice[0] = 1;
  slice[1] = n;
  
  H5::DataSpace mem_space(2, slice, NULL);
  H5::DataSpace data_space = dataset->getSpace();
  
  hsize_t offset[2];
  offset[1] = 0;
  
  GetRNGstate();
  
  for (hsize_t i=0; i<m; i++)
  {
    offset[0] = i;
    
    for (hsize_t j=0; j<n; j++)
      x[j] = (T) min + (max - min)*unif_rand();
    
    data_space.selectHyperslab(H5S_SELECT_SET, slice, offset);
    dataset->write(x, h5type, mem_space, data_space);
  }
  
  PutRNGstate();
  
  std::free(x);
}

extern "C" SEXP R_hdfmat_fill_runif(SEXP m_, SEXP n_, SEXP ds, SEXP min_, SEXP max_, SEXP type)
{
  // H5::Exception::dontPrint();
  H5::DataSet *dataset = (H5::DataSet*) getRptr(ds);
  
  const hsize_t m = (hsize_t) REAL(m_)[0];
  const hsize_t n = (hsize_t) REAL(n_)[0];
  
  const double min = DBL(min_);
  const double max = DBL(max_);
  
  if (INT(type) == TYPE_DOUBLE)
    fill_runif(min, max, m, n, dataset, H5::PredType::IEEE_F64LE);
  else // if (INT(type) == TYPE_FLOAT)
    fill_runif((float)min, (float)max, m, n, dataset, H5::PredType::IEEE_F32LE);
  
  return R_NilValue;
}



template <typename T>
static inline void fill_diag(const int vlen, const T *v, const hsize_t m,
  const hsize_t n, H5::DataSet *dataset, H5::PredType h5type)
{
  const hsize_t len = m<n?m:n;
  
  hsize_t slice[2];
  slice[0] = 1;
  slice[1] = 1;
  
  H5::DataSpace mem_space(2, slice, NULL);
  H5::DataSpace data_space = dataset->getSpace();
  
  hsize_t offset[2];
  
  for (hsize_t i=0; i<len; i++)
  {
    offset[0] = i;
    offset[1] = i;
    data_space.selectHyperslab(H5S_SELECT_SET, slice, offset);
    
    dataset->write(v+(i % vlen), h5type, mem_space, data_space);
  }
}

extern "C" SEXP R_hdfmat_fill_diag(SEXP m_, SEXP n_, SEXP ds, SEXP val_, SEXP type)
{
  // H5::Exception::dontPrint();
  H5::DataSet *dataset = (H5::DataSet*) getRptr(ds);
  
  const hsize_t m = (hsize_t) REAL(m_)[0];
  const hsize_t n = (hsize_t) REAL(n_)[0];
  const int vlen = LENGTH(val_);
  const double *v = REAL(val_);
  
  if (INT(type) == TYPE_DOUBLE)
    fill_diag(vlen, v, m, n, dataset, H5::PredType::IEEE_F64LE);
  else // if (INT(type) == TYPE_FLOAT)
  {
    float *v_f = (float*) std::malloc(vlen * sizeof(*v_f));
    for (int i=0; i<vlen; i++)
      v_f[i] = (float) v[i];
    
    fill_diag(vlen, v_f, m, n, dataset, H5::PredType::IEEE_F32LE);
    std::free(v_f);
  }
  
  return R_NilValue;
}
