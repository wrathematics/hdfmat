#include "lanczos.hh"

#include <fml/src/fml/cpu/cpumat.hh>
#include <fml/src/fml/cpu/cpuvec.hh>
#include <fml/src/fml/cpu/linalg/blas.hh>
#include <fml/src/fml/cpu/linalg/linalg_eigen.hh>

#include "hdfmat.h"
#include "extptr.h"
#include "types.h"


template <typename T>
static inline void lanczos(const hsize_t m, const hsize_t n, const int k,
  T *alpha, T *beta, T *q, H5::DataSet *dataset, H5::PredType h5type)
{
  // H5::Exception::dontPrint();
  T *A_j = (T*) std::malloc(m * sizeof(*A_j));
  T *v = (T*) std::malloc((m+n) * sizeof(*v));
  
  hsize_t slice[2];
  slice[0] = 1;
  slice[1] = n;
  
  H5::DataSpace mem_space(2, slice, NULL);
  H5::DataSpace data_space = dataset->getSpace();
  
  hsize_t offset[2];
  offset[1] = 0;
  
  for (int i=0; i<k; i++)
  {
    std::memset(v+m, 0, n*sizeof(*v));
    
    // v = A*q[, i] iterate over *row j* of A - A_j is 1xn
    for (hsize_t j=0; j<m; j++)
    {
      offset[0] = j;
      data_space.selectHyperslab(H5S_SELECT_SET, slice, offset);
      dataset->read(A_j, h5type, mem_space, data_space);
      
      #pragma omp for simd if(n > fml::omp::OMP_MIN_SIZE)
      for (hsize_t jj=0; jj<n; jj++)
        v[m+jj] += A_j[jj] * q[j + (m+n)*i];
      
      v[j] = dot(n, A_j, q + m+(m+n)*i);
    }
    
    alpha[i] = dot(m+n, q + (m+n)*i, v);
    
    if (i == 0)
    {
      // v = v - alpha[i]*q[, i]
      #pragma omp for simd if(m+n > fml::omp::OMP_MIN_SIZE)
      for (hsize_t j=0; j<m+n; j++)
        v[j] -= alpha[i] * q[j + (m+n)*i];
    }
    else
    {
      // v = v - beta[i-1]*q[, i-1] - alpha[i]*q[, i]
      #pragma omp for simd if(m+n > fml::omp::OMP_MIN_SIZE)
      for (hsize_t j=0; j<m+n; j++)
        v[j] -= beta[i-1] * q[j + (m+n)*(i-1)];
      
      #pragma omp for simd if(m+n > fml::omp::OMP_MIN_SIZE)
      for (hsize_t j=0; j<m+n; j++)
        v[j] -= alpha[i] * q[j + (m+n)*i];
    }
    
    beta[i] = l2norm(m+n, v);
    
    if (i < k-1)
    {
      #pragma omp for simd if(m+n > fml::omp::OMP_MIN_SIZE)
      for (hsize_t j=0; j<m+n; j++)
        q[j + (m+n)*(i+1)] = v[j] / beta[i];
    }
  }
  
  std::free(A_j);
  std::free(v);
}



template <typename T>
static inline void svd(const hsize_t m, const hsize_t n, const int k,
  T *values, H5::DataSet *dataset, H5::PredType h5type)
{
  T *alpha, *beta, *q;
  alloc(m+n, k, &alpha, &beta, &q);
  initialize(m+n, k, q);
  
  lanczos(m, n, k, alpha, beta, q, dataset, h5type);
  std::free(q);
  
  T *td = (T *) std::malloc(k*k * sizeof(*td));
  tridiagonal(k, alpha, beta, td);
  std::free(alpha);
  std::free(beta);
  
  fml::cpumat<T> td_mat(td, k, k, false);
  fml::cpuvec<T> values_vec(values, k, false);
  fml::linalg::eigen_sym(td_mat, values_vec);
  std::free(td);
  
  values_vec.rev();
  for (int i=0; i<k; i++)
  {
    if (values[i] < 0)
      values[i] = 0;
  }
}



extern "C" SEXP R_hdfmat_svd(SEXP k_, SEXP m_, SEXP n_, SEXP ds, SEXP type)
{
  SEXP values;
  H5::DataSet *dataset = (H5::DataSet*) getRptr(ds);
  
  const int k = INT(k_);
  const hsize_t m = (hsize_t) DBL(m_);
  const hsize_t n = (hsize_t) DBL(n_);
  
  if (INT(type) == TYPE_DOUBLE)
  {
    PROTECT(values = allocVector(REALSXP, k));
    svd(m, n, k, REAL(values), dataset, H5::PredType::IEEE_F64LE);
  }
  else // if (INT(type) == TYPE_FLOAT)
  {
    PROTECT(values = allocVector(INTSXP, k));
    svd(m, n, k, FLOAT(values), dataset, H5::PredType::IEEE_F32LE);
  }
  
  UNPROTECT(1);
  return values;
}
