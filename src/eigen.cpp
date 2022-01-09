#include "lanczos.hh"
#include "omp.h"

#include <fml/src/fml/cpu/cpumat.hh>
#include <fml/src/fml/cpu/cpuvec.hh>
#include <fml/src/fml/cpu/linalg/eigen.hh>

#include "hdfmat.h"
#include "extptr.h"
#include "types.h"


template <typename T>
static inline void lanczos(const hsize_t n, const int k,
  T *alpha, T *beta, T *q, H5::DataSet *dataset, H5::PredType h5type)
{
  T *A_j = (T*) std::malloc(n * sizeof(*A_j));
  T *v = (T*) std::malloc(n * sizeof(*v));
  
  hsize_t slice[2];
  slice[0] = 1;
  slice[1] = n;
  
  H5::DataSpace mem_space(2, slice, NULL);
  H5::DataSpace data_space = dataset->getSpace();
  
  hsize_t offset[2];
  offset[1] = 0;
  
  for (int i=0; i<k; i++)
  {
    // v = A*q[, i] iterate over *row j* of A - A_j is 1xn
    for (hsize_t j=0; j<n; j++)
    {
      offset[0] = j;
      data_space.selectHyperslab(H5S_SELECT_SET, slice, offset);
      dataset->read(A_j, h5type, mem_space, data_space);
      
      v[j] = dot(n, A_j, q+n*i);
    }
    
    alpha[i] = dot(n, q+n*i, v);
    
    if (i == 0)
    {
      // v = v - alpha[i]*q[, i]
      #pragma omp for simd
      for (hsize_t j=0; j<n; j++)
        v[j] -= alpha[i] * q[j + n*i];
    }
    else
    {
      // v = v - beta[i-1]*q[, i-1] - alpha[i]*q[, i]
      #pragma omp for simd
      for (hsize_t j=0; j<n; j++)
        v[j] -= beta[i-1] * q[j + n*(i-1)];
      
      #pragma omp for simd
      for (hsize_t j=0; j<n; j++)
        v[j] -= alpha[i] * q[j + n*i];
    }
    
    beta[i] = l2norm(n, v);
    
    if (i < k-1)
    {
      #pragma omp for simd
      for (hsize_t j=0; j<n; j++)
        q[j + n*(i+1)] = v[j] / beta[i];
    }
  }
  
  std::free(A_j);
  std::free(v);
}



template <typename T>
static inline void eigen_sym(const hsize_t n, const int k,
  T *values, H5::DataSet *dataset, H5::PredType h5type)
{
  T *alpha, *beta, *q;
  alloc(n, k, &alpha, &beta, &q);
  initialize(n, k, q);
  
  lanczos(n, k, alpha, beta, q, dataset, h5type);
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



extern "C" SEXP R_hdfmat_eigen_sym(SEXP k_, SEXP n_, SEXP ds, SEXP type)
{
  SEXP values;
  H5::DataSet *dataset = (H5::DataSet*) getRptr(ds);
  
  const int k = INT(k_);
  const hsize_t n = (hsize_t) DBL(n_);
  
  if (INT(type) == TYPE_DOUBLE)
  {
    PROTECT(values = allocVector(REALSXP, k));
    TRY_CATCH( eigen_sym(n, k, REAL(values), dataset, H5::PredType::IEEE_F64LE) );
  }
  else // if (INT(type) == TYPE_FLOAT)
  {
    PROTECT(values = allocVector(INTSXP, k));
    TRY_CATCH( eigen_sym(n, k, FLOAT(values), dataset, H5::PredType::IEEE_F32LE) );
  }
  
  UNPROTECT(1);
  return values;
}
