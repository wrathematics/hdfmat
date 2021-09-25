#include <cmath>
#include <cstdlib>

#include <float/float32.h>

#include <fml/src/fml/cpu/cpumat.hh>
#include <fml/src/fml/cpu/cpuvec.hh>
#include <fml/src/fml/cpu/linalg/blas.hh>
#include <fml/src/fml/cpu/linalg/linalg_eigen.hh>

#include "hdfmat.h"
#include "extptr.h"
#include "types.h"

#define OMP_MIN_SIZE 2500


template <typename T>
static inline T dot(const hsize_t n, const T *x, const T *y)
{
  T d = 0;
  #pragma omp parallel for simd reduction(+:d) if(n > OMP_MIN_SIZE)
  for (hsize_t i=0; i<n; i++)
    d += x[i] * y[i];
  
  return d;
}



template <typename T>
static inline T l2norm(const hsize_t n, const T *x)
{
  return sqrt(dot(n, x, x));
}



template <typename T>
static inline void alloc(const hsize_t n, const int k, T **alpha, T **beta, T **q)
{
  *alpha = (T *) malloc(k * sizeof(**alpha));
  *beta = (T *) malloc(k * sizeof(**beta));
  *q = (T *) malloc(n*k * sizeof(**q));
}



template <typename T>
static inline void initialize(const hsize_t n, const int k, T *q)
{
  memset(q, 0, n*k*sizeof(*q));
  GetRNGstate();
  
  for (hsize_t i=0; i<n; i++)
    q[i] = (T)unif_rand();
  
  T l2 = l2norm(n, q);
  #pragma omp parallel for simd if(n > OMP_MIN_SIZE)
  for (hsize_t i=0; i<n; i++)
    q[i] /= l2;
  
  PutRNGstate();
}



template <typename T>
static inline void tridiagonal(const int k, const T *alpha, const T *beta,
  T *td)
{
  memset(td, 0, k*k*sizeof(T));
  
  for (int i=0; i<k; i++)
    td[i + k*i] = alpha[i];
  
  for (int i=0; i<k-1; i++)
  {
    td[i + k*(i+1)] = beta[i];
    td[i+1 + k*i] = beta[i];
  }
}



template <typename T>
static inline void lanczos(const hsize_t n, const int k,
  T *alpha, T *beta, T *q, H5::DataSet *dataset, H5::PredType h5type)
{
  // H5::Exception::dontPrint();
  T *A_j = (T*) malloc(n * sizeof(*A_j));
  T *v = (T*) malloc(n * sizeof(*v));
  
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
      for (hsize_t j=0; j<n; j++)
        v[j] -= alpha[i] * q[j + n*i];
    }
    else
    {
      // v = v - beta[i-1]*q[, i-1] - alpha[i]*q[, i]
      for (hsize_t j=0; j<n; j++)
        v[j] -= beta[i-1] * q[j + n*(i-1)];
      for (hsize_t j=0; j<n; j++)
        v[j] -= alpha[i] * q[j + n*i];
    }
  
    beta[i] = l2norm(n, v);
  
    if (i < k-1)
    {
      for (hsize_t j=0; j<n; j++)
        q[j + n*(i+1)] = v[j] / beta[i];
    }
  }
  
  free(A_j);
  free(v);
}



template <typename T>
static inline void eigen_sym(const hsize_t n, const int k,
  T *values, H5::DataSet *dataset, H5::PredType h5type)
{
  T *alpha, *beta, *q;
  alloc(n, k, &alpha, &beta, &q);
  initialize(n, k, q);
  
  lanczos(n, k, alpha, beta, q, dataset, h5type);
  free(q);
  
  T *td = (T *) malloc(k*k * sizeof(*td));
  tridiagonal(k, alpha, beta, td);
  free(alpha);
  free(beta);
  
  fml::cpumat<T> td_mat(td, k, k, false);
  fml::cpuvec<T> values_vec(values, k, false);
  fml::linalg::eigen_sym(td_mat, values_vec);
  free(td);
  
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
    eigen_sym(n, k, REAL(values), dataset, H5::PredType::IEEE_F64LE);
  }
  else // if (INT(type) == TYPE_FLOAT)
  {
    PROTECT(values = allocVector(INTSXP, k));
    eigen_sym(n, k, FLOAT(values), dataset, H5::PredType::IEEE_F64LE);
  }
  
  UNPROTECT(1);
  return values;
}
