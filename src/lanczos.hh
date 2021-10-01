#ifndef LANCZOS_H
#define LANCZOS_H
#pragma once


#include <cmath>
#include <cstdlib>
#include <cstring>

#include <fml/src/fml/_internals/omp.hh>

#include <H5Cpp.h>

#include <R.h>
#include <Rinternals.h>


template <typename T>
static inline T dot(const hsize_t n, const T *x, const T *y)
{
  T d = 0;
  #pragma omp parallel for simd reduction(+:d) if(n > fml::omp::OMP_MIN_SIZE)
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
  std::memset(q, 0, n*k*sizeof(*q));
  GetRNGstate();
  
  for (hsize_t i=0; i<n; i++)
    q[i] = (T)unif_rand();
  
  T l2 = l2norm(n, q);
  #pragma omp parallel for simd if(n > fml::omp::OMP_MIN_SIZE)
  for (hsize_t i=0; i<n; i++)
    q[i] /= l2;
  
  PutRNGstate();
}



template <typename T>
static inline void tridiagonal(const int k, const T *alpha, const T *beta,
  T *td)
{
  std::memset(td, 0, k*k*sizeof(T));
  
  #pragma omp for simd
  for (int i=0; i<k; i++)
    td[i + k*i] = alpha[i];
  
  for (int i=0; i<k-1; i++)
  {
    td[i + k*(i+1)] = beta[i];
    td[i+1 + k*i] = beta[i];
  }
}


#endif
