#pragma once

/**
 * Optimized 4th order FD stencils
 *
 * @author Hari Sundar hari@cs.utah.edu
 * @date   14 May 2017
 *
 **/ 



#include <immintrin.h>
#define AVX_SIMD_LENGTH 4
#define ALIGNMENT 32


inline int dendro_dx(const double* in, double* out, int *sz, double h) {
  double fac1 = 1.0/12.0*h;
  double fac2 = 8.0/12.0*h;

  __m256d _f1 = _mm256_broadcast_sd(&fac1);
  __m256d _f2 = _mm256_broadcast_sd(&fac2);
  fac1 *= -1.0;
  fac2 *= -1.0;
  __m256d _mf1 = _mm256_broadcast_sd(&fac1);
  __m256d _mf2 = _mm256_broadcast_sd(&fac2);

  __m256d _o, _im1, _im2, _ip1, _ip2;
  
  int idx=2*sz[1]*sz[0];
  for (int k=2; k<sz[2]-2; ++k) {
    idx += 2*sz[0];
    for (int j=2; j<sz[1]-2; ++j) {
      idx += 2;
      for (int i=2; i<sz[0]-2; i+=AVX_SIMD_LENGTH) {
        _o = _mm256_setzero_pd(); 
  
        _im2 = _mm256_loadu_pd( in + idx - 2);
        _im1 = _mm256_loadu_pd( in + idx - 1);
        _ip1 = _mm256_loadu_pd( in + idx + 1);
        _ip2 = _mm256_loadu_pd( in + idx + 2);
        
        _o = _mm256_setzero_pd(); 
        _mm256_fmadd_pd(_im2, _f1, _o);
        _mm256_fmadd_pd(_im1, _mf2, _o);
        _mm256_fmadd_pd(_ip1, _f2, _o);
        _mm256_fmadd_pd(_ip2, _f1, _o);
        
        _mm256_store_pd(out+idx, _o);

        idx+=AVX_SIMD_LENGTH;
      } // i
      idx += 2;
    } // j 
    idx += 2*sz[0];
  } // k
}

int dendro_dy(const double* in, double* out, int *sz, double h) {
  double fac1 = 1.0/12.0*h;
  double fac2 = 8.0/12.0*h;
  int n = sz[0];
  int n2 = 2*sz[0];

  __m256d _f1 = _mm256_broadcast_sd(&fac1);
  __m256d _f2 = _mm256_broadcast_sd(&fac2);
  fac1 *= -1.0;
  fac2 *= -1.0;
  __m256d _mf1 = _mm256_broadcast_sd(&fac1);
  __m256d _mf2 = _mm256_broadcast_sd(&fac2);

  __m256d _o, _im1, _im2, _ip1, _ip2;
  int idx=2*sz[1]*sz[0];
  for (int k=2; k<sz[2]-2; ++k) {
    idx += 2*sz[0];
    for (int j=2; j<sz[1]-2; ++j) {
      idx += 2;
      for (int i=2; i<sz[0]-2; i+=AVX_SIMD_LENGTH) {
        // load 4
        _im2 = _mm256_load_pd( in + idx - n2);
        _im1 = _mm256_load_pd( in + idx - n);
        _ip1 = _mm256_load_pd( in + idx + n);
        _ip2 = _mm256_load_pd( in + idx + n2);
        
        _o = _mm256_setzero_pd(); 
        _mm256_fmadd_pd(_im2, _f1, _o);
        _mm256_fmadd_pd(_im1, _mf2, _o);
        _mm256_fmadd_pd(_ip1, _f2, _o);
        _mm256_fmadd_pd(_ip2, _f1, _o);

        _mm256_store_pd(out + idx, _o);

        idx += AVX_SIMD_LENGTH;
      } // i
      idx += 2;
    } // j 
    idx += 2*sz[0];
  } // k
}
// #endif 

int dendro_dz(const double* in, double* out, int *sz, double h) {
  double fac1 = 1.0/12.0*h;
  double fac2 = 8.0/12.0*h;
  int n = sz[0]*sz[1];
  int n2 = 2*n;

  __m256d _f1 = _mm256_broadcast_sd(&fac1);
  __m256d _f2 = _mm256_broadcast_sd(&fac2);
  fac1 *= -1.0;
  fac2 *= -1.0;
  __m256d _mf1 = _mm256_broadcast_sd(&fac1);
  __m256d _mf2 = _mm256_broadcast_sd(&fac2);

  __m256d _o, _im1, _im2, _ip1, _ip2;
  int idx=n2;
  for (int k=2; k<sz[2]-2; ++k) {
    idx += 2*sz[0];
    for (int j=2; j<sz[1]-2; ++j) {
      idx += 2;
      for (int i=2; i<sz[0]-2; i+=AVX_SIMD_LENGTH) {
        // load 4
        _im2 = _mm256_load_pd( in + idx - n2);
        _im1 = _mm256_load_pd( in + idx - n);
        _ip1 = _mm256_load_pd( in + idx + n);
        _ip2 = _mm256_load_pd( in + idx + n2);
        
        _o = _mm256_setzero_pd(); 
        _mm256_fmadd_pd(_im2, _f1, _o);
        _mm256_fmadd_pd(_im1, _mf2, _o);
        _mm256_fmadd_pd(_ip1, _f2, _o);
        _mm256_fmadd_pd(_ip2, _f1, _o);

        _mm256_store_pd(out + idx, _o);

        idx += AVX_SIMD_LENGTH;
      } // i
      idx += 2;
    } // j 
    idx += 2*sz[0];
  } // k
}

inline int dendro_adx(const double* in, double* out, int *sz, double h) {
  
  double fac1 = 1.0/12.0*h;
  double fac2 = 8.0/12.0*h;

  __m256d _f1 = _mm256_broadcast_sd(&fac1);
  __m256d _f2 = _mm256_broadcast_sd(&fac2);
  fac1 *= -1.0;
  fac2 *= -1.0;
  __m256d _mf1 = _mm256_broadcast_sd(&fac1);
  __m256d _mf2 = _mm256_broadcast_sd(&fac2);

  __m256d _o, _im1, _im2, _ip1, _ip2;
  
  int idx=2*sz[1]*sz[0];
  for (int k=2; k<sz[2]-2; ++k) {
    idx += 2*sz[0];
    for (int j=2; j<sz[1]-2; ++j) {
      idx += 2;
      for (int i=2; i<sz[0]-2; i+=AVX_SIMD_LENGTH) {
        _o = _mm256_setzero_pd(); 
  
        _im2 = _mm256_loadu_pd( in + idx - 2);
        _im1 = _mm256_loadu_pd( in + idx - 1);
        _ip1 = _mm256_loadu_pd( in + idx + 1);
        _ip2 = _mm256_loadu_pd( in + idx + 2);
        
        _o = _mm256_setzero_pd(); 
        _mm256_fmadd_pd(_im2, _f1, _o);
        _mm256_fmadd_pd(_im1, _mf2, _o);
        _mm256_fmadd_pd(_ip1, _f2, _o);
        _mm256_fmadd_pd(_ip2, _f1, _o);
        
        _mm256_store_pd(out+idx, _o);

        idx+=AVX_SIMD_LENGTH;
      } // i
      idx += 2;
    } // j 
    idx += 2*sz[0];
  } // k
}

int dendro_ady(const double* in, const double* beta, double* out, int *sz, double h) {
  const __m256 zeros = _mm256_setzero_ps(); 
  // __m256d _mm256_cmp_pd(__m256d m1, __m256d m2, _CMP_GE_OS);
  // int _mm256_movemask_pd(__m256d a); // check sign
  
  // load beta and check ... 
  

  double fac1 = 1.0/12.0*h;
  double fac2 = 8.0/12.0*h;
  int n = sz[0];
  int n2 = 2*sz[0];

  __m256d _f1 = _mm256_broadcast_sd(&fac1);
  __m256d _f2 = _mm256_broadcast_sd(&fac2);
  fac1 *= -1.0;
  fac2 *= -1.0;
  __m256d _mf1 = _mm256_broadcast_sd(&fac1);
  __m256d _mf2 = _mm256_broadcast_sd(&fac2);

  __m256d _o, _im1, _im2, _ip1, _ip2;
  int idx=2*sz[1]*sz[0];
  for (int k=2; k<sz[2]-2; ++k) {
    idx += 2*sz[0];
    for (int j=2; j<sz[1]-2; ++j) {
      idx += 2;
      for (int i=2; i<sz[0]-2; i+=AVX_SIMD_LENGTH) {
        // load 4
        _im2 = _mm256_load_pd( in + idx - n2);
        _im1 = _mm256_load_pd( in + idx - n);
        _ip1 = _mm256_load_pd( in + idx + n);
        _ip2 = _mm256_load_pd( in + idx + n2);
        
        _o = _mm256_setzero_pd(); 
        _mm256_fmadd_pd(_im2, _f1, _o);
        _mm256_fmadd_pd(_im1, _mf2, _o);
        _mm256_fmadd_pd(_ip1, _f2, _o);
        _mm256_fmadd_pd(_ip2, _f1, _o);

        _mm256_store_pd(out + idx, _o);

        idx += AVX_SIMD_LENGTH;
      } // i
      idx += 2;
    } // j 
    idx += 2*sz[0];
  } // k
}
// #endif 

int dendro_adz(const double* in, double* out, int *sz, double h) {
  double fac1 = 1.0/12.0*h;
  double fac2 = 8.0/12.0*h;
  int n = sz[0]*sz[1];
  int n2 = 2*n;

  __m256d _f1 = _mm256_broadcast_sd(&fac1);
  __m256d _f2 = _mm256_broadcast_sd(&fac2);
  fac1 *= -1.0;
  fac2 *= -1.0;
  __m256d _mf1 = _mm256_broadcast_sd(&fac1);
  __m256d _mf2 = _mm256_broadcast_sd(&fac2);

  __m256d _o, _im1, _im2, _ip1, _ip2;
  int idx=n2;
  for (int k=2; k<sz[2]-2; ++k) {
    idx += 2*sz[0];
    for (int j=2; j<sz[1]-2; ++j) {
      idx += 2;
      for (int i=2; i<sz[0]-2; i+=AVX_SIMD_LENGTH) {
        // load 4
        _im2 = _mm256_load_pd( in + idx - n2);
        _im1 = _mm256_load_pd( in + idx - n);
        _ip1 = _mm256_load_pd( in + idx + n);
        _ip2 = _mm256_load_pd( in + idx + n2);
        
        _o = _mm256_setzero_pd(); 
        _mm256_fmadd_pd(_im2, _f1, _o);
        _mm256_fmadd_pd(_im1, _mf2, _o);
        _mm256_fmadd_pd(_ip1, _f2, _o);
        _mm256_fmadd_pd(_ip2, _f1, _o);

        _mm256_store_pd(out + idx, _o);

        idx += AVX_SIMD_LENGTH;
      } // i
      idx += 2;
    } // j 
    idx += 2*sz[0];
  } // k
}

