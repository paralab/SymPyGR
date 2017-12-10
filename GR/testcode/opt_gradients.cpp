/*
#include "eqtest.h"
*/
#include <iostream>
#include <dollar.h>
#include <cstdlib>
#include <random>

// #if defined SSE3 || defined AVX
#include <immintrin.h>

#define AVX_SIMD_LENGTH 4
#define VECTOR_LENGTH 8
#define ALIGNMENT 32
// #endif

int apply_stencil_x(const double* in, double* out, int *sz, double h) { $
  // du = (u[pp-2] - 8.0 * u[pp-1] + 8.0 * u[pp+1] - u[pp+2]) / (12.0*h);
  double fac1 = 1.0/(12.0*h);
  double fac2 = 8.0/(12.0*h);

  int idx=2*sz[1]*sz[0];
  for (int k=2; k<sz[2]-2; ++k) {
    idx += 2*sz[0];
    for (int j=2; j<sz[1]-2; ++j) {
      idx += 2;
      for (int i=2; i<sz[0]-2; ++i) {
        out[idx] = fac1*(in[idx-2] - in[idx+2]) + fac2*(in[idx+1] - in[idx-1]) ;
        idx++;
      } // i
      idx += 2;
    } // j 
    idx += 2*sz[0];
  } // k
}

int apply_stencil_z(const double* in, double* out, int *sz, double h) { $
  // du = (u[pp-2*nx*ny] - 8.0 * u[pp-1*nx*ny] + 8.0 * u[pp+1*nx*ny] - u[pp+2*nx*ny]) / (12.0*h);
  double fac1 = 1.0/(12.0*h);
  double fac2 = 8.0/(12.0*h);
  int n = sz[0]*sz[1];

  int idx=2*sz[1]*sz[0];
  for (int k=2; k<sz[2]-2; ++k) {
    idx += 2*sz[0];
    for (int j=2; j<sz[1]-2; ++j) {
      idx += 2;
      for (int i=2; i<sz[0]-2; ++i) {
        out[idx] = fac1*(in[idx - (n<<1)] - in[idx+(n<<1)]) + fac2*(in[idx+n] - in[idx-n]) ; 
        idx++;
      } // i
      idx += 2;
    } // j 
    idx += 2*sz[0];
  } // k
}


int apply_stencil_y(const double* in, double* out, int *sz, double h) { $
  // du = (u[pp-2*nx] - 8.0 * u[pp-1*nx] + 8.0 * u[pp+1*nx] - u[pp+2*nx]) / (12.0*h);
  double fac1 = 1.0/(12.0*h);
  double fac2 = 8.0/(12.0*h);
  int n = sz[0];

  int idx=2*sz[1]*sz[0];
  for (int k=2; k<sz[2]-2; ++k) {
    idx += 2*sz[0];
    for (int j=2; j<sz[1]-2; ++j) {
      idx += 2;
      for (int i=2; i<sz[0]-2; ++i) {
        out[idx] = fac1*(in[idx - (n<<1)] - in[idx+(n<<1)]) + fac2*(in[idx+n] - in[idx-n]) ; 
        idx++;
      } // i
      idx += 2;
    } // j 
    idx += 2*sz[0];
  } // k
}

// #if defined SSE3 || defined AVX

int avx_apply_stencil_x(const double* in, double* out, int *sz, double h) { $
  double fac1 = 1.0/(12.0*h);
  double fac2 = 8.0/(12.0*h);

  __m256d _f1 = _mm256_broadcast_sd(&fac1);
  __m256d _f2 = _mm256_broadcast_sd(&fac2);

  __m256d _o, _im1, _im2, _ip1, _ip2;
  
  int idx=2*sz[1]*sz[0];
  for (int k=2; k<sz[2]-2; ++k) {
    idx += 2*sz[0];
    for (int j=2; j<sz[1]-2; ++j) {
      idx += 2;
      for (int i=2; i<sz[0]-2; i+=AVX_SIMD_LENGTH) {
  
        _im2 = _mm256_loadu_pd( in + idx - 2);
        _im1 = _mm256_loadu_pd( in + idx - 1);
        _ip1 = _mm256_loadu_pd( in + idx + 1);
        _ip2 = _mm256_loadu_pd( in + idx + 2);
        
        _o = _mm256_setzero_pd(); 
        _o = _mm256_fmadd_pd(_im2, _f1, _o);
        _o = _mm256_fnmadd_pd(_im1, _f2, _o);
        _o = _mm256_fmadd_pd(_ip1, _f2, _o);
        _o = _mm256_fnmadd_pd(_ip2, _f1, _o);
        
        _mm256_stream_pd(out+idx, _o);
        // _mm256_store_pd(out+idx, _o);

        idx+=AVX_SIMD_LENGTH;
      } // i
      idx += 2;
    } // j 
    idx += 2*sz[0];
  } // k
}

int avx_apply_stencil_y(const double* in, double* out, int *sz, double h) { $
  double fac1 = 1.0/(12.0*h);
  double fac2 = 8.0/(12.0*h);
  int n = sz[0];
  int n2 = 2*sz[0];
  
  std::cout << "f1: " << fac1 << ", f2: " << fac2 << std::endl;

  __m256d _f1 = _mm256_broadcast_sd(&fac1);
  __m256d _f2 = _mm256_broadcast_sd(&fac2);

  __m256d _o, _im1, _im2, _ip1, _ip2;
  int idx=2*sz[1]*sz[0];
  for (int k=2; k<sz[2]-2; ++k) {
    idx += 2*sz[0];
    for (int j=2; j<sz[1]-2; ++j) {
      idx += 2;
      for (int i=2; i<sz[0]-2; i+=AVX_SIMD_LENGTH) {
        /*
        _im2 = _mm256_loadu_pd( in + idx - n2);
        _im1 = _mm256_loadu_pd( in + idx - n);
        _ip1 = _mm256_loadu_pd( in + idx + n);
        _ip2 = _mm256_loadu_pd( in + idx + n2);
        
        _o = _mm256_setzero_pd(); 
        _mm256_fmadd_pd(_im2, _f1, _o);
        _mm256_fnmadd_pd(_im1, _f2, _o);
        _mm256_fmadd_pd(_ip1, _f2, _o);
        _mm256_fnmadd_pd(_ip2, _f1, _o);
        */
        _o = _mm256_load_pd(in + idx);
        _mm256_stream_pd(out + idx, _o);

        idx += AVX_SIMD_LENGTH;
      } // i
      idx += 2;
    } // j 
    idx += 2*sz[0];
  } // k
}
// #endif 

int avx_apply_stencil_z(const double* in, double* out, int *sz, double h) { $
  double fac1 = 1.0/(12.0*h);
  double fac2 = 8.0/(12.0*h);
  int n = sz[0]*sz[1];
  int n2 = 2*n;

  __m256d _f1 = _mm256_broadcast_sd(&fac1);
  __m256d _f2 = _mm256_broadcast_sd(&fac2);

  __m256d _o, _im1, _im2, _ip1, _ip2;
  int idx=n2;
  for (int k=2; k<sz[2]-2; ++k) {
    idx += 2*sz[0];
    for (int j=2; j<sz[1]-2; ++j) {
      idx += 2;
      for (int i=2; i<sz[0]-2; i+=AVX_SIMD_LENGTH) {
        _im2 = _mm256_loadu_pd( in + idx - n2);
        _im1 = _mm256_loadu_pd( in + idx - n);
        _ip1 = _mm256_loadu_pd( in + idx + n);
        _ip2 = _mm256_loadu_pd( in + idx + n2);
        
        _o = _mm256_setzero_pd(); 
        _mm256_fmadd_pd(_im2, _f1, _o);
        _mm256_fnmadd_pd(_im1, _f2, _o);
        _mm256_fmadd_pd(_ip1, _f2, _o);
        _mm256_fnmadd_pd(_ip2, _f1, _o);

        _mm256_stream_pd(out + idx, _o);

        idx += AVX_SIMD_LENGTH;
      } // i
      idx += 2;
    } // j 
    idx += 2*sz[0];
  } // k
}

int avx_apply_stencil_ddx(const double* in, double* out, int *sz, double h) { $
  double fac3 = -1.0/(12.0*h*h);
  double fac4 = 16.0/(12.0*h*h);
  double fac5 = -30.0/(12.0*h*h);
  __m256d _f3 = _mm256_broadcast_sd(&fac3);
  __m256d _f4 = _mm256_broadcast_sd(&fac4);
  __m256d _f5 = _mm256_broadcast_sd(&fac5);

  __m256d _ddo, _im1, _im2, _ip0, _ip1, _ip2;
  
  int idx=2*sz[1]*sz[0];
  for (int k=2; k<sz[2]-2; ++k) {
    idx += 2*sz[0];
    for (int j=2; j<sz[1]-2; ++j) {
      idx += 2;
      for (int i=2; i<sz[0]-2; i+=AVX_SIMD_LENGTH) {
  
        _im2 = _mm256_loadu_pd( in + idx - 2);
        _im1 = _mm256_loadu_pd( in + idx - 1);
        _ip0 = _mm256_load_pd( in + idx);
        _ip1 = _mm256_loadu_pd( in + idx + 1);
        _ip2 = _mm256_loadu_pd( in + idx + 2);
        
        _ddo = _mm256_setzero_pd(); 
        _mm256_fmadd_pd(_im2, _f3, _ddo);
        _mm256_fmadd_pd(_im1, _f4, _ddo);
        _mm256_fmadd_pd(_ip0, _f5, _ddo);
        _mm256_fmadd_pd(_ip1, _f4, _ddo);
        _mm256_fmadd_pd(_ip2, _f3, _ddo);
        
        _mm256_stream_pd(out+idx, _ddo);

        idx+=AVX_SIMD_LENGTH;
      } // i
      idx += 2;
    } // j 
    idx += 2*sz[0];
  } // k
}

int avx_apply_stencil_d_ddx(const double* in, double* out, double *const ddout, int *sz, double h) { $
  double fac1 = 1.0/(12.0*h);
  double fac2 = 8.0/(12.0*h);

  __m256d _f1 = _mm256_broadcast_sd(&fac1);
  __m256d _f2 = _mm256_broadcast_sd(&fac2);
  
  double fac3 = -1.0/(12.0*h*h);
  double fac4 = 16.0/(12.0*h*h);
  double fac5 = -30.0/(12.0*h*h);
  __m256d _f3 = _mm256_broadcast_sd(&fac3);
  __m256d _f4 = _mm256_broadcast_sd(&fac4);
  __m256d _f5 = _mm256_broadcast_sd(&fac5);

  __m256d _o, _ddo, _im1, _im2, _ip0, _ip1, _ip2;
  
  int idx=2*sz[1]*sz[0];
  for (int k=2; k<sz[2]-2; ++k) {
    idx += 2*sz[0];
    for (int j=2; j<sz[1]-2; ++j) {
      idx += 2;
      for (int i=2; i<sz[0]-2; i+=AVX_SIMD_LENGTH) {
  
        _im2 = _mm256_loadu_pd( in + idx - 2);
        _im1 = _mm256_loadu_pd( in + idx - 1);
        _ip0 = _mm256_load_pd( in + idx);
        _ip1 = _mm256_loadu_pd( in + idx + 1);
        _ip2 = _mm256_loadu_pd( in + idx + 2);
       
        _o = _mm256_setzero_pd(); 
        _mm256_fmadd_pd(_im2, _f1, _o);
        _mm256_fnmadd_pd(_im1, _f2, _o);
        _mm256_fmadd_pd(_ip1, _f2, _o);
        _mm256_fnmadd_pd(_ip2, _f1, _o);
        
        _mm256_stream_pd(out + idx, _o);
        
        _ddo = _mm256_setzero_pd(); 
        _mm256_fmadd_pd(_im2, _f3, _ddo);
        _mm256_fmadd_pd(_im1, _f4, _ddo);
        _mm256_fmadd_pd(_ip0, _f5, _ddo);
        _mm256_fmadd_pd(_ip1, _f4, _ddo);
        _mm256_fmadd_pd(_ip2, _f3, _ddo);
        
        _mm256_stream_pd(ddout+idx, _ddo);

        idx+=AVX_SIMD_LENGTH;
      } // i
      idx += 2;
    } // j 
    idx += 2*sz[0];
  } // k
}

int avx_apply_stencil_d_ddy(const double* in, double* const out, 
                          double *const ddout, int *sz, double h) { $
  // du = (u[pp-2*nx] - 8.0 * u[pp-1*nx] + 8.0 * u[pp+1*nx] - u[pp+2*nx]) / (12.0*h);
  double fac1 = 1.0/(12.0*h);
  double fac2 = 8.0/(12.0*h);
  int n = sz[0];
  int n2 = 2*sz[0];

  __m256d _f1 = _mm256_broadcast_sd(&fac1);
  __m256d _f2 = _mm256_broadcast_sd(&fac2);

  double fac3 = -1.0/(12.0*h*h);
  double fac4 = 16.0/(12.0*h*h);
  double fac5 = -30.0/(12.0*h*h);
  __m256d _f3 = _mm256_broadcast_sd(&fac3);
  __m256d _f4 = _mm256_broadcast_sd(&fac4);
  __m256d _f5 = _mm256_broadcast_sd(&fac5);

  __m256d _o, _ddo, _im1, _im2, _ip1, _ip2, _ip0;
  int idx=2*sz[1]*sz[0];
  for (int k=2; k<sz[2]-2; ++k) {
    idx += 2*sz[0];
    for (int j=2; j<sz[1]-2; ++j) {
      idx += 2;
      for (int i=2; i<sz[0]-2; i+=AVX_SIMD_LENGTH) {
        // out[idx] = fac1*(in[idx - (n<<1)] - in[idx+(n<<1)]) + fac2*(in[idx+n] - fac2*in[idx-n]) ; 
        // idx++;
        // load 4
        _im2 = _mm256_load_pd( in + idx - n2);
        _im1 = _mm256_load_pd( in + idx - n);
        _ip0 = _mm256_load_pd( in + idx);
        _ip1 = _mm256_load_pd( in + idx + n);
        _ip2 = _mm256_load_pd( in + idx + n2);
        
        _o = _mm256_setzero_pd(); 
        _mm256_fmadd_pd(_im2, _f1, _o);
        _mm256_fnmadd_pd(_im1, _f2, _o);
        _mm256_fmadd_pd(_ip1, _f2, _o);
        _mm256_fnmadd_pd(_ip2, _f1, _o);

        _mm256_stream_pd(out + idx, _o);

        _ddo = _mm256_setzero_pd(); 
        _mm256_fmadd_pd(_im2, _f3, _ddo);
        _mm256_fmadd_pd(_im1, _f4, _ddo);
        _mm256_fmadd_pd(_ip0, _f5, _ddo);
        _mm256_fmadd_pd(_ip1, _f4, _ddo);
        _mm256_fmadd_pd(_ip2, _f3, _ddo);
        _mm256_stream_pd(ddout + idx, _ddo);

        idx += AVX_SIMD_LENGTH;
      } // i
      idx += 2;
    } // j 
    idx += 2*sz[0];
  } // k
}

int avx_apply_stencil_ddy(const double* in, double* const out, 
                          int *sz, double h) { $
  int n = sz[0];
  int n2 = 2*sz[0];

  double fac3 = -1.0/(12.0*h*h);
  double fac4 = 16.0/(12.0*h*h);
  double fac5 = -30.0/(12.0*h*h);
  __m256d _f3 = _mm256_broadcast_sd(&fac3);
  __m256d _f4 = _mm256_broadcast_sd(&fac4);
  __m256d _f5 = _mm256_broadcast_sd(&fac5);

  __m256d _ddo, _im1, _im2, _ip1, _ip2, _ip0;
  int idx=2*sz[1]*sz[0];
  for (int k=2; k<sz[2]-2; ++k) {
    idx += 2*sz[0];
    for (int j=2; j<sz[1]-2; ++j) {
      idx += 2;
      for (int i=2; i<sz[0]-2; i+=AVX_SIMD_LENGTH) {
        // out[idx] = fac1*(in[idx - (n<<1)] - in[idx+(n<<1)]) + fac2*(in[idx+n] - fac2*in[idx-n]) ; 
        // idx++;
        // load 4
        _im2 = _mm256_load_pd( in + idx - n2);
        _im1 = _mm256_load_pd( in + idx - n);
        _ip0 = _mm256_load_pd( in + idx);
        _ip1 = _mm256_load_pd( in + idx + n);
        _ip2 = _mm256_load_pd( in + idx + n2);
        
        _ddo = _mm256_setzero_pd(); 
        _mm256_fmadd_pd(_im2, _f3, _ddo);
        _mm256_fmadd_pd(_im1, _f4, _ddo);
        _mm256_fmadd_pd(_ip0, _f5, _ddo);
        _mm256_fmadd_pd(_ip1, _f4, _ddo);
        _mm256_fmadd_pd(_ip2, _f3, _ddo);
        _mm256_stream_pd(out + idx, _ddo);

        idx += AVX_SIMD_LENGTH;
      } // i
      idx += 2;
    } // j 
    idx += 2*sz[0];
  } // k
}


int avx_apply_stencil_d_ddz(const double* in, double* const out, 
                           double* const ddout, int *sz, double h) { $
  double fac1 = 1.0/(12.0*h);
  double fac2 = 8.0/(12.0*h);
  int n = sz[0]*sz[1];
  int n2 = 2*n;

  __m256d _f1 = _mm256_broadcast_sd(&fac1);
  __m256d _f2 = _mm256_broadcast_sd(&fac2);
  fac1 *= -1.0;
  fac2 *= -1.0;
  __m256d _mf1 = _mm256_broadcast_sd(&fac1);
  __m256d _mf2 = _mm256_broadcast_sd(&fac2);

  double fac3 = -1.0/12.0*h*h;
  double fac4 = 16.0/12.0*h*h;
  double fac5 = -30.0/12.0*h*h;
  __m256d _f3 = _mm256_broadcast_sd(&fac3);
  __m256d _f4 = _mm256_broadcast_sd(&fac4);
  __m256d _f5 = _mm256_broadcast_sd(&fac5);

  __m256d _o, _ddo, _data; // _im1, _im2, _ip1, _ip2, _ip0;
  int idx=n2;
  for (int k=2; k<sz[2]-2; ++k) {
    idx += 2*sz[0];
    for (int j=2; j<sz[1]-2; ++j) {
      idx += 2;
      for (int i=2; i<sz[0]-2; i+=AVX_SIMD_LENGTH) {
        // load 4
        _o = _mm256_setzero_pd(); 

#define ABC
#ifdef ABC
        _ddo = _mm256_setzero_pd();  // faster
#else
        _ddo = _mm256_load_pd( in + idx ); // slower
        _ddo = _mm256_mul_pd( _ddo, _f5 );
#endif
        
				_data = _mm256_load_pd( in + idx - n2);
				_mm256_fmadd_pd(_data, _f1, _o);
        _mm256_fmadd_pd(_data, _f3, _ddo);

        _data = _mm256_load_pd( in + idx - n);
        _mm256_fmadd_pd(_data, _mf2, _o);
        _mm256_fmadd_pd(_data, _f4, _ddo);

#ifdef ABC	
        _mm256_fmadd_pd(_data, _f5, _ddo);
#endif

        _data = _mm256_load_pd( in + idx + n);
        _mm256_fmadd_pd(_data, _f2, _o);
        _mm256_fmadd_pd(_data, _f4, _ddo);

        _data = _mm256_load_pd( in + idx + n2);
        _mm256_fmadd_pd(_data, _mf1, _o);
        _mm256_fmadd_pd(_data, _f3, _ddo);

        _mm256_stream_pd(out + idx, _o);
        _mm256_stream_pd(ddout + idx, _ddo);

        idx += AVX_SIMD_LENGTH;
      } // i
      idx += 2;
    } // j 
    idx += 2*sz[0];
  } // k
}

int avx_apply_stencil_ddz(const double* in, double* const out, 
                           int *sz, double h) { $
  int n = sz[0]*sz[1];
  int n2 = 2*n;

  double fac3 = -1.0/(12.0*h*h);
  double fac4 = 16.0/(12.0*h*h);
  double fac5 = -30.0/(12.0*h*h);
  __m256d _f3 = _mm256_broadcast_sd(&fac3);
  __m256d _f4 = _mm256_broadcast_sd(&fac4);
  __m256d _f5 = _mm256_broadcast_sd(&fac5);

  __m256d _o, _im1, _im2, _ip1, _ip2, _ip0;
  int idx=n2;
  for (int k=2; k<sz[2]-2; ++k) {
    idx += 2*sz[0];
    for (int j=2; j<sz[1]-2; ++j) {
      idx += 2;
      for (int i=2; i<sz[0]-2; i+=AVX_SIMD_LENGTH) {
        // load 4
        _im2 = _mm256_load_pd( in + idx - n2);
        _im1 = _mm256_load_pd( in + idx - n);
        _ip0 = _mm256_load_pd( in + idx);
        _ip1 = _mm256_load_pd( in + idx + n);
        _ip2 = _mm256_load_pd( in + idx + n2);
        
        _o = _mm256_setzero_pd(); 
        _mm256_fmadd_pd(_im2, _f3, _o);
        _mm256_fmadd_pd(_im1, _f4, _o);
        _mm256_fmadd_pd(_ip0, _f5, _o);
        _mm256_fmadd_pd(_ip1, _f4, _o);
        _mm256_fmadd_pd(_ip2, _f3, _o);
        _mm256_stream_pd(out + idx, _o);


        idx += AVX_SIMD_LENGTH;
      } // i
      idx += 2;
    } // j 
    idx += 2*sz[0];
  } // k
}

// advective derivatives 


int avx_apply_adv_stencil_x(const double* in, double* out, double* beta, int *sz, double h) { $
  int idx=2*sz[1]*sz[0];
  int mask, pp;
  double eka = 1.0;
  double dva = 2.0;
  __m256d zeros = _mm256_setzero_pd();
  __m256d one =  _mm256_broadcast_sd(&eka);
  __m256d two =  _mm256_broadcast_sd(&dva);
  __m256d _cmp;
  
  double fac1 = 0.25/h;
  double fac2 = 0.833333333333333/h;
  double fac3 = 1.5/h;
  double fac4 = 0.5/h;
  double fac5 = 0.0833333333333333/h;

  __m256d _f1 = _mm256_broadcast_sd(&fac1);
  __m256d _f2 = _mm256_broadcast_sd(&fac2);
  __m256d _f3 = _mm256_broadcast_sd(&fac3);
  __m256d _f4 = _mm256_broadcast_sd(&fac4);
  __m256d _f5 = _mm256_broadcast_sd(&fac5);


  __m256d _o, _b, _1mb, _2bm1, _bb;
  __m256d _im3, _im2, _im1, _i0, _ip1, _ip2, _ip3;
  for (int k=2; k<sz[2]-3; ++k) {
    idx += 2*sz[0];
    for (int j=2; j<sz[1]-3; ++j) {
      idx += 2;
      for (int i=2; i<sz[0]-3; i+=AVX_SIMD_LENGTH) {
        _b = _mm256_load_pd( beta + idx);
        _cmp = _mm256_cmp_pd(_b, zeros, _CMP_GE_OS);
        
				_1mb  = _mm256_sub_pd( one, _cmp);
        _2bm1 = _mm256_fmsub_pd(two, _cmp, one);

        _im3 = _mm256_loadu_pd( in + idx - 3);
        _im2 = _mm256_loadu_pd( in + idx - 2);
        _im1 = _mm256_loadu_pd( in + idx - 1);
        _i0 = _mm256_load_pd( in + idx );
        _ip1 = _mm256_loadu_pd( in + idx + 1);
        _ip2 = _mm256_loadu_pd( in + idx + 2);
        _ip3 = _mm256_loadu_pd( in + idx + 3);

        // the actual computation 
        _o = _mm256_setzero_pd();
        
        // i - 3 
        _im3 = _mm256_mul_pd(_im3, _1mb);
        _mm256_fmadd_pd(_im3, _f5, _o);
        // i - 2 
        _im2 = _mm256_mul_pd(_im2, _1mb);
        _mm256_fmadd_pd(_im2, _f4, _o);
        // i - 1 
        _bb = _mm256_mul_pd(_b, _f1);
        _mm256_fmadd_pd(_1mb, _f3, _bb);
        _mm256_fnmadd_pd(_im1, _bb, _o);
        
        // i
        _i0 = _mm256_mul_pd(_i0, _2bm1);
        _mm256_fmadd_pd(_i0, _f2, _o);
        // i + 1
        _bb = _mm256_mul_pd(_b, _f3);
        _mm256_fmadd_pd(_1mb, _f1, _bb);
        _mm256_fmadd_pd(_ip1, _bb, _o);
        
        // i + 2 
        _ip2 = _mm256_mul_pd(_ip2, _b);
        _mm256_fmadd_pd(_ip2, _f4, _o);
        // i + 3 
        _ip3 = _mm256_mul_pd(_ip3, _b);
        _mm256_fmadd_pd(_ip3, _f5, _o);
        
        _mm256_stream_pd(out+idx, _o);

        idx+=AVX_SIMD_LENGTH;
      } // i
      idx += 2;
    } // j
    idx += 2*sz[0];
  } // k
}

int avx_apply_adv_stencil_x_with_if(const double* in, double* out, double* beta, int *sz, double h) { $
  double fac1 = 0.25/h;
  double fac2 = 0.833333333333333/h;
  double fac3 = 1.5/h;
  double fac4 = 0.5/h;
  double fac5 = 0.0833333333333333/h;

  __m256d _f1 = _mm256_broadcast_sd(&fac1);
  __m256d _f2 = _mm256_broadcast_sd(&fac2);
  __m256d _f3 = _mm256_broadcast_sd(&fac3);
  __m256d _f4 = _mm256_broadcast_sd(&fac4);
  __m256d _f5 = _mm256_broadcast_sd(&fac5);

  __m256d _o, _i1, _i2, _i3, _i4, _i5, _b;
  
  int idx=2*sz[1]*sz[0];
  int mask, pp;
  for (int k=2; k<sz[2]-3; ++k) {
    idx += 2*sz[0];
    for (int j=2; j<sz[1]-3; ++j) {
      idx += 2;
      for (int i=2; i<sz[0]-3; i+=AVX_SIMD_LENGTH) {
        _b = _mm256_loadu_pd( beta + idx);
        mask = _mm256_movemask_pd( _b );
        if (!mask) {
          // all positive 
          // std::cout << idx << ": pos (" << i << ", " << j << ", " << k << ")" << std::endl; 
          // du = ( u[pp-1]*(-0.25) + u[pp  ]*(-0.833333333333333) + u[pp+1]*(1.5) + u[pp+2]*(-0.5) + u[pp+3]*(0.0833333333333333)) / h;
          _i1 = _mm256_loadu_pd( in + idx - 1);
          _i2 = _mm256_loadu_pd( in + idx);
          _i3 = _mm256_loadu_pd( in + idx + 1);
          _i4 = _mm256_loadu_pd( in + idx + 2);
          _i5 = _mm256_loadu_pd( in + idx + 3);

          _o = _mm256_setzero_pd(); 
          _mm256_fnmadd_pd(_i1, _f1, _o);
          _mm256_fnmadd_pd(_i2, _f2, _o);
          _mm256_fmadd_pd(_i3, _f3, _o);
          _mm256_fnmadd_pd(_i4, _f4, _o);
          _mm256_fmadd_pd(_i5, _f5, _o);
           
          _mm256_stream_pd(out+idx, _o);
        } else if (mask == 15) {
          // std::cout << idx << ": neg" << std::endl; 
          // all negative 
          // du = ( u[pp-3]*(-0.0833333333333333) + u[pp-2]*(0.5) + u[pp-1]*(-1.5) + u[pp  ]*(0.833333333333333) + u[pp+1]*(0.25)) / h;
          _i1 = _mm256_loadu_pd( in + idx + 1);
          _i2 = _mm256_loadu_pd( in + idx);
          _i3 = _mm256_loadu_pd( in + idx - 1);
          _i4 = _mm256_loadu_pd( in + idx - 2);
          _i5 = _mm256_loadu_pd( in + idx - 3);
          
          _o = _mm256_setzero_pd(); 
          _mm256_fnmadd_pd(_i5, _f5, _o);
          _mm256_fmadd_pd(_i4, _f4, _o);
          _mm256_fnmadd_pd(_i3, _f3, _o);
          _mm256_fmadd_pd(_i2, _f2, _o);
          _mm256_fmadd_pd(_i1, _f1, _o);
          
          _mm256_stream_pd(out+idx, _o);
        } else {
          // process all 4
          for (int pp=idx; pp<idx+4; ++pp) {
            if (beta[pp] >= 0.0)
              out[pp] = -in[pp-1]*fac1 - in[pp]*fac2 + in[pp+1]*fac3 - in[pp+2]*fac4 + in[pp+3]*fac4;
            else 
              out[pp] = -in[pp-3]*fac5 + in[pp-2]*fac4 - in[pp-1]*fac3 + in[pp]*fac2 + in[pp-1]*fac1;
          } // s
        } // else non-vec

        idx+=AVX_SIMD_LENGTH;
      } // i
      idx += 2;
    } // j 
    idx += 2*sz[0];
  } // k
}

int avx_apply_adv_stencil_y(const double* in, double* out, double* beta, int *sz, double h) { $
  int idx=2*sz[1]*sz[0];
  int mask, pp;
  double eka = 1.0;
  double dva = 2.0;
  __m256d zeros = _mm256_setzero_pd();
  __m256d one =  _mm256_broadcast_sd(&eka);
  __m256d two =  _mm256_broadcast_sd(&dva);
  __m256d _cmp;
  
  int n = sz[0];
  int n2 = 2*sz[0];
  int n3 = 3*sz[0];
  
  double fac1 = 0.25/h;
  double fac2 = 0.833333333333333/h;
  double fac3 = 1.5/h;
  double fac4 = 0.5/h;
  double fac5 = 0.0833333333333333/h;

  __m256d _f1 = _mm256_broadcast_sd(&fac1);
  __m256d _f2 = _mm256_broadcast_sd(&fac2);
  __m256d _f3 = _mm256_broadcast_sd(&fac3);
  __m256d _f4 = _mm256_broadcast_sd(&fac4);
  __m256d _f5 = _mm256_broadcast_sd(&fac5);

  __m256d _o, _b, _1mb, _2bm1, _bb;
  __m256d _im3, _im2, _im1, _i0, _ip1, _ip2, _ip3;
  for (int k=2; k<sz[2]-3; ++k) {
    idx += 2*sz[0];
    for (int j=2; j<sz[1]-3; ++j) {
      idx += 2;
      for (int i=2; i<sz[0]-3; i+=AVX_SIMD_LENGTH) {
        _b = _mm256_load_pd( beta + idx);
        _cmp = _mm256_cmp_pd(_b, zeros, _CMP_GE_OS);
        
				_1mb  = _mm256_sub_pd( one, _cmp);
        _2bm1 = _mm256_fmsub_pd(two, _cmp, one);

        _im3 = _mm256_load_pd( in + idx - n3);
        _im2 = _mm256_load_pd( in + idx - n2);
        _im1 = _mm256_load_pd( in + idx - n);
        _i0 = _mm256_load_pd( in + idx );
        _ip1 = _mm256_load_pd( in + idx + n);
        _ip2 = _mm256_load_pd( in + idx + n2);
        _ip3 = _mm256_load_pd( in + idx + n3);

        // the actual computation 
        _o = _mm256_setzero_pd();
        
        // i - 3 
        _im3 = _mm256_mul_pd(_im3, _1mb);
        _mm256_fmadd_pd(_im3, _f5, _o);
        // i - 2 
        _im2 = _mm256_mul_pd(_im2, _1mb);
        _mm256_fmadd_pd(_im2, _f4, _o);
        // i - 1 
        _bb = _mm256_mul_pd(_b, _f1);
        _mm256_fmadd_pd(_1mb, _f3, _bb);
        _mm256_fnmadd_pd(_im1, _bb, _o);
        
        // i
        _i0 = _mm256_mul_pd(_i0, _2bm1);
        _mm256_fmadd_pd(_i0, _f2, _o);
        // i + 1
        _bb = _mm256_mul_pd(_b, _f3);
        _mm256_fmadd_pd(_1mb, _f1, _bb);
        _mm256_fmadd_pd(_ip1, _bb, _o);
        
        // i + 2 
        _ip2 = _mm256_mul_pd(_ip2, _b);
        _mm256_fmadd_pd(_ip2, _f4, _o);
        // i + 3 
        _ip3 = _mm256_mul_pd(_ip3, _b);
        _mm256_fmadd_pd(_ip3, _f5, _o);
        
        _mm256_stream_pd(out+idx, _o);

        idx+=AVX_SIMD_LENGTH;
      } // i
      idx += 2;
    } // j
    idx += 2*sz[0];
  } // k
}

int avx_apply_adv_stencil_z(const double* in, double* out, double* beta, int *sz, double h) { $
  int idx=2*sz[1]*sz[0];
  int mask, pp;
  double eka = 1.0;
  double dva = 2.0;
  __m256d zeros = _mm256_setzero_pd();
  __m256d one =  _mm256_broadcast_sd(&eka);
  __m256d two =  _mm256_broadcast_sd(&dva);
  __m256d _cmp;
  
  int n = sz[0]*sz[1];
  int n2 = 2*n;
  int n3 = 3*n;
  
  double fac1 = 0.25/h;
  double fac2 = 0.833333333333333/h;
  double fac3 = 1.5/h;
  double fac4 = 0.5/h;
  double fac5 = 0.0833333333333333/h;

  __m256d _f1 = _mm256_broadcast_sd(&fac1);
  __m256d _f2 = _mm256_broadcast_sd(&fac2);
  __m256d _f3 = _mm256_broadcast_sd(&fac3);
  __m256d _f4 = _mm256_broadcast_sd(&fac4);
  __m256d _f5 = _mm256_broadcast_sd(&fac5);

  __m256d _o, _b, _1mb, _2bm1, _bb;
  __m256d _im3, _im2, _im1, _i0, _ip1, _ip2, _ip3;
  for (int k=2; k<sz[2]-3; ++k) {
    idx += 2*sz[0];
    for (int j=2; j<sz[1]-3; ++j) {
      idx += 2;
      for (int i=2; i<sz[0]-3; i+=AVX_SIMD_LENGTH) {
        _b = _mm256_load_pd( beta + idx);
        _cmp = _mm256_cmp_pd(_b, zeros, _CMP_GE_OS);
        
				_1mb  = _mm256_sub_pd( one, _cmp);
        _2bm1 = _mm256_fmsub_pd(two, _cmp, one);

        _im3 = _mm256_load_pd( in + idx - n3);
        _im2 = _mm256_load_pd( in + idx - n2);
        _im1 = _mm256_load_pd( in + idx - n);
        _i0 = _mm256_load_pd( in + idx );
        _ip1 = _mm256_load_pd( in + idx + n);
        _ip2 = _mm256_load_pd( in + idx + n2);
        _ip3 = _mm256_load_pd( in + idx + n3);

        // the actual computation 
        _o = _mm256_setzero_pd();
        
        // i - 3 
        _im3 = _mm256_mul_pd(_im3, _1mb);
        _mm256_fmadd_pd(_im3, _f5, _o);
        // i - 2 
        _im2 = _mm256_mul_pd(_im2, _1mb);
        _mm256_fmadd_pd(_im2, _f4, _o);
        // i - 1 
        _bb = _mm256_mul_pd(_b, _f1);
        _mm256_fmadd_pd(_1mb, _f3, _bb);
        _mm256_fnmadd_pd(_im1, _bb, _o);
        
        // i
        _i0 = _mm256_mul_pd(_i0, _2bm1);
        _mm256_fmadd_pd(_i0, _f2, _o);
        // i + 1
        _bb = _mm256_mul_pd(_b, _f3);
        _mm256_fmadd_pd(_1mb, _f1, _bb);
        _mm256_fmadd_pd(_ip1, _bb, _o);
        
        // i + 2 
        _ip2 = _mm256_mul_pd(_ip2, _b);
        _mm256_fmadd_pd(_ip2, _f4, _o);
        // i + 3 
        _ip3 = _mm256_mul_pd(_ip3, _b);
        _mm256_fmadd_pd(_ip3, _f5, _o);
        
        _mm256_stream_pd(out+idx, _o);

        idx+=AVX_SIMD_LENGTH;
      } // i
      idx += 2;
    } // j
    idx += 2*sz[0];
  } // k
}

// void compute_derivatives(const double* u, const int *shp, double h, const double *beta, double** dd) {
//   double *dx = dd[0];
//   double *dy = dd[1];
//   double *dz = dd[2];
//
//   double *dxx = dd[3];
//   double *dyy = dd[4];
//   double *dzz = dd[5];
//
//   double *dxy = dd[6];
//   double *dyz = dd[7];
//   double *dzx = dd[8];
//
//   double *adx = dd[9];
//   double *ady = dd[10];
//   double *adz = dd[11];
//
//   // now compute all 12 derivatives ...
//   // -- default ordering is along x ...
//   
//   // order 
//   #<{(|******* ORDER ***************  
//    *  u
//    *  |_ [dx] _ [dxx]                     0,1
//    *  |_ [adx]                            2
//    *  |
//    *  |_ T _ dY _ dYY _ T _ [dyy]         4
//    *     |   |
//    *     |   |_ T _ [dy] _ [dyx]          6
//    *     |   |_ T _ dYZ _ T _ [dyz]       7
//    *     |
//    *     |_ adY _ T _ [ady]               8
//    *     |
//    *     |_ T _ dZ _ dZZ _ T _ [dzz]      9
//    *        |   |
//    *        |   |_ T _ [dz] _ [dzx]       11
//    *        |_ adZ _T _ [adz]             12
//    *
//    *  10 Transpose operations
//    *  12 Stencils
//    |)}># 
//   
//    //! x deriv.
//
// }

int
main(int argc, char* argv[]) {
  int sz[3] = {128, 128, 128};

  if (argc > 3) {
    sz[0] = std::atoi(argv[1]);
    sz[1] = std::atoi(argv[2]);
    sz[2] = std::atoi(argv[3]);
  }

  int n = sz[0]*sz[1]*sz[2];

  // allocate 
  // double *u  = new double[n];
  // double *ux = new double[n];
  // double *uy = new double[n];
  // double *uz = new double[n];
  // double *uxx = new double[n];
  // double *uyy = new double[n];
  // double *uzz = new double[n];
  // double *uxy = new double[n];
  // double *uyz = new double[n];
  // double *uzx = new double[n];

  double *u, *ux, *uy, *uz;
  double *uxx, *uyy, *uzz;
  double *uxy, *uyz, *uzx;
  double *beta;
  double *adx, *ady, *adz;
  size_t align = 32;

  posix_memalign((void **)&u, align, (n+AVX_SIMD_LENGTH)*sizeof(double)); 
  posix_memalign((void **)&ux, align, (n+AVX_SIMD_LENGTH)*sizeof(double)); 
  posix_memalign((void **)&uy, align, (n+AVX_SIMD_LENGTH)*sizeof(double)); 
  posix_memalign((void **)&uz, align, (n+AVX_SIMD_LENGTH)*sizeof(double)); 
  posix_memalign((void **)&uxx, align, (n+AVX_SIMD_LENGTH)*sizeof(double)); 
  posix_memalign((void **)&uyy, align, (n+AVX_SIMD_LENGTH)*sizeof(double)); 
  posix_memalign((void **)&uzz, align, (n+AVX_SIMD_LENGTH)*sizeof(double)); 
  posix_memalign((void **)&uxy, align, (n+AVX_SIMD_LENGTH)*sizeof(double)); 
  posix_memalign((void **)&uyz, align, (n+AVX_SIMD_LENGTH)*sizeof(double)); 
  posix_memalign((void **)&uzx, align, (n+AVX_SIMD_LENGTH)*sizeof(double)); 
  posix_memalign((void **)&beta, align, (n+AVX_SIMD_LENGTH)*sizeof(double)); 
  posix_memalign((void **)&adx, align, (n+AVX_SIMD_LENGTH)*sizeof(double)); 
  posix_memalign((void **)&ady, align, (n+AVX_SIMD_LENGTH)*sizeof(double)); 
  posix_memalign((void **)&adz, align, (n+AVX_SIMD_LENGTH)*sizeof(double)); 

  u+=2;
  ux+=2; uy+=2; uz+=2;
  uxx+=2; uyy+=2; uzz+=2;
  uxy+=2; uyz+=2; uzx+=2;
  adx += 2; ady += 2; adz += 2;
  beta += 2;

  double h = 1.0/sz[0];

  std::random_device rd; 
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 2.0);

  for (int i=0; i<n; ++i) { 
    u[i] = dis(gen);
    beta[i] = 1.0;
  }
  for (int i=256*128; i<256*130+40; ++i) { 
    beta[i] = -1.0;
  }

  beta[256*130+5] = 1.0;
  beta[256*130+9] = 1.0;
  beta[256*130+10] = 1.0;
  beta[256*130+13] = 1.0;
  beta[256*130+14] = 1.0;
  beta[256*130+15] = 1.0;
  

  for (int i=0; i<100; ++i) {
    // transpose_xz(u, uz, sz);
    // now the stencils ...
    /*
    {
      apply_stencil_x(u, ux, sz, h);
      apply_stencil_y(u, uy, sz, h);
      apply_stencil_z(u, uz, sz, h);
    }
    */
    // {
    //   apply_stencil_x(ux, uxx, sz, h);
    //   apply_stencil_y(uy, uyy, sz, h);
    //   apply_stencil_z(uz, uzz, sz, h);
    // }
    // {
    //   apply_stencil_x(uz, uzx, sz, h);
    //   apply_stencil_y(ux, uxy, sz, h);
    //   apply_stencil_z(uy, uyz, sz, h);
    // }
    // vectorized 
    {
      avx_apply_stencil_x(u, ux, sz, h);
      avx_apply_stencil_y(u, uy, sz, h);
      avx_apply_stencil_z(u, uz, sz, h);
    }
    // {
    //   avx_apply_stencil_x(ux, uxx, sz, h);
    //   avx_apply_stencil_y(uy, uyy, sz, h);
    //   avx_apply_stencil_z(uz, uzz, sz, h);
    // }
    // {
    //   avx_apply_stencil_x(uz, uzx, sz, h);
    //   avx_apply_stencil_y(ux, uxy, sz, h);
    //   avx_apply_stencil_z(uy, uyz, sz, h);
    // }
    /*
    {
      avx_apply_stencil_d_ddx(u, ux, uxx, sz, h);
      avx_apply_stencil_d_ddy(u, uy, uyy, sz, h);
      avx_apply_stencil_d_ddz(u, uz, uzz, sz, h);
      avx_apply_stencil_ddx(u, uxx, sz, h);
      avx_apply_stencil_ddy(u, uyy, sz, h);
      avx_apply_stencil_ddz(u, uzz, sz, h);
    }
    {
      avx_apply_adv_stencil_x(u, adx, beta, sz, h);
      avx_apply_adv_stencil_y(u, ady, beta, sz, h);
      avx_apply_adv_stencil_z(u, ady, beta, sz, h);
    }
    */

  }
  

  std::cout << "u = ";
  for (int i=0; i<20; ++i) 
    std::cout << u[50000+i] << " ";
  std::cout << std::endl;


  std::cout << "ux = ";
  for (int i=0; i<50; ++i) 
    std::cout << ux[50000+i] << " ";
  std::cout << std::endl;

  std::cout << "uy = ";
  for (int i=0; i<50; ++i) 
    std::cout << uy[50000+i] << " ";
  std::cout << std::endl;


  dollar::text(std::cout);                      
  dollar::clear();   

  // clean up
  // delete [] u;
  // delete [] ux;
  // delete [] uy;
  // delete [] uz;
  // delete [] uxx;
  // delete [] uyy;
  // delete [] uzz;
  // delete [] uxy;
  // delete [] uyz;
  // delete [] uzx;
  
  u-=2;
  ux-=2; uy-=2; uz-=2;
  uxx-=2; uyy-=2; uzz-=2;
  uxy-=2; uyz-=2; uzx-=2;
  beta -= 2;
  adx -=2; ady -= 2; adz -= 2;

  free(u);
  free(ux);
  free(uy);
  free(uz);
  free(uxx);
  free(uyy);
  free(uzz);
  free(uxy);
  free(uyz);
  free(uzx);
  free(adx);
  free(ady);
  free(adz);
  free(beta);
  return 0;
}


