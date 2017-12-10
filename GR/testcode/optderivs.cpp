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
  // du = (u[pp-2] - 8.0 * u[pp-1] + 8.0 * u[pp+1] - u[pp+2]) / (12.0*h);
  double fac1 = 1.0/(12.0*h);
  double fac2 = 8.0/(12.0*h);

  __m256d _f1 = _mm256_broadcast_sd(&fac1);
  __m256d _f2 = _mm256_broadcast_sd(&fac2);
  // fac1 *= -1.0;
  // fac2 *= -1.0;
  // __m256d _mf1 = _mm256_broadcast_sd(&fac1);
  // __m256d _mf2 = _mm256_broadcast_sd(&fac2);

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

        idx+=AVX_SIMD_LENGTH;
      } // i
      idx += 2;
    } // j
    idx += 2*sz[0];
  } // k
}

int avx_apply_stencil_y(const double* in, double* out, int *sz, double h) { $
  // du = (u[pp-2*nx] - 8.0 * u[pp-1*nx] + 8.0 * u[pp+1*nx] - u[pp+2*nx]) / (12.0*h);
  double fac1 = 1.0/(12.0*h);
  double fac2 = 8.0/(12.0*h);
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
        // out[idx] = fac1*(in[idx - (n<<1)] - in[idx+(n<<1)]) + fac2*(in[idx+n] - fac2*in[idx-n]) ;
        // idx++;
        // load 4
        _im2 = _mm256_load_pd( in + idx - n2);
        _im1 = _mm256_load_pd( in + idx - n);
        _ip1 = _mm256_load_pd( in + idx + n);
        _ip2 = _mm256_load_pd( in + idx + n2);

        _o = _mm256_setzero_pd();
        _o = _mm256_fmadd_pd(_im2, _f1, _o);
        _o = _mm256_fmadd_pd(_im1, _mf2, _o);
        _o = _mm256_fmadd_pd(_ip1, _f2, _o);
        _o = _mm256_fmadd_pd(_ip2, _mf1, _o);

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
        _o = _mm256_fmadd_pd(_im2, _f1, _o);
        _o = _mm256_fmadd_pd(_im1, _mf2, _o);
        _o = _mm256_fmadd_pd(_ip1, _f2, _o);
        _o = _mm256_fmadd_pd(_ip2, _mf1, _o);

        _mm256_stream_pd(out + idx, _o);

        idx += AVX_SIMD_LENGTH;
      } // i
      idx += 2;
    } // j
    idx += 2*sz[0];
  } // k
}

int avx_apply_stencil_d_ddx(const double* const in, double* const out,
                             double* const ddout, int *sz, double h) { $
  double fac1 = 1.0/(12.0*h);
  double fac2 = 8.0/(12.0*h);

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

  __m256d _o, _ddo, _im1, _im2, _ip1, _ip2, _ip0;

  int idx=2*sz[1]*sz[0];
  for (int k=2; k<sz[2]-2; ++k) {
    idx += 2*sz[0];
    for (int j=2; j<sz[1]-2; ++j) {
      idx += 2;
      for (int i=2; i<sz[0]-2; i+=AVX_SIMD_LENGTH) {
        // out[idx] = fac1*(in[idx-2] - in[idx+2]) + fac2*(in[idx+1] - fac2*in[idx-1]) ;
        _o = _mm256_setzero_pd();

        _im2 = _mm256_loadu_pd( in + idx - 2);
        _im1 = _mm256_loadu_pd( in + idx - 1);
        _ip0 = _mm256_loadu_pd( in + idx);
        _ip1 = _mm256_loadu_pd( in + idx + 1);
        _ip2 = _mm256_loadu_pd( in + idx + 2);

        _o = _mm256_setzero_pd();
        _o = _mm256_fmadd_pd(_im2, _f1, _o);
        _o = _mm256_fmadd_pd(_im1, _mf2, _o);
        _o = _mm256_fmadd_pd(_ip1, _f2, _o);
        _o = _mm256_fmadd_pd(_ip2, _mf1, _o);
        _mm256_stream_pd(out+idx, _o);

        _ddo = _mm256_setzero_pd();
        _ddo = _mm256_fmadd_pd(_im2, _f3, _ddo);
        _ddo = _mm256_fmadd_pd(_im1, _f4, _ddo);
        _ddo = _mm256_fmadd_pd(_ip0, _f5, _ddo);
        _ddo = _mm256_fmadd_pd(_ip1, _f4, _ddo);
        _ddo = _mm256_fmadd_pd(_ip2, _f3, _ddo);
        _mm256_stream_pd(ddout + idx, _ddo);


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
        _im2 = _mm256_loadu_pd( in + idx - n2);
        _im1 = _mm256_loadu_pd( in + idx - n);
        _ip0 = _mm256_loadu_pd( in + idx);
        _ip1 = _mm256_loadu_pd( in + idx + n);
        _ip2 = _mm256_loadu_pd( in + idx + n2);

        _o = _mm256_setzero_pd();
        _o = _mm256_fmadd_pd(_im2, _f1, _o);
        _o = _mm256_fnmadd_pd(_im1, _f2, _o);
        _o = _mm256_fmadd_pd(_ip1, _f2, _o);
        _o = _mm256_fnmadd_pd(_ip2, _f1, _o);

        _mm256_stream_pd(out + idx, _o);

        _ddo = _mm256_setzero_pd();
        _ddo = _mm256_fmadd_pd(_im2, _f3, _ddo);
        _ddo = _mm256_fmadd_pd(_im1, _f4, _ddo);
        _ddo = _mm256_fmadd_pd(_ip0, _f5, _ddo);
        _ddo = _mm256_fmadd_pd(_ip1, _f4, _ddo);
        _ddo = _mm256_fmadd_pd(_ip2, _f3, _ddo);

        _mm256_stream_pd(ddout + idx, _ddo);

        idx += AVX_SIMD_LENGTH;
      } // i
      idx += 2;
    } // j
    idx += 2*sz[0];
  } // k
}

 /*

int avx_apply_stencil_d_ddy(const double* in, double* const out,
                          double *const ddout, int *sz, double h) { $
  double fac1 = 1.0/12.0*h;
  double fac2 = 8.0/12.0*h;
  int n = sz[0];
  int n2 = 2*sz[0];

  __m256d _f1 = _mm256_broadcast_sd(&fac1);
  __m256d _f2 = _mm256_broadcast_sd(&fac2);

  double fac3 = -1.0/12.0*h*h;
  double fac4 = 16.0/12.0*h*h;
  double fac5 = -30.0/12.0*h*h;
  __m256d _f3 = _mm256_broadcast_sd(&fac3);
  __m256d _f4 = _mm256_broadcast_sd(&fac4);
  __m256d _f5 = _mm256_broadcast_sd(&fac5);

  __m256d _o, _ddo, _im1, _im2, _ip1, _ip2, _ip0;
  int idx, idx2; // =2*sz[1]*sz[0];
  for (int k=2; k<sz[2]-2; ++k) {
    idx2 = k*sz[0]*sz[1];

    for (int i=2; i<sz[0]-2; i+=AVX_SIMD_LENGTH) {
      idx = idx2 + n2 + i;

      _im2 = _mm256_load_pd( in + idx - n2);
      _im1 = _mm256_load_pd( in + idx - n);
      _ip0 = _mm256_load_pd( in + idx);
      _ip1 = _mm256_load_pd( in + idx + n);
      for (int j=2; j<sz[1]-2; ++j) {

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

        // update registers
        _im2 = _im1;
        _im1 = _ip0;
        _ip0 = _ip1;
        _ip1 = _ip2;

        idx += n;
      } // j
    } // i
  } // k
}
*/

int avx_apply_stencil_ddy(const double* in, double* const out,
                          int *sz, double h) { $
  // du = (u[pp-2*nx] - 8.0 * u[pp-1*nx] + 8.0 * u[pp+1*nx] - u[pp+2*nx]) / (12.0*h);
  double fac1 = 1.0/(12.0*h);
  double fac2 = 8.0/(12.0*h);
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
        _ddo = _mm256_fmadd_pd(_im2, _f3, _ddo);
        _ddo = _mm256_fmadd_pd(_im1, _f4, _ddo);
        _ddo = _mm256_fmadd_pd(_ip0, _f5, _ddo);
        _ddo = _mm256_fmadd_pd(_ip1, _f4, _ddo);
        _ddo = _mm256_fmadd_pd(_ip2, _f3, _ddo);
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

  double fac3 = -1.0/(12.0*h*h);
  double fac4 = 16.0/(12.0*h*h);
  double fac5 = -30.0/(12.0*h*h);
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
        _ddo = _mm256_loadu_pd( in + idx ); // slower
        _ddo = _mm256_mul_pd( _ddo, _f5 );
#endif

				_data = _mm256_loadu_pd( in + idx - n2);
				_o = _mm256_fmadd_pd(_data, _f1, _o);
        _ddo = _mm256_fmadd_pd(_data, _f3, _ddo);

        _data = _mm256_loadu_pd( in + idx - n);
        _o = _mm256_fmadd_pd(_data, _mf2, _o);
        _ddo = _mm256_fmadd_pd(_data, _f4, _ddo);

#ifdef ABC
        _ddo = _mm256_fmadd_pd(_data, _f5, _ddo);
#endif

        _data = _mm256_loadu_pd( in + idx + n);
        _o = _mm256_fmadd_pd(_data, _f2, _o);
        _ddo = _mm256_fmadd_pd(_data, _f4, _ddo);

        _data = _mm256_loadu_pd( in + idx + n2);
        _o = _mm256_fmadd_pd(_data, _mf1, _o);
        _ddo = _mm256_fmadd_pd(_data, _f3, _ddo);

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
        _o = _mm256_fmadd_pd(_im2, _f3, _o);
        _o = _mm256_fmadd_pd(_im1, _f4, _o);
        _o = _mm256_fmadd_pd(_ip0, _f5, _o);
        _o = _mm256_fmadd_pd(_ip1, _f4, _o);
        _o = _mm256_fmadd_pd(_ip2, _f3, _o);
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
        _b = _mm256_loadu_pd( beta + idx);
        _cmp = _mm256_cmp_pd(_b, zeros, _CMP_GE_OS);

				_1mb  = _mm256_sub_pd( one, _cmp);
        _2bm1 = _mm256_fmsub_pd(two, _cmp, one);

        _im3 = _mm256_loadu_pd( in + idx - 3);
        _im2 = _mm256_loadu_pd( in + idx - 2);
        _im1 = _mm256_loadu_pd( in + idx - 1);
        _i0 = _mm256_loadu_pd( in + idx );
        _ip1 = _mm256_loadu_pd( in + idx + 1);
        _ip2 = _mm256_loadu_pd( in + idx + 2);
        _ip3 = _mm256_loadu_pd( in + idx + 3);

        // the actual computation
        _o = _mm256_setzero_pd();

        // i - 3
        _im3 = _mm256_mul_pd(_im3, _1mb);
        _o = _mm256_fmadd_pd(_im3, _f5, _o);
        // i - 2
        _im2 = _mm256_mul_pd(_im2, _1mb);
        _o = _mm256_fmadd_pd(_im2, _f4, _o);
        // i - 1
        _bb = _mm256_mul_pd(_b, _f1);
        _bb = _mm256_fmadd_pd(_1mb, _f3, _bb);
        _o = _mm256_fnmadd_pd(_im1, _bb, _o);

        // i
        _i0 = _mm256_mul_pd(_i0, _2bm1);
        _o = _mm256_fmadd_pd(_i0, _f2, _o);
        // i + 1
        _bb = _mm256_mul_pd(_b, _f3);
        _bb = _mm256_fmadd_pd(_1mb, _f1, _bb);
        _o = _mm256_fmadd_pd(_ip1, _bb, _o);

        // i + 2
        _ip2 = _mm256_mul_pd(_ip2, _b);
        _o = _mm256_fmadd_pd(_ip2, _f4, _o);
        // i + 3
        _ip3 = _mm256_mul_pd(_ip3, _b);
        _o = _mm256_fmadd_pd(_ip3, _f5, _o);

        _mm256_stream_pd(out+idx, _o);

        idx+=AVX_SIMD_LENGTH;
      } // i
      idx += 6;
    } // j
    idx += 3*sz[0];
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
        _b = _mm256_loadu_pd( beta + idx);
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
        _o = _mm256_fmadd_pd(_im3, _f5, _o);
        // i - 2
        _im2 = _mm256_mul_pd(_im2, _1mb);
        _o = _mm256_fmadd_pd(_im2, _f4, _o);
        // i - 1
        _bb = _mm256_mul_pd(_b, _f1);
        _bb = _mm256_fmadd_pd(_1mb, _f3, _bb);
        _o = _mm256_fnmadd_pd(_im1, _bb, _o);

        // i
        _i0 = _mm256_mul_pd(_i0, _2bm1);
        _o = _mm256_fmadd_pd(_i0, _f2, _o);
        // i + 1
        _bb = _mm256_mul_pd(_b, _f3);
        _bb = _mm256_fmadd_pd(_1mb, _f1, _bb);
        _o = _mm256_fmadd_pd(_ip1, _bb, _o);

        // i + 2
        _ip2 = _mm256_mul_pd(_ip2, _b);
        _o = _mm256_fmadd_pd(_ip2, _f4, _o);
        // i + 3
        _ip3 = _mm256_mul_pd(_ip3, _b);
        _o = _mm256_fmadd_pd(_ip3, _f5, _o);

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
        _b = _mm256_loadu_pd( beta + idx);
        _cmp = _mm256_cmp_pd(_b, zeros, _CMP_GE_OS);

				_1mb  = _mm256_sub_pd( one, _cmp);
        _2bm1 = _mm256_fmsub_pd(two, _cmp, one);

        _im3 = _mm256_loadu_pd( in + idx - n3);
        _im2 = _mm256_loadu_pd( in + idx - n2);
        _im1 = _mm256_loadu_pd( in + idx - n);
        _i0 = _mm256_loadu_pd( in + idx );
        _ip1 = _mm256_loadu_pd( in + idx + n);
        _ip2 = _mm256_loadu_pd( in + idx + n2);
        _ip3 = _mm256_loadu_pd( in + idx + n3);

        // the actual computation
        _o = _mm256_setzero_pd();

        // i - 3
        _im3 = _mm256_mul_pd(_im3, _1mb);
        _o = _mm256_fmadd_pd(_im3, _f5, _o);
        // i - 2
        _im2 = _mm256_mul_pd(_im2, _1mb);
        _o = _mm256_fmadd_pd(_im2, _f4, _o);
        // i - 1
        _bb = _mm256_mul_pd(_b, _f1);
        _bb = _mm256_fmadd_pd(_1mb, _f3, _bb);
        _o = _mm256_fnmadd_pd(_im1, _bb, _o);

        // i
        _i0 = _mm256_mul_pd(_i0, _2bm1);
        _o = _mm256_fmadd_pd(_i0, _f2, _o);
        // i + 1
        _bb = _mm256_mul_pd(_b, _f3);
        _bb = _mm256_fmadd_pd(_1mb, _f1, _bb);
        _o = _mm256_fmadd_pd(_ip1, _bb, _o);

        // i + 2
        _ip2 = _mm256_mul_pd(_ip2, _b);
        _o = _mm256_fmadd_pd(_ip2, _f4, _o);
        // i + 3
        _ip3 = _mm256_mul_pd(_ip3, _b);
        _o = _mm256_fmadd_pd(_ip3, _f5, _o);

        _mm256_stream_pd(out+idx, _o);

        idx+=AVX_SIMD_LENGTH;
      } // i
      idx += 2;
    } // j
    idx += 2*sz[0];
  } // k
}
