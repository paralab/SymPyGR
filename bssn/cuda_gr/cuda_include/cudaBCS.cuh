#ifndef BCS_CUH
#define BCS_CUH

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream>
#include "def.h"

#define CHECK_ERROR( err, msg ) if( err != cudaSuccess ) { std::cerr << "ERROR:" << cudaGetErrorName ( err ) << "  |  " << "ERROR DES: " << cudaGetErrorString( err ) << "  |  " << "User msg: " << msg << std::endl; exit( 0 ); }

#define IDX(i,j,k) ( (i) + nx * ( (j) + ny * (k) ) )

void bssn_bcs(double * dev_var_out, double * dev_var_in, 
    int u_offset, double * dxf, double * dyf, double * dzf,
    const double * pmin, const double * pmax, const double f_falloff, const double f_asymptotic,
    const unsigned int * host_sz, int bflag, cudaStream_t stream);

#endif