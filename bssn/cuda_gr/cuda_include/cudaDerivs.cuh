/**
 * Created on: Feb 12, 2018
 * 		Author: Akila, Eranga, Eminda, Ruwan
 **/

#ifndef DERIVS_CUDA_CUH
#define DERIVS_CUDA_CUH

#include <iostream>
#include <stdio.h>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "deviceDerivs.cuh"

#define CHECK_ERROR( err, msg ) if( err != cudaSuccess ) { std::cerr << "\033[1;31mERROR:" << cudaGetErrorName ( err ) << "  |  " << "ERROR DES: " << cudaGetErrorString( err ) << "  |  " << "User msg: " << msg << "\033[0m" << std::endl; exit( 0 ); }

void calc_deriv_kernel_wrapper(double * dev_var_out, double * dev_var_in, double hx, double hy, double hz, const unsigned int * host_sz, int bflag, cudaStream_t stream,
    #include "para_derivs_offsets.h"
);

void calc_ko_deriv_kernel_wrapper(double * dev_var_out, double * dev_var_in, double hx, double hy, double hz, const unsigned int * host_sz, int bflag, cudaStream_t stream,
    #include "para_derivs_offsets.h"
);

#endif
