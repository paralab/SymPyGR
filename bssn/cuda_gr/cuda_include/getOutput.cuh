/**
 * Created on: Sep 21, 2018
 * 		Author: Akila, Eranga, Eminda, Ruwan
 **/

#ifndef GET_OUTPUT_CUH
#define GET_OUTPUT_CUH

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream>

#define CHECK_ERROR( err, msg ) if( err != cudaSuccess ) { std::cerr << "\033[1;31mERROR:" << cudaGetErrorName ( err ) << "  |  " << "ERROR DES: " << cudaGetErrorString( err ) << "  |  " << "User msg: " << msg << "\033[0m" << std::endl; exit( 0 ); }

void get_output_kernel_wrapper(double * dev_var_out, const unsigned int * host_sz, cudaStream_t stream,
    #include "para_derivs_offsets.h"
);

#endif
