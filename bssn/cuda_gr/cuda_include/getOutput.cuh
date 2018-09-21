/**
 * Created on: Sep 21, 2018
 * 		Author: Akila, Eranga, Eminda, Ruwan
 **/

#ifndef GET_OUTPUT_CUH
#define GET_OUTPUT_CUH

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream>

#define CHECK_ERROR( err, msg ) if( err != cudaSuccess ) { std::cerr << "ERROR:" << cudaGetErrorName ( err ) << "  |  " << "ERROR DES: " << cudaGetErrorString( err ) << "  |  " << "User msg: " << msg << std::endl; exit( 0 ); }

void get_output (double * dev_var_out, const unsigned int * host_sz, cudaStream_t stream,
    #include "list_of_offset_para.h"
    ,
    #include "list_of_para.h"
);

#endif
