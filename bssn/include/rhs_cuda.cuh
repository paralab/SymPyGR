/**
 * Created on: Feb 12, 2018
 * 		Author: Eminda, Akila, Eranga, Ruwan
 **/

#ifndef RHS_CUDA_H_
#define RHS_CUDA_H_
#include <cmath>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "derivs_cuda.cuh"
#include "profile_param.h"

// Please select only one of below
#define parallelized 1
#define parallel_async_hybrid 0
#define pure_async 0
#define pure_async_htod_dtoH_overlap 0

#define CHECK_ERROR( err, msg ) if( err != cudaSuccess ) { std::cerr << "ERROR:" << cudaGetErrorName ( err ) << "  |  " << "ERROR DES: " << cudaGetErrorString( err ) << "  |  " << "User msg: " << msg << std::endl; exit( 0 ); }

void cuda_bssnrhs(double * dev_var_out, double * dev_var_in, 
const unsigned int unzip_dof, 
const double * pmin, const double * pmax, const unsigned int *sz, 
const unsigned int& bflag, cudaStream_t stream,
#include "list_of_para.h"
);

void bssn_bcs(double * dev_var_out, double * dev_var_in, 
    int u_offset, double * dxf, double * dyf, double * dzf,
    const double * pmin, const double * pmax, const double f_falloff, const double f_asymptotic,
    const unsigned int * host_sz, int bflag, cudaStream_t stream);


void get_output (double * dev_var_out, const unsigned int * host_sz, cudaStream_t stream,
        #include "list_of_offset_para.h"
        ,
        #include "list_of_para.h"
    );

#endif
