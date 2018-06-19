/**
 * Created on: Feb 12, 2018
 * 		Author: Eminda, Akila, Eranga, Ruwan
 **/

#ifndef RHS_CUDA_H_
#define RHS_CUDA_H_
#include <cmath>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "derivs_cuda.h"
#include "profile_param.h"
#include "GPUConfig.h"

#define deriv_x cuda_deriv42_x
#define deriv_y cuda_deriv42_y
#define deriv_z cuda_deriv42_z

#define deriv_xx cuda_deriv42_xx
#define deriv_yy cuda_deriv42_yy
#define deriv_zz cuda_deriv42_zz

#define adv_deriv_x cuda_deriv42_adv_x
#define adv_deriv_y cuda_deriv42_adv_y
#define adv_deriv_z cuda_deriv42_adv_z

#define ko_deriv_x cuda_ko_deriv42_x
#define ko_deriv_y cuda_ko_deriv42_y
#define ko_deriv_z cuda_ko_deriv42_z

#define CHECK_ERROR( err, msg ) if( err != cudaSuccess ) { std::cerr << "ERROR:" << cudaGetErrorName ( err ) << "  |  " << "ERROR DES: " << cudaGetErrorString( err ) << "  |  " << "User msg: " << msg << std::endl; exit( 0 ); }

void cuda_bssnrhs(double * dev_var_out, double * dev_var_in, 
const unsigned int unzip_dof, const unsigned int& offset, 
const double * pmin,const double * pmax, const unsigned int *sz, 
const unsigned int& bflag, cudaStream_t stream, cudaStream_t streamAlt,
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
