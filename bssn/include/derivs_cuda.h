/**
 * Created on: Feb 12, 2018
 * 		Author: Eminda, Akila, Eranga, Ruwan
 **/

#ifndef DERIVS_CUDA_H_
#define DERIVS_CUDA_H_

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <iostream>
#include "def.h"

#define CHECK_ERROR( err, msg ) if( err != cudaSuccess ) { std::cerr << "ERROR:" << cudaGetErrorName ( err ) << "  |  " << "ERROR DES: " << cudaGetErrorString( err ) << "  |  " << "User msg: " << msg << std::endl; exit( 0 ); }

#define IDX(i,j,k) ( (i) + nx * ( (j) + ny * (k) ) )

// void cuda_deriv42_x(double * output, double * dev_var_in, int u_offset, double dx, int bflag, const unsigned int * host_sz, cudaStream_t stream);
// void cuda_deriv42_y(double * output, double * dev_var_in, int u_offset, double dy, int bflag, const unsigned int * host_sz, cudaStream_t stream);
// void cuda_deriv42_z(double * output, double * dev_var_in, int u_offset, double dz, int bflag, const unsigned int * host_sz, cudaStream_t stream);

// void cuda_deriv42_xx(double * output, double * dev_var_in, int u_offset, double dx, int bflag, const unsigned int * host_sz, cudaStream_t stream);
// void cuda_deriv42_yy(double * output, double * dev_var_in, int u_offset, double dy, int bflag, const unsigned int * host_sz, cudaStream_t stream);
// void cuda_deriv42_zz(double * output, double * dev_var_in, int u_offset, double dz, int bflag, const unsigned int * host_sz, cudaStream_t stream);

// void cuda_deriv42_adv_x(double * output, double * dev_var_in, int u_offset, double dx, int betax, int bflag, const unsigned int * host_sz, cudaStream_t stream);
// void cuda_deriv42_adv_y(double * output, double * dev_var_in, int u_offset, double dy, int betay, int bflag, const unsigned int * host_sz, cudaStream_t stream);
// void cuda_deriv42_adv_z(double * output, double * dev_var_in, int u_offset, double dz, int betaz, int bflag, const unsigned int * host_sz, cudaStream_t stream);

// void cuda_ko_deriv42_x(double * output, double * dev_var_in, int u_offset, double dx, int bflag, const unsigned int * host_sz, cudaStream_t stream);
// void cuda_ko_deriv42_y(double * output, double * dev_var_in, int u_offset, double dy, int bflag, const unsigned int * host_sz, cudaStream_t stream);
// void cuda_ko_deriv42_z(double * output, double * dev_var_in, int u_offset, double dz, int bflag, const unsigned int * host_sz, cudaStream_t stream);
   
void cuda_calc_all(double * dev_var_in, double hx, double hy, double hz, int sz_x, 
int sz_y, int sz_z, int bflag, cudaStream_t stream,
#include "list_of_para.h"
,
#include "list_of_offset_para.h"
);
void cuda_deriv_calc_all_adv(double * dev_var_in, double hx, double hy, double hz, int sz_x, 
int sz_y, int sz_z, int bflag, cudaStream_t stream,
#include "list_of_para.h"
,
#include "list_of_offset_para.h"
);

// void calc_ko_deriv_all(double * dev_var_in,double * dev_dy_hx,double * dev_dy_hy, double * dev_dy_hz, int * dev_sz,
//                        int* dev_bflag, const unsigned int * host_sz,
// #include "list_of_para.h"
// ,
// #include "list_of_offset_para.h"
// );

#endif
