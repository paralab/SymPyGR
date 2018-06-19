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
#include "GPUConfig.h"

#define CHECK_ERROR( err, msg ) if( err != cudaSuccess ) { std::cerr << "ERROR:" << cudaGetErrorName ( err ) << "  |  " << "ERROR DES: " << cudaGetErrorString( err ) << "  |  " << "User msg: " << msg << std::endl; exit( 0 ); }

#define IDX(i,j,k) ( (i) + nx * ( (j) + ny * (k) ) )

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

void calc_ko_deriv_all(double * dev_var_in, double hx, double hy, double hz, int sz_x, 
int sz_y, int sz_z, int bflag, cudaStream_t stream,
#include "list_of_para.h"
,
#include "list_of_offset_para.h"
);

#endif
