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

#define IDX(i,j,k) ( (i) + nx * ( (j) + ny * (k) ) )
void cuda_calc_all(double * dev_var_in, double * dev_dy_hx, double * dev_dy_hy, double *dev_dy_hz, int *dev_zero, int * dev_sz, int* dev_bflag,
                   const unsigned int * host_sz,
                    #include "list_of_para.h"
);
void cuda_deriv_calc_all_adv(double * dev_var_in, int * dev_sz, int* dev_bflag, const unsigned int * host_sz,
                        double *dev_dy_hx,double * dev_dy_hy,double *dev_dy_hz,
#include "list_of_para.h"
);

void calc_ko_deriv_all(double * dev_var_in,double * dev_dy_hx,double * dev_dy_hy, double * dev_dy_hz, int * dev_sz,
                       int* dev_bflag, const unsigned int * host_sz,
#include "list_of_para.h"
);

   
#endif
