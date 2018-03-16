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

void cuda_deriv42_y(double * output, double * dev_var_in, int * dev_u_offset, double * dev_dy, int * dev_sz, unsigned bflag, const unsigned int * host_sz);

#endif
