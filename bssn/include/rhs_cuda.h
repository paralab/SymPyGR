/**
 * Created on: Feb 12, 2018
 * 		Author: Eminda, Akila, Eranga, Ruwan
 **/

#ifndef RHS_CUDA_H_
#define RHS_CUDA_H_

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "profile_param.h"

#include "derivs_cuda.h"

#define deriv_x cuda_deriv42_x
#define deriv_y cuda_deriv42_y
#define deriv_z cuda_deriv42_z

void cuda_bssnrhs(double * dev_var_out, double * dev_var_in, 
const unsigned int unzip_dof, const unsigned int& offset, 
const double *pmin,const double *pmax, const unsigned int *sz, 
const unsigned int& bflag);

#endif
