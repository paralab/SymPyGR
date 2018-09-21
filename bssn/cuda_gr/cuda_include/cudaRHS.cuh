/**
 * Created on: Sep 21, 2018
 * 		Author: Akila, Eranga, Eminda, Ruwan
 **/

#ifndef CUDA_RHS_CUH
#define CUDA_RHS_CUH

#include "cuda_runtime.h"

#include "cudaDerivs.cuh"
#include "bssnEqns.cuh"
#include "getOutput.cuh"
#include "cudaBCS.cuh"

void cuda_bssnrhs(double * dev_var_out, double * dev_var_in, 
const unsigned int unzip_dof, 
const double * pmin, const double * pmax, const unsigned int *sz, 
const unsigned int& bflag, cudaStream_t stream,
#include "list_of_para.h"
);

#endif