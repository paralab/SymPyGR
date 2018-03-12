/**
 * kernal.h
 * 
 * Created on: Feb 12, 2018
 * 		Author: Eminda, Akila, Eranga, Ruwan
 **/

#ifndef KERNAL_Y_
#define KERNAL_Y_

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "def.h"


#define IDX(i,j,k) ( (i) + nx * ( (j) + ny * (k) ) )

__global__ void firstThreeForLoops(int *c, const int *a, const int *b, double ** dev_var_in, int offset);

void deriv42_yWithCuda(double * dev_var_in, int u_offset, double dy, const unsigned int *sz, unsigned bflag);

#endif /* KERNAL_H_ */
