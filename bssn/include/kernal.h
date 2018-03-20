/**
 * kernal.h
 * 
 * Created on: Feb 12, 2018
 * 		Author: Eminda, Akila, Eranga, Ruwan
 **/

#ifndef KERNAL_H_
#define KERNAL_H_

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

__global__ void addKernel(int *c, const int *a, const int *b);
void addWithCuda(int *c, const int *a, const int *b, unsigned int size);

#endif /* KERNAL_H_ */
