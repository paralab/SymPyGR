#ifndef KERNAL_DERIVE42_X_H
#define KERNAL_DERIVE42_X_H

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

    __inline __device__ int IDX_R(int i,int j,int k,int *nx,int *ny);
    __global__ void calc_deriv42_x(double *dev_Dxu, double *dev_u,
                                int* dev_idx_by_2, int* dev_idx_by_12,
                                 int* dev_ie, int* dev_je, int* dev_ke, int* dev_flag );
    void kernal_calc_deriv42_x(double *const Dxu, const double *const u,
                           const double dx, const unsigned int *sz, unsigned bflag);
#endif
