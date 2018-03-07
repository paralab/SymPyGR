#include <stdio.h>
#include "kernal_deriv42_x.h"
#include <math.h>
#include "def.h"


__inline __device__ int IDX_R(int i,int j,int k,int *nx,int *ny) 
{
    return (i) + (*nx) * ( (j) + (*ny) * (k) );
}

__global__ void calc_deriv42_x(double *dev_Dxu, double *dev_u,
                                double* dev_idx_by_2, double* dev_idx_by_12,
                                 int* dev_ie, int* dev_je, int* dev_ke, int* dev_flag, int* dev_nx, int* dev_ny ) {

    //ib, jb, kb values are accumulated to the x, y, z
    int i = 3 + threadIdx.x + blockIdx.x * blockDim.x;
    int j = 1 + threadIdx.y + blockIdx.x * blockDim.y;
    int k = 1 + threadIdx.z + blockIdx.x * blockDim.z;

    if(i <= *dev_ie && j <= *dev_je && k <= *dev_ke) {
        
        int pp = IDX_R(i, j, k, dev_nx, dev_ny);
        dev_Dxu[pp] = (dev_u[pp - 2] - 8.0 * dev_u[pp - 1] + 
                        8.0 * dev_u[pp + 1] - dev_u[pp + 2]) * *dev_idx_by_12;
        
        if (*dev_flag == 0 && (i == 3 || j == 4)) {
            dev_Dxu[IDX_R(3, j, k, dev_nx, dev_ny)] = (-3.0 * dev_u[IDX_R(3, j, k, dev_nx, dev_ny)]
                                        + 4.0 * dev_u[IDX_R(4, j, k, dev_nx, dev_ny)]
                                        - dev_u[IDX_R(5, j, k, dev_nx, dev_ny)]
                                        ) * *dev_idx_by_2;
            dev_Dxu[IDX_R(4, j, k, dev_nx, dev_ny)] = (-dev_u[IDX_R(3, j, k, dev_nx, dev_ny)]
                                        + dev_u[IDX_R(5, j, k, dev_nx, dev_ny)]
                                        ) * *dev_idx_by_2;
        }

        if (*dev_flag == 1 && (i == *dev_ie - 2 || i == *dev_ie - 1)) {
            dev_Dxu[IDX_R(*dev_ie - 2, j, k, dev_nx, dev_ny)] = (-dev_u[IDX_R(*dev_ie - 3, j, k, dev_nx, dev_ny)]
                                            + dev_u[IDX_R(*dev_ie - 1, j, k, dev_nx, dev_ny)]
                                            ) * *dev_idx_by_2;

            dev_Dxu[IDX_R(*dev_ie - 1, j, k, dev_nx, dev_ny)] = (dev_u[IDX_R(*dev_ie - 3, j, k, dev_nx, dev_ny)]
                                            - 4.0 * dev_u[IDX_R(*dev_ie - 2, j, k, dev_nx, dev_ny)]
                                            + 3.0 * dev_u[IDX_R(*dev_ie - 1, j, k, dev_nx, dev_ny)]
                                            ) * *dev_idx_by_2;
        }
    }
    
}



void kernal_calc_deriv42_x(double *const Dxu, const double *const u,
                           const double dx, const unsigned int *sz, unsigned bflag) {

    const double idx = 1.0 / dx;
    const double idx_by_2 = 0.5 * idx;
    const double idx_by_12 = idx / 12.0;

    const int nx = sz[0];
    const int ny = sz[1];
    const int nz = sz[2];
    // const int ib = 3;
    // const int jb = 1;
    // const int kb = 1;
    const int ie = nx - 3;
    const int je = ny - 1;
    const int ke = nz - 1;
    int flag = 10;
    int *dev_flag, *dev_ie, *dev_je, *dev_ke, *dev_nx, *dev_ny;
    double *dev_idx, *dev_idx_by_2, *dev_idx_by_12;
    double *dev_Dxu, *dev_u;

    cudaMalloc((void**) &dev_idx, sizeof(double));
    cudaMalloc((void**) &dev_idx_by_2, sizeof(double));
    cudaMalloc((void**) &dev_idx_by_12, sizeof(double));
    cudaMalloc((void**) &dev_ie, sizeof(int));
    cudaMalloc((void**) &dev_je, sizeof(int));
    cudaMalloc((void**) &dev_ke, sizeof(int));
    cudaMalloc((void**) &dev_u, sizeof(double)*sizeof(u));
    cudaMalloc((void**) &dev_Dxu, sizeof(double)*sizeof(Dxu));
    cudaMalloc((void**) &dev_nx, sizeof(int));
    cudaMalloc((void**) &dev_ny, sizeof(int));

    if (bflag & (1u<<OCT_DIR_DOWN)) {
        flag = 0;
    }
    
    if (bflag & (1u<<OCT_DIR_UP)) {
        flag = 1;
    }

    cudaMemcpy(dev_idx, &idx, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_idx_by_2, &idx_by_2, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_idx_by_12, &idx_by_12, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_flag, &flag, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_ie, &ie, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_je, &je, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_ke, &ke, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_Dxu, Dxu, sizeof(double) * sizeof(Dxu), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_u, u, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_nx, &nx, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_ny, &ny, sizeof(int), cudaMemcpyHostToDevice);

    int maximumIterations = max(max(ie, je), ke);
    int requiredBlocks = maximumIterations / 10;
    if (ie % 10 == 0 || je % 10 == 0 || ke % 10 == 0) {
        requiredBlocks++;
    }
    
    int threads_x = ie / requiredBlocks;
    int threads_y = je / requiredBlocks;
    int threads_z = ke / requiredBlocks;

    calc_deriv42_x <<< 1000, dim3(threads_x,threads_y,threads_z) >>> (dev_Dxu, dev_u, dev_idx_by_2,
                         dev_idx_by_12, dev_ie, dev_je, dev_ke, dev_flag, dev_nx, dev_ny);
                    
    cudaMemcpy(Dxu, dev_Dxu, sizeof(double)*sizeof(Dxu), cudaMemcpyDeviceToHost);

    
}