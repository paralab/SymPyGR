/**
 * kernal.cu
 * 
 * Created on: Feb 12, 2018
 * 		Author: Akila
 **/

 #include "derivs_cuda.h"
 #include "rhs.h"
 #include "cuda_runtime.h"
 #include "device_launch_parameters.h"
 #include <stdio.h>
 
 __global__ void firstThreeForLoops(double * dev_var_in, const int * dev_u_offset, double * dev_dy, int * dev_sz)
 {
    int x = threadIdx.x + blockIdx.x*10;
    int y = threadIdx.y + blockIdx.x*10;
    int z = threadIdx.z + blockIdx.x*10;

    int i;
    int j;
    int k;

    if( (dev_sz[0]-3-3)<=x ){ return; } else { i = x+3; } //i handler
    if( (dev_sz[1]-3-3)<=y ){ return; } else { j = x+3; } //j handler
    if( (dev_sz[2]-1-1)<=z ){ return; } else { k = x+3; } //k handler

    

    int nx = dev_sz[0];
    int ny = dev_sz[1];
    int pp = IDX(i, j, k);

    // printf("%f\n", (dev_var_in[(*dev_u_offset) + pp - 2*dev_sz[0]]));
    // double f = dev_var_in[(*dev_u_offset) + pp - 2*dev_sz[0]] - 8.0*dev_var_in[(*dev_u_offset) + pp - dev_sz[0]] + 8.0*dev_var_in[(*dev_u_offset) + pp + dev_sz[0]] - dev_var_in[(*dev_u_offset) + pp + 2*dev_sz[0]]; //*((1.0/dev_dy[0])/12.0);
    double f = (dev_var_in[(*dev_u_offset) + pp - 2*dev_sz[0]] - 8.0*dev_var_in[(*dev_u_offset) + pp - dev_sz[0]] + 8.0*dev_var_in[(*dev_u_offset) + pp + dev_sz[0]] - dev_var_in[(*dev_u_offset) + pp + 2*dev_sz[0]] )*((1.0/dev_dy[0])/12.0);
    printf("%f\n", f);

 }
 
 void deriv42_yWithCuda(double * dev_var_in, int u_offset, double dy, const unsigned int *sz, unsigned bflag)
 {

    int * dev_sz;
    double * dev_dy;
    int * dev_u_offset;

    // std::cout << sz[0] << std::endl;

    cudaMalloc((void **) &dev_dy, sizeof(double));
    cudaMalloc((void **) &dev_sz, 3*sizeof(int));
    cudaMalloc((void **) &dev_u_offset, sizeof(int));

    cudaMemcpy(dev_dy, &dy, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_sz, sz, 3*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_u_offset, &u_offset, sizeof(int), cudaMemcpyHostToDevice);

    int zblocks = ((sz[2]-1)/10)+1;
    int yblocks = ((sz[0]-3)/10)+1;
    int xblocks = ((sz[1]-3)/10)+1;
    int max1 = ( zblocks < yblocks ) ? yblocks : zblocks;
    int max = ( ( max1 < xblocks ) ? xblocks : max1 );

    firstThreeForLoops<<< max, dim3(10, 10, 10) >>>(dev_var_in, dev_u_offset, dev_dy, dev_sz);

    // for (int k = kb; k < ke; k++) {
    //     for (int i = ib; i < ie; i++) {
    //       for (int j = jb; j < je; j++) {
    //         int pp = IDX(i,j,k); //(i) + nx * ( (j) + ny * (k) )
    //         Dyu[pp] = (u[pp-2*nx] - 8.0*u[pp-nx] + 8.0*u[pp+nx] - u[pp+2*nx])*idy_by_12;
    //         // printf("%f\n", u[0]);
    //       }
    //     }
    //   }
    
    //   if (bflag & (1u<<OCT_DIR_DOWN)) {
    //     for (int k = kb; k < ke; k++) {
    //       for (int i = ib; i < ie; i++) {
    //         Dyu[IDX(i, 3,k)] = ( - 3.0 * u[IDX(i,3,k)]
    //                             +  4.0 * u[IDX(i,4,k)]
    //                             -        u[IDX(i,5,k)]
    //                           ) * idy_by_2;
    
    //         Dyu[IDX(i,4,k)] = ( - u[IDX(i,3,k)]
    //                             + u[IDX(i,5,k)]
    //                           ) * idy_by_2;
    //       }
    //     }
    //   }
    
    //   if (bflag & (1u<<OCT_DIR_UP)) {
    //     for (int k = kb; k < ke; k++) {
    //       for (int i = ib; i < ie; i++) {
    //         Dyu[IDX(i,je-2,k)] = ( - u[IDX(i,je-3,k)]
    //                                + u[IDX(i,je-1,k)]
    //                              ) * idy_by_2;
    
    //         Dyu[IDX(i,je-1,k)] = (        u[IDX(i,je-3,k)]
    //                               - 4.0 * u[IDX(i,je-2,k)]
    //                               + 3.0 * u[IDX(i,je-1,k)]
    //                           ) * idy_by_2;
    //       }
    //     }
    //   } 
    cudaFree(&dev_sz);
    cudaFree(&dev_dy);
    cudaFree(&dev_u_offset);
    
    // #ifdef DEBUG_DERIVS_COMP
    //   for (int k = 3; k < sz[2]-3; k++) {
    //     for (int j = 3; j < sz[1]-3; j++) {
    //       for (int i = 3; i < sz[0]-3; i++) {
    //         int pp = IDX(i,j,k);
    //         if(std::isnan(Dyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
    //       }
    //     }
    //   }
    // #endif
 }
 
 