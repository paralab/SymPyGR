/**
 * kernal.cu
 * 
 * Created on: Feb 12, 2018
 * 		Author: Eminda, Akila, Eranga, Ruwan
 **/

 #include "derivs_cuda.h"
 #include "rhs.h"
 #include "cuda_runtime.h"
 #include "device_launch_parameters.h"
 #include <stdio.h>
 
 __global__ void firstThreeFor(int *c, const int *a, const int *b)
 {
     int i = threadIdx.x;
     printf ("dff %d",i);
     c[i] = a[i] + b[i];
 }
 
 
 void deriv42_yWithCuda(double * const  Dyu, const double * const  u,
    const double dy, const unsigned int *sz, unsigned bflag)
 {
    const double idy = 1.0/dy;
    const double idy_by_2 = 0.50 * idy;
    const double idy_by_12 = idy / 12.0;
  
    const int nx = sz[0];
    const int ny = sz[1];
    const int nz = sz[2];
    const int ib = 3;
    const int jb = 3;
    const int kb = 1;
    const int ie = sz[0]-3;
    const int je = sz[1]-3;
    const int ke = sz[2]-1;
  
    const int n=nx;
    std::cout << ie <<std::endl;
    for (int k = kb; k < ke; k++) {
        for (int i = ib; i < ie; i++) {
          for (int j = jb; j < je; j++) {
            int pp = IDX(i,j,k); //(i) + nx * ( (j) + ny * (k) )
            Dyu[pp] = (u[pp-2*nx] - 8.0*u[pp-nx] + 8.0*u[pp+nx] - u[pp+2*nx])*idy_by_12;
          }
        }
      }
    
      if (bflag & (1u<<OCT_DIR_DOWN)) {
        for (int k = kb; k < ke; k++) {
          for (int i = ib; i < ie; i++) {
            Dyu[IDX(i, 3,k)] = ( - 3.0 * u[IDX(i,3,k)]
                                +  4.0 * u[IDX(i,4,k)]
                                -        u[IDX(i,5,k)]
                              ) * idy_by_2;
    
            Dyu[IDX(i,4,k)] = ( - u[IDX(i,3,k)]
                                + u[IDX(i,5,k)]
                              ) * idy_by_2;
          }
        }
      }
    
      if (bflag & (1u<<OCT_DIR_UP)) {
        for (int k = kb; k < ke; k++) {
          for (int i = ib; i < ie; i++) {
            Dyu[IDX(i,je-2,k)] = ( - u[IDX(i,je-3,k)]
                                   + u[IDX(i,je-1,k)]
                                 ) * idy_by_2;
    
            Dyu[IDX(i,je-1,k)] = (        u[IDX(i,je-3,k)]
                                  - 4.0 * u[IDX(i,je-2,k)]
                                  + 3.0 * u[IDX(i,je-1,k)]
                              ) * idy_by_2;
          }
        }
      }
    
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
    
    
 
//      int *dev_a = 0;
//      int *dev_b = 0;
//      int *dev_c = 0;
//      cudaError_t cudaStatus;
 
//      // Choose which GPU to run on, change this on a multi-GPU system.
//      cudaStatus = cudaSetDevice(0);
//      if (cudaStatus != cudaSuccess) {
//          fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
//          goto Error;
//      }
 
//      // Allocate GPU buffers for three vectors (two input, one output)    .
//      cudaStatus = cudaMalloc((void**)&dev_c, size * sizeof(int));
//      if (cudaStatus != cudaSuccess) {
//          fprintf(stderr, "cudaMalloc failed!");
//          goto Error;
//      }
 
//      cudaStatus = cudaMalloc((void**)&dev_a, size * sizeof(int));
//      if (cudaStatus != cudaSuccess) {
//          fprintf(stderr, "cudaMalloc failed!");
//          goto Error;
//      }
 
//      cudaStatus = cudaMalloc((void**)&dev_b, size * sizeof(int));
//      if (cudaStatus != cudaSuccess) {
//          fprintf(stderr, "cudaMalloc failed!");
//          goto Error;
//      }
 
//      // Copy input vectors from host memory to GPU buffers.
//      cudaStatus = cudaMemcpy(dev_a, a, size * sizeof(int), cudaMemcpyHostToDevice);
//      if (cudaStatus != cudaSuccess) {
//          fprintf(stderr, "cudaMemcpy failed!");
//          goto Error;
//      }
 
//      cudaStatus = cudaMemcpy(dev_b, b, size * sizeof(int), cudaMemcpyHostToDevice);
//      if (cudaStatus != cudaSuccess) {
//          fprintf(stderr, "cudaMemcpy failed!");
//          goto Error;
//      }
 
//      // Launch a kernel on the GPU with one thread for each element.
//      addKernel<<<1, size>>>(dev_c, dev_a, dev_b);
 
//      // Check for any errors launching the kernel
//      cudaStatus = cudaGetLastError();
//      if (cudaStatus != cudaSuccess) {
//          fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
//          goto Error;
//      }
 
//      // cudaDeviceSynchronize waits for the kernel to finish, and returns
//      // any errors encountered during the launch.
//      cudaStatus = cudaDeviceSynchronize();
//      if (cudaStatus != cudaSuccess) {
//          fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
//          goto Error;
//      }
 
//      // Copy output vector from GPU buffer to host memory.
//      cudaStatus = cudaMemcpy(c, dev_c, size * sizeof(int), cudaMemcpyDeviceToHost);
//      if (cudaStatus != cudaSuccess) {
//          fprintf(stderr, "cudaMemcpy failed!");
//          goto Error;
//      }
 
//  Error:
//      cudaFree(dev_c);
//      cudaFree(dev_a);
//      cudaFree(dev_b);
 
//      //return 0;
 }
 
 