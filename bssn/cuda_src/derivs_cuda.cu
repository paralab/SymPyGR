/**
 * Created on: Feb 12, 2018
 * 		Author: Akila
 **/

 #include "derivs_cuda.h"
 
 __global__ void firstThreeForLoops(double * output, double * dev_var_in, const int * dev_u_offset, double * dev_dy, int * dev_sz)
 {
    int x = threadIdx.x + blockIdx.x*10;
    int y = threadIdx.y + blockIdx.x*10;
    int z = threadIdx.z + blockIdx.x*10;

    int i;
    int j;
    int k;

    if( (dev_sz[0]-3-3)<=x ){ return; } else { i = x+3; } //i handler
    if( (dev_sz[1]-3-3)<=y ){ return; } else { j = y+3; } //j handler
    if( (dev_sz[2]-1-1)<=z ){ return; } else { k = z+1; } //k handler

    int nx = dev_sz[0]; 
    int ny = dev_sz[1]; 
    int pp = IDX(i, j, k);

    output[pp] = (dev_var_in[(*dev_u_offset) + pp - 2*dev_sz[0]] - 8.0*dev_var_in[(*dev_u_offset) + pp - dev_sz[0]] + 8.0*dev_var_in[(*dev_u_offset) + pp + dev_sz[0]] - dev_var_in[(*dev_u_offset) + pp + 2*dev_sz[0]] )*((1.0/dev_dy[0])/12.0);
    // printf("%f\n", output[pp]);
}

 __global__ void secondTwoForLoops(double * output, double * dev_var_in, const int * dev_u_offset, double * dev_dy, int * dev_sz)
 {
    int x = threadIdx.x + blockIdx.x*30;
    int z = threadIdx.y + blockIdx.x*30;

    int i;
    int k;

    if( (dev_sz[0]-3-3)<=x ){ return; } else { i = x+3; } //i handler
    if( (dev_sz[2]-1-1)<=z ){ return; } else { k = z+1; } //k handler

    int nx = dev_sz[0];
    int ny = dev_sz[1];

    int pp3 = IDX(i, 3, k);
    int pp4 = IDX(i, 4, k);
    int pp5 = IDX(i, 5, k);

    output[pp3] = ((-3)*dev_var_in[(*dev_u_offset) + pp3] +  4*dev_var_in[(*dev_u_offset) + pp4] - dev_var_in[(*dev_u_offset) + pp5]) * 0.5 / dev_dy[0];
    output[pp4] = (dev_var_in[(*dev_u_offset) + pp5] - dev_var_in[(*dev_u_offset) + pp3]) * (0.50/dev_dy[0]);
    // printf("%f\n", output[pp3]);
 }

 __global__ void thirdTwoForLoops(double * output, double * dev_var_in, const int * dev_u_offset, double * dev_dy, int * dev_sz)
 {
    int x = threadIdx.x + blockIdx.x*30;
    int z = threadIdx.y + blockIdx.x*30;

    int i;
    int k;

    if( (dev_sz[0]-3-3)<=x ){ return; } else { i = x+3; } //i handler
    if( (dev_sz[2]-1-1)<=z ){ return; } else { k = z+1; } //k handler

    int nx = dev_sz[0];
    int ny = dev_sz[1];

    int pp2 = IDX(i, dev_sz[1]-5, k); // IDX(i,je-2,k)
    int pp3 = IDX(i, dev_sz[1]-6, k); // IDX(i,je-3,k)
    int pp1 = IDX(i, dev_sz[1]-4, k); // IDX(i,je-1,k)

    output[pp2] = (dev_var_in[(*dev_u_offset) + pp1] - dev_var_in[(*dev_u_offset) + pp3]) * 0.50 / dev_dy[0];
    output[pp1] = (dev_var_in[(*dev_u_offset) + pp3]- 4.0 * dev_var_in[(*dev_u_offset) + pp2]+ 3.0 * dev_var_in[(*dev_u_offset) + pp1]) * 0.50 / dev_dy[0];
    printf("%f\n", output[pp1]);
 }
 
void cuda_deriv42_y(double * output, double * dev_var_in, int * dev_u_offset, double * dev_dy, int * dev_sz, unsigned bflag, const unsigned int * host_sz)
 {
    int zblocks = ((host_sz[2]-1)/10)+1;
    int yblocks = ((host_sz[0]-3)/10)+1;
    int xblocks = ((host_sz[1]-3)/10)+1;
    int max1 = ( zblocks < yblocks ) ? yblocks : zblocks;
    int max = ( ( max1 < xblocks ) ? xblocks : max1 );

    firstThreeForLoops<<< max, dim3(10, 10, 10) >>>(output, dev_var_in, dev_u_offset, dev_dy, dev_sz);

    // Check for any errors launching the kernel
    cudaError_t cudaStatus;
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "firstThreeForLoops Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        return;
    }
    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching firstThreeForLoops kernal!\n", cudaStatus);
        return;
    }


    if (bflag & (1u<<OCT_DIR_DOWN)) {
        // Not tested yet-----------------------------------------------------------------------------
        int yblocks = ((host_sz[2]-1)/30)+1;
        int xblocks = ((host_sz[0]-3)/30)+1;
        int max = ( xblocks < yblocks ) ? yblocks : xblocks;

        secondTwoForLoops<<< max, dim3(30, 30) >>>(output, dev_var_in, dev_u_offset, dev_dy, dev_sz);

        // Check for any errors launching the kernel
        cudaError_t cudaStatus;
        cudaStatus = cudaGetLastError();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "secondTwoForLoops Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
            return;
        }
        // cudaDeviceSynchronize waits for the kernel to finish, and returns
        // any errors encountered during the launch.
        cudaStatus = cudaDeviceSynchronize();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching secondTwoForLoops kernal!\n", cudaStatus);
            return;
        }
    }
    
    if (bflag & (1u<<OCT_DIR_UP)) {
        // Not tested yet-----------------------------------------------------------------------------
        int yblocks = ((host_sz[2]-1)/30)+1;
        int xblocks = ((host_sz[0]-3)/30)+1;
        int max = ( xblocks < yblocks ) ? yblocks : xblocks;

        thirdTwoForLoops<<< max, dim3(30, 30) >>>(output, dev_var_in, dev_u_offset, dev_dy, dev_sz);

        // Check for any errors launching the kernel
        cudaError_t cudaStatus;
        cudaStatus = cudaGetLastError();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "thirdTwoForLoops Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
            return;
        }
        // cudaDeviceSynchronize waits for the kernel to finish, and returns
        // any errors encountered during the launch.
        cudaStatus = cudaDeviceSynchronize();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching thirdTwoForLoops kernal!\n", cudaStatus);
            return;
        }
    } 

    // No GPU code for the following part
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
 
 