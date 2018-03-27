/**
 * Created on: March 15, 2018
 * 		Author: Akila
 **/

 #include "derivs_cuda.h"

__global__ void cuda_deriv42_y_firstThreeForLoops(double* output, double * dev_var_in, const int * dev_u_offset, double * dev_dy, int * dev_sz)
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

__global__ void cuda_deriv42_y_secondTwoForLoops(double * output, double * dev_var_in, const int * dev_u_offset, double * dev_dy, int * dev_sz)
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

 __global__ void cuda_deriv42_y_thirdTwoForLoops(double * output, double * dev_var_in, const int * dev_u_offset, double * dev_dy, int * dev_sz)
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
    // printf("%f\n", output[pp1]);
 }
 
 // Please some one verify the below kernals carefully -------------------------------------------------------------------------------
 __global__ void cuda_deriv42_x_firstThreeForLoops(double * output, double * dev_var_in, const int * dev_u_offset, double * dev_dy, int * dev_sz)
{
    int x = threadIdx.x + blockIdx.x*10;
    int y = threadIdx.y + blockIdx.x*10;
    int z = threadIdx.z + blockIdx.x*10;

    int i;
    int j;
    int k;

    if( (dev_sz[0]-3-3)<=x ){ return; } else { i = x+3; } //i handler
    if( (dev_sz[1]-1-1)<=y ){ return; } else { j = y+1; } //j handler
    if( (dev_sz[2]-1-1)<=z ){ return; } else { k = z+1; } //k handler

    int nx = dev_sz[0]; 
    int ny = dev_sz[1]; 
    int pp = IDX(i, j, k);

    output[pp] = (dev_var_in[(*dev_u_offset) + pp - 2] - 8.0*dev_var_in[(*dev_u_offset) + pp - 1] + 8.0*dev_var_in[(*dev_u_offset) + pp + 1] - dev_var_in[(*dev_u_offset) + pp + 2] )*((1.0/dev_dy[0])/12.0);
    // printf("%f\n", output[pp]);
}

__global__ void cuda_deriv42_x_secondTwoForLoops(double * output, double * dev_var_in, const int * dev_u_offset, double * dev_dy, int * dev_sz)
{
   int x = threadIdx.x + blockIdx.x*30;
   int z = threadIdx.y + blockIdx.x*30;

   int j;
   int k;

   if( (dev_sz[1]-1-1)<=x ){ return; } else { j = x+1; } 
   if( (dev_sz[2]-1-1)<=z ){ return; } else { k = z+1; }

   int nx = dev_sz[0];
   int ny = dev_sz[1];

   int pp3 = IDX(3, j, k);
   int pp4 = IDX(4, j, k);
   int pp5 = IDX(5, j, k);

   output[pp3] = ((-3)*dev_var_in[(*dev_u_offset) + pp3] + 4*dev_var_in[(*dev_u_offset) + pp4] - dev_var_in[(*dev_u_offset) + pp5]) * 0.5 / dev_dy[0];
   output[pp4] = (dev_var_in[(*dev_u_offset) + pp5] - dev_var_in[(*dev_u_offset) + pp3]) * (0.50/dev_dy[0]);
   // printf("%f\n", output[pp3]);
}

__global__ void cuda_deriv42_x_thirdTwoForLoops(double * output, double * dev_var_in, const int * dev_u_offset, double * dev_dy, int * dev_sz)
{
   int x = threadIdx.x + blockIdx.x*30;
   int z = threadIdx.y + blockIdx.x*30;

   int j;
   int k;

   if( (dev_sz[1]-1-1)<=x ){ return; } else { j = x+1; } 
   if( (dev_sz[2]-1-1)<=z ){ return; } else { k = z+1; }

   int nx = dev_sz[0];
   int ny = dev_sz[1];

   int pp2 = IDX(dev_sz[0]-5, j, k); // IDX(ie-2,j,k)
   int pp3 = IDX(dev_sz[0]-6, j, k); // IDX(ie-3,j,k)
   int pp1 = IDX(dev_sz[0]-4,j,k); // IDX(ie-1,j,k)

   output[pp2] = (dev_var_in[(*dev_u_offset) + pp1] - dev_var_in[(*dev_u_offset) + pp3]) * 0.50 / dev_dy[0];
   output[pp1] = (dev_var_in[(*dev_u_offset) + pp3]- 4.0 * dev_var_in[(*dev_u_offset) + pp2]+ 3.0 * dev_var_in[(*dev_u_offset) + pp1]) * 0.50 / dev_dy[0];
   // printf("%f\n", output[pp1]);
}

__global__ void cuda_deriv42_z_firstThreeForLoops(double * output, double * dev_var_in, const int * dev_u_offset, double * dev_dy, int * dev_sz)
{
   int x = threadIdx.x + blockIdx.x*10;
   int y = threadIdx.y + blockIdx.x*10;
   int z = threadIdx.z + blockIdx.x*10;

   int i;
   int j;
   int k;

   if( (dev_sz[2]-3-3)<=x ){ return; } else { k = x+3; }
   if( (dev_sz[0]-3-3)<=y ){ return; } else { i = y+3; }
   if( (dev_sz[1]-3-3)<=z ){ return; } else { j = z+3; }

   int nx = dev_sz[0]; 
   int ny = dev_sz[1]; 
   int n = nx * ny;
   int pp = IDX(i, j, k);
   
   output[pp] = (dev_var_in[(*dev_u_offset) + pp - 2*n] - 8.0*dev_var_in[(*dev_u_offset) + pp - n] + 8.0*dev_var_in[(*dev_u_offset) + pp + n] - dev_var_in[(*dev_u_offset) + pp + 2*n]) * ((1.0/dev_dy[0])/12);
   // printf("%f\n", output[pp]);
}

__global__ void cuda_deriv42_z_secondTwoForLoops(double * output, double * dev_var_in, const int * dev_u_offset, double * dev_dy, int * dev_sz)
{
   int x = threadIdx.x + blockIdx.x*30;
   int z = threadIdx.y + blockIdx.x*30;

   int j;
   int i;

   if( (dev_sz[0]-3-3)<=x ){ return; } else { i = x+3; } 
   if( (dev_sz[1]-3-3)<=z ){ return; } else { j = z+3; }

   int nx = dev_sz[0];
   int ny = dev_sz[1];

   int pp3 = IDX(i, j, 3); // IDX(i, j, 3)
   int pp4 = IDX(i, j, 4); // IDX(i,j,4)
   int pp5 = IDX(i, j, 5); // IDX(i,j,5)

   output[pp3] = ((-3)*dev_var_in[(*dev_u_offset) + pp3] + 4*dev_var_in[(*dev_u_offset) + pp4] - dev_var_in[(*dev_u_offset) + pp5]) * 0.5 / dev_dy[0];
   output[pp4] = (dev_var_in[(*dev_u_offset) + pp5] - dev_var_in[(*dev_u_offset) + pp3]) * (0.50/dev_dy[0]);
   // printf("%f\n", output[pp3]);
}

__global__ void cuda_deriv42_z_thirdTwoForLoops(double * output, double * dev_var_in, const int * dev_u_offset, double * dev_dy, int * dev_sz)
{
   int x = threadIdx.x + blockIdx.x*30;
   int z = threadIdx.y + blockIdx.x*30;

   int i;
   int j;

   if( (dev_sz[0]-3-3)<=x ){ return; } else { i = x+3; } 
   if( (dev_sz[1]-3-3)<=z ){ return; } else { j = z+3; }

   int nx = dev_sz[0];
   int ny = dev_sz[1];

   int pp2 = IDX(i, j, dev_sz[2]-5); // IDX(i,j,ke-2)
   int pp3 = IDX(i, j, dev_sz[2]-6); // IDX(i,j,ke-3)
   int pp1 = IDX(i, j, dev_sz[2]-4); // IDX(i,j,ke-1)

   output[pp2] = (dev_var_in[(*dev_u_offset) + pp1] - dev_var_in[(*dev_u_offset) + pp3]) * 0.50 / dev_dy[0];
   output[pp1] = (dev_var_in[(*dev_u_offset) + pp3]- 4.0 * dev_var_in[(*dev_u_offset) + pp2]+ 3.0 * dev_var_in[(*dev_u_offset) + pp1]) * 0.50 / dev_dy[0];
   // printf("%f\n", output[pp1]);
}

void cuda_deriv42_y(double * output, double * dev_var_in, int * dev_u_offset, double * dev_dy, int * dev_sz, unsigned bflag, const unsigned int * host_sz)
 {
    int zblocks = ((host_sz[2]-1)/10)+1;
    int yblocks = ((host_sz[0]-3)/10)+1;
    int xblocks = ((host_sz[1]-3)/10)+1;
    int max1 = ( zblocks < yblocks ) ? yblocks : zblocks;
    int max = ( ( max1 < xblocks ) ? xblocks : max1 );

    cuda_deriv42_y_firstThreeForLoops<<< max, dim3(10, 10, 10) >>>(output, dev_var_in, dev_u_offset, dev_dy, dev_sz);

    // Check for any errors launching the kernel
    cudaError_t cudaStatus;
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cuda_deriv42_y_firstThreeForLoops Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        return;
    }
    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching cuda_deriv42_y_firstThreeForLoops kernal!\n", cudaStatus);
        return;
    }


    if (bflag & (1u<<OCT_DIR_DOWN)) {
        // Not tested yet-----------------------------------------------------------------------------
        int yblocks = ((host_sz[2]-1)/30)+1;
        int xblocks = ((host_sz[0]-3)/30)+1;
        int max = ( xblocks < yblocks ) ? yblocks : xblocks;

        cuda_deriv42_y_secondTwoForLoops<<< max, dim3(30, 30) >>>(output, dev_var_in, dev_u_offset, dev_dy, dev_sz);

        // Check for any errors launching the kernel
        cudaError_t cudaStatus;
        cudaStatus = cudaGetLastError();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cuda_deriv42_y_secondTwoForLoops Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
            return;
        }
        // cudaDeviceSynchronize waits for the kernel to finish, and returns
        // any errors encountered during the launch.
        cudaStatus = cudaDeviceSynchronize();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching cuda_deriv42_y_secondTwoForLoops kernal!\n", cudaStatus);
            return;
        }
    }
    
    if (bflag & (1u<<OCT_DIR_UP)) {
        // Not tested yet-----------------------------------------------------------------------------
        int yblocks = ((host_sz[2]-1)/30)+1;
        int xblocks = ((host_sz[0]-3)/30)+1;
        int max = ( xblocks < yblocks ) ? yblocks : xblocks;

        cuda_deriv42_y_thirdTwoForLoops<<< max, dim3(30, 30) >>>(output, dev_var_in, dev_u_offset, dev_dy, dev_sz);

        // Check for any errors launching the kernel
        cudaError_t cudaStatus;
        cudaStatus = cudaGetLastError();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cuda_deriv42_y_thirdTwoForLoops Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
            return;
        }
        // cudaDeviceSynchronize waits for the kernel to finish, and returns
        // any errors encountered during the launch.
        cudaStatus = cudaDeviceSynchronize();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching cuda_deriv42_y_thirdTwoForLoops kernal!\n", cudaStatus);
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
 
void cuda_deriv42_x(double * output, double * dev_var_in, int * dev_u_offset, double * dev_dy, int * dev_sz, unsigned bflag, const unsigned int * host_sz)
{
    int zblocks = ((host_sz[2]-1)/10)+1; // k
    int yblocks = ((host_sz[1]-1)/10)+1; // j
    int xblocks = ((host_sz[0]-3)/10)+1; // i
    int max1 = ( zblocks < yblocks ) ? yblocks : zblocks;
    int max = ( ( max1 < xblocks ) ? xblocks : max1 );

    cuda_deriv42_x_firstThreeForLoops<<< max, dim3(10, 10, 10) >>>(output, dev_var_in, dev_u_offset, dev_dy, dev_sz);

    // Check for any errors launching the kernel
    cudaError_t cudaStatus;
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cuda_deriv42_x_firstThreeForLoops Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        return;
    }
    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching cuda_deriv42_x_firstThreeForLoops kernal!\n", cudaStatus);
        return;
    }

    if (bflag & (1u<<OCT_DIR_LEFT)) {
        // Not tested yet-----------------------------------------------------------------------------
        int yblocks = ((host_sz[2]-1)/30)+1;
        int xblocks = ((host_sz[1]-1)/30)+1;
        int max = ( xblocks < yblocks ) ? yblocks : xblocks;

        cuda_deriv42_x_secondTwoForLoops<<< max, dim3(30, 30) >>>(output, dev_var_in, dev_u_offset, dev_dy, dev_sz);

        // Check for any errors launching the kernel
        cudaError_t cudaStatus;
        cudaStatus = cudaGetLastError();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cuda_deriv42_x_secondTwoForLoops Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
            return;
        }
        // cudaDeviceSynchronize waits for the kernel to finish, and returns
        // any errors encountered during the launch.
        cudaStatus = cudaDeviceSynchronize();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching cuda_deriv42_x_secondTwoForLoops kernal!\n", cudaStatus);
            return;
        }
    }

    if (bflag & (1u<<OCT_DIR_RIGHT)) {
        // Not tested yet-----------------------------------------------------------------------------
        int yblocks = ((host_sz[2]-1)/30)+1;
        int xblocks = ((host_sz[1]-1)/30)+1;
        int max = ( xblocks < yblocks ) ? yblocks : xblocks;

        cuda_deriv42_x_thirdTwoForLoops<<< max, dim3(30, 30) >>>(output, dev_var_in, dev_u_offset, dev_dy, dev_sz);

        // Check for any errors launching the kernel
        cudaError_t cudaStatus;
        cudaStatus = cudaGetLastError();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cuda_deriv42_x_thirdTwoForLoops Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
            return;
        }
        // cudaDeviceSynchronize waits for the kernel to finish, and returns
        // any errors encountered during the launch.
        cudaStatus = cudaDeviceSynchronize();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching cuda_deriv42_x_thirdTwoForLoops kernal!\n", cudaStatus);
            return;
        }
    } 

    // No GPU code for the following part
    // #ifdef DEBUG_DERIVS_COMP
    //     #pragma message("DEBUG_DERIVS_COMP: ON")
    //     for (int k = 3; k < sz[2]-3; k++) {
    //         for (int j = 3; j < sz[1]-3; j++) {
    //             for (int i = 3; i < sz[0]-3; i++) {
    //                 int pp = IDX(i,j,k);
    //                 if(isnan(Dxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
    //                 }
    //             }
    //         }
    // #endif
}

void cuda_deriv42_z(double * output, double * dev_var_in, int * dev_u_offset, 
    double * dev_dy, int * dev_sz, unsigned bflag, const unsigned int * host_sz)
{
    int zblocks = ((host_sz[2]-1)/10)+1; // k
    int yblocks = ((host_sz[1]-1)/10)+1; // j
    int xblocks = ((host_sz[0]-3)/10)+1; // i
    int max1 = ( zblocks < yblocks ) ? yblocks : zblocks;
    int max = ( ( max1 < xblocks ) ? xblocks : max1 );

    cuda_deriv42_z_firstThreeForLoops<<< max, dim3(10, 10, 10) >>>(output, dev_var_in, dev_u_offset, dev_dy, dev_sz);

    // Check for any errors launching the kernel
    cudaError_t cudaStatus;
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cuda_deriv42_z_firstThreeForLoops Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        return;
    }
    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching cuda_deriv42_z_firstThreeForLoops kernal!\n", cudaStatus);
        return;
    }

    if (bflag & (1u<<OCT_DIR_BACK)) {
        // Not tested yet-----------------------------------------------------------------------------
        int yblocks = ((host_sz[0]-3)/30)+1;
        int xblocks = ((host_sz[0]-3)/30)+1;
        int max = ( xblocks < yblocks ) ? yblocks : xblocks;

        cuda_deriv42_z_secondTwoForLoops<<< max, dim3(30, 30) >>>(output, dev_var_in, dev_u_offset, dev_dy, dev_sz);

        // Check for any errors launching the kernel
        cudaError_t cudaStatus;
        cudaStatus = cudaGetLastError();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cuda_deriv42_z_secondTwoForLoops Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
            return;
        }
        // cudaDeviceSynchronize waits for the kernel to finish, and returns
        // any errors encountered during the launch.
        cudaStatus = cudaDeviceSynchronize();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching cuda_deriv42_z_secondTwoForLoops kernal!\n", cudaStatus);
            return;
        }
    }

    if (bflag & (1u<<OCT_DIR_FRONT)) {
        // Not tested yet-----------------------------------------------------------------------------
        int yblocks = ((host_sz[0]-3)/30)+1;
        int xblocks = ((host_sz[0]-3)/30)+1;
        int max = ( xblocks < yblocks ) ? yblocks : xblocks;

        cuda_deriv42_z_thirdTwoForLoops<<< max, dim3(30, 30) >>>(output, dev_var_in, dev_u_offset, dev_dy, dev_sz);

        // Check for any errors launching the kernel
        cudaError_t cudaStatus;
        cudaStatus = cudaGetLastError();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cuda_deriv42_z_thirdTwoForLoops Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
            return;
        }
        // cudaDeviceSynchronize waits for the kernel to finish, and returns
        // any errors encountered during the launch.
        cudaStatus = cudaDeviceSynchronize();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching cuda_deriv42_z_thirdTwoForLoops kernal!\n", cudaStatus);
            return;
        }
    } 

    //   #ifdef DEBUG_DERIVS_COMP
    //     for (int k = kb; k < ke; k++) {
    //       for (int j = jb; j < je; j++) {
    //         for (int i = ib; i < ie; i++) {
    //           int pp = IDX(i,j,k);
    //           if(std::isnan(Dzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
    //         }
    //       }
    //     }
    //   #endif
}

__global__ void calc_deriv42_adv_x(double * output, double * dev_var_in, int * dev_betax,
     double *dev_dx, int* dev_lflag, int* dev_rflag, int* dev_sz, int* dev_u_offset) {
    
    //ib, jb, kb values are accumulated to the x, y, z
    int i = 3 + threadIdx.x + blockIdx.x * blockDim.x;
    int j = 3 + threadIdx.y + blockIdx.x * blockDim.y;
    int k = 3 + threadIdx.z + blockIdx.x * blockDim.z;

    int idx_by_2 = 0.50 * (1.0 / dev_dx[0]);
    int idx_by_12 = (1.0 / dev_dx[0])/12;
    
    if(i >= dev_sz[0]-3 || j >= dev_sz[1]-3 || k >= dev_sz[2]-3) return;
   
    int pp = IDX(i, j, k);
    //printf("pp = %f\n", dev_var_in[*dev_betax + pp]);
    if (dev_var_in[*dev_betax + pp] > 0.0 ) {
        output[pp] = ( -  3.0 * dev_var_in[*dev_u_offset + pp - 1]
                    - 10.0 * dev_var_in[*dev_u_offset + pp]
                    + 18.0 * dev_var_in[*dev_u_offset + pp + 1]
                    -  6.0 * dev_var_in[*dev_u_offset + pp + 2]
                    +        dev_var_in[*dev_u_offset + pp + 3]
                  ) * idx_by_12;
    }
    else {
        output[pp] = ( -        dev_var_in[*dev_u_offset + pp - 3]
                    +  6.0 * dev_var_in[*dev_u_offset + pp - 2]
                    - 18.0 * dev_var_in[*dev_u_offset + pp - 1]
                    + 10.0 * dev_var_in[*dev_u_offset + pp]
                    +  3.0 * dev_var_in[*dev_u_offset + pp +1]
                  ) * idx_by_12;
    }
    
    if (*dev_lflag && (i == 0)) {
        
        output[IDX(3,j,k)] = ( -  3.0 * dev_var_in[*dev_u_offset + IDX(3,j,k)]
                +  4.0 * dev_var_in[*dev_u_offset + IDX(4,j,k)]
                -        dev_var_in[*dev_u_offset + IDX(5,j,k)]
                ) * idx_by_2;

        if (dev_var_in[*dev_betax + IDX(4,j,k)] > 0.0) {
            output[IDX(4,j,k)] = ( -  3.0 * dev_var_in[*dev_u_offset + IDX(4,j,k)]
                            +  4.0 * dev_var_in[*dev_u_offset + IDX(5,j,k)]
                            -        dev_var_in[*dev_u_offset + IDX(6,j,k)]
                        ) * idx_by_2;
        }
        else {
            output[IDX(4,j,k)] = ( -         dev_var_in[*dev_u_offset + IDX(3,j,k)]
                            +        dev_var_in[*dev_u_offset + IDX(5,j,k)]
                        ) * idx_by_2;
        }

        if (dev_var_in[*dev_betax + IDX(5,j,k)] > 0.0 ) {
            output[IDX(5,j,k)] = (-  3.0 * dev_var_in[*dev_u_offset + IDX(4,j,k)]
                        - 10.0 * dev_var_in[*dev_u_offset + IDX(5,j,k)]
                        + 18.0 * dev_var_in[*dev_u_offset + IDX(6,j,k)]
                        -  6.0 * dev_var_in[*dev_u_offset + IDX(7,j,k)]
                        +        dev_var_in[*dev_u_offset + IDX(8,j,k)]
                        ) * idx_by_12;
        }
        else {
            output[IDX(5,j,k)] = (           dev_var_in[*dev_u_offset + IDX(3,j,k)]
                            -  4.0 * dev_var_in[*dev_u_offset + IDX(4,j,k)]
                            +  3.0 * dev_var_in[*dev_u_offset + IDX(5,j,k)]
                        ) * idx_by_2;
        }
    }

    if (*dev_rflag && (i == 1)) {
        
        const int ie = dev_sz[0] - 3;
        
        if ( dev_var_in[*dev_betax + IDX(ie-3,j,k)] < 0.0 ) {
            output[IDX(ie-3,j,k)] = (  - 3.0 * dev_var_in[*dev_u_offset + IDX(ie-3,j,k)]
                                    + 4.0 * dev_var_in[*dev_u_offset + IDX(ie-2,j,k)]
                                    -       dev_var_in[*dev_u_offset + IDX(ie-1,j,k)]
                                 ) * idx_by_2;
        }
        else {
            output[IDX(ie-3,j,k)] = ( -   dev_var_in[*dev_u_offset + IDX(ie-6,j,k)]
                              +  6.0 * dev_var_in[*dev_u_offset + IDX(ie-5,j,k)]
                              - 18.0 * dev_var_in[*dev_u_offset + IDX(ie-4,j,k)]
                              + 10.0 * dev_var_in[*dev_u_offset + IDX(ie-3  ,j,k)]
                              +  3.0 * dev_var_in[*dev_u_offset + IDX(ie-2,j,k)]
                            ) * idx_by_12;
        }
  
          if (dev_var_in[*dev_betax + IDX(ie-2,j,k)] > 0.0 ) {
            output[IDX(ie-2,j,k)] = (  -  dev_var_in[*dev_u_offset + IDX(ie-3,j,k)]
                                    +  dev_var_in[*dev_u_offset + IDX(ie-1,j,k)]
                                 ) * idx_by_2;
          }
          else {
            output[IDX(ie-2,j,k)] = (     dev_var_in[*dev_u_offset + IDX(ie-4,j,k)]
                               - 4.0 * dev_var_in[*dev_u_offset + IDX(ie-3,j,k)]
                               + 3.0 * dev_var_in[*dev_u_offset + IDX(ie-2,j,k)]
                                 ) * idx_by_2;
          }
  
          output[IDX(ie-1,j,k)] = (          dev_var_in[*dev_u_offset + IDX(ie-3,j,k)]
                                  - 4.0 * dev_var_in[*dev_u_offset + IDX(ie-2,j,k)]
                                  + 3.0 * dev_var_in[*dev_u_offset + IDX(ie-1,j,k)]
                               ) * idx_by_2;
    }
}

void cuda_deriv42_adv_x(double * output, double * dev_var_in, 
    int * dev_u_offset, double * dev_dx, int * dev_sz,
    int * dev_betax, int* dev_lbflag, int* dev_rbflag, const unsigned int * host_sz)
{
    cudaError_t cudaStatus;
    const int ie = host_sz[0] - 3;//x direction
    const int je = host_sz[1] - 3;//y direction
    const int ke = host_sz[2] - 3;//z direction
    //printf("i = %d, j = %d, k = %d\n", ie, je, ke);

    int temp_max = (ie>je)? ie : je;
    int maximumIterations = (temp_max>ke) ? temp_max: ke;
    
    int requiredBlocks = maximumIterations / 10;
    if (ie % 10 != 0 || je % 10 != 0 || ke % 10 != 0) {
        requiredBlocks++;
    }
    
    int threads_x = ie / requiredBlocks;
    int threads_y = je / requiredBlocks;
    int threads_z = ke / requiredBlocks;
   
    calc_deriv42_adv_x <<< requiredBlocks, dim3(threads_x,threads_y,threads_z) >>> (output, dev_var_in, dev_betax,
        dev_dx, dev_lbflag, dev_rbflag, dev_sz, dev_u_offset);
    
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching cuda_deriv42_adv_x kernal!\n", cudaStatus);
            return;
    }
                    
    // cudaMemcpy(Dxu, dev_Dxu, sizeof(double)*sizeof(Dxu), cudaMemcpyDeviceToHost);

}

__global__ void calc_deriv42_adv_y(double * output, double * dev_var_in, int * dev_betay,
    double *dev_dy, int* dev_up_bflag, int* dev_down_bflag, int* dev_sz, int* dev_u_offset) {
   
   //ib, jb, kb values are accumulated to the x, y, z
   int i = 3 + threadIdx.x + blockIdx.x * blockDim.x;
   int j = 3 + threadIdx.y + blockIdx.x * blockDim.y;
   int k = 3 + threadIdx.z + blockIdx.x * blockDim.z;

   int idy_by_2 = 0.50 * (1.0 / dev_dy[0]);
   int idy_by_12 = (1.0 / dev_dy[0])/12.0;
   int nx = dev_sz[0];
   
   if(i >= nx-3 || j >= dev_sz[1]-3 || k >= dev_sz[2]-3) return;
  
   int pp = IDX(i, j, k);
   //printf("pp = %f\n", dev_var_in[*dev_betax + pp]);
   if (dev_var_in[*dev_betay + pp] > 0.0 ) {
        output[pp] = ( -  3.0 * dev_var_in[*dev_u_offset + pp - nx]
                    - 10.0 * dev_var_in[*dev_u_offset + pp]
                    + 18.0 * dev_var_in[*dev_u_offset + pp + nx]
                    -  6.0 * dev_var_in[*dev_u_offset + pp + 2*nx]
                    +        dev_var_in[*dev_u_offset + pp + 3*nx]
                  ) * idy_by_12;
   }
   else {
       output[pp] = ( -        dev_var_in[*dev_u_offset + pp - 3*nx]
                   +  6.0 * dev_var_in[*dev_u_offset + pp - 2*nx]
                   - 18.0 * dev_var_in[*dev_u_offset + pp - nx]
                   + 10.0 * dev_var_in[*dev_u_offset + pp]
                   +  3.0 * dev_var_in[*dev_u_offset + pp +nx]
                 ) * idy_by_12;
               
   }
   
   if (*dev_down_bflag && (j == 0)) {
       
       output[IDX(i,3,k)] = ( -  3.0 * dev_var_in[*dev_u_offset + IDX(i,3,k)]
               +  4.0 * dev_var_in[*dev_u_offset + IDX(i,4,k)]
               -        dev_var_in[*dev_u_offset + IDX(i,5,k)]
               ) * idy_by_2;
               
       if (dev_var_in[*dev_betay + IDX(i,4,k)] > 0.0) {
           output[IDX(i,4,k)] = ( -  3.0 * dev_var_in[*dev_u_offset + IDX(i,4,k)]
                           +  4.0 * dev_var_in[*dev_u_offset + IDX(i,5,k)]
                           -        dev_var_in[*dev_u_offset + IDX(i,6,k)]
                       ) * idy_by_2;

       }
       else {
           output[IDX(i,4,k)] = ( -         dev_var_in[*dev_u_offset + IDX(i,3,k)]
                           +        dev_var_in[*dev_u_offset + IDX(i,5,k)]
                       ) * idy_by_2;
                       
       }

       if (dev_var_in[*dev_betay + IDX(i,5,k)] > 0.0 ) {
           output[IDX(i,5,k)] = (-  3.0 * dev_var_in[*dev_u_offset + IDX(i,4,k)]
                       - 10.0 * dev_var_in[*dev_u_offset + IDX(i,5,k)]
                       + 18.0 * dev_var_in[*dev_u_offset + IDX(i,6,k)]
                       -  6.0 * dev_var_in[*dev_u_offset + IDX(i,7,k)]
                       +        dev_var_in[*dev_u_offset + IDX(i,8,k)]
                       ) * idy_by_12;
       }
       else {
           output[IDX(i,5,k)] = (           dev_var_in[*dev_u_offset + IDX(i,3,k)]
                           -  4.0 * dev_var_in[*dev_u_offset + IDX(i,4,k)]
                           +  3.0 * dev_var_in[*dev_u_offset + IDX(i,5,k)]
                       ) * idy_by_2;
       }
   }

   if (*dev_up_bflag && (j == 1)) {
       
       const int je = dev_sz[1] - 3;
       
       if ( dev_var_in[*dev_betay + IDX(i,je-3,k)] < 0.0 ) {
           output[IDX(i,je-3,k)] = (  - 3.0 * dev_var_in[*dev_u_offset + IDX(i,je-3,k)]
                                   + 4.0 * dev_var_in[*dev_u_offset + IDX(i,je-2,k)]
                                   -       dev_var_in[*dev_u_offset + IDX(i,je-1,k)]
                                ) * idy_by_2;
       }
       else {
           output[IDX(i,je-3,k)] = ( -   dev_var_in[*dev_u_offset + IDX(i,je-6,k)]
                             +  6.0 * dev_var_in[*dev_u_offset + IDX(i,je-5,k)]
                             - 18.0 * dev_var_in[*dev_u_offset + IDX(i,je-4,k)]
                             + 10.0 * dev_var_in[*dev_u_offset + IDX(i,je-3,k)]
                             +  3.0 * dev_var_in[*dev_u_offset + IDX(i,je-2,k)]
                           ) * idy_by_12;
       }
 
         if (dev_var_in[*dev_betay + IDX(i,je-2,k)] > 0.0 ) {
           output[IDX(i,je-2,k)] = (  -  dev_var_in[*dev_u_offset + IDX(i,je-3,k)]
                                   +  dev_var_in[*dev_u_offset + IDX(i,je-1,k)]
                                ) * idy_by_2;
         }
         else {
           output[IDX(i,je-2,k)] = (     dev_var_in[*dev_u_offset + IDX(i,je-4,k)]
                              - 4.0 * dev_var_in[*dev_u_offset + IDX(i,je-3,k)]
                              + 3.0 * dev_var_in[*dev_u_offset + IDX(i,je-2,k)]
                                ) * idy_by_2;
         }
 
         output[IDX(i,je-1,k)]  = (          dev_var_in[*dev_u_offset + IDX(i,je-3,k)]
                                 - 4.0 * dev_var_in[*dev_u_offset + IDX(i,je-2,k)]
                                 + 3.0 * dev_var_in[*dev_u_offset + IDX(i,je-1,k)]
                              ) * idy_by_2;
   }
}

void cuda_deriv42_adv_y(double * output, double * dev_var_in, 
    int * dev_u_offset, double * dev_dy, int * dev_sz,
    int * dev_betay, int* dev_up_bflag, int* dev_down_bflag, const unsigned int * host_sz)
{
    cudaError_t cudaStatus;
    const int ie = host_sz[0] - 3;//x direction
    const int je = host_sz[1] - 3;//y direction
    const int ke = host_sz[2] - 3;//z direction
    //printf("i = %d, j = %d, k = %d\n", ie, je, ke);

    int temp_max = (ie>je)? ie : je;
    int maximumIterations = (temp_max>ke) ? temp_max: ke;
    
    int requiredBlocks = maximumIterations / 10;
    if (ie % 10 != 0 || je % 10 != 0 || ke % 10 != 0) {
        requiredBlocks++;
    }
    
    int threads_x = ie / requiredBlocks;
    int threads_y = je / requiredBlocks;
    int threads_z = ke / requiredBlocks;
    
    calc_deriv42_adv_x <<< requiredBlocks, dim3(threads_x,threads_y,threads_z) >>> (output, dev_var_in, dev_betay,
        dev_dy, dev_up_bflag, dev_down_bflag, dev_sz, dev_u_offset);
    
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching cuda_deriv42_adv_x kernal!\n", cudaStatus);
            return;
    }
                    
    // cudaMemcpy(Dxu, dev_Dxu, sizeof(double)*sizeof(Dxu), cudaMemcpyDeviceToHost);

}