/**
 * Created on: March 15, 2018
 * 		Author: Akila
 **/

 #include "derivs_cuda.h"
 

 __global__ void calc_deriv42_x(double * output, double * dev_var_in, 
        const int * dev_u_offset, double * dev_dy, int * dev_sz, int* dev_bflag)
 {

    int i = 3 + threadIdx.x + blockIdx.x * blockDim.x;
    int j = 1 + threadIdx.y + blockIdx.y * blockDim.y;
    int k = 1 + threadIdx.z + blockIdx.z * blockDim.z;

    int nx = dev_sz[0]; 
    int ny = dev_sz[1]; 
    if(i >= nx-3 || j >= ny-1 || k >= dev_sz[2]-1) return;
    
    int pp = IDX(i, j, k);
 
    output[pp] = (dev_var_in[(*dev_u_offset) + pp - 2] - 8.0*dev_var_in[(*dev_u_offset)
                     + pp - 1] + 8.0*dev_var_in[(*dev_u_offset) + pp + 1] 
                     - dev_var_in[(*dev_u_offset) + pp + 2] )*((1.0/dev_dy[0])/12.0);

    if ((*dev_bflag & (1u<<OCT_DIR_LEFT)) && i==3)  {
        int pp3 = IDX(3, j, k);
        int pp4 = IDX(4, j, k);
        int pp5 = IDX(5, j, k);
        output[pp3] = ((-3)*dev_var_in[(*dev_u_offset) + pp3] + 4*dev_var_in[(*dev_u_offset) 
                    + pp4] - dev_var_in[(*dev_u_offset) + pp5]) * 0.5 / dev_dy[0];
        output[pp4] = (dev_var_in[(*dev_u_offset) + pp5] - dev_var_in[(*dev_u_offset) 
                    + pp3]) * (0.50/dev_dy[0]);
    }

    if ((*dev_bflag & (1u<<OCT_DIR_RIGHT)) && i==4)  {
        int pp2 = IDX(nx-5, j, k); // IDX(ie-2,j,k)
        int pp3 = IDX(nx-6, j, k); // IDX(ie-3,j,k)
        int pp1 = IDX(nx-4,j,k); // IDX(ie-1,j,k)
        output[pp2] = (dev_var_in[(*dev_u_offset) + pp1] - dev_var_in[(*dev_u_offset) + pp3]) 
                    * 0.50 / dev_dy[0];
        output[pp1] = (dev_var_in[(*dev_u_offset) + pp3]- 4.0 * dev_var_in[(*dev_u_offset) + pp2]
                    + 3.0 * dev_var_in[(*dev_u_offset) + pp1]) * 0.50 / dev_dy[0];

    }
    
 }

 void cuda_deriv42_x(double * output, double * dev_var_in, int * dev_u_offset, double * dev_dy, 
        int * dev_sz, int* dev_bflag, const unsigned int * host_sz)
 {
 
    const int ie = host_sz[0] - 3;//x direction
    const int je = host_sz[1] - 1;//y direction
    const int ke = host_sz[2] - 1;//z direction
  
    int temp_max = (ie>je)? ie : je;
    int maximumIterations = (temp_max>ke) ? temp_max: ke;
     
    int requiredBlocks = (9+maximumIterations) / 10;
    
    calc_deriv42_x <<< dim3(requiredBlocks, requiredBlocks, requiredBlocks),
                      dim3((ie + requiredBlocks -1)/requiredBlocks,
                      (je + requiredBlocks -1)/requiredBlocks, 
                      (ke + requiredBlocks -1)/requiredBlocks) >>> (output, dev_var_in,
                        dev_u_offset, dev_dy, dev_sz, dev_bflag);
 
    CHECK_ERROR(cudaGetLastError(), "calc_deriv42_x Kernel launch failed");
 
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

__global__ void calc_deriv42_y(double* output, double * dev_var_in, const int * dev_u_offset,
                         double * dev_dy, int * dev_sz, int* dev_bflag)
{
    int i = 3 + threadIdx.x + blockIdx.x * blockDim.x;
    int j = 3 + threadIdx.y + blockIdx.y * blockDim.y;
    int k = 1 + threadIdx.z + blockIdx.z * blockDim.z;

    int nx = dev_sz[0]; 
    int ny = dev_sz[1]; 
    if(i >= nx-3 || j >= ny-3 || k >= dev_sz[2]-1) return;
    
    int pp = IDX(i, j, k);

    output[pp] = (dev_var_in[*dev_u_offset + pp - 2*nx] 
                - 8.0*dev_var_in[*dev_u_offset + pp - nx] 
                + 8.0*dev_var_in[*dev_u_offset + pp + nx] 
                - dev_var_in[*dev_u_offset + pp + 2*nx] )*((1.0/dev_dy[0])/12.0);
    
            
    if ((*dev_bflag & (1u<<OCT_DIR_DOWN)) && j==3)  {
        int pp3 = IDX(i, 3, k);
        int pp4 = IDX(i, 4, k);
        int pp5 = IDX(i, 5, k);

        output[pp3] = ((-3)*dev_var_in[(*dev_u_offset) + pp3] +  4*dev_var_in[(*dev_u_offset) + pp4] 
                    - dev_var_in[(*dev_u_offset) + pp5]) * 0.5 / dev_dy[0];
        output[pp4] = (dev_var_in[(*dev_u_offset) + pp5] - dev_var_in[(*dev_u_offset) + pp3]) 
                    * (0.50/dev_dy[0]);
        
    }

    if ((*dev_bflag & (1u<<OCT_DIR_UP)) && j==4)  {
        int pp2 = IDX(i, ny-5, k); // IDX(i,je-2,k)
        int pp3 = IDX(i, ny-6, k); // IDX(i,je-3,k)
        int pp1 = IDX(i, ny-4, k); // IDX(i,je-1,k)
    
        output[pp2] = (dev_var_in[(*dev_u_offset) + pp1] - dev_var_in[(*dev_u_offset) + pp3]) 
                    * 0.50 / dev_dy[0];
        output[pp1] = (dev_var_in[(*dev_u_offset) + pp3]- 4.0 * dev_var_in[(*dev_u_offset) + pp2]
                    + 3.0 * dev_var_in[(*dev_u_offset) + pp1]) * 0.50 / dev_dy[0];
        
    }
}

void cuda_deriv42_y(double * output, double * dev_var_in, int * dev_u_offset, double * dev_dy, 
                int * dev_sz, int* dev_bflag, const unsigned int * host_sz)
 {
    const int ie = host_sz[0] - 3;//x direction
    const int je = host_sz[1] - 3;//y direction
    const int ke = host_sz[2] - 1;//z direction
  
    int temp_max = (ie>je)? ie : je;
    int maximumIterations = (temp_max>ke) ? temp_max: ke;
     
    int requiredBlocks = (9+maximumIterations) / 10;

    calc_deriv42_y<<<dim3(requiredBlocks, requiredBlocks, requiredBlocks),
                    dim3((ie + requiredBlocks -1)/requiredBlocks,
                    (je + requiredBlocks -1)/requiredBlocks, 
                    (ke + requiredBlocks -1)/requiredBlocks)>>> 
                    (output, dev_var_in, dev_u_offset, dev_dy, dev_sz, dev_bflag);

    CHECK_ERROR(cudaGetLastError(), "calc_deriv42_y Kernel launch failed");
    

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

__global__ void calc_deriv42_z(double * output, double * dev_var_in, const int * dev_u_offset,
                             double * dev_dy, int * dev_sz, int* dev_bflag)
{
    int i = 3 + threadIdx.x + blockIdx.x * blockDim.x;
    int j = 3 + threadIdx.y + blockIdx.y * blockDim.y;
    int k = 3 + threadIdx.z + blockIdx.z * blockDim.z;

    int nx = dev_sz[0]; 
    int ny = dev_sz[1]; 
    if(i >= nx-3 || j >= ny-3 || k >= dev_sz[2]-3) return;

    int n = nx * ny;
    int pp = IDX(i, j, k);
   
   output[pp] = (dev_var_in[(*dev_u_offset) + pp - 2*n] - 8.0*dev_var_in[(*dev_u_offset) + pp - n] 
                + 8.0*dev_var_in[(*dev_u_offset) + pp + n] - dev_var_in[(*dev_u_offset) + pp + 2*n]) 
                * ((1.0/dev_dy[0])/12);
    
    if ((*dev_bflag & (1u<<OCT_DIR_BACK)) && k==3)  {
        int pp3 = IDX(i, j, 3); // IDX(i, j, 3)
        int pp4 = IDX(i, j, 4); // IDX(i,j,4)
        int pp5 = IDX(i, j, 5); // IDX(i,j,5)

        output[pp3] = ((-3)*dev_var_in[(*dev_u_offset) + pp3] + 4*dev_var_in[(*dev_u_offset) + pp4] 
                    - dev_var_in[(*dev_u_offset) + pp5]) * 0.5 / dev_dy[0];
        output[pp4] = (dev_var_in[(*dev_u_offset) + pp5] - dev_var_in[(*dev_u_offset) + pp3])
                     * (0.50/dev_dy[0]);
    }
            
    if ((*dev_bflag & (1u<<OCT_DIR_FRONT)) && k==4)  {
        int pp2 = IDX(i, j, dev_sz[2]-5); // IDX(i,j,ke-2)
        int pp3 = IDX(i, j, dev_sz[2]-6); // IDX(i,j,ke-3)
        int pp1 = IDX(i, j, dev_sz[2]-4); // IDX(i,j,ke-1)

        output[pp2] = (dev_var_in[(*dev_u_offset) + pp1] - dev_var_in[(*dev_u_offset) + pp3]) 
                    * 0.50 / dev_dy[0];
        output[pp1] = (dev_var_in[(*dev_u_offset) + pp3]- 4.0 * dev_var_in[(*dev_u_offset) + pp2]
                    + 3.0 * dev_var_in[(*dev_u_offset) + pp1]) * 0.50 / dev_dy[0];
    }
  
}

void cuda_deriv42_z(double * output, double * dev_var_in, int * dev_u_offset, 
    double * dev_dy, int * dev_sz, int* dev_bflag, const unsigned int * host_sz)
{
     const int ie = host_sz[0] - 3;//x direction
    const int je = host_sz[1] - 3;//y direction
    const int ke = host_sz[2] - 3;//z direction
  
    int temp_max = (ie>je)? ie : je;
    int maximumIterations = (temp_max>ke) ? temp_max: ke;
     
    int requiredBlocks = (9+maximumIterations) / 10;

    calc_deriv42_z<<<dim3(requiredBlocks, requiredBlocks, requiredBlocks),
                    dim3((ie + requiredBlocks -1)/requiredBlocks,
                    (je + requiredBlocks -1)/requiredBlocks, 
                    (ke + requiredBlocks -1)/requiredBlocks)>>> 
                    (output, dev_var_in, dev_u_offset, dev_dy, dev_sz, dev_bflag);

    CHECK_ERROR(cudaGetLastError(), "calc_deriv42_z Kernel launch failed");

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

__global__ void calc_deriv42_xx(double * output, double * dev_var_in, const int * dev_u_offset,
             double * dev_dy, int * dev_sz, int* dev_bflag)
 {
    int i = 3 + threadIdx.x + blockIdx.x * blockDim.x;
    int j = 3 + threadIdx.y + blockIdx.y * blockDim.y;
    int k = 3 + threadIdx.z + blockIdx.z * blockDim.z;

    int nx = dev_sz[0]; 
    int ny = dev_sz[1]; 
    if(i >= nx-3 || j >= ny-3 || k >= dev_sz[2]-3) return;

    int pp = IDX(i, j, k);

    output[pp] = ((-1)*dev_var_in[(*dev_u_offset) + pp - 2] 
                + 16.0*dev_var_in[(*dev_u_offset) + pp - 1] 
                - 30.0*dev_var_in[(*dev_u_offset) + pp] 
                + 16.0*dev_var_in[(*dev_u_offset) + pp + 1] 
                - dev_var_in[(*dev_u_offset) + pp + 2] 
            )*(1.0/(dev_dy[0]*dev_dy[0]))/12.0;
            
    if ((*dev_bflag & (1u<<OCT_DIR_LEFT)) && i==3)  {
        int pp3 = IDX(3, j, k); 
        int pp4 = IDX(4, j, k); 
        int pp5 = IDX(5, j, k); 
        int pp6 = IDX(6, j, k); 
     
        output[pp3] = (
                 2.0 *   dev_var_in[(*dev_u_offset) + pp3] 
             -   5.0 *   dev_var_in[(*dev_u_offset) + pp4] 
             +   4.0 *   dev_var_in[(*dev_u_offset) + pp5] 
             -           dev_var_in[(*dev_u_offset) + pp6]
            ) * 1.0/(dev_dy[0]*dev_dy[0]);
     
        output[pp4] = (
                         dev_var_in[(*dev_u_offset) + pp3]
             -   2.0 *   dev_var_in[(*dev_u_offset) + pp4]
             +           dev_var_in[(*dev_u_offset) + pp5]
         ) * 1.0/(dev_dy[0]*dev_dy[0]);
     
    }
                    
    if ((*dev_bflag & (1u<<OCT_DIR_RIGHT)) && i==4)  {
        int pp1 = IDX(dev_sz[0] - 4, j, k); // IDX(ie-1,j,k)
        int pp2 = IDX(dev_sz[0] - 5, j, k); // IDX(ie-2,j,k)
        int pp3 = IDX(dev_sz[0] - 6, j, k); // IDX(ie-3,j,k)
        int pp4 = IDX(dev_sz[0] - 7, j, k); // IDX(ie-4,j,k)

        output[pp2] = (
                            dev_var_in[(*dev_u_offset) + pp3] 
                -   2.0 *   dev_var_in[(*dev_u_offset) + pp2] 
                +           dev_var_in[(*dev_u_offset) + pp1] 
                ) * 1.0/(dev_dy[0]*dev_dy[0]);


            output[pp1] = (
                -   1.0 *   dev_var_in[(*dev_u_offset) + pp4] 
                +   4.0 *   dev_var_in[(*dev_u_offset) + pp3] 
                -   5.0 *   dev_var_in[(*dev_u_offset) + pp2] 
                +   2.0 *   dev_var_in[(*dev_u_offset) + pp1]
                ) * 1.0/(dev_dy[0]*dev_dy[0]);
    }
}

void cuda_deriv42_xx(double * output, double * dev_var_in, int * dev_u_offset, 
                double * dev_dy, int * dev_sz, int* dev_bflag, const unsigned int * host_sz)
{
    const int ie = host_sz[0] - 3;//x direction
    const int je = host_sz[1] - 3;//y direction
    const int ke = host_sz[2] - 3;//z direction
  
    int temp_max = (ie>je)? ie : je;
    int maximumIterations = (temp_max>ke) ? temp_max: ke;
     
    int requiredBlocks = (9+maximumIterations) / 10;

    calc_deriv42_xx<<<dim3(requiredBlocks, requiredBlocks, requiredBlocks),
                    dim3((ie + requiredBlocks -1)/requiredBlocks,
                    (je + requiredBlocks -1)/requiredBlocks, 
                    (ke + requiredBlocks -1)/requiredBlocks)>>> 
                    (output, dev_var_in, dev_u_offset, dev_dy, dev_sz, dev_bflag);

    // Check for any errors launching the kernel
    cudaError_t cudaStatus;
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "calc_deriv42_xx Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        exit(0);
    }
    
    // No GPU code for the following part
    // #ifdef DEBUG_DERIVS_COMP
    // for (int k = kb; k < ke; k++) {
    //   for (int j = jb; j < je; j++) {
    //     for (int i = ib; i < ie; i++) {
    //       int pp = IDX(i,j,k);
    //       if(std::isnan(DxDxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
    //     }
    //   }
    // }
    // #endif    
}



__global__ void calc_deriv42_yy(double * output, double * dev_var_in, const int * dev_u_offset, 
                double * dev_dy, int * dev_sz, int* dev_bflag)
 {
    int i = 3 + threadIdx.x + blockIdx.x * blockDim.x;
    int j = 3 + threadIdx.y + blockIdx.y * blockDim.y;
    int k = 3 + threadIdx.z + blockIdx.z * blockDim.z;

    int nx = dev_sz[0]; 
    int ny = dev_sz[1]; 
    if(i >= nx-3 || j >= ny-3 || k >= dev_sz[2]-3) return;

    int pp = IDX(i, j, k);

    output[pp] = ((-1)*dev_var_in[(*dev_u_offset) + pp - 2*nx] 
                + 16.0*dev_var_in[(*dev_u_offset) + pp - nx] 
                - 30.0*dev_var_in[(*dev_u_offset) + pp] 
                + 16.0*dev_var_in[(*dev_u_offset) + pp + nx] 
                - dev_var_in[(*dev_u_offset) + pp + 2*nx] 
            )*(1.0/(dev_dy[0]*dev_dy[0]))/12.0;
            
    if ((*dev_bflag & (1u<<OCT_DIR_DOWN)) && j==3)  {
        int pp3 = IDX(i, 3, k); 
        int pp4 = IDX(i, 4, k); 
        int pp5 = IDX(i, 5, k); 
        int pp6 = IDX(i, 6, k); 
     
        output[pp3] = (
                 2.0 *   dev_var_in[(*dev_u_offset) + pp3] 
             -   5.0 *   dev_var_in[(*dev_u_offset) + pp4] 
             +   4.0 *   dev_var_in[(*dev_u_offset) + pp5] 
             -           dev_var_in[(*dev_u_offset) + pp6]
            ) * 1.0/(dev_dy[0]*dev_dy[0]);
     
        output[pp4] = (
                         dev_var_in[(*dev_u_offset) + pp3]
             -   2.0 *   dev_var_in[(*dev_u_offset) + pp4]
             +           dev_var_in[(*dev_u_offset) + pp5]
         ) * 1.0/(dev_dy[0]*dev_dy[0]);
    }
                            
    if ((*dev_bflag & (1u<<OCT_DIR_UP)) && j==4)  {
        int pp1 = IDX(i, dev_sz[1] - 4, k); 
        int pp2 = IDX(i, dev_sz[1] - 5, k); 
        int pp3 = IDX(i, dev_sz[1] - 6, k); 
        int pp4 = IDX(i, dev_sz[1] - 7, k); 
     
        output[pp2] = (
                         dev_var_in[(*dev_u_offset) + pp3] 
             -   2.0 *   dev_var_in[(*dev_u_offset) + pp2] 
             +           dev_var_in[(*dev_u_offset) + pp1] 
             ) * 1.0/(dev_dy[0]*dev_dy[0]);
     
     
         output[pp1] = (
             -   1.0 *   dev_var_in[(*dev_u_offset) + pp4] 
             +   4.0 *   dev_var_in[(*dev_u_offset) + pp3] 
             -   5.0 *   dev_var_in[(*dev_u_offset) + pp2] 
             +   2.0 *   dev_var_in[(*dev_u_offset) + pp1]
             ) * 1.0/(dev_dy[0]*dev_dy[0]);
    
    }
    
}

void cuda_deriv42_yy(double * output, double * dev_var_in, int * dev_u_offset, double * dev_dy, 
                int * dev_sz, int* dev_bflag, const unsigned int * host_sz)
{
    const int ie = host_sz[0] - 3;//x direction
    const int je = host_sz[1] - 3;//y direction
    const int ke = host_sz[2] - 3;//z direction
  
    int temp_max = (ie>je)? ie : je;
    int maximumIterations = (temp_max>ke) ? temp_max: ke;
     
    int requiredBlocks = (9+maximumIterations) / 10;

    calc_deriv42_yy<<<dim3(requiredBlocks, requiredBlocks, requiredBlocks),
                    dim3((ie + requiredBlocks -1)/requiredBlocks,
                    (je + requiredBlocks -1)/requiredBlocks, 
                    (ke + requiredBlocks -1)/requiredBlocks)>>> 
                    (output, dev_var_in, dev_u_offset, dev_dy, dev_sz, dev_bflag);

    CHECK_ERROR(cudaGetLastError(), "calc_deriv42_yy Kernel launch failed");
    

    // No GPU code for the following part
    // #ifdef DEBUG_DERIVS_COMP
    // for (int k = kb; k < ke; k++) {
    // for (int j = jb; j < je; j++) {
    //     for (int i = ib; i < ie; i++) {
    //     int pp = IDX(i,j,k);
    //     if(std::isnan(DyDyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
    //     }
    // }
    // }
    // #endif 
}


__global__ void calc_deriv42_zz(double * output, double * dev_var_in, const int * dev_u_offset,
                 double * dev_dy, int * dev_sz, int* dev_bflag)
 {
    int i = 3 + threadIdx.x + blockIdx.x * blockDim.x;
    int j = 3 + threadIdx.y + blockIdx.y * blockDim.y;
    int k = 3 + threadIdx.z + blockIdx.z * blockDim.z;

    int nx = dev_sz[0]; 
    int ny = dev_sz[1]; 
    if(i >= nx-3 || j >= ny-3 || k >= dev_sz[2]-3) return;

    int pp = IDX(i, j, k);
    int n = nx * ny;

    output[pp] = ((-1)*dev_var_in[(*dev_u_offset) + pp - 2*n] 
                + 16.0*dev_var_in[(*dev_u_offset) + pp - n] 
                - 30.0*dev_var_in[(*dev_u_offset) + pp] 
                + 16.0*dev_var_in[(*dev_u_offset) + pp + n] 
                - dev_var_in[(*dev_u_offset) + pp + 2*n] 
            )*(1.0/(dev_dy[0]*dev_dy[0]))/12.0;

    if ((*dev_bflag & (1u<<OCT_DIR_BACK)) && k==3)  {
        int pp3 = IDX(i, j, 3); 
        int pp4 = IDX(i, j, 4); 
        int pp5 = IDX(i, j, 5); 
        int pp6 = IDX(i, j, 6); 
     
        output[pp3] = (
                 2.0 *   dev_var_in[(*dev_u_offset) + pp3] 
             -   5.0 *   dev_var_in[(*dev_u_offset) + pp4] 
             +   4.0 *   dev_var_in[(*dev_u_offset) + pp5] 
             -           dev_var_in[(*dev_u_offset) + pp6]
            ) * 1.0/(dev_dy[0]*dev_dy[0]);
     
        output[pp4] = (
                         dev_var_in[(*dev_u_offset) + pp3]
             -   2.0 *   dev_var_in[(*dev_u_offset) + pp4]
             +           dev_var_in[(*dev_u_offset) + pp5]
         ) * 1.0/(dev_dy[0]*dev_dy[0]);
    }
                                    
    if ((*dev_bflag & (1u<<OCT_DIR_FRONT)) && k==4)  {
        int pp1 = IDX(i, j, dev_sz[2] - 4); 
        int pp2 = IDX(i, j, dev_sz[2] - 5); 
        int pp3 = IDX(i, j, dev_sz[2] - 6); 
        int pp4 = IDX(i, j, dev_sz[2] - 7); 

        output[pp2] = (
                            dev_var_in[(*dev_u_offset) + pp3] 
                -   2.0 *   dev_var_in[(*dev_u_offset) + pp2] 
                +           dev_var_in[(*dev_u_offset) + pp1] 
                ) * 1.0/(dev_dy[0]*dev_dy[0]);


            output[pp1] = (
                -   1.0 *   dev_var_in[(*dev_u_offset) + pp4] 
                +   4.0 *   dev_var_in[(*dev_u_offset) + pp3] 
                -   5.0 *   dev_var_in[(*dev_u_offset) + pp2] 
                +   2.0 *   dev_var_in[(*dev_u_offset) + pp1]
                ) * 1.0/(dev_dy[0]*dev_dy[0]);
    }
}


void cuda_deriv42_zz(double * output, double * dev_var_in, int * dev_u_offset, 
        double * dev_dy, int * dev_sz, int* dev_bflag, const unsigned int * host_sz)
{
    const int ie = host_sz[0] - 3;//x direction
    const int je = host_sz[1] - 3;//y direction
    const int ke = host_sz[2] - 3;//z direction
  
    int temp_max = (ie>je)? ie : je;
    int maximumIterations = (temp_max>ke) ? temp_max: ke;
     
    int requiredBlocks = (9+maximumIterations) / 10;

    calc_deriv42_zz<<<dim3(requiredBlocks, requiredBlocks, requiredBlocks),
                    dim3((ie + requiredBlocks -1)/requiredBlocks,
                    (je + requiredBlocks -1)/requiredBlocks, 
                    (ke + requiredBlocks -1)/requiredBlocks)>>> 
                    (output, dev_var_in, dev_u_offset, dev_dy, dev_sz, dev_bflag);

    CHECK_ERROR(cudaGetLastError(), "calc_deriv42_zz Kernel launch failed");
    

    // No GPU code for the following part
    // #ifdef DEBUG_DERIVS_COMP
    // for (int k = kb; k < ke; k++) {
    //   for (int j = jb; j < je; j++) {
    //     for (int i = ib; i < ie; i++) {
    //       int pp = IDX(i,j,k);
    //       if(std::isnan(DzDzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
    //     }
    //   }
    // }
    // #endif
}


__global__ void calc_deriv42_adv_x(double * output, double * dev_var_in, int * dev_betax,
     double *dev_dx, int* dev_bflag, int* dev_sz, int* dev_u_offset) {
    
    //ib, jb, kb values are accumulated to the x, y, z
    int i = 3 + threadIdx.x + blockIdx.x * blockDim.x;
    int j = 3 + threadIdx.y + blockIdx.y * blockDim.y;
    int k = 3 + threadIdx.z + blockIdx.z * blockDim.z;

    double idx_by_2 = 0.50 * (1.0 / dev_dx[0]);
    double idx_by_12 = (1.0 / dev_dx[0])/12;
    int nx = dev_sz[0];
    int ny = dev_sz[1];

    if(i >= nx-3 || j >= ny-3 || k >= dev_sz[2]-3) return;
    
    int pp = IDX(i, j, k);
    //printf("ie = %d, je = %d, ke = %d\n", i, j, k);
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
    
    if ((*dev_bflag & (1u<<OCT_DIR_LEFT)) && (i == 3)) {
        
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

    if ((*dev_bflag & (1u<<OCT_DIR_RIGHT)) && (i == 4)) {
        
        const int ie = nx - 3;
        
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
    int * dev_betax, int* dev_bflag, const unsigned int * host_sz)
{
    cudaError_t cudaStatus;
    const int ie = host_sz[0] - 3;//x direction
    const int je = host_sz[1] - 3;//y direction
    const int ke = host_sz[2] - 3;//z direction
    //printf("ie = %d, je = %d, ke = %d\n", ie, je, ke);

    int temp_max = (ie>je)? ie : je;
    int maximumIterations = (temp_max>ke) ? temp_max: ke;
    
    int requiredBlocks = (9+maximumIterations) / 10;
  
    calc_deriv42_adv_x <<< dim3(requiredBlocks, requiredBlocks, requiredBlocks),
                        dim3((ie + requiredBlocks -1)/requiredBlocks,
                        (je + requiredBlocks -1)/requiredBlocks, 
                        (ke + requiredBlocks -1)/requiredBlocks) >>> (output, dev_var_in, dev_betax,
                            dev_dx, dev_bflag, dev_sz, dev_u_offset);
    
    CHECK_ERROR(cudaGetLastError(), "calc_deriv42_adv_x Kernel launch failed Kernel launch failed");
}

__global__ void calc_deriv42_adv_y(double * output, double * dev_var_in, int * dev_betay,
    double *dev_dy, int* dev_bflag, int* dev_sz, int* dev_u_offset) {
   
   //ib, jb, kb values are accumulated to the x, y, z
   int i = 3 + threadIdx.x + blockIdx.x * blockDim.x;
   int j = 3 + threadIdx.y + blockIdx.y * blockDim.y;
   int k = 3 + threadIdx.z + blockIdx.z * blockDim.z;

   double idy_by_2 = 0.50 * (1.0 / dev_dy[0]);
   double idy_by_12 = (1.0 / dev_dy[0])/12.0;
   int nx = dev_sz[0];
   int ny = dev_sz[1];
   
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
   
   if ((*dev_bflag & (1u<<OCT_DIR_DOWN)) && (j == 3)) {
       
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

   if ((*dev_bflag & (1u<<OCT_DIR_UP)) && (j == 4)) {
       
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
    int * dev_betay, int* dev_bflag, const unsigned int * host_sz)
{
    cudaError_t cudaStatus;
    const int ie = host_sz[0] - 3;//x direction
    const int je = host_sz[1] - 3;//y direction
    const int ke = host_sz[2] - 3;//z direction

    int temp_max = (ie>je)? ie : je;
    int maximumIterations = (temp_max>ke) ? temp_max: ke;
    
    int requiredBlocks = (9+maximumIterations) / 10;
  
    calc_deriv42_adv_y <<< dim3(requiredBlocks, requiredBlocks, requiredBlocks),
                        dim3((ie + requiredBlocks -1)/requiredBlocks,
                        (je + requiredBlocks -1)/requiredBlocks, 
                        (ke + requiredBlocks -1)/requiredBlocks) >>> (output, dev_var_in, dev_betay,
                            dev_dy, dev_bflag, dev_sz, dev_u_offset);
        
    CHECK_ERROR(cudaGetLastError(), "calc_deriv42_adv_y Kernel launch failed Kernel launch failed");
}

__global__ void calc_deriv42_adv_z(double * output, double * dev_var_in, int * dev_betaz,
    double *dev_dz, int* dev_bflag, int* dev_sz, int* dev_u_offset) {
   
   //ib, jb, kb values are accumulated to the x, y, z
   int i = 3 + threadIdx.x + blockIdx.x * blockDim.x;
   int j = 3 + threadIdx.y + blockIdx.y * blockDim.y;
   int k = 3 + threadIdx.z + blockIdx.z * blockDim.z;

   double idz_by_2 = 0.50 * (1.0 / dev_dz[0]);
   double idz_by_12 = (1.0 / dev_dz[0])/12.0;
   int nx = dev_sz[0];
   int ny = dev_sz[1];
   
   if(i >= nx-3 || j >= ny-3 || k >= dev_sz[2]-3) return;

   int n = nx * ny;
   int pp = IDX(i, j, k);
   
   if (dev_var_in[*dev_betaz + pp] > 0.0 ) {
        output[pp] = ( -  3.0 * dev_var_in[*dev_u_offset + pp - n]
                    - 10.0 * dev_var_in[*dev_u_offset + pp]
                    + 18.0 * dev_var_in[*dev_u_offset + pp + n]
                    -  6.0 * dev_var_in[*dev_u_offset + pp + 2*n]
                    +        dev_var_in[*dev_u_offset + pp + 3*n]
                  ) * idz_by_12;
   }
   else {
       output[pp] = ( -        dev_var_in[*dev_u_offset + pp - 3*n]
                   +  6.0 * dev_var_in[*dev_u_offset + pp - 2*n]
                   - 18.0 * dev_var_in[*dev_u_offset + pp - n]
                   + 10.0 * dev_var_in[*dev_u_offset + pp]
                   +  3.0 * dev_var_in[*dev_u_offset + pp +n]
                 ) * idz_by_12;
               
   }
   
   if ((*dev_bflag & (1u<<OCT_DIR_BACK)) && (k == 3)) {
       
       output[IDX(i,j,3)] = ( -  3.0 * dev_var_in[*dev_u_offset + IDX(i,j,3)]
               +  4.0 * dev_var_in[*dev_u_offset + IDX(i,j,4)]
               -        dev_var_in[*dev_u_offset + IDX(i,j,5)]
               ) * idz_by_2;
               
       if (dev_var_in[*dev_betaz + IDX(i,j,4)] > 0.0) {
           output[IDX(i,j,4)] = ( -  3.0 * dev_var_in[*dev_u_offset + IDX(i,j,4)]
                           +  4.0 * dev_var_in[*dev_u_offset + IDX(i,j,5)]
                           -        dev_var_in[*dev_u_offset + IDX(i,j,6)]
                       ) * idz_by_2;

       }
       else {
           output[IDX(i,j,4)] = ( -         dev_var_in[*dev_u_offset + IDX(i,j,3)]
                           +        dev_var_in[*dev_u_offset + IDX(i,j,5)]
                       ) * idz_by_2;
                       
       }

       if (dev_var_in[*dev_betaz + IDX(i,j,5)] > 0.0 ) {
           output[IDX(i,j,5)] = (-  3.0 * dev_var_in[*dev_u_offset + IDX(i,j,4)]
                       - 10.0 * dev_var_in[*dev_u_offset + IDX(i,j,5)]
                       + 18.0 * dev_var_in[*dev_u_offset + IDX(i,j,6)]
                       -  6.0 * dev_var_in[*dev_u_offset + IDX(i,j,7)]
                       +        dev_var_in[*dev_u_offset + IDX(i,j,8)]
                       ) * idz_by_12;
       }
       else {
           output[IDX(i,j,5)] = (           dev_var_in[*dev_u_offset + IDX(i,j,3)]
                           -  4.0 * dev_var_in[*dev_u_offset + IDX(i,j,4)]
                           +  3.0 * dev_var_in[*dev_u_offset + IDX(i,j,5)]
                       ) * idz_by_2;
       }
   }

   if ((*dev_bflag & (1u<<OCT_DIR_FRONT)) && (k == 4)) {
       
       const int ke = dev_sz[12] - 3;
       
       if ( dev_var_in[*dev_betaz + IDX(i,j,ke-3)] < 0.0 ) {
           output[IDX(i,j,ke-3)] = (  - 3.0 * dev_var_in[*dev_u_offset + IDX(i,j,ke-3)]
                                   + 4.0 * dev_var_in[*dev_u_offset + IDX(i,j,ke-2)]
                                   -       dev_var_in[*dev_u_offset + IDX(i,j,ke-1)]
                                ) * idz_by_2;
       }
       else {
           output[IDX(i,j,ke-3)] = ( -   dev_var_in[*dev_u_offset + IDX(i,j,ke-6)]
                             +  6.0 * dev_var_in[*dev_u_offset + IDX(i,j,ke-5)]
                             - 18.0 * dev_var_in[*dev_u_offset + IDX(i,j,ke-4)]
                             + 10.0 * dev_var_in[*dev_u_offset + IDX(i,j,ke-3)]
                             +  3.0 * dev_var_in[*dev_u_offset + IDX(i,j,ke-2)]
                           ) * idz_by_12;
       }
 
         if (dev_var_in[*dev_betaz + IDX(i,j,ke-2)] > 0.0 ) {
           output[IDX(i,j,ke-2)] = (  -  dev_var_in[*dev_u_offset + IDX(i,j,ke-3)]
                                   +  dev_var_in[*dev_u_offset + IDX(i,j,ke-1)]
                                ) * idz_by_2;
         }
         else {
           output[IDX(i,j,ke-2)] = (     dev_var_in[*dev_u_offset + IDX(i,j,ke-4)]
                              - 4.0 * dev_var_in[*dev_u_offset + IDX(i,j,ke-3)]
                              + 3.0 * dev_var_in[*dev_u_offset + IDX(i,j,ke-2)]
                                ) * idz_by_2;
         }
 
         output[IDX(i,j,ke-1)]  = (          dev_var_in[*dev_u_offset + IDX(i,j,ke-3)]
                                 - 4.0 * dev_var_in[*dev_u_offset + IDX(i,j,ke-2)]
                                 + 3.0 * dev_var_in[*dev_u_offset + IDX(i,j,ke-1)]
                              ) * idz_by_2;
   }
}

void cuda_deriv42_adv_z(double * output, double * dev_var_in, 
    int * dev_u_offset, double * dev_dz, int * dev_sz,
    int * dev_betaz, int* dev_bflag, const unsigned int * host_sz)
{
    cudaError_t cudaStatus;
    const int ie = host_sz[0] - 3;//x direction
    const int je = host_sz[1] - 3;//y direction
    const int ke = host_sz[2] - 3;//z direction

    int temp_max = (ie>je)? ie : je;
    int maximumIterations = (temp_max>ke) ? temp_max: ke;
    
    int requiredBlocks = (9+maximumIterations) / 10;
  
    calc_deriv42_adv_z <<< dim3(requiredBlocks, requiredBlocks, requiredBlocks),
                        dim3((ie + requiredBlocks -1)/requiredBlocks,
                        (je + requiredBlocks -1)/requiredBlocks, 
                        (ke + requiredBlocks -1)/requiredBlocks) >>> (output, dev_var_in, dev_betaz,
                            dev_dz, dev_bflag, dev_sz, dev_u_offset);
    
    CHECK_ERROR(cudaGetLastError(), "calc_deriv42_adv_z Kernel launch failed Kernel launch failed");
}

__global__ void calc_ko_deriv42_x(double * output, double * dev_var_in,
    double *dev_dx, int* dev_bflag, int* dev_sz, int* dev_u_offset) {
   
   //ib, jb, kb values are accumulated to the x, y, z
   int i = 4 + threadIdx.x + blockIdx.x * blockDim.x;
   int j = 3 + threadIdx.y + blockIdx.y * blockDim.y;
   int k = 3 + threadIdx.z + blockIdx.z * blockDim.z;

   int nx = dev_sz[0];
   int ny = dev_sz[1];

   if(i >= nx-4 || j >= ny-3 || k >= dev_sz[2]-3) return;

    if(i==4) {
        int ib=3;
        output[IDX(3, j, k)] = (-1.0 / 64.0 / dev_dx[0]) *
                         (
                         -      dev_var_in[*dev_u_offset + IDX(ib+4,j,k)]
                         +  6.0*dev_var_in[*dev_u_offset + IDX(ib+3,j,k)]
                         - 15.0*dev_var_in[*dev_u_offset + IDX(ib+2,j,k)]
                         + 20.0*dev_var_in[*dev_u_offset + IDX(ib+1,j,k)]
                         - 15.0*dev_var_in[*dev_u_offset + IDX(ib,j,k)]
                         +  6.0*dev_var_in[*dev_u_offset + IDX(ib-1,j,k)]
                         -      dev_var_in[*dev_u_offset + IDX(ib-2,j,k)]
                         );
    }

   int pp = IDX(i, j, k);
   
   output[pp] = (-1.0 / 64.0 / dev_dx[0]) *
                         (
                         -      dev_var_in[*dev_u_offset + pp - 3]
                         +  6.0*dev_var_in[*dev_u_offset + pp - 2]
                         - 15.0*dev_var_in[*dev_u_offset + pp - 1]
                         + 20.0*dev_var_in[*dev_u_offset + pp ]
                         - 15.0*dev_var_in[*dev_u_offset + pp + 1]
                         +  6.0*dev_var_in[*dev_u_offset + pp + 2]
                         -      dev_var_in[*dev_u_offset + pp + 3]
                         );

    if(i==5) {
        int ie = nx-3;
        output[IDX(ie-1, j, k)] = (-1.0 / 64.0 / dev_dx[0]) *
                         (
                         -      dev_var_in[*dev_u_offset + IDX(ie+1,j,k)]
                         +  6.0*dev_var_in[*dev_u_offset + IDX(ie,j,k)]
                         - 15.0*dev_var_in[*dev_u_offset + IDX(ie-1,j,k)]
                         + 20.0*dev_var_in[*dev_u_offset + IDX(ie-2,j,k)]
                         - 15.0*dev_var_in[*dev_u_offset + IDX(ie-3,j,k)]
                         +  6.0*dev_var_in[*dev_u_offset + IDX(ie-4,j,k)]
                         -      dev_var_in[*dev_u_offset + IDX(ie-5,j,k)]
                         );
    }
   
   if ((*dev_bflag & (1u<<OCT_DIR_LEFT)) && (i == 4)) {

    output[IDX(3,j,k)] =  (      dev_var_in[*dev_u_offset + IDX(6,j,k)]
                                - 3.0*dev_var_in[*dev_u_offset + IDX(5,j,k)]
                                + 3.0*dev_var_in[*dev_u_offset + IDX(4,j,k)]
                                -     dev_var_in[*dev_u_offset + IDX(3,j,k)]
                            )/59.0/48.0*64*dev_dx[0];
    output[IDX(4,j,k)] =  (     dev_var_in[*dev_u_offset + IDX(7,j,k)]
                                -  6.0*dev_var_in[*dev_u_offset + IDX(6,j,k)]
                                + 12.0*dev_var_in[*dev_u_offset + IDX(5,j,k)]
                                - 10.0*dev_var_in[*dev_u_offset + IDX(4,j,k)]
                                +  3.0*dev_var_in[*dev_u_offset + IDX(3,j,k)]
                                )/43.0/48.0*64*dev_dx[0];
    output[IDX(5,j,k)] =  (     dev_var_in[*dev_u_offset + IDX(8,j,k)]
                                -  6.0*dev_var_in[*dev_u_offset + IDX(7,j,k)]
                                + 15.0*dev_var_in[*dev_u_offset + IDX(6,j,k)]
                                - 19.0*dev_var_in[*dev_u_offset + IDX(5,j,k)]
                                + 12.0*dev_var_in[*dev_u_offset + IDX(4,j,k)]
                                -  3.0*dev_var_in[*dev_u_offset + IDX(3,j,k)]
                                )/49.0/48.0*64*dev_dx[0];
    }

   if ((*dev_bflag & (1u<<OCT_DIR_RIGHT)) && (i == 5)) {
       
       const int ie = nx - 3;
       output[IDX(ie-3,j,k)] = ( dev_var_in[*dev_u_offset + IDX(ie-6,j,k)]
                                - 6.0*dev_var_in[*dev_u_offset + IDX(ie-5,j,k)]
                                + 15.0*dev_var_in[*dev_u_offset + IDX(ie-4,j,k)]
                                - 19.0*dev_var_in[*dev_u_offset + IDX(ie-3,j,k)]
                                + 12.0*dev_var_in[*dev_u_offset + IDX(ie-2,j,k)]
                                -  3.0*dev_var_in[*dev_u_offset + IDX(ie-1,j,k)]
                                )/49.0/48.0*64*dev_dx[0];
        
        output[IDX(ie-2,j,k)] =  ( dev_var_in[*dev_u_offset + IDX(ie-5,j,k)]
                                -  6.0*dev_var_in[*dev_u_offset + IDX(ie-4,j,k)]
                                + 12.0*dev_var_in[*dev_u_offset + IDX(ie-3,j,k)]
                                - 10.0*dev_var_in[*dev_u_offset + IDX(ie-2,j,k)]
                                +  3.0*dev_var_in[*dev_u_offset + IDX(ie-1,j,k)]
                                )/43.0/48.0*64*dev_dx[0];
       
 
        output[IDX(ie-1,j,k)] = ( dev_var_in[*dev_u_offset + IDX(ie-4,j,k)]
                                -  3.0*dev_var_in[*dev_u_offset + IDX(ie-3,j,k)]
                                +  3.0*dev_var_in[*dev_u_offset + IDX(ie-2,j,k)]
                                -      dev_var_in[*dev_u_offset + IDX(ie-1,j,k)]
                                )/59.0/48.0*64*dev_dx[0];
   }
}

void cuda_ko_deriv42_x(double * output, double * dev_var_in, 
   int * dev_u_offset, double * dev_dx, int * dev_sz,
   int* dev_bflag, const unsigned int * host_sz)
{
    cudaError_t cudaStatus;
    const int ie = host_sz[0] - 4;//x direction
    const int je = host_sz[1] - 3;//y direction
    const int ke = host_sz[2] - 3;//z direction

    int temp_max = (ie>je)? ie : je;
    int maximumIterations = (temp_max>ke) ? temp_max: ke;

    int requiredBlocks = (9+maximumIterations) / 10;
    
    calc_ko_deriv42_x <<< dim3(requiredBlocks, requiredBlocks, requiredBlocks),
                            dim3((ie + requiredBlocks -1)/requiredBlocks,
                            (je + requiredBlocks -1)/requiredBlocks, 
                            (ke + requiredBlocks -1)/requiredBlocks) >>> (output, dev_var_in,
                                dev_dx, dev_bflag, dev_sz, dev_u_offset);
    
    CHECK_ERROR(cudaGetLastError(), "calc_ko_deriv42_x Kernel launch failed Kernel launch failed");
}

__global__ void calc_ko_deriv42_y(double * output, double * dev_var_in,
    double *dev_dy, int* dev_bflag, int* dev_sz, int* dev_u_offset) {
   
   //ib, jb, kb values are accumulated to the x, y, z
   int i = 3 + threadIdx.x + blockIdx.x * blockDim.x;
   int j = 4 + threadIdx.y + blockIdx.y * blockDim.y;
   int k = 3 + threadIdx.z + blockIdx.z * blockDim.z;

   int nx = dev_sz[0];
   int ny = dev_sz[1];

   if(i >= nx-3 || j >= ny-4 || k >= dev_sz[2]-3) return;

   if(j==4) {
    int jb=3;
    output[IDX(i,jb,k)] = (-1.0 / 64.0 / dev_dy[0]) *
                (
                    -      dev_var_in[*dev_u_offset + IDX(i,jb+4,k)]
                    +  6.0*dev_var_in[*dev_u_offset + IDX(i,jb+3,k)]
                    - 15.0*dev_var_in[*dev_u_offset + IDX(i,jb+2,k)]
                    + 20.0*dev_var_in[*dev_u_offset + IDX(i,jb+1,k)]
                    - 15.0*dev_var_in[*dev_u_offset + IDX(i,jb,k)]
                    +  6.0*dev_var_in[*dev_u_offset + IDX(i,jb-1,k)]
                    -      dev_var_in[*dev_u_offset + IDX(i,jb-2,k)]
                    );
    }

   int pp = IDX(i, j, k);
   
   output[pp] = (-1.0 / 64.0 / dev_dy[0]) *
                (
                    -      dev_var_in[*dev_u_offset + pp-3*nx]
                    +  6.0*dev_var_in[*dev_u_offset + pp-2*nx]
                    - 15.0*dev_var_in[*dev_u_offset + pp-nx]
                    + 20.0*dev_var_in[*dev_u_offset + pp]
                    - 15.0*dev_var_in[*dev_u_offset + pp+nx]
                    +  6.0*dev_var_in[*dev_u_offset + pp+2*nx]
                    -      dev_var_in[*dev_u_offset + pp+3*nx]
                    );

    if(j==5) {
        int je = ny - 3;
        output[IDX(i,je-1,k)] = (-1.0 / 64.0 / dev_dy[0]) *
                (
                    -      dev_var_in[*dev_u_offset + IDX(i,je+1,k)]
                    +  6.0*dev_var_in[*dev_u_offset + IDX(i,je,k)]
                    - 15.0*dev_var_in[*dev_u_offset + IDX(i,je-1,k)]
                    + 20.0*dev_var_in[*dev_u_offset + IDX(i,je-2,k)]
                    - 15.0*dev_var_in[*dev_u_offset + IDX(i,je-3,k)]
                    +  6.0*dev_var_in[*dev_u_offset + IDX(i,je-4,k)]
                    -      dev_var_in[*dev_u_offset + IDX(i,je-5,k)]
                    );                   
    }
   if ((*dev_bflag & (1u<<OCT_DIR_DOWN)) && (j == 4)) {

    output[IDX(i,3,k)] =  (      dev_var_in[*dev_u_offset +IDX(i,6,k)]
                                - 3.0*dev_var_in[*dev_u_offset +IDX(i,5,k)]
                                + 3.0*dev_var_in[*dev_u_offset + IDX(i,4,k)]
                                -     dev_var_in[*dev_u_offset + IDX(i,3,k)]
                            )/59.0/48.0*64*dev_dy[0];
    output[IDX(i,4,k)] =  (     dev_var_in[*dev_u_offset + IDX(i,7,k)]
                                -  6.0*dev_var_in[*dev_u_offset + IDX(i,6,k)]
                                + 12.0*dev_var_in[*dev_u_offset + IDX(i,5,k)]
                                - 10.0*dev_var_in[*dev_u_offset + IDX(i,4,k)]
                                +  3.0*dev_var_in[*dev_u_offset + IDX(i,3,k)]
                                )/43.0/48.0*64*dev_dy[0];
    output[IDX(i,5,k)] =  (     dev_var_in[*dev_u_offset + IDX(i,8,k)]
                                -  6.0*dev_var_in[*dev_u_offset + IDX(i,7,k)]
                                + 15.0*dev_var_in[*dev_u_offset + IDX(i,6,k)]
                                - 19.0*dev_var_in[*dev_u_offset + IDX(i,5,k)]
                                + 12.0*dev_var_in[*dev_u_offset + IDX(i,4,k)]
                                -  3.0*dev_var_in[*dev_u_offset + IDX(i,3,k)]
                                )/49.0/48.0*64*dev_dy[0];
    }

   if ((*dev_bflag & (1u<<OCT_DIR_UP)) && (j == 5)) {
       
       const int je = ny - 3;
       output[IDX(i,je-3,k)] = (dev_var_in[*dev_u_offset + IDX(i,je-6,k)]
                                -  6.0*dev_var_in[*dev_u_offset + IDX(i,je-5,k)]
                                + 15.0*dev_var_in[*dev_u_offset + IDX(i,je-4,k)]
                                - 19.0*dev_var_in[*dev_u_offset + IDX(i,je-3,k)]
                                + 12.0*dev_var_in[*dev_u_offset + IDX(i,je-2,k)]
                                -  3.0*dev_var_in[*dev_u_offset + IDX(i,je-1,k)]
                                )/49.0/48.0*64*dev_dy[0];
        
        output[IDX(i,je-2,k)] = (dev_var_in[*dev_u_offset + IDX(i,je-5,k)]
                                -  6.0*dev_var_in[*dev_u_offset + IDX(i,je-4,k)]
                                + 12.0*dev_var_in[*dev_u_offset + IDX(i,je-3,k)]
                                - 10.0*dev_var_in[*dev_u_offset + IDX(i,je-2,k)]
                                +  3.0*dev_var_in[*dev_u_offset + IDX(i,je-1,k)]
                                )/43.0/48.0*64*dev_dy[0];
       
 
        output[IDX(i,je-1,k)] = ( dev_var_in[*dev_u_offset + IDX(i,je-4,k)]
                                -  3.0*dev_var_in[*dev_u_offset + IDX(i,je-3,k)]
                                +  3.0*dev_var_in[*dev_u_offset + IDX(i,je-2,k)]
                                -      dev_var_in[*dev_u_offset + IDX(i,je-1,k)]
                                )/59.0/48.0*64*dev_dy[0];
   }
}

void cuda_ko_deriv42_y(double * output, double * dev_var_in, 
   int * dev_u_offset, double * dev_dy, int * dev_sz,
   int* dev_bflag, const unsigned int * host_sz)
{
    cudaError_t cudaStatus;
    const int ie = host_sz[0] - 3;//x direction
    const int je = host_sz[1] - 4;//y direction
    const int ke = host_sz[2] - 3;//z direction

    int temp_max = (ie>je)? ie : je;
    int maximumIterations = (temp_max>ke) ? temp_max: ke;

    int requiredBlocks = (9+maximumIterations) / 10;
    
    calc_ko_deriv42_y <<< dim3(requiredBlocks, requiredBlocks, requiredBlocks),
                        dim3((ie + requiredBlocks -1)/requiredBlocks,
                        (je + requiredBlocks -1)/requiredBlocks, 
                        (ke + requiredBlocks -1)/requiredBlocks) >>> (output, dev_var_in,
                            dev_dy, dev_bflag, dev_sz, dev_u_offset);

    CHECK_ERROR(cudaGetLastError(), "calc_ko_deriv42_y Kernel launch failed Kernel launch failed");
}

__global__ void calc_ko_deriv42_z(double * output, double * dev_var_in,
    double *dev_dz, int* dev_bflag, int* dev_sz, int* dev_u_offset) {
   
   //ib, jb, kb values are accumulated to the x, y, z
   int i = 3 + threadIdx.x + blockIdx.x * blockDim.x;
   int j = 3 + threadIdx.y + blockIdx.y * blockDim.y;
   int k = 4 + threadIdx.z + blockIdx.z * blockDim.z;

   int nx = dev_sz[0];
   int ny = dev_sz[1];

   if(i >= nx-3 || j >= ny-3 || k >= dev_sz[2]-4) return;

   if(k==4) {
    int kb=3;
    output[IDX(i,j,kb)] = (-1.0 / 64.0 / dev_dz[0]) *
                (
                    -      dev_var_in[*dev_u_offset + IDX(i,j,kb+4)]
                    +  6.0*dev_var_in[*dev_u_offset + IDX(i,j,kb+3)]
                    - 15.0*dev_var_in[*dev_u_offset + IDX(i,j,kb+2)]
                    + 20.0*dev_var_in[*dev_u_offset + IDX(i,j,kb+1)]
                    - 15.0*dev_var_in[*dev_u_offset + IDX(i,j,kb)]
                    +  6.0*dev_var_in[*dev_u_offset + IDX(i,j,kb-1)]
                    -      dev_var_in[*dev_u_offset + IDX(i,j,kb-2)]
                    );
    }

    int pp = IDX(i, j, k);
    int n = nx * ny;
    output[pp] = (-1.0 / 64.0 / dev_dz[0]) *
                (
                    -      dev_var_in[*dev_u_offset + pp-3*n]
                    +  6.0*dev_var_in[*dev_u_offset + pp-2*n]
                    - 15.0*dev_var_in[*dev_u_offset + pp-n]
                    + 20.0*dev_var_in[*dev_u_offset + pp]
                    - 15.0*dev_var_in[*dev_u_offset + pp+n]
                    +  6.0*dev_var_in[*dev_u_offset + pp+2*n]
                    -      dev_var_in[*dev_u_offset + pp+3*n]
                    );

    if(k==5) {
        int ke = dev_sz[2] - 3;
        output[IDX(i,j,ke-1)] = (-1.0 / 64.0 / dev_dz[0]) *
        (
            -      dev_var_in[*dev_u_offset + IDX(i,j,ke+1)]
            +  6.0*dev_var_in[*dev_u_offset + IDX(i,j,ke)]
            - 15.0*dev_var_in[*dev_u_offset + IDX(i,j,ke-1)]
            + 20.0*dev_var_in[*dev_u_offset + IDX(i,j,ke-2)]
            - 15.0*dev_var_in[*dev_u_offset + IDX(i,j,ke-3)]
            +  6.0*dev_var_in[*dev_u_offset + IDX(i,j,ke-4)]
            -      dev_var_in[*dev_u_offset + IDX(i,j,ke-5)]
            );               
    }
   
   

   
   if ((*dev_bflag & (1u<<OCT_DIR_BACK)) && (k == 4)) {

    output[IDX(i,3,k)] =  (      dev_var_in[*dev_u_offset +IDX(i,k,6)]
                                - 3.0*dev_var_in[*dev_u_offset +IDX(i,k,5)]
                                + 3.0*dev_var_in[*dev_u_offset + IDX(i,k,4)]
                                -     dev_var_in[*dev_u_offset + IDX(i,k,3)]
                            )/59.0/48.0*64*dev_dz[0];
    output[IDX(i,j,4)] =  (     dev_var_in[*dev_u_offset + IDX(i,j,7)]
                                -  6.0*dev_var_in[*dev_u_offset + IDX(i,j,6)]
                                + 12.0*dev_var_in[*dev_u_offset + IDX(i,j,5)]
                                - 10.0*dev_var_in[*dev_u_offset + IDX(i,j,4)]
                                +  3.0*dev_var_in[*dev_u_offset + IDX(i,j,3)]
                                )/43.0/48.0*64*dev_dz[0];
    output[IDX(i,j,5)] =  (     dev_var_in[*dev_u_offset + IDX(i,j,8)]
                                -  6.0*dev_var_in[*dev_u_offset + IDX(i,j,7)]
                                + 15.0*dev_var_in[*dev_u_offset + IDX(i,j,6)]
                                - 19.0*dev_var_in[*dev_u_offset + IDX(i,j,5)]
                                + 12.0*dev_var_in[*dev_u_offset + IDX(i,j,4)]
                                -  3.0*dev_var_in[*dev_u_offset + IDX(i,j,3)]
                                )/49.0/48.0*64*dev_dz[0];
    }

   if ((*dev_bflag & (1u<<OCT_DIR_FRONT)) && (k == 5)) {
       
       const int ke = dev_sz[2] - 3;
       output[IDX(i,j,ke-3)] = (    dev_var_in[*dev_u_offset + IDX(i,j,ke-6)]
                                    -  6.0*dev_var_in[*dev_u_offset + IDX(i,j,ke-5)]
                                    + 15.0*dev_var_in[*dev_u_offset + IDX(i,j,ke-4)]
                                    - 19.0*dev_var_in[*dev_u_offset + IDX(i,j,ke-3)]
                                    + 12.0*dev_var_in[*dev_u_offset + IDX(i,j,ke-2)]
                                    -  3.0*dev_var_in[*dev_u_offset + IDX(i,j,ke-1)]
                                    )/49.0/48.0*64*dev_dz[0];
        
        output[IDX(i,j,ke-2)] = (   dev_var_in[*dev_u_offset + IDX(i,j,ke-5)]
                                    -  6.0*dev_var_in[*dev_u_offset + IDX(i,j,ke-4)]
                                    + 12.0*dev_var_in[*dev_u_offset + IDX(i,j,ke-3)]
                                    - 10.0*dev_var_in[*dev_u_offset + IDX(i,j,ke-2)]
                                    +  3.0*dev_var_in[*dev_u_offset + IDX(i,j,ke-1)]
                                    )/43.0/48.0*64*dev_dz[0];
       
 
        output[IDX(i,j,ke-1)] = (   dev_var_in[*dev_u_offset + IDX(i,j,ke-4)]
                                    -  3.0*dev_var_in[*dev_u_offset + IDX(i,j,ke-3)]
                                    +  3.0*dev_var_in[*dev_u_offset + IDX(i,j,ke-2)]
                                    -      dev_var_in[*dev_u_offset + IDX(i,j,ke-1)]
                                    )/59.0/48.0*64*dev_dz[0];
   }
}

void cuda_ko_deriv42_z(double * output, double * dev_var_in, 
   int * dev_u_offset, double * dev_dz, int * dev_sz,
   int* dev_bflag, const unsigned int * host_sz)
{
    cudaError_t cudaStatus;
    const int ie = host_sz[0] - 3;//x direction
    const int je = host_sz[1] - 3;//y direction
    const int ke = host_sz[2] - 4;//z direction

    int temp_max = (ie>je)? ie : je;
    int maximumIterations = (temp_max>ke) ? temp_max: ke;
    
    int requiredBlocks = (9+maximumIterations) / 10;
    
    calc_ko_deriv42_z <<< dim3(requiredBlocks, requiredBlocks, requiredBlocks),
                        dim3((ie + requiredBlocks -1)/requiredBlocks,
                        (je + requiredBlocks -1)/requiredBlocks, 
                        (ke + requiredBlocks -1)/requiredBlocks) >>> (output, dev_var_in,
                            dev_dz, dev_bflag, dev_sz, dev_u_offset);
   
    CHECK_ERROR(cudaGetLastError(), "calc_ko_deriv42_z Kernel launch failed Kernel launch failed");
}