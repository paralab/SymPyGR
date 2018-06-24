#include "derivs_cuda.h"

 
__device__ void device_calc_deriv_x(double * output, double * dev_var_in,
    const int offset, double hx, int bflag,
    int nx,int ny,int nz, int sz_x, int sz_y, int sz_z){
    
    int thread_id = blockIdx.x*threads_per_block + threadIdx.x;

    for (int id = thread_id*thread_load_deriv; id<(thread_id+1)*thread_load_deriv; id++){

        int i = id%(sz_x-6) + 3;
        int j = ((id/(sz_x-6))%(sz_y-2)) + 1;
        int k = (id/(sz_y-2)/(sz_x-6)) + 1;

        if(i >= nx-3 || j >= ny-1 || k >= nz-1) return;

        int pp = IDX(i, j, k);

        output[pp] = (dev_var_in[offset + pp - 2] - 8.0*dev_var_in[offset
                                                + pp - 1] + 8.0*dev_var_in[offset + pp + 1]
                                                - dev_var_in[offset + pp + 2] )*((1.0/hx)/12.0);

        if ((bflag & (1u<<OCT_DIR_LEFT)) && i==3)  {
            int pp3 = IDX(3, j, k);
            int pp4 = IDX(4, j, k);
            int pp5 = IDX(5, j, k);
            output[pp3] = ((-3)*dev_var_in[offset + pp3] + 4*dev_var_in[offset
                                                        + pp4] - dev_var_in[offset + pp5]) * 0.5 / hx;
            output[pp4] = (dev_var_in[offset + pp5] - dev_var_in[offset
                                                + pp3]) * (0.50/hx);
        }

        if ((bflag & (1u<<OCT_DIR_RIGHT)) && i==4)  {
            int pp2 = IDX(nx-5, j, k); // IDX(ie-2,j,k)
            int pp3 = IDX(nx-6, j, k); // IDX(ie-3,j,k)
            int pp1 = IDX(nx-4,j,k); // IDX(ie-1,j,k)
            output[pp2] = (dev_var_in[offset + pp1] - dev_var_in[offset + pp3])
            * 0.50 / hx;
            output[pp1] = (dev_var_in[offset + pp3]- 4.0 * dev_var_in[offset + pp2]
            + 3.0 * dev_var_in[offset + pp1]) * 0.50 / hx;

        }
    }

} 
 
 __device__ void device_calc_deriv_y(double * output, double * dev_var_in,
    const int offset, double hx, int bflag,
    int nx,int ny,int nz, int sz_x, int sz_y, int sz_z){
        
    int thread_id = blockIdx.x*threads_per_block + threadIdx.x;

    for (int id = thread_id*thread_load_deriv; id<(thread_id+1)*thread_load_deriv; id++){

        int i = id%(sz_x-6) + 3;
        int j = ((id/(sz_x-6))%(sz_y-6)) + 3;
        int k = (id/(sz_y-6)/(sz_x-6)) + 1;

        if(i >= nx-3 || j >= ny-3 || k >= nz-1) return;

        int pp = IDX(i, j, k);

        output[pp] = (dev_var_in[offset + pp - 2*nx]
        - 8.0*dev_var_in[offset + pp - nx]
        + 8.0*dev_var_in[offset + pp + nx]
        - dev_var_in[offset + pp + 2*nx] )*((1.0/hx)/12.0);


        if ((bflag & (1u<<OCT_DIR_DOWN)) && j==3)  {
            int pp3 = IDX(i, 3, k);
            int pp4 = IDX(i, 4, k);
            int pp5 = IDX(i, 5, k);

            output[pp3] = ((-3)*dev_var_in[(offset) + pp3] +  4*dev_var_in[(offset) + pp4]
            - dev_var_in[(offset) + pp5]) * 0.5 / hx;
            output[pp4] = (dev_var_in[(offset) + pp5] - dev_var_in[(offset) + pp3])
            * (0.50/hx);

        }

        if ((bflag & (1u<<OCT_DIR_UP)) && j==4)  {
            int pp2 = IDX(i, ny-5, k); // IDX(i,je-2,k)
            int pp3 = IDX(i, ny-6, k); // IDX(i,je-3,k)
            int pp1 = IDX(i, ny-4, k); // IDX(i,je-1,k)

            output[pp2] = (dev_var_in[(offset) + pp1] - dev_var_in[(offset) + pp3])
            * 0.50 / hx;
            output[pp1] = (dev_var_in[(offset) + pp3]- 4.0 * dev_var_in[(offset) + pp2]
            + 3.0 * dev_var_in[(offset) + pp1]) * 0.50 / hx;

        }
    }
}

 
 __device__ void device_calc_deriv_z(double * output, double * dev_var_in,
    const int offset, double hx, int bflag,
    int nx,int ny,int nz, int sz_x, int sz_y, int sz_z){
        
    int thread_id = blockIdx.x*threads_per_block + threadIdx.x;

    for (int id = thread_id*thread_load_deriv; id<(thread_id+1)*thread_load_deriv; id++){
    
        int i = id%(sz_x-6) + 3;
        int j = ((id/(sz_x-6))%(sz_y-6)) + 3;
        int k = (id/(sz_y-6)/(sz_x-6)) + 3;
        
        if(i >= nx-3 || j >= ny-3 || k >= nz-3) return;

        int pp = IDX(i, j, k);

        int n = nx * ny;

        output[pp] = (dev_var_in[offset + pp - 2*n] - 8.0*dev_var_in[offset + pp - n]
        + 8.0*dev_var_in[offset + pp + n] - dev_var_in[offset + pp + 2*n])
        * ((1.0/hx)/12);

        if ((bflag & (1u<<OCT_DIR_BACK)) && k==3)  {
            int pp3 = IDX(i, j, 3); // IDX(i, j, 3)
            int pp4 = IDX(i, j, 4); // IDX(i,j,4)
            int pp5 = IDX(i, j, 5); // IDX(i,j,5)

            output[pp3] = ((-3)*dev_var_in[offset + pp3] + 4*dev_var_in[offset + pp4]
            - dev_var_in[offset + pp5]) * 0.5 / hx;
            output[pp4] = (dev_var_in[offset + pp5] - dev_var_in[offset + pp3])
            * (0.50/hx);
        }

        if ((bflag & (1u<<OCT_DIR_FRONT)) && k==4)  {
            int pp2 = IDX(i, j, nz-5); // IDX(i,j,ke-2)
            int pp3 = IDX(i, j, nz-6); // IDX(i,j,ke-3)
            int pp1 = IDX(i, j, nz-4); // IDX(i,j,ke-1)

            output[pp2] = (dev_var_in[offset + pp1] - dev_var_in[offset + pp3])
            * 0.50 / hx;
            output[pp1] = (dev_var_in[offset + pp3]- 4.0 * dev_var_in[offset + pp2]
            + 3.0 * dev_var_in[offset + pp1]) * 0.50 / hx;
        }
    }
}

__device__ void device_calc_deriv_xx(double * output, double * dev_var_in,
    const int offset, double hx, int bflag,
    int nx,int ny,int nz, int sz_x, int sz_y, int sz_z){
        
    int thread_id = blockIdx.x*threads_per_block + threadIdx.x;

    for (int id = thread_id*thread_load_deriv; id<(thread_id+1)*thread_load_deriv; id++){

        int i = id%(sz_x-6) + 3;
        int j = ((id/(sz_x-6))%(sz_y-6)) + 3;
        int k = (id/(sz_y-6)/(sz_x-6)) + 3;

        if(i >= nx-3 || j >= ny-3 || k >= nz-3) return;

        int pp = IDX(i, j, k);

        output[pp] = ((-1)*dev_var_in[offset + pp - 2]
        + 16.0*dev_var_in[offset + pp - 1]
        - 30.0*dev_var_in[offset + pp]
        + 16.0*dev_var_in[offset + pp + 1]
        - dev_var_in[offset + pp + 2]
        )*(1.0/(hx*hx))/12.0;

        if ((bflag & (1u<<OCT_DIR_LEFT)) && i==3)  {
            int pp3 = IDX(3, j, k);
            int pp4 = IDX(4, j, k);
            int pp5 = IDX(5, j, k);
            int pp6 = IDX(6, j, k);

            output[pp3] = (
            2.0 *   dev_var_in[offset + pp3]
            -   5.0 *   dev_var_in[offset + pp4]
            +   4.0 *   dev_var_in[offset + pp5]
            -           dev_var_in[offset + pp6]
            ) * 1.0/(hx*hx);

            output[pp4] = (
            dev_var_in[offset + pp3]
            -   2.0 *   dev_var_in[offset + pp4]
            +           dev_var_in[offset + pp5]
            ) * 1.0/(hx*hx);

        }

        if ((bflag & (1u<<OCT_DIR_RIGHT)) && i==4)  {
            int pp1 = IDX(nx - 4, j, k); // IDX(ie-1,j,k)
            int pp2 = IDX(nx - 5, j, k); // IDX(ie-2,j,k)
            int pp3 = IDX(nx - 6, j, k); // IDX(ie-3,j,k)
            int pp4 = IDX(nx - 7, j, k); // IDX(ie-4,j,k)

            output[pp2] = (
            dev_var_in[offset + pp3]
            -   2.0 *   dev_var_in[offset + pp2]
            +           dev_var_in[offset + pp1]
            ) * 1.0/(hx*hx);


            output[pp1] = (
            -   1.0 *   dev_var_in[offset + pp4]
            +   4.0 *   dev_var_in[offset + pp3]
            -   5.0 *   dev_var_in[offset + pp2]
            +   2.0 *   dev_var_in[offset + pp1]
            ) * 1.0/(hx*hx);
        }
    }
}
 
__device__ void device_calc_deriv_yy(double * output, double * dev_var_in,
    const int offset, double hy, int bflag,
    int nx,int ny,int nz, int sz_x, int sz_y, int sz_z){
        
    int thread_id = blockIdx.x*threads_per_block + threadIdx.x;

    for (int id = thread_id*thread_load_deriv; id<(thread_id+1)*thread_load_deriv; id++){
        
        int i = id%(sz_x-6) + 3;
        int j = ((id/(sz_x-6))%(sz_y-6)) + 3;
        int k = (id/(sz_y-6)/(sz_x-6)) + 3;
        
        if(i >= nx-3 || j >= ny-3 || k >= nz-3) return;

        int pp = IDX(i, j, k);

        output[pp] = ((-1)*dev_var_in[offset + pp - 2*nx]
        + 16.0*dev_var_in[offset + pp - nx]
        - 30.0*dev_var_in[offset + pp]
        + 16.0*dev_var_in[offset + pp + nx]
        - dev_var_in[offset + pp + 2*nx]
        )*(1.0/(hy*hy))/12.0;

        if ((bflag & (1u<<OCT_DIR_DOWN)) && j==3)  {
            int pp3 = IDX(i, 3, k);
            int pp4 = IDX(i, 4, k);
            int pp5 = IDX(i, 5, k);
            int pp6 = IDX(i, 6, k);

            output[pp3] = (
            2.0 *   dev_var_in[offset + pp3]
            -   5.0 *   dev_var_in[offset + pp4]
            +   4.0 *   dev_var_in[offset + pp5]
            -           dev_var_in[offset + pp6]
            ) * 1.0/(hy*hy);

            output[pp4] = (
            dev_var_in[offset + pp3]
            -   2.0 *   dev_var_in[offset + pp4]
            +           dev_var_in[offset + pp5]
            ) * 1.0/(hy*hy);
        }

        if ((bflag & (1u<<OCT_DIR_UP)) && j==4)  {
            int pp1 = IDX(i, ny - 4, k);
            int pp2 = IDX(i, ny - 5, k);
            int pp3 = IDX(i, ny - 6, k);
            int pp4 = IDX(i, ny - 7, k);

            output[pp2] = (
            dev_var_in[offset + pp3]
            -   2.0 *   dev_var_in[offset + pp2]
            +           dev_var_in[offset + pp1]
            ) * 1.0/(hy*hy);


            output[pp1] = (
            -   1.0 *   dev_var_in[offset + pp4]
            +   4.0 *   dev_var_in[offset + pp3]
            -   5.0 *   dev_var_in[offset + pp2]
            +   2.0 *   dev_var_in[offset + pp1]
            ) * 1.0/(hy*hy);

        }
    }
}

__device__ void device_calc_deriv_zz(double * output, double * dev_var_in,
    const int offset, double hz, int bflag,
    int nx,int ny,int nz, int sz_x, int sz_y, int sz_z){
        
    int thread_id = blockIdx.x*threads_per_block + threadIdx.x;

    for (int id = thread_id*thread_load_deriv; id<(thread_id+1)*thread_load_deriv; id++){
        
        int i = id%(sz_x-6) + 3;
        int j = ((id/(sz_x-6))%(sz_y-6)) + 3;
        int k = (id/(sz_y-6)/(sz_x-6)) + 3;
        
            if(i >= nx-3 || j >= ny-3 || k >= nz-3) return;

            int pp = IDX(i, j, k);

            int n = nx * ny;

            output[pp] = ((-1)*dev_var_in[offset + pp - 2*n]
            + 16.0*dev_var_in[offset + pp - n]
            - 30.0*dev_var_in[offset + pp]
            + 16.0*dev_var_in[offset + pp + n]
            - dev_var_in[offset + pp + 2*n]
            )*(1.0/(hz*hz))/12.0;

            if ((bflag & (1u<<OCT_DIR_BACK)) && k==3)  {
                int pp3 = IDX(i, j, 3);
                int pp4 = IDX(i, j, 4);
                int pp5 = IDX(i, j, 5);
                int pp6 = IDX(i, j, 6);

                output[pp3] = (
                2.0 *   dev_var_in[offset + pp3]
                -   5.0 *   dev_var_in[offset + pp4]
                +   4.0 *   dev_var_in[offset + pp5]
                -           dev_var_in[offset + pp6]
                ) * 1.0/(hz*hz);

                output[pp4] = (
                dev_var_in[offset + pp3]
                -   2.0 *   dev_var_in[offset + pp4]
                +           dev_var_in[offset + pp5]
                ) * 1.0/(hz*hz);
            }

            if ((bflag & (1u<<OCT_DIR_FRONT)) && k==4)  {
                int pp1 = IDX(i, j, nz - 4);
                int pp2 = IDX(i, j, nz - 5);
                int pp3 = IDX(i, j, nz - 6);
                int pp4 = IDX(i, j, nz - 7);

                output[pp2] = (
                dev_var_in[offset + pp3]
                -   2.0 *   dev_var_in[offset + pp2]
                +           dev_var_in[offset + pp1]
                ) * 1.0/(hz*hz);


                output[pp1] = (
                -   1.0 *   dev_var_in[offset + pp4]
                +   4.0 *   dev_var_in[offset + pp3]
                -   5.0 *   dev_var_in[offset + pp2]
                +   2.0 *   dev_var_in[offset + pp1]
                ) * 1.0/(hz*hz);
            }
    }
}
__global__ void calc_deriv42_first_part(double * dev_var_in, double hx, double hy, double hz, 
    int sz_x, int sz_y, int sz_z, int bflag,
    #include "list_of_para.h"
    ,
    #include "list_of_offset_para.h"
 ){
 
     int nx = sz_x;
     int ny = sz_y;
     int nz = sz_z;

 
 #include "bssnrhs_cuda_derivs_first_part.h"
 
 }
 
__global__ void calc_deriv42_second_part(double * dev_var_in, double hx, double hy, 
    double hz, int sz_x, int sz_y, int sz_z, int bflag,
    #include "list_of_para.h"
    ,
    #include "list_of_offset_para.h"
){
 
     int nx = sz_x;
     int ny = sz_y;
     int nz = sz_z;
 
 #include "bssnrhs_cuda_derivs_secondd_part.h"
 
 }
 
void cuda_calc_all(double * dev_var_in, double hx, double hy, double hz, int sz_x, 
    int sz_y, int sz_z, int bflag, cudaStream_t stream,
    #include "list_of_para.h"
    ,
    #include "list_of_offset_para.h"
    ){
 
    const int ie = sz_x - 1;//x direction
    const int je = sz_y - 1;//y direction
    const int ke = sz_z - 1;//z direction
 
    int total_points = ceil(1.0*ie*je*ke/thread_load_deriv);
    int blocks = ceil(1.0*total_points/threads_per_block);

    calc_deriv42_first_part <<< blocks, threads_per_block, 0, stream >>> (
                      dev_var_in, hx, hy, hz, sz_x, sz_y, sz_z, bflag,
                    #include "list_of_args.h"
                    ,
                    #include "list_of_offset_args.h"
             );

    calc_deriv42_second_part <<< blocks, threads_per_block, 0, stream >>> (
                      dev_var_in, hx, hy, hz, sz_x, sz_y, sz_z, bflag,
                    #include "list_of_args.h"
                    ,
                    #include "list_of_offset_args.h"
             );
}

