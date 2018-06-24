#include "derivs_cuda.h"

__device__ void device_calc_ko_deriv_x(double * output, double * dev_var_in,
    const int offset, double hx, int bflag,
    int nx,int ny,int nz, int sz_x, int sz_y, int sz_z){

    int thread_id = blockIdx.x*threads_per_block + threadIdx.x;

    for (int id = thread_id*thread_load_ko_deriv; id<(thread_id+1)*thread_load_ko_deriv; id++){
            
        int i = id%(sz_x-8) + 4;
        int j = ((id/(sz_x-8))%(sz_y-6)) + 3;
        int k = (id/(sz_y-6)/(sz_x-8)) + 3;
        
        if(i >= nx-4 || j >= ny-3 || k >= nz-3) return;

        if(i==4) {
            int ib=3;
            output[IDX(3, j, k)] = (-1.0 / 64.0 / hx) *
            (
                    -      dev_var_in[offset + IDX(ib+4,j,k)]
                    +  6.0*dev_var_in[offset + IDX(ib+3,j,k)]
                    - 15.0*dev_var_in[offset + IDX(ib+2,j,k)]
                    + 20.0*dev_var_in[offset + IDX(ib+1,j,k)]
                    - 15.0*dev_var_in[offset + IDX(ib,j,k)]
                    +  6.0*dev_var_in[offset + IDX(ib-1,j,k)]
                    -      dev_var_in[offset + IDX(ib-2,j,k)]
            );
        }

        int pp = IDX(i, j, k);

        output[pp] = (-1.0 / 64.0 / hx) *
        (
        -      dev_var_in[offset + pp - 3]
        +  6.0*dev_var_in[offset + pp - 2]
        - 15.0*dev_var_in[offset + pp - 1]
        + 20.0*dev_var_in[offset + pp ]
        - 15.0*dev_var_in[offset + pp + 1]
        +  6.0*dev_var_in[offset + pp + 2]
        -      dev_var_in[offset + pp + 3]
        );

        if(i==5) {
            int ie = nx-3;
            output[IDX(ie-1, j, k)] = (-1.0 / 64.0 / hx) *
                (
                        -      dev_var_in[offset + IDX(ie+1,j,k)]
                        +  6.0*dev_var_in[offset + IDX(ie,j,k)]
                        - 15.0*dev_var_in[offset + IDX(ie-1,j,k)]
                        + 20.0*dev_var_in[offset + IDX(ie-2,j,k)]
                        - 15.0*dev_var_in[offset + IDX(ie-3,j,k)]
                        +  6.0*dev_var_in[offset + IDX(ie-4,j,k)]
                        -      dev_var_in[offset + IDX(ie-5,j,k)]
                );
        }

            if ((bflag & (1u<<OCT_DIR_LEFT)) && (i == 4)) {

            output[IDX(3,j,k)] =  (      dev_var_in[offset + IDX(6,j,k)]
                - 3.0*dev_var_in[offset + IDX(5,j,k)]
                + 3.0*dev_var_in[offset + IDX(4,j,k)]
                -     dev_var_in[offset + IDX(3,j,k)]
            )/59.0/48.0*64*hx;
            output[IDX(4,j,k)] =  (     dev_var_in[offset + IDX(7,j,k)]
                -  6.0*dev_var_in[offset + IDX(6,j,k)]
                + 12.0*dev_var_in[offset + IDX(5,j,k)]
                - 10.0*dev_var_in[offset + IDX(4,j,k)]
                +  3.0*dev_var_in[offset + IDX(3,j,k)]
            )/43.0/48.0*64*hx;
            output[IDX(5,j,k)] =  (     dev_var_in[offset + IDX(8,j,k)]
                -  6.0*dev_var_in[offset + IDX(7,j,k)]
                + 15.0*dev_var_in[offset + IDX(6,j,k)]
                - 19.0*dev_var_in[offset + IDX(5,j,k)]
                + 12.0*dev_var_in[offset + IDX(4,j,k)]
                -  3.0*dev_var_in[offset + IDX(3,j,k)]
            )/49.0/48.0*64*hx;
        }

            if ((bflag & (1u<<OCT_DIR_RIGHT)) && (i == 5)) {

            const int ie = nx - 3;
            output[IDX(ie-3,j,k)] = ( dev_var_in[offset + IDX(ie-6,j,k)]
                - 6.0*dev_var_in[offset + IDX(ie-5,j,k)]
                + 15.0*dev_var_in[offset + IDX(ie-4,j,k)]
                - 19.0*dev_var_in[offset + IDX(ie-3,j,k)]
                + 12.0*dev_var_in[offset + IDX(ie-2,j,k)]
                -  3.0*dev_var_in[offset + IDX(ie-1,j,k)]
            )/49.0/48.0*64*hx;

            output[IDX(ie-2,j,k)] =  ( dev_var_in[offset + IDX(ie-5,j,k)]
                -  6.0*dev_var_in[offset + IDX(ie-4,j,k)]
                + 12.0*dev_var_in[offset + IDX(ie-3,j,k)]
                - 10.0*dev_var_in[offset + IDX(ie-2,j,k)]
                +  3.0*dev_var_in[offset + IDX(ie-1,j,k)]
            )/43.0/48.0*64*hx;


            output[IDX(ie-1,j,k)] = ( dev_var_in[offset + IDX(ie-4,j,k)]
                -  3.0*dev_var_in[offset + IDX(ie-3,j,k)]
                +  3.0*dev_var_in[offset + IDX(ie-2,j,k)]
                -      dev_var_in[offset + IDX(ie-1,j,k)]
            )/59.0/48.0*64*hx;
        }
    }
} 

 
__device__ void device_calc_ko_deriv_y(double * output, double * dev_var_in,
    const int offset, double hy, int bflag,
    int nx,int ny,int nz, int sz_x, int sz_y, int sz_z){

    int thread_id = blockIdx.x*threads_per_block + threadIdx.x;

    for (int id = thread_id*thread_load_ko_deriv; id<(thread_id+1)*thread_load_ko_deriv; id++){
            
        int i = id%(sz_x-6) + 3;
        int j = ((id/(sz_x-6))%(sz_y-8)) + 4;
        int k = (id/(sz_y-8)/(sz_x-6)) + 3;

        if(i >= nx-3 || j >= ny-4 || k >= nz-3) return;

        if(j==4) {
            int jb=3;
            output[IDX(i,jb,k)] = (-1.0 / 64.0 / hy) *
            (
                    -      dev_var_in[offset + IDX(i,jb+4,k)]
                    +  6.0*dev_var_in[offset + IDX(i,jb+3,k)]
                    - 15.0*dev_var_in[offset + IDX(i,jb+2,k)]
                    + 20.0*dev_var_in[offset + IDX(i,jb+1,k)]
                    - 15.0*dev_var_in[offset + IDX(i,jb,k)]
                    +  6.0*dev_var_in[offset + IDX(i,jb-1,k)]
                    -      dev_var_in[offset + IDX(i,jb-2,k)]
            );
        }

        int pp = IDX(i, j, k);

        output[pp] = (-1.0 / 64.0 / hy) *
        (
        -      dev_var_in[offset + pp-3*nx]
        +  6.0*dev_var_in[offset + pp-2*nx]
        - 15.0*dev_var_in[offset + pp-nx]
        + 20.0*dev_var_in[offset + pp]
        - 15.0*dev_var_in[offset + pp+nx]
        +  6.0*dev_var_in[offset + pp+2*nx]
        -      dev_var_in[offset + pp+3*nx]
        );

        if(j==5) {
            int je = ny - 3;
            output[IDX(i,je-1,k)] = (-1.0 / 64.0 / hy) *
            (
                    -      dev_var_in[offset + IDX(i,je+1,k)]
                    +  6.0*dev_var_in[offset + IDX(i,je,k)]
                    - 15.0*dev_var_in[offset + IDX(i,je-1,k)]
                    + 20.0*dev_var_in[offset + IDX(i,je-2,k)]
                    - 15.0*dev_var_in[offset + IDX(i,je-3,k)]
                    +  6.0*dev_var_in[offset + IDX(i,je-4,k)]
                    -      dev_var_in[offset + IDX(i,je-5,k)]
            );
        }
            if ((bflag & (1u<<OCT_DIR_DOWN)) && (j == 4)) {

            output[IDX(i,3,k)] =  (      dev_var_in[offset +IDX(i,6,k)]
                - 3.0*dev_var_in[offset +IDX(i,5,k)]
                + 3.0*dev_var_in[offset + IDX(i,4,k)]
                -     dev_var_in[offset + IDX(i,3,k)]
            )/59.0/48.0*64*hy;
            output[IDX(i,4,k)] =  (     dev_var_in[offset + IDX(i,7,k)]
                -  6.0*dev_var_in[offset + IDX(i,6,k)]
                + 12.0*dev_var_in[offset + IDX(i,5,k)]
                - 10.0*dev_var_in[offset + IDX(i,4,k)]
                +  3.0*dev_var_in[offset + IDX(i,3,k)]
            )/43.0/48.0*64*hy;
            output[IDX(i,5,k)] =  (     dev_var_in[offset + IDX(i,8,k)]
                -  6.0*dev_var_in[offset + IDX(i,7,k)]
                + 15.0*dev_var_in[offset + IDX(i,6,k)]
                - 19.0*dev_var_in[offset + IDX(i,5,k)]
                + 12.0*dev_var_in[offset + IDX(i,4,k)]
                -  3.0*dev_var_in[offset + IDX(i,3,k)]
            )/49.0/48.0*64*hy;
        }

            if ((bflag & (1u<<OCT_DIR_UP)) && (j == 5)) {

            const int je = ny - 3;
            output[IDX(i,je-3,k)] = (dev_var_in[offset + IDX(i,je-6,k)]
            -  6.0*dev_var_in[offset + IDX(i,je-5,k)]
            + 15.0*dev_var_in[offset + IDX(i,je-4,k)]
            - 19.0*dev_var_in[offset + IDX(i,je-3,k)]
            + 12.0*dev_var_in[offset + IDX(i,je-2,k)]
            -  3.0*dev_var_in[offset + IDX(i,je-1,k)]
            )/49.0/48.0*64*hy;

            output[IDX(i,je-2,k)] = (dev_var_in[offset + IDX(i,je-5,k)]
            -  6.0*dev_var_in[offset + IDX(i,je-4,k)]
            + 12.0*dev_var_in[offset + IDX(i,je-3,k)]
            - 10.0*dev_var_in[offset + IDX(i,je-2,k)]
            +  3.0*dev_var_in[offset + IDX(i,je-1,k)]
            )/43.0/48.0*64*hy;


            output[IDX(i,je-1,k)] = ( dev_var_in[offset + IDX(i,je-4,k)]
                -  3.0*dev_var_in[offset + IDX(i,je-3,k)]
                +  3.0*dev_var_in[offset + IDX(i,je-2,k)]
                -      dev_var_in[offset + IDX(i,je-1,k)]
            )/59.0/48.0*64*hy;
        }
    }
}



__device__ void device_calc_ko_deriv_z(double * output, double * dev_var_in,
    const int offset, double hz, int bflag,
    int nx,int ny,int nz, int sz_x, int sz_y, int sz_z){

    int thread_id = blockIdx.x*threads_per_block + threadIdx.x;

    for (int id = thread_id*thread_load_ko_deriv; id<(thread_id+1)*thread_load_ko_deriv; id++){
            
        int i = id%(sz_x-6) + 3;
        int j = ((id/(sz_x-6))%(sz_y-6)) + 3;
        int k = (id/(sz_y-6)/(sz_x-6)) + 4;

        if(i >= nx-3 || j >= ny-3 || k >= nz-4) return;

        if(k==4) {
            int kb=3;
            output[IDX(i,j,kb)] = (-1.0 / 64.0 / hz) *
            (
                    -      dev_var_in[offset + IDX(i,j,kb+4)]
                    +  6.0*dev_var_in[offset + IDX(i,j,kb+3)]
                    - 15.0*dev_var_in[offset + IDX(i,j,kb+2)]
                    + 20.0*dev_var_in[offset + IDX(i,j,kb+1)]
                    - 15.0*dev_var_in[offset + IDX(i,j,kb)]
                    +  6.0*dev_var_in[offset + IDX(i,j,kb-1)]
                    -      dev_var_in[offset + IDX(i,j,kb-2)]
            );
        }

        int pp = IDX(i, j, k);
        int n = nx * ny;
        output[pp] = (-1.0 / 64.0 / hz) *
        (
        -      dev_var_in[offset + pp-3*n]
        +  6.0*dev_var_in[offset + pp-2*n]
        - 15.0*dev_var_in[offset + pp-n]
        + 20.0*dev_var_in[offset + pp]
        - 15.0*dev_var_in[offset + pp+n]
        +  6.0*dev_var_in[offset + pp+2*n]
        -      dev_var_in[offset + pp+3*n]
        );

        if(k==5) {
            int ke = nz - 3;
            output[IDX(i,j,ke-1)] = (-1.0 / 64.0 / hz) *
            (
                    -      dev_var_in[offset + IDX(i,j,ke+1)]
                    +  6.0*dev_var_in[offset + IDX(i,j,ke)]
                    - 15.0*dev_var_in[offset + IDX(i,j,ke-1)]
                    + 20.0*dev_var_in[offset + IDX(i,j,ke-2)]
                    - 15.0*dev_var_in[offset + IDX(i,j,ke-3)]
                    +  6.0*dev_var_in[offset + IDX(i,j,ke-4)]
                    -      dev_var_in[offset + IDX(i,j,ke-5)]
            );
        }

        if ((bflag & (1u<<OCT_DIR_BACK)) && (k == 4)) {

            output[IDX(i,3,k)] =  (      dev_var_in[offset +IDX(i,k,6)]
                - 3.0*dev_var_in[offset +IDX(i,k,5)]
                + 3.0*dev_var_in[offset + IDX(i,k,4)]
                -     dev_var_in[offset + IDX(i,k,3)]
            )/59.0/48.0*64*hz;
            output[IDX(i,j,4)] =  (     dev_var_in[offset + IDX(i,j,7)]
                -  6.0*dev_var_in[offset + IDX(i,j,6)]
                + 12.0*dev_var_in[offset + IDX(i,j,5)]
                - 10.0*dev_var_in[offset + IDX(i,j,4)]
                +  3.0*dev_var_in[offset + IDX(i,j,3)]
            )/43.0/48.0*64*hz;
            output[IDX(i,j,5)] =  (     dev_var_in[offset + IDX(i,j,8)]
                -  6.0*dev_var_in[offset + IDX(i,j,7)]
                + 15.0*dev_var_in[offset + IDX(i,j,6)]
                - 19.0*dev_var_in[offset + IDX(i,j,5)]
                + 12.0*dev_var_in[offset + IDX(i,j,4)]
                -  3.0*dev_var_in[offset + IDX(i,j,3)]
            )/49.0/48.0*64*hz;
        }

        if ((bflag & (1u<<OCT_DIR_FRONT)) && (k == 5)) {

            const int ke = nz - 3;
            output[IDX(i,j,ke-3)] = (    dev_var_in[offset + IDX(i,j,ke-6)]
                -  6.0*dev_var_in[offset + IDX(i,j,ke-5)]
                + 15.0*dev_var_in[offset + IDX(i,j,ke-4)]
                - 19.0*dev_var_in[offset + IDX(i,j,ke-3)]
                + 12.0*dev_var_in[offset + IDX(i,j,ke-2)]
                -  3.0*dev_var_in[offset + IDX(i,j,ke-1)]
            )/49.0/48.0*64*hz;

            output[IDX(i,j,ke-2)] = (   dev_var_in[offset + IDX(i,j,ke-5)]
                -  6.0*dev_var_in[offset + IDX(i,j,ke-4)]
                + 12.0*dev_var_in[offset + IDX(i,j,ke-3)]
                - 10.0*dev_var_in[offset + IDX(i,j,ke-2)]
                +  3.0*dev_var_in[offset + IDX(i,j,ke-1)]
            )/43.0/48.0*64*hz;


            output[IDX(i,j,ke-1)] = (   dev_var_in[offset + IDX(i,j,ke-4)]
                -  3.0*dev_var_in[offset + IDX(i,j,ke-3)]
                +  3.0*dev_var_in[offset + IDX(i,j,ke-2)]
                -      dev_var_in[offset + IDX(i,j,ke-1)]
            )/59.0/48.0*64*hz;
        }
    }
}

 
 __global__ void cuda_calc_ko_deriv_all(double * dev_var_in, double hx, double hy, double hz, 
    int sz_x, int sz_y, int sz_z, int bflag,
    #include "list_of_para.h"
    ,
    #include "list_of_offset_para.h"
 ) {
    int nx = sz_x;
    int ny = sz_y;
    int nz = sz_z;

    #include "bssnrhs_cuda_ko_derivs.h"
 }

 void calc_ko_deriv_all( double * dev_var_in, double hx, double hy, double hz, int sz_x, 
    int sz_y, int sz_z, int bflag, cudaStream_t stream,
    #include "list_of_para.h"
    ,
    #include "list_of_offset_para.h"
 )
 {
    const int ie = sz_x - 3;//x direction
    const int je = sz_y - 3;//y direction
    const int ke = sz_z - 3;//z direction
 
    int total_points = ceil(1.0*ie*je*ke/thread_load_ko_deriv);
    int blocks = ceil(1.0*total_points/threads_per_block);

    cuda_calc_ko_deriv_all <<< blocks, threads_per_block, 0, stream >>> (
                    dev_var_in, hx, hy, hz, sz_x, sz_y, sz_z, bflag,
                  #include "list_of_args.h"
                  ,
                  #include "list_of_offset_args.h"
                );
 }