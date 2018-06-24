#include "derivs_cuda.h"

__device__ void device_calc_adv_x(double * output, double * dev_var_in,
    const int offset, double hx, int bflag,
    int nx,int ny,int nz, int betax, int sz_x, int sz_y, int sz_z){

    int thread_id = blockIdx.x*threads_per_block + threadIdx.x;

    for (int id = thread_id*thread_load_adv_deriv; id<(thread_id+1)*thread_load_adv_deriv; id++){
            
        int i = id%(sz_x-6) + 3;
        int j = ((id/(sz_x-6))%(sz_y-6)) + 3;
        int k = (id/(sz_y-6)/(sz_x-6)) + 3;

        double idx_by_2 = 0.50 * (1.0 / hx);
        double idx_by_12 = (1.0 / hx)/12;

        if(i >= nx-3 || j >= ny-3 || k >= nz-3) return;

        int pp = IDX(i, j, k);

        if (dev_var_in[betax + pp] > 0.0 ) {
            output[pp] = ( -  3.0 * dev_var_in[offset + pp - 1]
            - 10.0 * dev_var_in[offset + pp]
            + 18.0 * dev_var_in[offset + pp + 1]
            -  6.0 * dev_var_in[offset + pp + 2]
            +        dev_var_in[offset + pp + 3]
            ) * idx_by_12;
        }
        else {
            output[pp] = ( -        dev_var_in[offset + pp - 3]
            +  6.0 * dev_var_in[offset + pp - 2]
            - 18.0 * dev_var_in[offset + pp - 1]
            + 10.0 * dev_var_in[offset + pp]
            +  3.0 * dev_var_in[offset + pp +1]
            ) * idx_by_12;
        }

        if ((bflag & (1u<<OCT_DIR_LEFT)) && (i == 3)) {

            output[IDX(3,j,k)] = ( -  3.0 * dev_var_in[offset + IDX(3,j,k)]
            +  4.0 * dev_var_in[offset + IDX(4,j,k)]
            -        dev_var_in[offset + IDX(5,j,k)]
            ) * idx_by_2;

            if (dev_var_in[betax + IDX(4,j,k)] > 0.0) {
            output[IDX(4,j,k)] = ( -  3.0 * dev_var_in[offset + IDX(4,j,k)]
                +  4.0 * dev_var_in[offset + IDX(5,j,k)]
                -        dev_var_in[offset + IDX(6,j,k)]
            ) * idx_by_2;
            }
            else {
            output[IDX(4,j,k)] = ( -         dev_var_in[offset + IDX(3,j,k)]
                +        dev_var_in[offset + IDX(5,j,k)]
            ) * idx_by_2;
            }

            if (dev_var_in[betax + IDX(5,j,k)] > 0.0 ) {
            output[IDX(5,j,k)] = (-  3.0 * dev_var_in[offset + IDX(4,j,k)]
                - 10.0 * dev_var_in[offset + IDX(5,j,k)]
                + 18.0 * dev_var_in[offset + IDX(6,j,k)]
                -  6.0 * dev_var_in[offset + IDX(7,j,k)]
                +        dev_var_in[offset + IDX(8,j,k)]
            ) * idx_by_12;
            }
            else {
            output[IDX(5,j,k)] = (           dev_var_in[offset + IDX(3,j,k)]
                        -  4.0 * dev_var_in[offset + IDX(4,j,k)]
                        +  3.0 * dev_var_in[offset + IDX(5,j,k)]
            ) * idx_by_2;
            }
        }

        if ((bflag & (1u<<OCT_DIR_RIGHT)) && (i == 4)) {

            const int ie = nx - 3;

            if ( dev_var_in[betax + IDX(ie-3,j,k)] < 0.0 ) {
                output[IDX(ie-3,j,k)] = (  - 3.0 * dev_var_in[offset + IDX(ie-3,j,k)]
                        + 4.0 * dev_var_in[offset + IDX(ie-2,j,k)]
                        -       dev_var_in[offset + IDX(ie-1,j,k)]
                    ) * idx_by_2;
            }
            else {
                output[IDX(ie-3,j,k)] = ( -   dev_var_in[offset + IDX(ie-6,j,k)]
                        +  6.0 * dev_var_in[offset + IDX(ie-5,j,k)]
                        - 18.0 * dev_var_in[offset + IDX(ie-4,j,k)]
                        + 10.0 * dev_var_in[offset + IDX(ie-3  ,j,k)]
                        +  3.0 * dev_var_in[offset + IDX(ie-2,j,k)]
                    ) * idx_by_12;
            }

            if (dev_var_in[betax + IDX(ie-2,j,k)] > 0.0 ) {
                output[IDX(ie-2,j,k)] = (  -  dev_var_in[offset + IDX(ie-3,j,k)]
                        +  dev_var_in[offset + IDX(ie-1,j,k)]
                    ) * idx_by_2;
            }
            else {
                output[IDX(ie-2,j,k)] = (     dev_var_in[offset + IDX(ie-4,j,k)]
                            - 4.0 * dev_var_in[offset + IDX(ie-3,j,k)]
                            + 3.0 * dev_var_in[offset + IDX(ie-2,j,k)]
                    ) * idx_by_2;
            }

            output[IDX(ie-1,j,k)] = (          dev_var_in[offset + IDX(ie-3,j,k)]
                        - 4.0 * dev_var_in[offset + IDX(ie-2,j,k)]
                        + 3.0 * dev_var_in[offset + IDX(ie-1,j,k)]
            ) * idx_by_2;
        }
    }
}
 
__device__ void device_calc_adv_y(double * output, double * dev_var_in,
    const int offset, double hy, int bflag,
    int nx,int ny,int nz, int betay, int sz_x, int sz_y, int sz_z){

    int thread_id = blockIdx.x*threads_per_block + threadIdx.x;

    for (int id = thread_id*thread_load_adv_deriv; id<(thread_id+1)*thread_load_adv_deriv; id++){
            
        int i = id%(sz_x-6) + 3;
        int j = ((id/(sz_x-6))%(sz_y-6)) + 3;
        int k = (id/(sz_y-6)/(sz_x-6)) + 3;

        double idy_by_2 = 0.50 * (1.0 / hy);
        double idy_by_12 = (1.0 / hy)/12.0;

        if(i >= nx-3 || j >= ny-3 || k >= nz-3) return;

        int pp = IDX(i, j, k);

        if (dev_var_in[betay + pp] > 0.0 ) {
            output[pp] = ( -  3.0 * dev_var_in[offset + pp - nx]
            - 10.0 * dev_var_in[offset + pp]
            + 18.0 * dev_var_in[offset + pp + nx]
            -  6.0 * dev_var_in[offset + pp + 2*nx]
            +        dev_var_in[offset + pp + 3*nx]
            ) * idy_by_12;
        }
        else {
            output[pp] = ( -        dev_var_in[offset + pp - 3*nx]
            +  6.0 * dev_var_in[offset + pp - 2*nx]
            - 18.0 * dev_var_in[offset + pp - nx]
            + 10.0 * dev_var_in[offset + pp]
            +  3.0 * dev_var_in[offset + pp +nx]
            ) * idy_by_12;

        }

        if ((bflag & (1u<<OCT_DIR_DOWN)) && (j == 3)) {

            output[IDX(i,3,k)] = ( -  3.0 * dev_var_in[offset + IDX(i,3,k)]
            +  4.0 * dev_var_in[offset + IDX(i,4,k)]
            -        dev_var_in[offset + IDX(i,5,k)]
            ) * idy_by_2;

            if (dev_var_in[betay + IDX(i,4,k)] > 0.0) {
                output[IDX(i,4,k)] = ( -  3.0 * dev_var_in[offset + IDX(i,4,k)]
                    +  4.0 * dev_var_in[offset + IDX(i,5,k)]
                    -        dev_var_in[offset + IDX(i,6,k)]
                ) * idy_by_2;

            }
            else {
                output[IDX(i,4,k)] = ( -         dev_var_in[offset + IDX(i,3,k)]
                    +        dev_var_in[offset + IDX(i,5,k)]
                ) * idy_by_2;

            }

            if (dev_var_in[betay + IDX(i,5,k)] > 0.0 ) {
                output[IDX(i,5,k)] = (-  3.0 * dev_var_in[offset + IDX(i,4,k)]
                    - 10.0 * dev_var_in[offset + IDX(i,5,k)]
                    + 18.0 * dev_var_in[offset + IDX(i,6,k)]
                    -  6.0 * dev_var_in[offset + IDX(i,7,k)]
                    +        dev_var_in[offset + IDX(i,8,k)]
                ) * idy_by_12;
            }
            else {
                output[IDX(i,5,k)] = (           dev_var_in[offset + IDX(i,3,k)]
                            -  4.0 * dev_var_in[offset + IDX(i,4,k)]
                            +  3.0 * dev_var_in[offset + IDX(i,5,k)]
                ) * idy_by_2;
            }
        }

        if ((bflag & (1u<<OCT_DIR_UP)) && (j == 4)) {

            const int je = ny - 3;

            if ( dev_var_in[betay + IDX(i,je-3,k)] < 0.0 ) {
            output[IDX(i,je-3,k)] = (  - 3.0 * dev_var_in[offset + IDX(i,je-3,k)]
                    + 4.0 * dev_var_in[offset + IDX(i,je-2,k)]
                    -       dev_var_in[offset + IDX(i,je-1,k)]
                ) * idy_by_2;
            }
            else {
            output[IDX(i,je-3,k)] = ( -   dev_var_in[offset + IDX(i,je-6,k)]
                    +  6.0 * dev_var_in[offset + IDX(i,je-5,k)]
                    - 18.0 * dev_var_in[offset + IDX(i,je-4,k)]
                    + 10.0 * dev_var_in[offset + IDX(i,je-3,k)]
                    +  3.0 * dev_var_in[offset + IDX(i,je-2,k)]
                ) * idy_by_12;
            }

            if (dev_var_in[betay + IDX(i,je-2,k)] > 0.0 ) {
            output[IDX(i,je-2,k)] = (  -  dev_var_in[offset + IDX(i,je-3,k)]
                    +  dev_var_in[offset + IDX(i,je-1,k)]
                ) * idy_by_2;
            }
            else {
            output[IDX(i,je-2,k)] = (     dev_var_in[offset + IDX(i,je-4,k)]
                        - 4.0 * dev_var_in[offset + IDX(i,je-3,k)]
                        + 3.0 * dev_var_in[offset + IDX(i,je-2,k)]
                ) * idy_by_2;
            }

            output[IDX(i,je-1,k)]  = (          dev_var_in[offset + IDX(i,je-3,k)]
                        - 4.0 * dev_var_in[offset + IDX(i,je-2,k)]
                        + 3.0 * dev_var_in[offset + IDX(i,je-1,k)]
            ) * idy_by_2;
        }
    }
}

 
__device__ void device_calc_adv_z(double * output, double * dev_var_in,
    const int offset, double hz, int bflag,
    int nx,int ny,int nz, int betaz, int sz_x, int sz_y, int sz_z){

    int thread_id = blockIdx.x*threads_per_block + threadIdx.x;

    for (int id = thread_id*thread_load_adv_deriv; id<(thread_id+1)*thread_load_adv_deriv; id++){
        
        int i = id%(sz_x-6) + 3;
        int j = ((id/(sz_x-6))%(sz_y-6)) + 3;
        int k = (id/(sz_y-6)/(sz_x-6)) + 3;
        

        double idz_by_2 = 0.50 * (1.0 / hz);
        double idz_by_12 = (1.0 / hz)/12.0;

        if(i >= nx-3 || j >= ny-3 || k >= nz-3) return;

        int n = nx * ny;
        int pp = IDX(i, j, k);

        if (dev_var_in[betaz + pp] > 0.0 ) {
            output[pp] = ( -  3.0 * dev_var_in[offset + pp - n]
            - 10.0 * dev_var_in[offset + pp]
            + 18.0 * dev_var_in[offset + pp + n]
            -  6.0 * dev_var_in[offset + pp + 2*n]
            +        dev_var_in[offset + pp + 3*n]
            ) * idz_by_12;
        }
        else {
            output[pp] = ( -        dev_var_in[offset + pp - 3*n]
            +  6.0 * dev_var_in[offset + pp - 2*n]
            - 18.0 * dev_var_in[offset + pp - n]
            + 10.0 * dev_var_in[offset + pp]
            +  3.0 * dev_var_in[offset + pp +n]
            ) * idz_by_12;

        }

        if ((bflag & (1u<<OCT_DIR_BACK)) && (k == 3)) {

            output[IDX(i,j,3)] = ( -  3.0 * dev_var_in[offset + IDX(i,j,3)]
            +  4.0 * dev_var_in[offset + IDX(i,j,4)]
            -        dev_var_in[offset + IDX(i,j,5)]
            ) * idz_by_2;

            if (dev_var_in[betaz + IDX(i,j,4)] > 0.0) {
            output[IDX(i,j,4)] = ( -  3.0 * dev_var_in[offset + IDX(i,j,4)]
                +  4.0 * dev_var_in[offset + IDX(i,j,5)]
                -        dev_var_in[offset + IDX(i,j,6)]
            ) * idz_by_2;

            }
            else {
            output[IDX(i,j,4)] = ( -         dev_var_in[offset + IDX(i,j,3)]
                +        dev_var_in[offset + IDX(i,j,5)]
            ) * idz_by_2;

            }

            if (dev_var_in[betaz + IDX(i,j,5)] > 0.0 ) {
            output[IDX(i,j,5)] = (-  3.0 * dev_var_in[offset + IDX(i,j,4)]
                - 10.0 * dev_var_in[offset + IDX(i,j,5)]
                + 18.0 * dev_var_in[offset + IDX(i,j,6)]
                -  6.0 * dev_var_in[offset + IDX(i,j,7)]
                +        dev_var_in[offset + IDX(i,j,8)]
            ) * idz_by_12;
            }
            else {
            output[IDX(i,j,5)] = (           dev_var_in[offset + IDX(i,j,3)]
                        -  4.0 * dev_var_in[offset + IDX(i,j,4)]
                        +  3.0 * dev_var_in[offset + IDX(i,j,5)]
            ) * idz_by_2;
            }
        }

        if ((bflag & (1u<<OCT_DIR_FRONT)) && (k == 4)) {

            const int ke = nz - 3;

            if ( dev_var_in[betaz + IDX(i,j,ke-3)] < 0.0 ) {
            output[IDX(i,j,ke-3)] = (  - 3.0 * dev_var_in[offset + IDX(i,j,ke-3)]
                    + 4.0 * dev_var_in[offset + IDX(i,j,ke-2)]
                    -       dev_var_in[offset + IDX(i,j,ke-1)]
                ) * idz_by_2;
            }
            else {
            output[IDX(i,j,ke-3)] = ( -   dev_var_in[offset + IDX(i,j,ke-6)]
                    +  6.0 * dev_var_in[offset + IDX(i,j,ke-5)]
                    - 18.0 * dev_var_in[offset + IDX(i,j,ke-4)]
                    + 10.0 * dev_var_in[offset + IDX(i,j,ke-3)]
                    +  3.0 * dev_var_in[offset + IDX(i,j,ke-2)]
                ) * idz_by_12;
            }

            if (dev_var_in[betaz + IDX(i,j,ke-2)] > 0.0 ) {
            output[IDX(i,j,ke-2)] = (  -  dev_var_in[offset + IDX(i,j,ke-3)]
                    +  dev_var_in[offset + IDX(i,j,ke-1)]
                ) * idz_by_2;
            }
            else {
            output[IDX(i,j,ke-2)] = (     dev_var_in[offset + IDX(i,j,ke-4)]
                        - 4.0 * dev_var_in[offset + IDX(i,j,ke-3)]
                        + 3.0 * dev_var_in[offset + IDX(i,j,ke-2)]
                ) * idz_by_2;
            }

            output[IDX(i,j,ke-1)]  = (          dev_var_in[offset + IDX(i,j,ke-3)]
                        - 4.0 * dev_var_in[offset + IDX(i,j,ke-2)]
                        + 3.0 * dev_var_in[offset + IDX(i,j,ke-1)]
            ) * idz_by_2;
        }
    }
}
 

 
 __global__ void calc_all_adv(double * dev_var_in, double hx, double hy, double hz, 
    int sz_x, int sz_y, int sz_z, int bflag,
    #include "list_of_para.h"
    ,
    #include "list_of_offset_para.h"
 ) {
    int nx = sz_x;
    int ny = sz_y;
    int nz = sz_z;
 
     #include "bssnrhs_cuda_derivs_adv.h"
     //ib, jb, kb values are accumulated to the x, y, z
 
 }

 void cuda_deriv_calc_all_adv(double * dev_var_in, double hx, double hy, double hz, int sz_x, 
    int sz_y, int sz_z, int bflag, cudaStream_t stream,
    #include "list_of_para.h"
    ,
    #include "list_of_offset_para.h"
    ){
     const int ie = sz_x - 3;//x direction
     const int je = sz_y - 3;//y direction
     const int ke = sz_z - 3;//z direction
     
     int total_points = ceil(1.0*ie*je*ke/thread_load_adv_deriv);
     int blocks = ceil(1.0*total_points/threads_per_block);

     calc_all_adv <<< blocks, threads_per_block, 0, stream >>> (
                    dev_var_in, hx, hy, hz, sz_x, sz_y, sz_z, bflag,
                  #include "list_of_args.h"
                  ,
                  #include "list_of_offset_args.h"
           );
}
