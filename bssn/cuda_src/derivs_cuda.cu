#include "derivs_cuda.h"

 
__device__ void device_calc_deriv_x(double * output, double * dev_var_in,
    const int offset, double hx, int bflag,
    int nx,int ny,int nz, int sz_x, int sz_y, int sz_z){
    
    int id = blockIdx.x*threads_per_block + threadIdx.x;

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
 
 __device__ void device_calc_deriv_y(double * output, double * dev_var_in,
    const int offset, double hx, int bflag,
    int nx,int ny,int nz, int sz_x, int sz_y, int sz_z){
        
    int id = blockIdx.x*threads_per_block + threadIdx.x;

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

 
 __device__ void device_calc_deriv_z(double * output, double * dev_var_in,
    const int offset, double hx, int bflag,
    int nx,int ny,int nz, int sz_x, int sz_y, int sz_z){
        
    int id = blockIdx.x*threads_per_block + threadIdx.x;

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

__device__ void device_calc_deriv_xx(double * output, double * dev_var_in,
    const int offset, double hx, int bflag,
    int nx,int ny,int nz, int sz_x, int sz_y, int sz_z){
        
    int id = blockIdx.x*threads_per_block + threadIdx.x;

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
 
__device__ void device_calc_deriv_yy(double * output, double * dev_var_in,
    const int offset, double hy, int bflag,
    int nx,int ny,int nz, int sz_x, int sz_y, int sz_z){
        
    int id = blockIdx.x*threads_per_block + threadIdx.x;

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

__device__ void device_calc_deriv_zz(double * output, double * dev_var_in,
    const int offset, double hz, int bflag,
    int nx,int ny,int nz, int sz_x, int sz_y, int sz_z){
        
    int id = blockIdx.x*threads_per_block + threadIdx.x;

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
 
    int total_points = ie*je*ke;
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



__device__ void device_calc_adv_x(double * output, double * dev_var_in,
    const int offset, double hx, int bflag,
    int nx,int ny,int nz, int betax, int sz_x, int sz_y, int sz_z){

    int id = blockIdx.x*threads_per_block + threadIdx.x;
    
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
 
__device__ void device_calc_adv_y(double * output, double * dev_var_in,
    const int offset, double hy, int bflag,
    int nx,int ny,int nz, int betay, int sz_x, int sz_y, int sz_z){

    int id = blockIdx.x*threads_per_block + threadIdx.x;
    
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

 
__device__ void device_calc_adv_z(double * output, double * dev_var_in,
    const int offset, double hz, int bflag,
    int nx,int ny,int nz, int betaz, int sz_x, int sz_y, int sz_z){

    int id = blockIdx.x*threads_per_block + threadIdx.x;

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
     
     int total_points = ie*je*ke;
     int blocks = ceil(1.0*total_points/threads_per_block);

     calc_all_adv <<< blocks, threads_per_block, 0, stream >>> (
                    dev_var_in, hx, hy, hz, sz_x, sz_y, sz_z, bflag,
                  #include "list_of_args.h"
                  ,
                  #include "list_of_offset_args.h"
           );
}

__device__ void device_calc_ko_deriv_x(double * output, double * dev_var_in,
    const int offset, double hx, int bflag,
    int nx,int ny,int nz, int sz_x, int sz_y, int sz_z){

    int id = blockIdx.x*threads_per_block + threadIdx.x;

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

 
__device__ void device_calc_ko_deriv_y(double * output, double * dev_var_in,
    const int offset, double hy, int bflag,
    int nx,int ny,int nz, int sz_x, int sz_y, int sz_z){

    int id = blockIdx.x*threads_per_block + threadIdx.x;

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



__device__ void device_calc_ko_deriv_z(double * output, double * dev_var_in,
    const int offset, double hz, int bflag,
    int nx,int ny,int nz, int sz_x, int sz_y, int sz_z){

    int id = blockIdx.x*threads_per_block + threadIdx.x;

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
 
    int total_points = ie*je*ke;
    int blocks = ceil(1.0*total_points/threads_per_block);

    cuda_calc_ko_deriv_all <<< blocks, threads_per_block, 0, stream >>> (
                    dev_var_in, hx, hy, hz, sz_x, sz_y, sz_z, bflag,
                  #include "list_of_args.h"
                  ,
                  #include "list_of_offset_args.h"
                );
 }