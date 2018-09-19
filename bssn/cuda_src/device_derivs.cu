/**
 * Created on: March 15, 2018
 * 		Author: Akila
 **/

 #include "derivs_cuda.cuh"
 
 __device__ void calc_deriv42_x(int id, double * output, double * dev_var_in, const int u_offset, double dx, 
    const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag){

    int nx = host_sz_x; 
    int ny = host_sz_y; 

    const double idx = 1.0/dx;
    const double idx_by_12 = idx / 12.0;

    int i = id%(host_sz_x-6) + 3;
    int j = ((id/(host_sz_x-6))%(host_sz_y-2)) + 1;
    int k = (id/(host_sz_z-2)/(host_sz_x-6)) + 1; 
    if (k>=host_sz_z-1) return;
    int pp = IDX(i, j, k);

    output[pp] = (dev_var_in[(u_offset) + pp - 2] - 8.0*dev_var_in[(u_offset)
                    + pp - 1] + 8.0*dev_var_in[(u_offset) + pp + 1] 
                    - dev_var_in[(u_offset) + pp + 2] )*idx_by_12;
 }

 __device__ void calc_deriv42_y(int id, double* output, double * dev_var_in, const int u_offset, double dy, 
    const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag){

    int nx = host_sz_x; 
    int ny = host_sz_y; 

    const double idy = 1.0/dy;
    const double idy_by_12 = idy / 12.0;
    
    int i = id%(host_sz_x-6) + 3;
    int j = ((id/(host_sz_x-6))%(host_sz_y-6)) + 3;
    int k = (id/(host_sz_z-6)/(host_sz_x-6)) + 1;
    if (k>=host_sz_z-1) return;
    int pp = IDX(i, j, k);
    
    output[pp] = (dev_var_in[u_offset + pp - 2*nx] 
                - 8.0*dev_var_in[u_offset + pp - nx] 
                + 8.0*dev_var_in[u_offset + pp + nx] 
                - dev_var_in[u_offset + pp + 2*nx] )*idy_by_12;
    
}

__device__ void calc_deriv42_z(int id, double* output, double * dev_var_in, const int u_offset, double dz, 
    const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag){
    int nx = host_sz_x; 
    int ny = host_sz_y; 
    int n = nx * ny;

    const double idz = 1.0/dz;
    const double idz_by_12 = idz / 12.0;

    int i = id%(host_sz_x-6) + 3;
    int j = ((id/(host_sz_x-6))%(host_sz_y-6)) + 3;
    int k = (id/(host_sz_z-6)/(host_sz_x-6)) + 3;
    if (k>=host_sz_z-3) return;
    int pp = IDX(i, j, k);
    output[pp] = (dev_var_in[(u_offset) + pp - 2*n] - 8.0*dev_var_in[(u_offset) + pp - n] 
                    + 8.0*dev_var_in[(u_offset) + pp + n] - dev_var_in[(u_offset) + pp + 2*n]) 
                    * idz_by_12;
}

__device__ void calc_deriv42_xx(int id, double* output, double * dev_var_in, const int u_offset, double dx, 
    const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag){
    int nx = host_sz_x; 
    int ny = host_sz_y; 

    const double idx_sqrd = 1.0/(dx*dx);
    const double idx_sqrd_by_12 = idx_sqrd / 12.0;
    
    int i = id%(host_sz_x-6) + 3;
    int j = ((id/(host_sz_x-6))%(host_sz_y-6)) + 3;
    int k = (id/(host_sz_z-6)/(host_sz_x-6)) + 3;
    if (k>=host_sz_z-3) return;
    int pp = IDX(i, j, k);
    output[pp] = ((-1)*dev_var_in[(u_offset) + pp - 2] 
                + 16.0*dev_var_in[(u_offset) + pp - 1] 
                - 30.0*dev_var_in[(u_offset) + pp] 
                + 16.0*dev_var_in[(u_offset) + pp + 1] 
                - dev_var_in[(u_offset) + pp + 2] 
            )*idx_sqrd_by_12;

}

__device__ void calc_deriv42_yy(int id, double* output, double * dev_var_in, const int u_offset, double dy, 
    const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag){
    int nx = host_sz_x; 
    int ny = host_sz_y; 

    const double idy_sqrd = 1.0/(dy*dy);
    const double idy_sqrd_by_12 = idy_sqrd / 12.0;

    int i = id%(host_sz_x-6) + 3;
    int j = ((id/(host_sz_x-6))%(host_sz_y-6)) + 3;
    int k = (id/(host_sz_z-6)/(host_sz_x-6)) + 3;
    if (k>=host_sz_z-3) return;
    int pp = IDX(i, j, k);
    output[pp] = ((-1)*dev_var_in[(u_offset) + pp - 2*nx] 
                + 16.0*dev_var_in[(u_offset) + pp - nx] 
                - 30.0*dev_var_in[(u_offset) + pp] 
                + 16.0*dev_var_in[(u_offset) + pp + nx] 
                - dev_var_in[(u_offset) + pp + 2*nx] 
            )*idy_sqrd_by_12;
            
}

__device__ void calc_deriv42_zz(int id, double* output, double * dev_var_in, const int u_offset, double dz, 
    const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag){
    int nx = host_sz_x; 
    int ny = host_sz_y; 
    int n = nx * ny;

    const double idz_sqrd = 1.0/(dz*dz);
    const double idz_sqrd_by_12 = idz_sqrd / 12.0;

    int i = id%(host_sz_x-6) + 3;
    int j = ((id/(host_sz_x-6))%(host_sz_y-6)) + 3;
    int k = (id/(host_sz_z-6)/(host_sz_x-6)) + 3;
    if (k>=host_sz_z-3) return;
    int pp = IDX(i, j, k);
    output[pp] = ((-1)*dev_var_in[(u_offset) + pp - 2*n] 
                + 16.0*dev_var_in[(u_offset) + pp - n] 
                - 30.0*dev_var_in[(u_offset) + pp] 
                + 16.0*dev_var_in[(u_offset) + pp + n] 
                - dev_var_in[(u_offset) + pp + 2*n] 
            )*idz_sqrd_by_12;
}

__device__ void calc_deriv42_adv_x(int id, double * output, double * dev_var_in, int u_offset, double dx, int betax,
const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag) {
    int nx = host_sz_x;
    int ny = host_sz_y;

    const double idx = 1.0/dx;
    // const double idx_by_2 = 0.50 * idx;
    const double idx_by_12 = idx / 12.0;

    int i = id%(host_sz_x-6) + 3;
    int j = ((id/(host_sz_x-6))%(host_sz_y-6)) + 3;
    int k = (id/(host_sz_z-6)/(host_sz_x-6)) + 3;
    if (k>=host_sz_z-3) return;
    int pp = IDX(i, j, k);

    if (dev_var_in[betax + pp] > 0.0 ) {
        output[pp] = ( -  3.0 * dev_var_in[u_offset + pp - 1]
                    - 10.0 * dev_var_in[u_offset + pp]
                    + 18.0 * dev_var_in[u_offset + pp + 1]
                    -  6.0 * dev_var_in[u_offset + pp + 2]
                    +        dev_var_in[u_offset + pp + 3]
                ) * idx_by_12;
    }
    else {
        output[pp] = ( -        dev_var_in[u_offset + pp - 3]
                    +  6.0 * dev_var_in[u_offset + pp - 2]
                    - 18.0 * dev_var_in[u_offset + pp - 1]
                    + 10.0 * dev_var_in[u_offset + pp]
                    +  3.0 * dev_var_in[u_offset + pp +1]
                ) * idx_by_12;
    }
}

__device__ void calc_deriv42_adv_y(int id, double * output, double * dev_var_in, int u_offset, double dy, int betay,
const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag) {
    
    int nx = host_sz_x;
    int ny = host_sz_y;

    const double idy = 1.0/dy;
    // const double idy_by_2 = 0.50 * idy;
    const double idy_by_12 = idy / 12.0;

    int i = id%(host_sz_x-6) + 3;
    int j = ((id/(host_sz_x-6))%(host_sz_y-6)) + 3;
    int k = (id/(host_sz_z-6)/(host_sz_x-6)) + 3;
    if (k>=host_sz_z-3) return;
    int pp = IDX(i, j, k);

    if (dev_var_in[betay + pp] > 0.0 ) {
            output[pp] = ( -  3.0 * dev_var_in[u_offset + pp - nx]
                        - 10.0 * dev_var_in[u_offset + pp]
                        + 18.0 * dev_var_in[u_offset + pp + nx]
                        -  6.0 * dev_var_in[u_offset + pp + 2*nx]
                        +        dev_var_in[u_offset + pp + 3*nx]
                    ) * idy_by_12;
    }
    else {
        output[pp] = ( -        dev_var_in[u_offset + pp - 3*nx]
                    +  6.0 * dev_var_in[u_offset + pp - 2*nx]
                    - 18.0 * dev_var_in[u_offset + pp - nx]
                    + 10.0 * dev_var_in[u_offset + pp]
                    +  3.0 * dev_var_in[u_offset + pp +nx]
                    ) * idy_by_12;
                
    }
}

__device__ void calc_deriv42_adv_z(int id, double * output, double * dev_var_in, int u_offset, double dz, int betaz, 
const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag) {

    int nx = host_sz_x;
    int ny = host_sz_y;
    int n = nx * ny;
    
    const double idz = 1.0/dz;
    // const double idz_by_2 = 0.50 * idz;
    const double idz_by_12 = idz / 12.0;

    int i = id%(host_sz_x-6) + 3;
    int j = ((id/(host_sz_x-6))%(host_sz_y-6)) + 3;
    int k = (id/(host_sz_z-6)/(host_sz_x-6)) + 3;
    if (k>=host_sz_z-3) return;
    int pp = IDX(i, j, k);

    if (dev_var_in[betaz + pp] > 0.0 ) {
            output[pp] = ( -  3.0 * dev_var_in[u_offset + pp - n]
                        - 10.0 * dev_var_in[u_offset + pp]
                        + 18.0 * dev_var_in[u_offset + pp + n]
                        -  6.0 * dev_var_in[u_offset + pp + 2*n]
                        +        dev_var_in[u_offset + pp + 3*n]
                    ) * idz_by_12;
    }
    else {
        output[pp] = ( -        dev_var_in[u_offset + pp - 3*n]
                    +  6.0 * dev_var_in[u_offset + pp - 2*n]
                    - 18.0 * dev_var_in[u_offset + pp - n]
                    + 10.0 * dev_var_in[u_offset + pp]
                    +  3.0 * dev_var_in[u_offset + pp +n]
                    ) * idz_by_12;
                
    }
}

__global__ void calc_derivs1(
    double * dev_var_out, double * dev_var_in, double hx, double hy, double hz, 
    const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, 
    int bflag,
    #include "list_of_offset_para.h"
    ,
    #include "list_of_para.h"
){
    int tid = blockIdx.x*64 + threadIdx.x;

    calc_deriv42_x(tid, grad_0_alpha, dev_var_in, alphaInt, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_xx(tid, grad2_0_0_alpha, dev_var_in, alphaInt, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad_1_alpha, dev_var_in, alphaInt, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_yy(tid, grad2_1_1_alpha, dev_var_in, alphaInt, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad_2_alpha, dev_var_in, alphaInt, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_zz(tid, grad2_2_2_alpha, dev_var_in, alphaInt, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x(tid, grad_0_beta0, dev_var_in, beta0Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_xx(tid, grad2_0_0_beta0, dev_var_in, beta0Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad_1_beta0, dev_var_in, beta0Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_yy(tid, grad2_1_1_beta0, dev_var_in, beta0Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad_2_beta0, dev_var_in, beta0Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_zz(tid, grad2_2_2_beta0, dev_var_in, beta0Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x(tid, grad_0_beta1, dev_var_in, beta1Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_xx(tid, grad2_0_0_beta1, dev_var_in, beta1Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad_1_beta1, dev_var_in, beta1Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_yy(tid, grad2_1_1_beta1, dev_var_in, beta1Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad_2_beta1, dev_var_in, beta1Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_zz(tid, grad2_2_2_beta1, dev_var_in, beta1Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x(tid, grad_0_beta2, dev_var_in, beta2Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_xx(tid, grad2_0_0_beta2, dev_var_in, beta2Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad_1_beta2, dev_var_in, beta2Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_yy(tid, grad2_1_1_beta2, dev_var_in, beta2Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad_2_beta2, dev_var_in, beta2Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_zz(tid, grad2_2_2_beta2, dev_var_in, beta2Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x(tid, grad_0_B0, dev_var_in, B0Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad_1_B0, dev_var_in, B0Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad_2_B0, dev_var_in, B0Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x(tid, grad_0_B1, dev_var_in, B1Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad_1_B1, dev_var_in, B1Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad_2_B1, dev_var_in, B1Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x(tid, grad_0_B2, dev_var_in, B2Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad_1_B2, dev_var_in, B2Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad_2_B2, dev_var_in, B2Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x(tid, grad_0_chi, dev_var_in, chiInt, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_xx(tid, grad2_0_0_chi, dev_var_in, chiInt, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad_1_chi, dev_var_in, chiInt, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_yy(tid, grad2_1_1_chi, dev_var_in, chiInt, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad_2_chi, dev_var_in, chiInt, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_zz(tid, grad2_2_2_chi, dev_var_in, chiInt, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x(tid, grad_0_Gt0, dev_var_in, Gt0Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad_1_Gt0, dev_var_in, Gt0Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad_2_Gt0, dev_var_in, Gt0Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x(tid, grad_0_gt0, dev_var_in, gt0Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_xx(tid, grad2_0_0_gt0, dev_var_in, gt0Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad_1_gt0, dev_var_in, gt0Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_yy(tid, grad2_1_1_gt0, dev_var_in, gt0Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad_2_gt0, dev_var_in, gt0Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_zz(tid, grad2_2_2_gt0, dev_var_in, gt0Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x(tid, grad_0_Gt1, dev_var_in, Gt1Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad_1_Gt1, dev_var_in, Gt1Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad_2_Gt1, dev_var_in, Gt1Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x(tid, grad_0_Gt2, dev_var_in, Gt2Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad_1_Gt2, dev_var_in, Gt2Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad_2_Gt2, dev_var_in, Gt2Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x(tid, grad_0_K, dev_var_in, KInt, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad_1_K, dev_var_in, KInt, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad_2_K, dev_var_in, KInt, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x(tid, grad_0_gt1, dev_var_in, gt1Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_xx(tid, grad2_0_0_gt1, dev_var_in, gt1Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad_1_gt1, dev_var_in, gt1Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_yy(tid, grad2_1_1_gt1, dev_var_in, gt1Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad_2_gt1, dev_var_in, gt1Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_zz(tid, grad2_2_2_gt1, dev_var_in, gt1Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x(tid, grad_0_gt2, dev_var_in, gt2Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_xx(tid, grad2_0_0_gt2, dev_var_in, gt2Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad_1_gt2, dev_var_in, gt2Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_yy(tid, grad2_1_1_gt2, dev_var_in, gt2Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad_2_gt2, dev_var_in, gt2Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_zz(tid, grad2_2_2_gt2, dev_var_in, gt2Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x(tid, grad_0_gt3, dev_var_in, gt3Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_xx(tid, grad2_0_0_gt3, dev_var_in, gt3Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad_1_gt3, dev_var_in, gt3Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_yy(tid, grad2_1_1_gt3, dev_var_in, gt3Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad_2_gt3, dev_var_in, gt3Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_zz(tid, grad2_2_2_gt3, dev_var_in, gt3Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x(tid, grad_0_gt4, dev_var_in, gt4Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_xx(tid, grad2_0_0_gt4, dev_var_in, gt4Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad_1_gt4, dev_var_in, gt4Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_yy(tid, grad2_1_1_gt4, dev_var_in, gt4Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad_2_gt4, dev_var_in, gt4Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_zz(tid, grad2_2_2_gt4, dev_var_in, gt4Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x(tid, grad_0_gt5, dev_var_in, gt5Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_xx(tid, grad2_0_0_gt5, dev_var_in, gt5Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad_1_gt5, dev_var_in, gt5Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_yy(tid, grad2_1_1_gt5, dev_var_in, gt5Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad_2_gt5, dev_var_in, gt5Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_zz(tid, grad2_2_2_gt5, dev_var_in, gt5Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x(tid, grad_0_At0, dev_var_in, At0Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad_1_At0, dev_var_in, At0Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad_2_At0, dev_var_in, At0Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x(tid, grad_0_At1, dev_var_in, At1Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad_1_At1, dev_var_in, At1Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad_2_At1, dev_var_in, At1Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x(tid, grad_0_At2, dev_var_in, At2Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad_1_At2, dev_var_in, At2Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad_2_At2, dev_var_in, At2Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x(tid, grad_0_At3, dev_var_in, At3Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad_1_At3, dev_var_in, At3Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad_2_At3, dev_var_in, At3Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x(tid, grad_0_At4, dev_var_in, At4Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad_1_At4, dev_var_in, At4Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad_2_At4, dev_var_in, At4Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x(tid, grad_0_At5, dev_var_in, At5Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad_1_At5, dev_var_in, At5Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad_2_At5, dev_var_in, At5Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
}

__global__ void calc_derivs2(
    double * dev_var_out, double * dev_var_in, double hx, double hy, double hz, 
    const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, 
    int bflag,
    #include "list_of_offset_para.h"
    ,
    #include "list_of_para.h"
){
    int tid = blockIdx.x*64 + threadIdx.x;

    calc_deriv42_y(tid, grad2_0_1_gt0, grad_0_gt0, 0, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad2_0_2_gt0, grad_0_gt0, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad2_1_2_gt0, grad_1_gt0, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad2_0_1_gt1, grad_0_gt1, 0, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad2_0_2_gt1, grad_0_gt1, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad2_1_2_gt1, grad_1_gt1, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad2_0_1_gt2, grad_0_gt2, 0, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad2_0_2_gt2, grad_0_gt2, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad2_1_2_gt2, grad_1_gt2, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad2_0_1_gt3, grad_0_gt3, 0, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad2_0_2_gt3, grad_0_gt3, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad2_1_2_gt3, grad_1_gt3, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad2_0_1_gt4, grad_0_gt4, 0, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad2_0_2_gt4, grad_0_gt4, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad2_1_2_gt4, grad_1_gt4, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad2_0_1_gt5, grad_0_gt5, 0, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad2_0_2_gt5, grad_0_gt5, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad2_1_2_gt5, grad_1_gt5, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad2_0_1_chi, grad_0_chi, 0, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad2_0_2_chi, grad_0_chi, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad2_1_2_chi, grad_1_chi, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad2_0_1_alpha, grad_0_alpha, 0, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad2_0_2_alpha, grad_0_alpha, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad2_1_2_alpha, grad_1_alpha, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad2_0_1_beta0, grad_0_beta0, 0, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad2_0_2_beta0, grad_0_beta0, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad2_1_2_beta0, grad_1_beta0, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad2_0_1_beta1, grad_0_beta1, 0, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad2_0_2_beta1, grad_0_beta1, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad2_1_2_beta1, grad_1_beta1, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y(tid, grad2_0_1_beta2, grad_0_beta2, 0, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad2_0_2_beta2, grad_0_beta2, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z(tid, grad2_1_2_beta2, grad_1_beta2, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_adv_x(tid, agrad_0_gt0, dev_var_in, gt0Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y(tid, agrad_1_gt0, dev_var_in, gt0Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_adv_x(tid, agrad_0_gt1, dev_var_in, gt1Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y(tid, agrad_1_gt1, dev_var_in, gt1Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_adv_x(tid, agrad_0_gt2, dev_var_in, gt2Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y(tid, agrad_1_gt2, dev_var_in, gt2Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_adv_x(tid, agrad_0_gt3, dev_var_in, gt3Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y(tid, agrad_1_gt3, dev_var_in, gt3Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_adv_x(tid, agrad_0_gt4, dev_var_in, gt4Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y(tid, agrad_1_gt4, dev_var_in, gt4Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_adv_x(tid, agrad_0_gt5, dev_var_in, gt5Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y(tid, agrad_1_gt5, dev_var_in, gt5Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);


    calc_deriv42_adv_x(tid, agrad_0_At0, dev_var_in, At0Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y(tid, agrad_1_At0, dev_var_in, At0Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_adv_x(tid, agrad_0_At1, dev_var_in, At1Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y(tid, agrad_1_At1, dev_var_in, At1Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_adv_x(tid, agrad_0_At2, dev_var_in, At2Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y(tid, agrad_1_At2, dev_var_in, At2Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_adv_x(tid, agrad_0_At3, dev_var_in, At3Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y(tid, agrad_1_At3, dev_var_in, At3Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_adv_x(tid, agrad_0_At4, dev_var_in, At4Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y(tid, agrad_1_At4, dev_var_in, At4Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_adv_x(tid, agrad_0_At5, dev_var_in, At5Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y(tid, agrad_1_At5, dev_var_in, At5Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_adv_x(tid, agrad_0_alpha, dev_var_in, alphaInt, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y(tid, agrad_1_alpha, dev_var_in, alphaInt, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_x(tid, agrad_0_beta0, dev_var_in, beta0Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_x(tid, agrad_0_beta1, dev_var_in, beta1Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_x(tid, agrad_0_beta2, dev_var_in, beta2Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_x(tid, agrad_0_chi, dev_var_in, chiInt, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_x(tid, agrad_0_Gt0, dev_var_in, Gt0Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_x(tid, agrad_0_Gt1, dev_var_in, Gt1Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_x(tid, agrad_0_Gt2, dev_var_in, Gt2Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_x(tid, agrad_0_K, dev_var_in, KInt, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_x(tid, agrad_0_B0, dev_var_in, B0Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_x(tid, agrad_0_B1, dev_var_in, B1Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_x(tid, agrad_0_B2, dev_var_in, B2Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y(tid, agrad_1_beta0, dev_var_in, beta0Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y(tid, agrad_1_beta1, dev_var_in, beta1Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y(tid, agrad_1_beta2, dev_var_in, beta2Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y(tid, agrad_1_chi, dev_var_in, chiInt, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y(tid, agrad_1_Gt0, dev_var_in, Gt0Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y(tid, agrad_1_Gt1, dev_var_in, Gt1Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y(tid, agrad_1_Gt2, dev_var_in, Gt2Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y(tid, agrad_1_K, dev_var_in, KInt, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y(tid, agrad_1_B0, dev_var_in, B0Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y(tid, agrad_1_B1, dev_var_in, B1Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y(tid, agrad_1_B2, dev_var_in, B2Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z(tid, agrad_2_gt0, dev_var_in, gt0Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z(tid, agrad_2_gt1, dev_var_in, gt1Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z(tid, agrad_2_gt2, dev_var_in, gt2Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z(tid, agrad_2_gt3, dev_var_in, gt3Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z(tid, agrad_2_gt4, dev_var_in, gt4Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z(tid, agrad_2_gt5, dev_var_in, gt5Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z(tid, agrad_2_At0, dev_var_in, At0Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z(tid, agrad_2_At1, dev_var_in, At1Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z(tid, agrad_2_At2, dev_var_in, At2Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z(tid, agrad_2_At3, dev_var_in, At3Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z(tid, agrad_2_At4, dev_var_in, At4Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z(tid, agrad_2_At5, dev_var_in, At5Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z(tid, agrad_2_alpha, dev_var_in, alphaInt, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z(tid, agrad_2_beta0, dev_var_in, beta0Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z(tid, agrad_2_beta1, dev_var_in, beta1Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z(tid, agrad_2_beta2, dev_var_in, beta2Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z(tid, agrad_2_chi, dev_var_in, chiInt, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z(tid, agrad_2_Gt0, dev_var_in, Gt0Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z(tid, agrad_2_Gt1, dev_var_in, Gt1Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z(tid, agrad_2_Gt2, dev_var_in, Gt2Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z(tid, agrad_2_K, dev_var_in, KInt, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z(tid, agrad_2_B0, dev_var_in, B0Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z(tid, agrad_2_B1, dev_var_in, B1Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z(tid, agrad_2_B2, dev_var_in, B2Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
}

void calc_deriv_wrapper(double * dev_var_out, double * dev_var_in, double hx, double hy, double hz, const unsigned int * host_sz, int bflag, cudaStream_t stream,
    #include "list_of_offset_para.h"
    ,
    #include "list_of_para.h"
){
    const int ib = 1;
    const int jb = 1;
    const int kb = 1;
    const int ie = host_sz[0] - 1;
    const int je = host_sz[1] - 1;
    const int ke = host_sz[2] - 1;
    const unsigned int host_sz_x = host_sz[0];
    const unsigned int host_sz_y = host_sz[1];
    const unsigned int host_sz_z = host_sz[2];
 
    int number_of_threads_required;
    int number_of_blocks;

    number_of_threads_required=ceil((ie-ib)*(je-jb)*(ke-kb));
    number_of_blocks=ceil(1.0*number_of_threads_required/64);
    calc_derivs1 <<< number_of_blocks, 64, 0, stream>>> (dev_var_out, dev_var_in, hx, hy, hz, host_sz_x, host_sz_y, host_sz_z, bflag,
        #include "list_of_offset_args.h"
        ,
        #include "list_of_args.h"
    );

    number_of_threads_required=ceil((ie-ib)*(je-jb)*(ke-kb));
    number_of_blocks=ceil(1.0*number_of_threads_required/64);
    calc_derivs2 <<< number_of_blocks, 64, 0, stream>>> (dev_var_out, dev_var_in, hx, hy, hz, host_sz_x, host_sz_y, host_sz_z, bflag,
        #include "list_of_offset_args.h"
        ,
        #include "list_of_args.h"
    );
    CHECK_ERROR(cudaGetLastError(), "deriv Kernel launch failed");
}

__device__ void calc_ko_deriv42_x(int id, double * output, double * dev_var_in, int u_offset, double dx, 
const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag)
{
    int nx = host_sz_x;
    int ny = host_sz_y;

    int i = id%(host_sz_x-6) + 3;
    int j = ((id/(host_sz_x-6))%(host_sz_y-6)) + 3;
    int k = (id/(host_sz_z-6)/(host_sz_x-6)) + 3;
    if (k>=host_sz_z-3) return;
    int pp = IDX(i, j, k);

    
    if(i==4) {
        int ib=3;
        output[IDX(3, j, k)] = (-1.0 / 64.0 / dx) *
                        (
                        -      dev_var_in[u_offset + IDX(ib+4,j,k)]
                        +  6.0*dev_var_in[u_offset + IDX(ib+3,j,k)]
                        - 15.0*dev_var_in[u_offset + IDX(ib+2,j,k)]
                        + 20.0*dev_var_in[u_offset + IDX(ib+1,j,k)]
                        - 15.0*dev_var_in[u_offset + IDX(ib,j,k)]
                        +  6.0*dev_var_in[u_offset + IDX(ib-1,j,k)]
                        -      dev_var_in[u_offset + IDX(ib-2,j,k)]
                        );
    }

    output[pp] = (-1.0 / 64.0 / dx) *
                            (
                            -      dev_var_in[u_offset + pp - 3]
                            +  6.0*dev_var_in[u_offset + pp - 2]
                            - 15.0*dev_var_in[u_offset + pp - 1]
                            + 20.0*dev_var_in[u_offset + pp ]
                            - 15.0*dev_var_in[u_offset + pp + 1]
                            +  6.0*dev_var_in[u_offset + pp + 2]
                            -      dev_var_in[u_offset + pp + 3]
                            );

    if(i==5) {
        int ie = nx-3;
        output[IDX(ie-1, j, k)] = (-1.0 / 64.0 / dx) *
                        (
                        -      dev_var_in[u_offset + IDX(ie+1,j,k)]
                        +  6.0*dev_var_in[u_offset + IDX(ie,j,k)]
                        - 15.0*dev_var_in[u_offset + IDX(ie-1,j,k)]
                        + 20.0*dev_var_in[u_offset + IDX(ie-2,j,k)]
                        - 15.0*dev_var_in[u_offset + IDX(ie-3,j,k)]
                        +  6.0*dev_var_in[u_offset + IDX(ie-4,j,k)]
                        -      dev_var_in[u_offset + IDX(ie-5,j,k)]
                        );
    }
}

__device__ void calc_ko_deriv42_y(int id, double * output, double * dev_var_in, int u_offset, double dy,
const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag)
{
    
    int nx = host_sz_x;
    int ny = host_sz_y;

    int i = id%(host_sz_x-6) + 3;
    int j = ((id/(host_sz_x-6))%(host_sz_y-6)) + 3;
    int k = (id/(host_sz_z-6)/(host_sz_x-6)) + 3;
    if (k>=host_sz_z-3) return;
    int pp = IDX(i, j, k);

    if(j==4) {
        int jb=3;
        output[IDX(i,jb,k)] = (-1.0 / 64.0 / dy) *
                    (
                        -      dev_var_in[u_offset + IDX(i,jb+4,k)]
                        +  6.0*dev_var_in[u_offset + IDX(i,jb+3,k)]
                        - 15.0*dev_var_in[u_offset + IDX(i,jb+2,k)]
                        + 20.0*dev_var_in[u_offset + IDX(i,jb+1,k)]
                        - 15.0*dev_var_in[u_offset + IDX(i,jb,k)]
                        +  6.0*dev_var_in[u_offset + IDX(i,jb-1,k)]
                        -      dev_var_in[u_offset + IDX(i,jb-2,k)]
                        );
        }
    output[pp] = (-1.0 / 64.0 / dy) *
                    (
                        -      dev_var_in[u_offset + pp-3*nx]
                        +  6.0*dev_var_in[u_offset + pp-2*nx]
                        - 15.0*dev_var_in[u_offset + pp-nx]
                        + 20.0*dev_var_in[u_offset + pp]
                        - 15.0*dev_var_in[u_offset + pp+nx]
                        +  6.0*dev_var_in[u_offset + pp+2*nx]
                        -      dev_var_in[u_offset + pp+3*nx]
                        );

    if(j==5) {
        int je = ny - 3;
        output[IDX(i,je-1,k)] = (-1.0 / 64.0 / dy) *
                (
                    -      dev_var_in[u_offset + IDX(i,je+1,k)]
                    +  6.0*dev_var_in[u_offset + IDX(i,je,k)]
                    - 15.0*dev_var_in[u_offset + IDX(i,je-1,k)]
                    + 20.0*dev_var_in[u_offset + IDX(i,je-2,k)]
                    - 15.0*dev_var_in[u_offset + IDX(i,je-3,k)]
                    +  6.0*dev_var_in[u_offset + IDX(i,je-4,k)]
                    -      dev_var_in[u_offset + IDX(i,je-5,k)]
                    );                   
    }
}

__device__ void calc_ko_deriv42_z(int id, double * output, double * dev_var_in, int u_offset, double dz, 
const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag)
{  
    int nx = host_sz_x;
    int ny = host_sz_y;
    int n = nx * ny;

    int i = id%(host_sz_x-6) + 3;
    int j = ((id/(host_sz_x-6))%(host_sz_y-6)) + 3;
    int k = (id/(host_sz_z-6)/(host_sz_x-6)) + 3;
    if (k>=host_sz_z-3) return;
    int pp = IDX(i, j, k);
    
    if(k==4) {
        int kb=3;
        output[IDX(i,j,kb)] = (-1.0 / 64.0 / dz) *
                    (
                        -      dev_var_in[u_offset + IDX(i,j,kb+4)]
                        +  6.0*dev_var_in[u_offset + IDX(i,j,kb+3)]
                        - 15.0*dev_var_in[u_offset + IDX(i,j,kb+2)]
                        + 20.0*dev_var_in[u_offset + IDX(i,j,kb+1)]
                        - 15.0*dev_var_in[u_offset + IDX(i,j,kb)]
                        +  6.0*dev_var_in[u_offset + IDX(i,j,kb-1)]
                        -      dev_var_in[u_offset + IDX(i,j,kb-2)]
                        );
        }

    output[pp] = (-1.0 / 64.0 / dz) *
                (
                    -      dev_var_in[u_offset + pp-3*n]
                    +  6.0*dev_var_in[u_offset + pp-2*n]
                    - 15.0*dev_var_in[u_offset + pp-n]
                    + 20.0*dev_var_in[u_offset + pp]
                    - 15.0*dev_var_in[u_offset + pp+n]
                    +  6.0*dev_var_in[u_offset + pp+2*n]
                    -      dev_var_in[u_offset + pp+3*n]
                    );
    if(k==5) {
        int ke = host_sz_z - 3;
        output[IDX(i,j,ke-1)] = (-1.0 / 64.0 / dz) *
        (
            -      dev_var_in[u_offset + IDX(i,j,ke+1)]
            +  6.0*dev_var_in[u_offset + IDX(i,j,ke)]
            - 15.0*dev_var_in[u_offset + IDX(i,j,ke-1)]
            + 20.0*dev_var_in[u_offset + IDX(i,j,ke-2)]
            - 15.0*dev_var_in[u_offset + IDX(i,j,ke-3)]
            +  6.0*dev_var_in[u_offset + IDX(i,j,ke-4)]
            -      dev_var_in[u_offset + IDX(i,j,ke-5)]
            );               
    }
        
}


__global__ void calc_ko_derivs(
    double * dev_var_out, double * dev_var_in, double hx, double hy, double hz, 
    const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, 
    int bflag,
    #include "list_of_offset_para.h"
    ,
    #include "list_of_para.h"
){
    int tid = blockIdx.x*64 + threadIdx.x;

    calc_ko_deriv42_x(tid, grad_0_gt0, dev_var_in, gt0Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y(tid, grad_1_gt0, dev_var_in, gt0Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z(tid, grad_2_gt0, dev_var_in, gt0Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x(tid, grad_0_gt1, dev_var_in, gt1Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y(tid, grad_1_gt1, dev_var_in, gt1Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z(tid, grad_2_gt1, dev_var_in, gt1Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x(tid, grad_0_gt2, dev_var_in, gt2Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y(tid, grad_1_gt2, dev_var_in, gt2Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z(tid, grad_2_gt2, dev_var_in, gt2Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x(tid, grad_0_gt3, dev_var_in, gt3Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y(tid, grad_1_gt3, dev_var_in, gt3Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z(tid, grad_2_gt3, dev_var_in, gt3Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x(tid, grad_0_gt4, dev_var_in, gt4Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y(tid, grad_1_gt4, dev_var_in, gt4Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z(tid, grad_2_gt4, dev_var_in, gt4Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x(tid, grad_0_gt5, dev_var_in, gt5Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y(tid, grad_1_gt5, dev_var_in, gt5Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z(tid, grad_2_gt5, dev_var_in, gt5Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x(tid, grad_0_At0, dev_var_in, At0Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y(tid, grad_1_At0, dev_var_in, At0Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z(tid, grad_2_At0, dev_var_in, At0Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x(tid, grad_0_At1, dev_var_in, At1Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y(tid, grad_1_At1, dev_var_in, At1Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z(tid, grad_2_At1, dev_var_in, At1Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x(tid, grad_0_At2, dev_var_in, At2Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y(tid, grad_1_At2, dev_var_in, At2Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z(tid, grad_2_At2, dev_var_in, At2Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x(tid, grad_0_At3, dev_var_in, At3Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y(tid, grad_1_At3, dev_var_in, At3Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z(tid, grad_2_At3, dev_var_in, At3Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x(tid, grad_0_At4, dev_var_in, At4Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y(tid, grad_1_At4, dev_var_in, At4Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z(tid, grad_2_At4, dev_var_in, At4Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x(tid, grad_0_At5, dev_var_in, At5Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y(tid, grad_1_At5, dev_var_in, At5Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z(tid, grad_2_At5, dev_var_in, At5Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x(tid, grad_0_alpha, dev_var_in, alphaInt, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y(tid, grad_1_alpha, dev_var_in, alphaInt, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z(tid, grad_2_alpha, dev_var_in, alphaInt, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x(tid, grad_0_beta0, dev_var_in, beta0Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y(tid, grad_1_beta0, dev_var_in, beta0Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z(tid, grad_2_beta0, dev_var_in, beta0Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x(tid, grad_0_beta1, dev_var_in, beta1Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y(tid, grad_1_beta1, dev_var_in, beta1Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z(tid, grad_2_beta1, dev_var_in, beta1Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x(tid, grad_0_beta2, dev_var_in, beta2Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y(tid, grad_1_beta2, dev_var_in, beta2Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z(tid, grad_2_beta2, dev_var_in, beta2Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x(tid, grad_0_chi, dev_var_in, chiInt, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y(tid, grad_1_chi, dev_var_in, chiInt, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z(tid, grad_2_chi, dev_var_in, chiInt, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x(tid, grad_0_Gt0, dev_var_in, Gt0Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y(tid, grad_1_Gt0, dev_var_in, Gt0Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z(tid, grad_2_Gt0, dev_var_in, Gt0Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x(tid, grad_0_Gt1, dev_var_in, Gt1Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y(tid, grad_1_Gt1, dev_var_in, Gt1Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z(tid, grad_2_Gt1, dev_var_in, Gt1Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x(tid, grad_0_Gt2, dev_var_in, Gt2Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y(tid, grad_1_Gt2, dev_var_in, Gt2Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z(tid, grad_2_Gt2, dev_var_in, Gt2Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x(tid, grad_0_K, dev_var_in, KInt, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y(tid, grad_1_K, dev_var_in, KInt, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z(tid, grad_2_K, dev_var_in, KInt, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x(tid, grad_0_B0, dev_var_in, B0Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y(tid, grad_1_B0, dev_var_in, B0Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z(tid, grad_2_B0, dev_var_in, B0Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x(tid, grad_0_B1, dev_var_in, B1Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y(tid, grad_1_B1, dev_var_in, B1Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z(tid, grad_2_B1, dev_var_in, B1Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x(tid, grad_0_B2, dev_var_in, B2Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y(tid, grad_1_B2, dev_var_in, B2Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z(tid, grad_2_B2, dev_var_in, B2Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);


}

void calc_ko_deriv_wrapper(double * dev_var_out, double * dev_var_in, double hx, double hy, double hz, const unsigned int * host_sz, int bflag, cudaStream_t stream,
    #include "list_of_offset_para.h"
    ,
    #include "list_of_para.h"
    )
{
    const int ib = 1;
    const int jb = 1;
    const int kb = 1;
    const int ie = host_sz[0] - 1;
    const int je = host_sz[1] - 1;
    const int ke = host_sz[2] - 1;
    const unsigned int host_sz_x = host_sz[0];
    const unsigned int host_sz_y = host_sz[1];
    const unsigned int host_sz_z = host_sz[2];

    int number_of_threads_required;
    int number_of_blocks;

    number_of_threads_required=ceil((ie-ib)*(je-jb)*(ke-kb));
    number_of_blocks=ceil(1.0*number_of_threads_required/64);

    calc_ko_derivs <<< number_of_blocks, 64, 0, stream>>> (dev_var_out, dev_var_in, hx, hy, hz, host_sz_x, host_sz_y, host_sz_z, bflag,
        #include "list_of_offset_args.h"
        ,
        #include "list_of_args.h"
    );
    CHECK_ERROR(cudaGetLastError(), "ko deriv Kernel launch failed");
}

// device methods with bflag

__device__ void calc_deriv42_x_bflag(int id, double * output, double * dev_var_in, const int u_offset, double dx, 
    const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag){

    int nx = host_sz_x; 
    int ny = host_sz_y; 

    const double idx = 1.0/dx;
    const double idx_by_2 = 0.50 * idx;
    const double idx_by_12 = idx / 12.0;

    int i = id%(host_sz_x-6) + 3;
    int j = ((id/(host_sz_x-6))%(host_sz_y-2)) + 1;
    int k = (id/(host_sz_z-2)/(host_sz_x-6)) + 1; 
    if (k>=host_sz_z-1) return;
    int pp = IDX(i, j, k);

    output[pp] = (dev_var_in[(u_offset) + pp - 2] - 8.0*dev_var_in[(u_offset)
                    + pp - 1] + 8.0*dev_var_in[(u_offset) + pp + 1] 
                    - dev_var_in[(u_offset) + pp + 2] )*idx_by_12;

    if ((bflag & (1u<<OCT_DIR_LEFT)) && i==3)  {
        int pp3 = IDX(3, j, k);
        int pp4 = IDX(4, j, k);
        int pp5 = IDX(5, j, k);
        output[pp3] = ((-3)*dev_var_in[(u_offset) + pp3] + 4*dev_var_in[(u_offset) + pp4] - dev_var_in[(u_offset) + pp5]) * idx_by_2;
        output[pp4] = (dev_var_in[(u_offset) + pp5] - dev_var_in[(u_offset) + pp3]) * idx_by_2;
    }

    if ((bflag & (1u<<OCT_DIR_RIGHT)) && i==4)  {
        int pp2 = IDX(nx-5, j, k);
        int pp3 = IDX(nx-6, j, k);
        int pp1 = IDX(nx-4,j,k);
        output[pp2] = (dev_var_in[(u_offset) + pp1] - dev_var_in[(u_offset) + pp3]) * idx_by_2;
        output[pp1] = (dev_var_in[(u_offset) + pp3]- 4.0 * dev_var_in[(u_offset) + pp2]+ 3.0 * dev_var_in[(u_offset) + pp1]) * idx_by_2;
    }
 }

 __device__ void calc_deriv42_y_bflag(int id, double* output, double * dev_var_in, const int u_offset, double dy, 
    const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag){

    int nx = host_sz_x; 
    int ny = host_sz_y; 

    const double idy = 1.0/dy;
    const double idy_by_2 = 0.50 * idy;
    const double idy_by_12 = idy / 12.0;
    
    int i = id%(host_sz_x-6) + 3;
    int j = ((id/(host_sz_x-6))%(host_sz_y-6)) + 3;
    int k = (id/(host_sz_z-6)/(host_sz_x-6)) + 1;
    if (k>=host_sz_z-1) return;
    int pp = IDX(i, j, k);
    
    output[pp] = (dev_var_in[u_offset + pp - 2*nx] 
                - 8.0*dev_var_in[u_offset + pp - nx] 
                + 8.0*dev_var_in[u_offset + pp + nx] 
                - dev_var_in[u_offset + pp + 2*nx] )*idy_by_12;
    
    if ((bflag & (1u<<OCT_DIR_DOWN)) && j==3)  {
        int pp3 = IDX(i, 3, k);
        int pp4 = IDX(i, 4, k);
        int pp5 = IDX(i, 5, k);
        output[pp3] = ((-3)*dev_var_in[(u_offset) + pp3] +  4*dev_var_in[(u_offset) + pp4] - dev_var_in[(u_offset) + pp5]) * idy_by_2;
        output[pp4] = (dev_var_in[(u_offset) + pp5] - dev_var_in[(u_offset) + pp3]) * idy_by_2;
    }

    if ((bflag & (1u<<OCT_DIR_UP)) && j==4)  {
        int pp2 = IDX(i, ny-5, k); // IDX(i,je-2,k)
        int pp3 = IDX(i, ny-6, k); // IDX(i,je-3,k)
        int pp1 = IDX(i, ny-4, k); // IDX(i,je-1,k)
        output[pp2] = (dev_var_in[(u_offset) + pp1] - dev_var_in[(u_offset) + pp3]) * idy_by_2;
        output[pp1] = (dev_var_in[(u_offset) + pp3]- 4.0 * dev_var_in[(u_offset) + pp2] + 3.0 * dev_var_in[(u_offset) + pp1]) * idy_by_2;
    }
    
}

__device__ void calc_deriv42_z_bflag(int id, double* output, double * dev_var_in, const int u_offset, double dz, 
    const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag){
    int nx = host_sz_x; 
    int ny = host_sz_y; 
    int n = nx * ny;

    const double idz = 1.0/dz;
    const double idz_by_2 = 0.50 * idz;
    const double idz_by_12 = idz / 12.0;

    int i = id%(host_sz_x-6) + 3;
    int j = ((id/(host_sz_x-6))%(host_sz_y-6)) + 3;
    int k = (id/(host_sz_z-6)/(host_sz_x-6)) + 3;
    if (k>=host_sz_z-3) return;
    int pp = IDX(i, j, k);
    output[pp] = (dev_var_in[(u_offset) + pp - 2*n] - 8.0*dev_var_in[(u_offset) + pp - n] 
                    + 8.0*dev_var_in[(u_offset) + pp + n] - dev_var_in[(u_offset) + pp + 2*n]) 
                    * idz_by_12;
            
    if ((bflag & (1u<<OCT_DIR_BACK)) && k==3)  {
        int pp3 = IDX(i, j, 3); 
        int pp4 = IDX(i, j, 4); 
        int pp5 = IDX(i, j, 5);
        output[pp3] = ((-3)*dev_var_in[(u_offset) + pp3] + 4*dev_var_in[(u_offset) + pp4] - dev_var_in[(u_offset) + pp5]) * idz_by_2;
        output[pp4] = (dev_var_in[(u_offset) + pp5] - dev_var_in[(u_offset) + pp3]) * idz_by_2;
    }
                
    if ((bflag & (1u<<OCT_DIR_FRONT)) && k==4)  {
        int pp2 = IDX(i, j, host_sz_z-5); 
        int pp3 = IDX(i, j, host_sz_z-6); 
        int pp1 = IDX(i, j, host_sz_z-4);
        output[pp2] = (dev_var_in[(u_offset) + pp1] - dev_var_in[(u_offset) + pp3]) * idz_by_2;
        output[pp1] = (dev_var_in[(u_offset) + pp3]- 4.0 * dev_var_in[(u_offset) + pp2] + 3.0 * dev_var_in[(u_offset) + pp1]) * idz_by_2;
    }
}

__device__ void calc_deriv42_xx_bflag(int id, double* output, double * dev_var_in, const int u_offset, double dx, 
    const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag){
    int nx = host_sz_x; 
    int ny = host_sz_y; 

    const double idx_sqrd = 1.0/(dx*dx);
    const double idx_sqrd_by_12 = idx_sqrd / 12.0;
    
    int i = id%(host_sz_x-6) + 3;
    int j = ((id/(host_sz_x-6))%(host_sz_y-6)) + 3;
    int k = (id/(host_sz_z-6)/(host_sz_x-6)) + 3;
    if (k>=host_sz_z-3) return;
    int pp = IDX(i, j, k);
    output[pp] = ((-1)*dev_var_in[(u_offset) + pp - 2] 
                + 16.0*dev_var_in[(u_offset) + pp - 1] 
                - 30.0*dev_var_in[(u_offset) + pp] 
                + 16.0*dev_var_in[(u_offset) + pp + 1] 
                - dev_var_in[(u_offset) + pp + 2] 
            )*idx_sqrd_by_12;

    if ((bflag & (1u<<OCT_DIR_LEFT)) && i==3)  {
        int pp3 = IDX(3, j, k); 
        int pp4 = IDX(4, j, k); 
        int pp5 = IDX(5, j, k); 
        int pp6 = IDX(6, j, k); 
    
        output[pp3] = (
                2.0     *       dev_var_in[(u_offset) + pp3] 
            -   5.0     *       dev_var_in[(u_offset) + pp4] 
            +   4.0     *       dev_var_in[(u_offset) + pp5] 
            -                   dev_var_in[(u_offset) + pp6]
            )*idx_sqrd;
    
        output[pp4] = (
                            dev_var_in[(u_offset) + pp3]
            -   2.0     *   dev_var_in[(u_offset) + pp4]
            +               dev_var_in[(u_offset) + pp5]
        )*idx_sqrd;
    }
                        
    if ((bflag & (1u<<OCT_DIR_RIGHT)) && i==4)  {
        int pp1 = IDX(host_sz_x - 4, j, k); // IDX(ie-1,j,k)
        int pp2 = IDX(host_sz_x - 5, j, k); // IDX(ie-2,j,k)
        int pp3 = IDX(host_sz_x - 6, j, k); // IDX(ie-3,j,k)
        int pp4 = IDX(host_sz_x - 7, j, k); // IDX(ie-4,j,k)

        output[pp2] = (
                                dev_var_in[(u_offset) + pp3] 
                -   2.0     *   dev_var_in[(u_offset) + pp2] 
                +               dev_var_in[(u_offset) + pp1] 
                )*idx_sqrd;


        output[pp1] = (
            -   1.0 *   dev_var_in[(u_offset) + pp4] 
            +   4.0 *   dev_var_in[(u_offset) + pp3] 
            -   5.0 *   dev_var_in[(u_offset) + pp2] 
            +   2.0 *   dev_var_in[(u_offset) + pp1]
            )*idx_sqrd;
    }

}

__device__ void calc_deriv42_yy_bflag(int id, double* output, double * dev_var_in, const int u_offset, double dy, 
    const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag){
    int nx = host_sz_x; 
    int ny = host_sz_y; 

    const double idy_sqrd = 1.0/(dy*dy);
    const double idy_sqrd_by_12 = idy_sqrd / 12.0;

    int i = id%(host_sz_x-6) + 3;
    int j = ((id/(host_sz_x-6))%(host_sz_y-6)) + 3;
    int k = (id/(host_sz_z-6)/(host_sz_x-6)) + 3;
    if (k>=host_sz_z-3) return;
    int pp = IDX(i, j, k);
    output[pp] = ((-1)*dev_var_in[(u_offset) + pp - 2*nx] 
                + 16.0*dev_var_in[(u_offset) + pp - nx] 
                - 30.0*dev_var_in[(u_offset) + pp] 
                + 16.0*dev_var_in[(u_offset) + pp + nx] 
                - dev_var_in[(u_offset) + pp + 2*nx] 
            )*idy_sqrd_by_12;

    if ((bflag & (1u<<OCT_DIR_DOWN)) && j==3)  {
        int pp3 = IDX(i, 3, k); 
        int pp4 = IDX(i, 4, k); 
        int pp5 = IDX(i, 5, k); 
        int pp6 = IDX(i, 6, k); 
    
        output[pp3] = (
                2.0 *   dev_var_in[(u_offset) + pp3] 
            -   5.0 *   dev_var_in[(u_offset) + pp4] 
            +   4.0 *   dev_var_in[(u_offset) + pp5] 
            -           dev_var_in[(u_offset) + pp6]
            ) * idy_sqrd;
    
        output[pp4] = (
                        dev_var_in[(u_offset) + pp3]
            -   2.0 *   dev_var_in[(u_offset) + pp4]
            +           dev_var_in[(u_offset) + pp5]
        ) * idy_sqrd;
    }
                            
    if ((bflag & (1u<<OCT_DIR_UP)) && j==4)  {
        int pp1 = IDX(i, host_sz_y - 4, k); 
        int pp2 = IDX(i, host_sz_y - 5, k); 
        int pp3 = IDX(i, host_sz_y - 6, k); 
        int pp4 = IDX(i, host_sz_y - 7, k); 
    
        output[pp2] = (
                        dev_var_in[(u_offset) + pp3] 
            -   2.0 *   dev_var_in[(u_offset) + pp2] 
            +           dev_var_in[(u_offset) + pp1] 
            ) * idy_sqrd;
    
    
        output[pp1] = (
            -   1.0 *   dev_var_in[(u_offset) + pp4] 
            +   4.0 *   dev_var_in[(u_offset) + pp3] 
            -   5.0 *   dev_var_in[(u_offset) + pp2] 
            +   2.0 *   dev_var_in[(u_offset) + pp1]
            ) * idy_sqrd;
    
    }
            
}

__device__ void calc_deriv42_zz_bflag(int id, double* output, double * dev_var_in, const int u_offset, double dz, 
    const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag){
    int nx = host_sz_x; 
    int ny = host_sz_y; 
    int n = nx * ny;

    const double idz_sqrd = 1.0/(dz*dz);
    const double idz_sqrd_by_12 = idz_sqrd / 12.0;

    int i = id%(host_sz_x-6) + 3;
    int j = ((id/(host_sz_x-6))%(host_sz_y-6)) + 3;
    int k = (id/(host_sz_z-6)/(host_sz_x-6)) + 3;
    if (k>=host_sz_z-3) return;
    int pp = IDX(i, j, k);
    output[pp] = ((-1)*dev_var_in[(u_offset) + pp - 2*n] 
                + 16.0*dev_var_in[(u_offset) + pp - n] 
                - 30.0*dev_var_in[(u_offset) + pp] 
                + 16.0*dev_var_in[(u_offset) + pp + n] 
                - dev_var_in[(u_offset) + pp + 2*n] 
            )*idz_sqrd_by_12;
    
    if ((bflag & (1u<<OCT_DIR_BACK)) && k==3)  {
        int pp3 = IDX(i, j, 3); 
        int pp4 = IDX(i, j, 4); 
        int pp5 = IDX(i, j, 5); 
        int pp6 = IDX(i, j, 6); 
    
        output[pp3] = (
                2.0 *   dev_var_in[(u_offset) + pp3] 
            -   5.0 *   dev_var_in[(u_offset) + pp4] 
            +   4.0 *   dev_var_in[(u_offset) + pp5] 
            -           dev_var_in[(u_offset) + pp6]
            ) * idz_sqrd;
    
        output[pp4] = (
                        dev_var_in[(u_offset) + pp3]
            -   2.0 *   dev_var_in[(u_offset) + pp4]
            +           dev_var_in[(u_offset) + pp5]
        ) * idz_sqrd;
    }
                                        
    if ((bflag & (1u<<OCT_DIR_FRONT)) && k==4)  {
        int pp1 = IDX(i, j, host_sz_z - 4); 
        int pp2 = IDX(i, j, host_sz_z - 5); 
        int pp3 = IDX(i, j, host_sz_z - 6); 
        int pp4 = IDX(i, j, host_sz_z - 7); 

        output[pp2] = (
                            dev_var_in[(u_offset) + pp3] 
                -   2.0 *   dev_var_in[(u_offset) + pp2] 
                +           dev_var_in[(u_offset) + pp1] 
                ) * idz_sqrd;


        output[pp1] = (
            -   1.0 *   dev_var_in[(u_offset) + pp4] 
            +   4.0 *   dev_var_in[(u_offset) + pp3] 
            -   5.0 *   dev_var_in[(u_offset) + pp2] 
            +   2.0 *   dev_var_in[(u_offset) + pp1]
            ) * idz_sqrd;
    }
}

__device__ void calc_deriv42_adv_x_bflag(int id, double * output, double * dev_var_in, int u_offset, double dx, int betax,
const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag) {
    int nx = host_sz_x;
    int ny = host_sz_y;

    const double idx = 1.0/dx;
    const double idx_by_2 = 0.50 * idx;
    const double idx_by_12 = idx / 12.0;

    int i = id%(host_sz_x-6) + 3;
    int j = ((id/(host_sz_x-6))%(host_sz_y-6)) + 3;
    int k = (id/(host_sz_z-6)/(host_sz_x-6)) + 3;
    if (k>=host_sz_z-3) return;
    int pp = IDX(i, j, k);

    if (dev_var_in[betax + pp] > 0.0 ) {
        output[pp] = ( -  3.0 * dev_var_in[u_offset + pp - 1]
                    - 10.0 * dev_var_in[u_offset + pp]
                    + 18.0 * dev_var_in[u_offset + pp + 1]
                    -  6.0 * dev_var_in[u_offset + pp + 2]
                    +        dev_var_in[u_offset + pp + 3]
                ) * idx_by_12;
    }
    else {
        output[pp] = ( -        dev_var_in[u_offset + pp - 3]
                    +  6.0 * dev_var_in[u_offset + pp - 2]
                    - 18.0 * dev_var_in[u_offset + pp - 1]
                    + 10.0 * dev_var_in[u_offset + pp]
                    +  3.0 * dev_var_in[u_offset + pp +1]
                ) * idx_by_12;
    }

    if ((bflag & (1u<<OCT_DIR_LEFT)) && (i == 3)) {
        output[IDX(3,j,k)] = ( -  3.0 * dev_var_in[u_offset + IDX(3,j,k)]
                +  4.0 * dev_var_in[u_offset + IDX(4,j,k)]
                -        dev_var_in[u_offset + IDX(5,j,k)]
                ) * idx_by_2;

        if (dev_var_in[betax + IDX(4,j,k)] > 0.0) {
            output[IDX(4,j,k)] = ( -  3.0 * dev_var_in[u_offset + IDX(4,j,k)]
                            +  4.0 * dev_var_in[u_offset + IDX(5,j,k)]
                            -        dev_var_in[u_offset + IDX(6,j,k)]
                        ) * idx_by_2;
        } else {
            output[IDX(4,j,k)] = ( -         dev_var_in[u_offset + IDX(3,j,k)]
                            +        dev_var_in[u_offset + IDX(5,j,k)]
                        ) * idx_by_2;
        }

        if (dev_var_in[betax + IDX(5,j,k)] > 0.0 ) {
            output[IDX(5,j,k)] = (-  3.0 * dev_var_in[u_offset + IDX(4,j,k)]
                        - 10.0 * dev_var_in[u_offset + IDX(5,j,k)]
                        + 18.0 * dev_var_in[u_offset + IDX(6,j,k)]
                        -  6.0 * dev_var_in[u_offset + IDX(7,j,k)]
                        +        dev_var_in[u_offset + IDX(8,j,k)]
                        ) * idx_by_12;
        } else {
            output[IDX(5,j,k)] = (           dev_var_in[u_offset + IDX(3,j,k)]
                            -  4.0 * dev_var_in[u_offset + IDX(4,j,k)]
                            +  3.0 * dev_var_in[u_offset + IDX(5,j,k)]
                        ) * idx_by_2;
        }
    }

    if ((bflag & (1u<<OCT_DIR_RIGHT)) && (i == 4)) {
        const int ie = nx - 3;
        if ( dev_var_in[betax + IDX(ie-3,j,k)] < 0.0 ) {
            output[IDX(ie-3,j,k)] = (  - 3.0 * dev_var_in[u_offset + IDX(ie-3,j,k)]
                                    + 4.0 * dev_var_in[u_offset + IDX(ie-2,j,k)]
                                    -       dev_var_in[u_offset + IDX(ie-1,j,k)]
                                ) * idx_by_2;
        } else {
            output[IDX(ie-3,j,k)] = ( -   dev_var_in[u_offset + IDX(ie-6,j,k)]
                            +  6.0 * dev_var_in[u_offset + IDX(ie-5,j,k)]
                            - 18.0 * dev_var_in[u_offset + IDX(ie-4,j,k)]
                            + 10.0 * dev_var_in[u_offset + IDX(ie-3  ,j,k)]
                            +  3.0 * dev_var_in[u_offset + IDX(ie-2,j,k)]
                            ) * idx_by_12;
        }

        if (dev_var_in[betax + IDX(ie-2,j,k)] > 0.0 ) {
            output[IDX(ie-2,j,k)] = (  -  dev_var_in[u_offset + IDX(ie-3,j,k)]
                                    +  dev_var_in[u_offset + IDX(ie-1,j,k)]
                                ) * idx_by_2;
        } else {
            output[IDX(ie-2,j,k)] = (     dev_var_in[u_offset + IDX(ie-4,j,k)]
                            - 4.0 * dev_var_in[u_offset + IDX(ie-3,j,k)]
                            + 3.0 * dev_var_in[u_offset + IDX(ie-2,j,k)]
                                ) * idx_by_2;
        }

        output[IDX(ie-1,j,k)] = (          dev_var_in[u_offset + IDX(ie-3,j,k)]
                                - 4.0 * dev_var_in[u_offset + IDX(ie-2,j,k)]
                                + 3.0 * dev_var_in[u_offset + IDX(ie-1,j,k)]
                            ) * idx_by_2;
    }
}

__device__ void calc_deriv42_adv_y_bflag(int id, double * output, double * dev_var_in, int u_offset, double dy, int betay,
const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag) {
    
    int nx = host_sz_x;
    int ny = host_sz_y;

    const double idy = 1.0/dy;
    const double idy_by_2 = 0.50 * idy;
    const double idy_by_12 = idy / 12.0;

    int i = id%(host_sz_x-6) + 3;
    int j = ((id/(host_sz_x-6))%(host_sz_y-6)) + 3;
    int k = (id/(host_sz_z-6)/(host_sz_x-6)) + 3;
    if (k>=host_sz_z-3) return;
    int pp = IDX(i, j, k);

    if (dev_var_in[betay + pp] > 0.0 ) {
            output[pp] = ( -  3.0 * dev_var_in[u_offset + pp - nx]
                        - 10.0 * dev_var_in[u_offset + pp]
                        + 18.0 * dev_var_in[u_offset + pp + nx]
                        -  6.0 * dev_var_in[u_offset + pp + 2*nx]
                        +        dev_var_in[u_offset + pp + 3*nx]
                    ) * idy_by_12;
    }
    else {
        output[pp] = ( -        dev_var_in[u_offset + pp - 3*nx]
                    +  6.0 * dev_var_in[u_offset + pp - 2*nx]
                    - 18.0 * dev_var_in[u_offset + pp - nx]
                    + 10.0 * dev_var_in[u_offset + pp]
                    +  3.0 * dev_var_in[u_offset + pp +nx]
                    ) * idy_by_12;
                
    }

    if ((bflag & (1u<<OCT_DIR_DOWN)) && (j == 3)) {
            
        output[IDX(i,3,k)] = ( -  3.0 * dev_var_in[u_offset + IDX(i,3,k)]
                +  4.0 * dev_var_in[u_offset + IDX(i,4,k)]
                -        dev_var_in[u_offset + IDX(i,5,k)]
                ) * idy_by_2;
                
        if (dev_var_in[betay + IDX(i,4,k)] > 0.0) {
            output[IDX(i,4,k)] = ( -  3.0 * dev_var_in[u_offset + IDX(i,4,k)]
                            +  4.0 * dev_var_in[u_offset + IDX(i,5,k)]
                            -        dev_var_in[u_offset + IDX(i,6,k)]
                        ) * idy_by_2;

        }
        else {
            output[IDX(i,4,k)] = ( -         dev_var_in[u_offset + IDX(i,3,k)]
                            +        dev_var_in[u_offset + IDX(i,5,k)]
                        ) * idy_by_2;
                        
        }

        if (dev_var_in[betay + IDX(i,5,k)] > 0.0 ) {
            output[IDX(i,5,k)] = (-  3.0 * dev_var_in[u_offset + IDX(i,4,k)]
                        - 10.0 * dev_var_in[u_offset + IDX(i,5,k)]
                        + 18.0 * dev_var_in[u_offset + IDX(i,6,k)]
                        -  6.0 * dev_var_in[u_offset + IDX(i,7,k)]
                        +        dev_var_in[u_offset + IDX(i,8,k)]
                        ) * idy_by_12;
        }
        else {
            output[IDX(i,5,k)] = (           dev_var_in[u_offset + IDX(i,3,k)]
                            -  4.0 * dev_var_in[u_offset + IDX(i,4,k)]
                            +  3.0 * dev_var_in[u_offset + IDX(i,5,k)]
                        ) * idy_by_2;
        }
    }

    if ((bflag & (1u<<OCT_DIR_UP)) && (j == 4)) {
        
        const int je = host_sz_y - 3;
        
        if ( dev_var_in[betay + IDX(i,je-3,k)] < 0.0 ) {
            output[IDX(i,je-3,k)] = (  - 3.0 * dev_var_in[u_offset + IDX(i,je-3,k)]
                                    + 4.0 * dev_var_in[u_offset + IDX(i,je-2,k)]
                                    -       dev_var_in[u_offset + IDX(i,je-1,k)]
                                    ) * idy_by_2;
        }
        else {
            output[IDX(i,je-3,k)] = ( -   dev_var_in[u_offset + IDX(i,je-6,k)]
                                +  6.0 * dev_var_in[u_offset + IDX(i,je-5,k)]
                                - 18.0 * dev_var_in[u_offset + IDX(i,je-4,k)]
                                + 10.0 * dev_var_in[u_offset + IDX(i,je-3,k)]
                                +  3.0 * dev_var_in[u_offset + IDX(i,je-2,k)]
                            ) * idy_by_12;
        }
    
            if (dev_var_in[betay + IDX(i,je-2,k)] > 0.0 ) {
            output[IDX(i,je-2,k)] = (  -  dev_var_in[u_offset + IDX(i,je-3,k)]
                                    +  dev_var_in[u_offset + IDX(i,je-1,k)]
                                    ) * idy_by_2;
            }
            else {
            output[IDX(i,je-2,k)] = (     dev_var_in[u_offset + IDX(i,je-4,k)]
                                - 4.0 * dev_var_in[u_offset + IDX(i,je-3,k)]
                                + 3.0 * dev_var_in[u_offset + IDX(i,je-2,k)]
                                    ) * idy_by_2;
            }
    
            output[IDX(i,je-1,k)]  = (          dev_var_in[u_offset + IDX(i,je-3,k)]
                                    - 4.0 * dev_var_in[u_offset + IDX(i,je-2,k)]
                                    + 3.0 * dev_var_in[u_offset + IDX(i,je-1,k)]
                                ) * idy_by_2;
    }
}

__device__ void calc_deriv42_adv_z_bflag(int id, double * output, double * dev_var_in, int u_offset, double dz, int betaz, 
const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag) {

    int nx = host_sz_x;
    int ny = host_sz_y;
    int n = nx * ny;
    
    const double idz = 1.0/dz;
    const double idz_by_2 = 0.50 * idz;
    const double idz_by_12 = idz / 12.0;

    int i = id%(host_sz_x-6) + 3;
    int j = ((id/(host_sz_x-6))%(host_sz_y-6)) + 3;
    int k = (id/(host_sz_z-6)/(host_sz_x-6)) + 3;
    if (k>=host_sz_z-3) return;
    int pp = IDX(i, j, k);

    if (dev_var_in[betaz + pp] > 0.0 ) {
            output[pp] = ( -  3.0 * dev_var_in[u_offset + pp - n]
                        - 10.0 * dev_var_in[u_offset + pp]
                        + 18.0 * dev_var_in[u_offset + pp + n]
                        -  6.0 * dev_var_in[u_offset + pp + 2*n]
                        +        dev_var_in[u_offset + pp + 3*n]
                    ) * idz_by_12;
    }
    else {
        output[pp] = ( -        dev_var_in[u_offset + pp - 3*n]
                    +  6.0 * dev_var_in[u_offset + pp - 2*n]
                    - 18.0 * dev_var_in[u_offset + pp - n]
                    + 10.0 * dev_var_in[u_offset + pp]
                    +  3.0 * dev_var_in[u_offset + pp +n]
                    ) * idz_by_12;
                
    }

    if ((bflag & (1u<<OCT_DIR_BACK)) && (k == 3)) {
            
        output[IDX(i,j,3)] = ( -  3.0 * dev_var_in[u_offset + IDX(i,j,3)]
                +  4.0 * dev_var_in[u_offset + IDX(i,j,4)]
                -        dev_var_in[u_offset + IDX(i,j,5)]
                ) * idz_by_2;
                
        if (dev_var_in[betaz + IDX(i,j,4)] > 0.0) {
            output[IDX(i,j,4)] = ( -  3.0 * dev_var_in[u_offset + IDX(i,j,4)]
                            +  4.0 * dev_var_in[u_offset + IDX(i,j,5)]
                            -        dev_var_in[u_offset + IDX(i,j,6)]
                        ) * idz_by_2;

        }
        else {
            output[IDX(i,j,4)] = ( -         dev_var_in[u_offset + IDX(i,j,3)]
                            +        dev_var_in[u_offset + IDX(i,j,5)]
                        ) * idz_by_2;
                        
        }

        if (dev_var_in[betaz + IDX(i,j,5)] > 0.0 ) {
            output[IDX(i,j,5)] = (-  3.0 * dev_var_in[u_offset + IDX(i,j,4)]
                        - 10.0 * dev_var_in[u_offset + IDX(i,j,5)]
                        + 18.0 * dev_var_in[u_offset + IDX(i,j,6)]
                        -  6.0 * dev_var_in[u_offset + IDX(i,j,7)]
                        +        dev_var_in[u_offset + IDX(i,j,8)]
                        ) * idz_by_12;
        }
        else {
            output[IDX(i,j,5)] = (           dev_var_in[u_offset + IDX(i,j,3)]
                            -  4.0 * dev_var_in[u_offset + IDX(i,j,4)]
                            +  3.0 * dev_var_in[u_offset + IDX(i,j,5)]
                        ) * idz_by_2;
        }
    }

    if ((bflag & (1u<<OCT_DIR_FRONT)) && (k == 4)) {
        
        const int ke = host_sz_z - 3; // Here I changed
        
        if ( dev_var_in[betaz + IDX(i,j,ke-3)] < 0.0 ) {
            output[IDX(i,j,ke-3)] = (  - 3.0 * dev_var_in[u_offset + IDX(i,j,ke-3)]
                                    + 4.0 * dev_var_in[u_offset + IDX(i,j,ke-2)]
                                    -       dev_var_in[u_offset + IDX(i,j,ke-1)]
                                    ) * idz_by_2;
        }
        else {
            output[IDX(i,j,ke-3)] = ( -   dev_var_in[u_offset + IDX(i,j,ke-6)]
                                +  6.0 * dev_var_in[u_offset + IDX(i,j,ke-5)]
                                - 18.0 * dev_var_in[u_offset + IDX(i,j,ke-4)]
                                + 10.0 * dev_var_in[u_offset + IDX(i,j,ke-3)]
                                +  3.0 * dev_var_in[u_offset + IDX(i,j,ke-2)]
                            ) * idz_by_12;
        }
    
            if (dev_var_in[betaz + IDX(i,j,ke-2)] > 0.0 ) {
            output[IDX(i,j,ke-2)] = (  -  dev_var_in[u_offset + IDX(i,j,ke-3)]
                                    +  dev_var_in[u_offset + IDX(i,j,ke-1)]
                                    ) * idz_by_2;
            }
            else {
            output[IDX(i,j,ke-2)] = (     dev_var_in[u_offset + IDX(i,j,ke-4)]
                                - 4.0 * dev_var_in[u_offset + IDX(i,j,ke-3)]
                                + 3.0 * dev_var_in[u_offset + IDX(i,j,ke-2)]
                                    ) * idz_by_2;
            }
    
            output[IDX(i,j,ke-1)]  = (          dev_var_in[u_offset + IDX(i,j,ke-3)]
                                    - 4.0 * dev_var_in[u_offset + IDX(i,j,ke-2)]
                                    + 3.0 * dev_var_in[u_offset + IDX(i,j,ke-1)]
                                ) * idz_by_2;
    }
}


__device__ void calc_ko_deriv42_x_bflag(int id, double * output, double * dev_var_in, int u_offset, double dx, 
const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag)
{
    int nx = host_sz_x;
    int ny = host_sz_y;

    int i = id%(host_sz_x-6) + 3;
    int j = ((id/(host_sz_x-6))%(host_sz_y-6)) + 3;
    int k = (id/(host_sz_z-6)/(host_sz_x-6)) + 3;
    if (k>=host_sz_z-3) return;
    int pp = IDX(i, j, k);

    
    if(i==4) {
        int ib=3;
        output[IDX(3, j, k)] = (-1.0 / 64.0 / dx) *
                        (
                        -      dev_var_in[u_offset + IDX(ib+4,j,k)]
                        +  6.0*dev_var_in[u_offset + IDX(ib+3,j,k)]
                        - 15.0*dev_var_in[u_offset + IDX(ib+2,j,k)]
                        + 20.0*dev_var_in[u_offset + IDX(ib+1,j,k)]
                        - 15.0*dev_var_in[u_offset + IDX(ib,j,k)]
                        +  6.0*dev_var_in[u_offset + IDX(ib-1,j,k)]
                        -      dev_var_in[u_offset + IDX(ib-2,j,k)]
                        );
    }

    output[pp] = (-1.0 / 64.0 / dx) *
                            (
                            -      dev_var_in[u_offset + pp - 3]
                            +  6.0*dev_var_in[u_offset + pp - 2]
                            - 15.0*dev_var_in[u_offset + pp - 1]
                            + 20.0*dev_var_in[u_offset + pp ]
                            - 15.0*dev_var_in[u_offset + pp + 1]
                            +  6.0*dev_var_in[u_offset + pp + 2]
                            -      dev_var_in[u_offset + pp + 3]
                            );

    if(i==5) {
        int ie = nx-3;
        output[IDX(ie-1, j, k)] = (-1.0 / 64.0 / dx) *
                        (
                        -      dev_var_in[u_offset + IDX(ie+1,j,k)]
                        +  6.0*dev_var_in[u_offset + IDX(ie,j,k)]
                        - 15.0*dev_var_in[u_offset + IDX(ie-1,j,k)]
                        + 20.0*dev_var_in[u_offset + IDX(ie-2,j,k)]
                        - 15.0*dev_var_in[u_offset + IDX(ie-3,j,k)]
                        +  6.0*dev_var_in[u_offset + IDX(ie-4,j,k)]
                        -      dev_var_in[u_offset + IDX(ie-5,j,k)]
                        );
    }

    if ((bflag & (1u<<OCT_DIR_LEFT)) && (i == 4)) {

        output[IDX(3,j,k)] =  (      dev_var_in[u_offset + IDX(6,j,k)]
                                    - 3.0*dev_var_in[u_offset + IDX(5,j,k)]
                                    + 3.0*dev_var_in[u_offset + IDX(4,j,k)]
                                    -     dev_var_in[u_offset + IDX(3,j,k)]
                                )/59.0/48.0*64*dx;
        output[IDX(4,j,k)] =  (     dev_var_in[u_offset + IDX(7,j,k)]
                                    -  6.0*dev_var_in[u_offset + IDX(6,j,k)]
                                    + 12.0*dev_var_in[u_offset + IDX(5,j,k)]
                                    - 10.0*dev_var_in[u_offset + IDX(4,j,k)]
                                    +  3.0*dev_var_in[u_offset + IDX(3,j,k)]
                                    )/43.0/48.0*64*dx;
        output[IDX(5,j,k)] =  (     dev_var_in[u_offset + IDX(8,j,k)]
                                    -  6.0*dev_var_in[u_offset + IDX(7,j,k)]
                                    + 15.0*dev_var_in[u_offset + IDX(6,j,k)]
                                    - 19.0*dev_var_in[u_offset + IDX(5,j,k)]
                                    + 12.0*dev_var_in[u_offset + IDX(4,j,k)]
                                    -  3.0*dev_var_in[u_offset + IDX(3,j,k)]
                                    )/49.0/48.0*64*dx;
        }

    if ((bflag & (1u<<OCT_DIR_RIGHT)) && (i == 5)) {
        
        const int ie = nx - 3;
        output[IDX(ie-3,j,k)] = ( dev_var_in[u_offset + IDX(ie-6,j,k)]
                                    - 6.0*dev_var_in[u_offset + IDX(ie-5,j,k)]
                                    + 15.0*dev_var_in[u_offset + IDX(ie-4,j,k)]
                                    - 19.0*dev_var_in[u_offset + IDX(ie-3,j,k)]
                                    + 12.0*dev_var_in[u_offset + IDX(ie-2,j,k)]
                                    -  3.0*dev_var_in[u_offset + IDX(ie-1,j,k)]
                                    )/49.0/48.0*64*dx;
            
            output[IDX(ie-2,j,k)] =  ( dev_var_in[u_offset + IDX(ie-5,j,k)]
                                    -  6.0*dev_var_in[u_offset + IDX(ie-4,j,k)]
                                    + 12.0*dev_var_in[u_offset + IDX(ie-3,j,k)]
                                    - 10.0*dev_var_in[u_offset + IDX(ie-2,j,k)]
                                    +  3.0*dev_var_in[u_offset + IDX(ie-1,j,k)]
                                    )/43.0/48.0*64*dx;
        
    
            output[IDX(ie-1,j,k)] = ( dev_var_in[u_offset + IDX(ie-4,j,k)]
                                    -  3.0*dev_var_in[u_offset + IDX(ie-3,j,k)]
                                    +  3.0*dev_var_in[u_offset + IDX(ie-2,j,k)]
                                    -      dev_var_in[u_offset + IDX(ie-1,j,k)]
                                    )/59.0/48.0*64*dx;
    }
}

__device__ void calc_ko_deriv42_y_bflag(int id, double * output, double * dev_var_in, int u_offset, double dy,
const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag)
{
    
    int nx = host_sz_x;
    int ny = host_sz_y;

    int i = id%(host_sz_x-6) + 3;
    int j = ((id/(host_sz_x-6))%(host_sz_y-6)) + 3;
    int k = (id/(host_sz_z-6)/(host_sz_x-6)) + 3;
    if (k>=host_sz_z-3) return;
    int pp = IDX(i, j, k);

    if(j==4) {
        int jb=3;
        output[IDX(i,jb,k)] = (-1.0 / 64.0 / dy) *
                    (
                        -      dev_var_in[u_offset + IDX(i,jb+4,k)]
                        +  6.0*dev_var_in[u_offset + IDX(i,jb+3,k)]
                        - 15.0*dev_var_in[u_offset + IDX(i,jb+2,k)]
                        + 20.0*dev_var_in[u_offset + IDX(i,jb+1,k)]
                        - 15.0*dev_var_in[u_offset + IDX(i,jb,k)]
                        +  6.0*dev_var_in[u_offset + IDX(i,jb-1,k)]
                        -      dev_var_in[u_offset + IDX(i,jb-2,k)]
                        );
        }
    output[pp] = (-1.0 / 64.0 / dy) *
                    (
                        -      dev_var_in[u_offset + pp-3*nx]
                        +  6.0*dev_var_in[u_offset + pp-2*nx]
                        - 15.0*dev_var_in[u_offset + pp-nx]
                        + 20.0*dev_var_in[u_offset + pp]
                        - 15.0*dev_var_in[u_offset + pp+nx]
                        +  6.0*dev_var_in[u_offset + pp+2*nx]
                        -      dev_var_in[u_offset + pp+3*nx]
                        );

    if(j==5) {
        int je = ny - 3;
        output[IDX(i,je-1,k)] = (-1.0 / 64.0 / dy) *
                (
                    -      dev_var_in[u_offset + IDX(i,je+1,k)]
                    +  6.0*dev_var_in[u_offset + IDX(i,je,k)]
                    - 15.0*dev_var_in[u_offset + IDX(i,je-1,k)]
                    + 20.0*dev_var_in[u_offset + IDX(i,je-2,k)]
                    - 15.0*dev_var_in[u_offset + IDX(i,je-3,k)]
                    +  6.0*dev_var_in[u_offset + IDX(i,je-4,k)]
                    -      dev_var_in[u_offset + IDX(i,je-5,k)]
                    );                   
    }

    if ((bflag & (1u<<OCT_DIR_DOWN)) && (j == 4)) {

        output[IDX(i,3,k)] =  (      dev_var_in[u_offset +IDX(i,6,k)]
                                    - 3.0*dev_var_in[u_offset +IDX(i,5,k)]
                                    + 3.0*dev_var_in[u_offset + IDX(i,4,k)]
                                    -     dev_var_in[u_offset + IDX(i,3,k)]
                                )/59.0/48.0*64*dy;
        output[IDX(i,4,k)] =  (     dev_var_in[u_offset + IDX(i,7,k)]
                                    -  6.0*dev_var_in[u_offset + IDX(i,6,k)]
                                    + 12.0*dev_var_in[u_offset + IDX(i,5,k)]
                                    - 10.0*dev_var_in[u_offset + IDX(i,4,k)]
                                    +  3.0*dev_var_in[u_offset + IDX(i,3,k)]
                                    )/43.0/48.0*64*dy;
        output[IDX(i,5,k)] =  (     dev_var_in[u_offset + IDX(i,8,k)]
                                    -  6.0*dev_var_in[u_offset + IDX(i,7,k)]
                                    + 15.0*dev_var_in[u_offset + IDX(i,6,k)]
                                    - 19.0*dev_var_in[u_offset + IDX(i,5,k)]
                                    + 12.0*dev_var_in[u_offset + IDX(i,4,k)]
                                    -  3.0*dev_var_in[u_offset + IDX(i,3,k)]
                                    )/49.0/48.0*64*dy;
        }

    if ((bflag & (1u<<OCT_DIR_UP)) && (j == 5)) {
        
        const int je = ny - 3;
        output[IDX(i,je-3,k)] = (dev_var_in[u_offset + IDX(i,je-6,k)]
                                    -  6.0*dev_var_in[u_offset + IDX(i,je-5,k)]
                                    + 15.0*dev_var_in[u_offset + IDX(i,je-4,k)]
                                    - 19.0*dev_var_in[u_offset + IDX(i,je-3,k)]
                                    + 12.0*dev_var_in[u_offset + IDX(i,je-2,k)]
                                    -  3.0*dev_var_in[u_offset + IDX(i,je-1,k)]
                                    )/49.0/48.0*64*dy;
            
            output[IDX(i,je-2,k)] = (dev_var_in[u_offset + IDX(i,je-5,k)]
                                    -  6.0*dev_var_in[u_offset + IDX(i,je-4,k)]
                                    + 12.0*dev_var_in[u_offset + IDX(i,je-3,k)]
                                    - 10.0*dev_var_in[u_offset + IDX(i,je-2,k)]
                                    +  3.0*dev_var_in[u_offset + IDX(i,je-1,k)]
                                    )/43.0/48.0*64*dy;
        
    
            output[IDX(i,je-1,k)] = ( dev_var_in[u_offset + IDX(i,je-4,k)]
                                    -  3.0*dev_var_in[u_offset + IDX(i,je-3,k)]
                                    +  3.0*dev_var_in[u_offset + IDX(i,je-2,k)]
                                    -      dev_var_in[u_offset + IDX(i,je-1,k)]
                                    )/59.0/48.0*64*dy;
    }
}

__device__ void calc_ko_deriv42_z_bflag(int id, double * output, double * dev_var_in, int u_offset, double dz, 
const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag)
{  
    int nx = host_sz_x;
    int ny = host_sz_y;
    int n = nx * ny;

    int i = id%(host_sz_x-6) + 3;
    int j = ((id/(host_sz_x-6))%(host_sz_y-6)) + 3;
    int k = (id/(host_sz_z-6)/(host_sz_x-6)) + 3;
    if (k>=host_sz_z-3) return;
    int pp = IDX(i, j, k);
    
    if(k==4) {
        int kb=3;
        output[IDX(i,j,kb)] = (-1.0 / 64.0 / dz) *
                    (
                        -      dev_var_in[u_offset + IDX(i,j,kb+4)]
                        +  6.0*dev_var_in[u_offset + IDX(i,j,kb+3)]
                        - 15.0*dev_var_in[u_offset + IDX(i,j,kb+2)]
                        + 20.0*dev_var_in[u_offset + IDX(i,j,kb+1)]
                        - 15.0*dev_var_in[u_offset + IDX(i,j,kb)]
                        +  6.0*dev_var_in[u_offset + IDX(i,j,kb-1)]
                        -      dev_var_in[u_offset + IDX(i,j,kb-2)]
                        );
        }

    output[pp] = (-1.0 / 64.0 / dz) *
                (
                    -      dev_var_in[u_offset + pp-3*n]
                    +  6.0*dev_var_in[u_offset + pp-2*n]
                    - 15.0*dev_var_in[u_offset + pp-n]
                    + 20.0*dev_var_in[u_offset + pp]
                    - 15.0*dev_var_in[u_offset + pp+n]
                    +  6.0*dev_var_in[u_offset + pp+2*n]
                    -      dev_var_in[u_offset + pp+3*n]
                    );
    if(k==5) {
        int ke = host_sz_z - 3;
        output[IDX(i,j,ke-1)] = (-1.0 / 64.0 / dz) *
        (
            -      dev_var_in[u_offset + IDX(i,j,ke+1)]
            +  6.0*dev_var_in[u_offset + IDX(i,j,ke)]
            - 15.0*dev_var_in[u_offset + IDX(i,j,ke-1)]
            + 20.0*dev_var_in[u_offset + IDX(i,j,ke-2)]
            - 15.0*dev_var_in[u_offset + IDX(i,j,ke-3)]
            +  6.0*dev_var_in[u_offset + IDX(i,j,ke-4)]
            -      dev_var_in[u_offset + IDX(i,j,ke-5)]
            );               
    }

    if ((bflag & (1u<<OCT_DIR_BACK)) && (k == 4)) {

        output[IDX(i,3,k)] =  (      dev_var_in[u_offset +IDX(i,k,6)]
                                    - 3.0*dev_var_in[u_offset +IDX(i,k,5)]
                                    + 3.0*dev_var_in[u_offset + IDX(i,k,4)]
                                    -     dev_var_in[u_offset + IDX(i,k,3)]
                                )/59.0/48.0*64*dz;

        output[IDX(i,j,4)] =  (     dev_var_in[u_offset + IDX(i,j,7)]
                                    -  6.0*dev_var_in[u_offset + IDX(i,j,6)]
                                    + 12.0*dev_var_in[u_offset + IDX(i,j,5)]
                                    - 10.0*dev_var_in[u_offset + IDX(i,j,4)]
                                    +  3.0*dev_var_in[u_offset + IDX(i,j,3)]
                                    )/43.0/48.0*64*dz;

        output[IDX(i,j,5)] =  (     dev_var_in[u_offset + IDX(i,j,8)]
                                    -  6.0*dev_var_in[u_offset + IDX(i,j,7)]
                                    + 15.0*dev_var_in[u_offset + IDX(i,j,6)]
                                    - 19.0*dev_var_in[u_offset + IDX(i,j,5)]
                                    + 12.0*dev_var_in[u_offset + IDX(i,j,4)]
                                    -  3.0*dev_var_in[u_offset + IDX(i,j,3)]
                                    )/49.0/48.0*64*dz;
        }

    if ((bflag & (1u<<OCT_DIR_FRONT)) && (k == 5)) {
        
        const int ke = host_sz_z - 3;
        output[IDX(i,j,ke-3)] = (    dev_var_in[u_offset + IDX(i,j,ke-6)]
                                        -  6.0*dev_var_in[u_offset + IDX(i,j,ke-5)]
                                        + 15.0*dev_var_in[u_offset + IDX(i,j,ke-4)]
                                        - 19.0*dev_var_in[u_offset + IDX(i,j,ke-3)]
                                        + 12.0*dev_var_in[u_offset + IDX(i,j,ke-2)]
                                        -  3.0*dev_var_in[u_offset + IDX(i,j,ke-1)]
                                        )/49.0/48.0*64*dz;
            
            output[IDX(i,j,ke-2)] = (   dev_var_in[u_offset + IDX(i,j,ke-5)]
                                        -  6.0*dev_var_in[u_offset + IDX(i,j,ke-4)]
                                        + 12.0*dev_var_in[u_offset + IDX(i,j,ke-3)]
                                        - 10.0*dev_var_in[u_offset + IDX(i,j,ke-2)]
                                        +  3.0*dev_var_in[u_offset + IDX(i,j,ke-1)]
                                        )/43.0/48.0*64*dz;
        
    
            output[IDX(i,j,ke-1)] = (   dev_var_in[u_offset + IDX(i,j,ke-4)]
                                        -  3.0*dev_var_in[u_offset + IDX(i,j,ke-3)]
                                        +  3.0*dev_var_in[u_offset + IDX(i,j,ke-2)]
                                        -      dev_var_in[u_offset + IDX(i,j,ke-1)]
                                        )/59.0/48.0*64*dz;
    }
        
}