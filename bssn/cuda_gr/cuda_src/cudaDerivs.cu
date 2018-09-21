/**
 * Created on: Sep 21, 2018
 * 		Author: Akila, Eranga, Eminda, Ruwan
 **/
 
#include "cudaDerivs.cuh"

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

__global__ void calc_derivs1_bflag(
    double * dev_var_out, double * dev_var_in, double hx, double hy, double hz, 
    const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, 
    int bflag,
    #include "list_of_offset_para.h"
    ,
    #include "list_of_para.h"
){
    int tid = blockIdx.x*64 + threadIdx.x;

    calc_deriv42_x_bflag(tid, grad_0_alpha, dev_var_in, alphaInt, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_xx_bflag(tid, grad2_0_0_alpha, dev_var_in, alphaInt, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad_1_alpha, dev_var_in, alphaInt, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_yy_bflag(tid, grad2_1_1_alpha, dev_var_in, alphaInt, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad_2_alpha, dev_var_in, alphaInt, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_zz_bflag(tid, grad2_2_2_alpha, dev_var_in, alphaInt, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x_bflag(tid, grad_0_beta0, dev_var_in, beta0Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_xx_bflag(tid, grad2_0_0_beta0, dev_var_in, beta0Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad_1_beta0, dev_var_in, beta0Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_yy_bflag(tid, grad2_1_1_beta0, dev_var_in, beta0Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad_2_beta0, dev_var_in, beta0Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_zz_bflag(tid, grad2_2_2_beta0, dev_var_in, beta0Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x_bflag(tid, grad_0_beta1, dev_var_in, beta1Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_xx_bflag(tid, grad2_0_0_beta1, dev_var_in, beta1Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad_1_beta1, dev_var_in, beta1Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_yy_bflag(tid, grad2_1_1_beta1, dev_var_in, beta1Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad_2_beta1, dev_var_in, beta1Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_zz_bflag(tid, grad2_2_2_beta1, dev_var_in, beta1Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x_bflag(tid, grad_0_beta2, dev_var_in, beta2Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_xx_bflag(tid, grad2_0_0_beta2, dev_var_in, beta2Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad_1_beta2, dev_var_in, beta2Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_yy_bflag(tid, grad2_1_1_beta2, dev_var_in, beta2Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad_2_beta2, dev_var_in, beta2Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_zz_bflag(tid, grad2_2_2_beta2, dev_var_in, beta2Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x_bflag(tid, grad_0_B0, dev_var_in, B0Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad_1_B0, dev_var_in, B0Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad_2_B0, dev_var_in, B0Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x_bflag(tid, grad_0_B1, dev_var_in, B1Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad_1_B1, dev_var_in, B1Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad_2_B1, dev_var_in, B1Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x_bflag(tid, grad_0_B2, dev_var_in, B2Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad_1_B2, dev_var_in, B2Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad_2_B2, dev_var_in, B2Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x_bflag(tid, grad_0_chi, dev_var_in, chiInt, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_xx_bflag(tid, grad2_0_0_chi, dev_var_in, chiInt, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad_1_chi, dev_var_in, chiInt, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_yy_bflag(tid, grad2_1_1_chi, dev_var_in, chiInt, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad_2_chi, dev_var_in, chiInt, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_zz_bflag(tid, grad2_2_2_chi, dev_var_in, chiInt, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x_bflag(tid, grad_0_Gt0, dev_var_in, Gt0Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad_1_Gt0, dev_var_in, Gt0Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad_2_Gt0, dev_var_in, Gt0Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x_bflag(tid, grad_0_gt0, dev_var_in, gt0Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_xx_bflag(tid, grad2_0_0_gt0, dev_var_in, gt0Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad_1_gt0, dev_var_in, gt0Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_yy_bflag(tid, grad2_1_1_gt0, dev_var_in, gt0Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad_2_gt0, dev_var_in, gt0Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_zz_bflag(tid, grad2_2_2_gt0, dev_var_in, gt0Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x_bflag(tid, grad_0_Gt1, dev_var_in, Gt1Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad_1_Gt1, dev_var_in, Gt1Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad_2_Gt1, dev_var_in, Gt1Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x_bflag(tid, grad_0_Gt2, dev_var_in, Gt2Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad_1_Gt2, dev_var_in, Gt2Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad_2_Gt2, dev_var_in, Gt2Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x_bflag(tid, grad_0_K, dev_var_in, KInt, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad_1_K, dev_var_in, KInt, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad_2_K, dev_var_in, KInt, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x_bflag(tid, grad_0_gt1, dev_var_in, gt1Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_xx_bflag(tid, grad2_0_0_gt1, dev_var_in, gt1Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad_1_gt1, dev_var_in, gt1Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_yy_bflag(tid, grad2_1_1_gt1, dev_var_in, gt1Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad_2_gt1, dev_var_in, gt1Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_zz_bflag(tid, grad2_2_2_gt1, dev_var_in, gt1Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x_bflag(tid, grad_0_gt2, dev_var_in, gt2Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_xx_bflag(tid, grad2_0_0_gt2, dev_var_in, gt2Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad_1_gt2, dev_var_in, gt2Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_yy_bflag(tid, grad2_1_1_gt2, dev_var_in, gt2Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad_2_gt2, dev_var_in, gt2Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_zz_bflag(tid, grad2_2_2_gt2, dev_var_in, gt2Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x_bflag(tid, grad_0_gt3, dev_var_in, gt3Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_xx_bflag(tid, grad2_0_0_gt3, dev_var_in, gt3Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad_1_gt3, dev_var_in, gt3Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_yy_bflag(tid, grad2_1_1_gt3, dev_var_in, gt3Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad_2_gt3, dev_var_in, gt3Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_zz_bflag(tid, grad2_2_2_gt3, dev_var_in, gt3Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x_bflag(tid, grad_0_gt4, dev_var_in, gt4Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_xx_bflag(tid, grad2_0_0_gt4, dev_var_in, gt4Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad_1_gt4, dev_var_in, gt4Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_yy_bflag(tid, grad2_1_1_gt4, dev_var_in, gt4Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad_2_gt4, dev_var_in, gt4Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_zz_bflag(tid, grad2_2_2_gt4, dev_var_in, gt4Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x_bflag(tid, grad_0_gt5, dev_var_in, gt5Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_xx_bflag(tid, grad2_0_0_gt5, dev_var_in, gt5Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad_1_gt5, dev_var_in, gt5Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_yy_bflag(tid, grad2_1_1_gt5, dev_var_in, gt5Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad_2_gt5, dev_var_in, gt5Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_zz_bflag(tid, grad2_2_2_gt5, dev_var_in, gt5Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x_bflag(tid, grad_0_At0, dev_var_in, At0Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad_1_At0, dev_var_in, At0Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad_2_At0, dev_var_in, At0Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x_bflag(tid, grad_0_At1, dev_var_in, At1Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad_1_At1, dev_var_in, At1Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad_2_At1, dev_var_in, At1Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x_bflag(tid, grad_0_At2, dev_var_in, At2Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad_1_At2, dev_var_in, At2Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad_2_At2, dev_var_in, At2Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x_bflag(tid, grad_0_At3, dev_var_in, At3Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad_1_At3, dev_var_in, At3Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad_2_At3, dev_var_in, At3Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x_bflag(tid, grad_0_At4, dev_var_in, At4Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad_1_At4, dev_var_in, At4Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad_2_At4, dev_var_in, At4Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_x_bflag(tid, grad_0_At5, dev_var_in, At5Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad_1_At5, dev_var_in, At5Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad_2_At5, dev_var_in, At5Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
}

__global__ void calc_derivs2_bflag(
    double * dev_var_out, double * dev_var_in, double hx, double hy, double hz, 
    const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, 
    int bflag,
    #include "list_of_offset_para.h"
    ,
    #include "list_of_para.h"
){
    int tid = blockIdx.x*64 + threadIdx.x;

    calc_deriv42_y_bflag(tid, grad2_0_1_gt0, grad_0_gt0, 0, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad2_0_2_gt0, grad_0_gt0, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad2_1_2_gt0, grad_1_gt0, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad2_0_1_gt1, grad_0_gt1, 0, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad2_0_2_gt1, grad_0_gt1, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad2_1_2_gt1, grad_1_gt1, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad2_0_1_gt2, grad_0_gt2, 0, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad2_0_2_gt2, grad_0_gt2, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad2_1_2_gt2, grad_1_gt2, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad2_0_1_gt3, grad_0_gt3, 0, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad2_0_2_gt3, grad_0_gt3, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad2_1_2_gt3, grad_1_gt3, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad2_0_1_gt4, grad_0_gt4, 0, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad2_0_2_gt4, grad_0_gt4, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad2_1_2_gt4, grad_1_gt4, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad2_0_1_gt5, grad_0_gt5, 0, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad2_0_2_gt5, grad_0_gt5, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad2_1_2_gt5, grad_1_gt5, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad2_0_1_chi, grad_0_chi, 0, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad2_0_2_chi, grad_0_chi, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad2_1_2_chi, grad_1_chi, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad2_0_1_alpha, grad_0_alpha, 0, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad2_0_2_alpha, grad_0_alpha, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad2_1_2_alpha, grad_1_alpha, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad2_0_1_beta0, grad_0_beta0, 0, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad2_0_2_beta0, grad_0_beta0, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad2_1_2_beta0, grad_1_beta0, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad2_0_1_beta1, grad_0_beta1, 0, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad2_0_2_beta1, grad_0_beta1, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad2_1_2_beta1, grad_1_beta1, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_y_bflag(tid, grad2_0_1_beta2, grad_0_beta2, 0, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad2_0_2_beta2, grad_0_beta2, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_z_bflag(tid, grad2_1_2_beta2, grad_1_beta2, 0, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_adv_x_bflag(tid, agrad_0_gt0, dev_var_in, gt0Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y_bflag(tid, agrad_1_gt0, dev_var_in, gt0Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_adv_x_bflag(tid, agrad_0_gt1, dev_var_in, gt1Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y_bflag(tid, agrad_1_gt1, dev_var_in, gt1Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_adv_x_bflag(tid, agrad_0_gt2, dev_var_in, gt2Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y_bflag(tid, agrad_1_gt2, dev_var_in, gt2Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_adv_x_bflag(tid, agrad_0_gt3, dev_var_in, gt3Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y_bflag(tid, agrad_1_gt3, dev_var_in, gt3Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_adv_x_bflag(tid, agrad_0_gt4, dev_var_in, gt4Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y_bflag(tid, agrad_1_gt4, dev_var_in, gt4Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_adv_x_bflag(tid, agrad_0_gt5, dev_var_in, gt5Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y_bflag(tid, agrad_1_gt5, dev_var_in, gt5Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);


    calc_deriv42_adv_x_bflag(tid, agrad_0_At0, dev_var_in, At0Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y_bflag(tid, agrad_1_At0, dev_var_in, At0Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_adv_x_bflag(tid, agrad_0_At1, dev_var_in, At1Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y_bflag(tid, agrad_1_At1, dev_var_in, At1Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_adv_x_bflag(tid, agrad_0_At2, dev_var_in, At2Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y_bflag(tid, agrad_1_At2, dev_var_in, At2Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_adv_x_bflag(tid, agrad_0_At3, dev_var_in, At3Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y_bflag(tid, agrad_1_At3, dev_var_in, At3Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_adv_x_bflag(tid, agrad_0_At4, dev_var_in, At4Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y_bflag(tid, agrad_1_At4, dev_var_in, At4Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_adv_x_bflag(tid, agrad_0_At5, dev_var_in, At5Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y_bflag(tid, agrad_1_At5, dev_var_in, At5Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_deriv42_adv_x_bflag(tid, agrad_0_alpha, dev_var_in, alphaInt, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y_bflag(tid, agrad_1_alpha, dev_var_in, alphaInt, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_x_bflag(tid, agrad_0_beta0, dev_var_in, beta0Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_x_bflag(tid, agrad_0_beta1, dev_var_in, beta1Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_x_bflag(tid, agrad_0_beta2, dev_var_in, beta2Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_x_bflag(tid, agrad_0_chi, dev_var_in, chiInt, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_x_bflag(tid, agrad_0_Gt0, dev_var_in, Gt0Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_x_bflag(tid, agrad_0_Gt1, dev_var_in, Gt1Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_x_bflag(tid, agrad_0_Gt2, dev_var_in, Gt2Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_x_bflag(tid, agrad_0_K, dev_var_in, KInt, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_x_bflag(tid, agrad_0_B0, dev_var_in, B0Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_x_bflag(tid, agrad_0_B1, dev_var_in, B1Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_x_bflag(tid, agrad_0_B2, dev_var_in, B2Int, hx, beta0Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y_bflag(tid, agrad_1_beta0, dev_var_in, beta0Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y_bflag(tid, agrad_1_beta1, dev_var_in, beta1Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y_bflag(tid, agrad_1_beta2, dev_var_in, beta2Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y_bflag(tid, agrad_1_chi, dev_var_in, chiInt, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y_bflag(tid, agrad_1_Gt0, dev_var_in, Gt0Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y_bflag(tid, agrad_1_Gt1, dev_var_in, Gt1Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y_bflag(tid, agrad_1_Gt2, dev_var_in, Gt2Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y_bflag(tid, agrad_1_K, dev_var_in, KInt, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y_bflag(tid, agrad_1_B0, dev_var_in, B0Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y_bflag(tid, agrad_1_B1, dev_var_in, B1Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_y_bflag(tid, agrad_1_B2, dev_var_in, B2Int, hy, beta1Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z_bflag(tid, agrad_2_gt0, dev_var_in, gt0Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z_bflag(tid, agrad_2_gt1, dev_var_in, gt1Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z_bflag(tid, agrad_2_gt2, dev_var_in, gt2Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z_bflag(tid, agrad_2_gt3, dev_var_in, gt3Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z_bflag(tid, agrad_2_gt4, dev_var_in, gt4Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z_bflag(tid, agrad_2_gt5, dev_var_in, gt5Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z_bflag(tid, agrad_2_At0, dev_var_in, At0Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z_bflag(tid, agrad_2_At1, dev_var_in, At1Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z_bflag(tid, agrad_2_At2, dev_var_in, At2Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z_bflag(tid, agrad_2_At3, dev_var_in, At3Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z_bflag(tid, agrad_2_At4, dev_var_in, At4Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z_bflag(tid, agrad_2_At5, dev_var_in, At5Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z_bflag(tid, agrad_2_alpha, dev_var_in, alphaInt, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z_bflag(tid, agrad_2_beta0, dev_var_in, beta0Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z_bflag(tid, agrad_2_beta1, dev_var_in, beta1Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z_bflag(tid, agrad_2_beta2, dev_var_in, beta2Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z_bflag(tid, agrad_2_chi, dev_var_in, chiInt, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z_bflag(tid, agrad_2_Gt0, dev_var_in, Gt0Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z_bflag(tid, agrad_2_Gt1, dev_var_in, Gt1Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z_bflag(tid, agrad_2_Gt2, dev_var_in, Gt2Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z_bflag(tid, agrad_2_K, dev_var_in, KInt, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z_bflag(tid, agrad_2_B0, dev_var_in, B0Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z_bflag(tid, agrad_2_B1, dev_var_in, B1Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_deriv42_adv_z_bflag(tid, agrad_2_B2, dev_var_in, B2Int, hz, beta2Int, host_sz_x, host_sz_y, host_sz_z, bflag);
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

    if (bflag!=0){
        number_of_threads_required=ceil((ie-ib)*(je-jb)*(ke-kb));
        number_of_blocks=ceil(1.0*number_of_threads_required/64);
        calc_derivs1_bflag <<< number_of_blocks, 64, 0, stream>>> (dev_var_out, dev_var_in, hx, hy, hz, host_sz_x, host_sz_y, host_sz_z, bflag,
            #include "list_of_offset_args.h"
            ,
            #include "list_of_args.h"
        );
    
        number_of_threads_required=ceil((ie-ib)*(je-jb)*(ke-kb));
        number_of_blocks=ceil(1.0*number_of_threads_required/64);
        calc_derivs2_bflag <<< number_of_blocks, 64, 0, stream>>> (dev_var_out, dev_var_in, hx, hy, hz, host_sz_x, host_sz_y, host_sz_z, bflag,
            #include "list_of_offset_args.h"
            ,
            #include "list_of_args.h"
        );
    }else{
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
    }
    CHECK_ERROR(cudaGetLastError(), "deriv Kernel launch failed");
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


__global__ void calc_ko_derivs_bflag(
    double * dev_var_out, double * dev_var_in, double hx, double hy, double hz, 
    const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, 
    int bflag,
    #include "list_of_offset_para.h"
    ,
    #include "list_of_para.h"
){
    int tid = blockIdx.x*64 + threadIdx.x;

    calc_ko_deriv42_x_bflag(tid, grad_0_gt0, dev_var_in, gt0Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y_bflag(tid, grad_1_gt0, dev_var_in, gt0Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z_bflag(tid, grad_2_gt0, dev_var_in, gt0Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x_bflag(tid, grad_0_gt1, dev_var_in, gt1Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y_bflag(tid, grad_1_gt1, dev_var_in, gt1Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z_bflag(tid, grad_2_gt1, dev_var_in, gt1Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x_bflag(tid, grad_0_gt2, dev_var_in, gt2Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y_bflag(tid, grad_1_gt2, dev_var_in, gt2Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z_bflag(tid, grad_2_gt2, dev_var_in, gt2Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x_bflag(tid, grad_0_gt3, dev_var_in, gt3Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y_bflag(tid, grad_1_gt3, dev_var_in, gt3Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z_bflag(tid, grad_2_gt3, dev_var_in, gt3Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x_bflag(tid, grad_0_gt4, dev_var_in, gt4Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y_bflag(tid, grad_1_gt4, dev_var_in, gt4Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z_bflag(tid, grad_2_gt4, dev_var_in, gt4Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x_bflag(tid, grad_0_gt5, dev_var_in, gt5Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y_bflag(tid, grad_1_gt5, dev_var_in, gt5Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z_bflag(tid, grad_2_gt5, dev_var_in, gt5Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x_bflag(tid, grad_0_At0, dev_var_in, At0Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y_bflag(tid, grad_1_At0, dev_var_in, At0Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z_bflag(tid, grad_2_At0, dev_var_in, At0Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x_bflag(tid, grad_0_At1, dev_var_in, At1Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y_bflag(tid, grad_1_At1, dev_var_in, At1Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z_bflag(tid, grad_2_At1, dev_var_in, At1Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x_bflag(tid, grad_0_At2, dev_var_in, At2Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y_bflag(tid, grad_1_At2, dev_var_in, At2Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z_bflag(tid, grad_2_At2, dev_var_in, At2Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x_bflag(tid, grad_0_At3, dev_var_in, At3Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y_bflag(tid, grad_1_At3, dev_var_in, At3Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z_bflag(tid, grad_2_At3, dev_var_in, At3Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x_bflag(tid, grad_0_At4, dev_var_in, At4Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y_bflag(tid, grad_1_At4, dev_var_in, At4Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z_bflag(tid, grad_2_At4, dev_var_in, At4Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x_bflag(tid, grad_0_At5, dev_var_in, At5Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y_bflag(tid, grad_1_At5, dev_var_in, At5Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z_bflag(tid, grad_2_At5, dev_var_in, At5Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x_bflag(tid, grad_0_alpha, dev_var_in, alphaInt, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y_bflag(tid, grad_1_alpha, dev_var_in, alphaInt, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z_bflag(tid, grad_2_alpha, dev_var_in, alphaInt, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x_bflag(tid, grad_0_beta0, dev_var_in, beta0Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y_bflag(tid, grad_1_beta0, dev_var_in, beta0Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z_bflag(tid, grad_2_beta0, dev_var_in, beta0Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x_bflag(tid, grad_0_beta1, dev_var_in, beta1Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y_bflag(tid, grad_1_beta1, dev_var_in, beta1Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z_bflag(tid, grad_2_beta1, dev_var_in, beta1Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x_bflag(tid, grad_0_beta2, dev_var_in, beta2Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y_bflag(tid, grad_1_beta2, dev_var_in, beta2Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z_bflag(tid, grad_2_beta2, dev_var_in, beta2Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x_bflag(tid, grad_0_chi, dev_var_in, chiInt, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y_bflag(tid, grad_1_chi, dev_var_in, chiInt, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z_bflag(tid, grad_2_chi, dev_var_in, chiInt, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x_bflag(tid, grad_0_Gt0, dev_var_in, Gt0Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y_bflag(tid, grad_1_Gt0, dev_var_in, Gt0Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z_bflag(tid, grad_2_Gt0, dev_var_in, Gt0Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x_bflag(tid, grad_0_Gt1, dev_var_in, Gt1Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y_bflag(tid, grad_1_Gt1, dev_var_in, Gt1Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z_bflag(tid, grad_2_Gt1, dev_var_in, Gt1Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x_bflag(tid, grad_0_Gt2, dev_var_in, Gt2Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y_bflag(tid, grad_1_Gt2, dev_var_in, Gt2Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z_bflag(tid, grad_2_Gt2, dev_var_in, Gt2Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x_bflag(tid, grad_0_K, dev_var_in, KInt, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y_bflag(tid, grad_1_K, dev_var_in, KInt, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z_bflag(tid, grad_2_K, dev_var_in, KInt, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x_bflag(tid, grad_0_B0, dev_var_in, B0Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y_bflag(tid, grad_1_B0, dev_var_in, B0Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z_bflag(tid, grad_2_B0, dev_var_in, B0Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x_bflag(tid, grad_0_B1, dev_var_in, B1Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y_bflag(tid, grad_1_B1, dev_var_in, B1Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z_bflag(tid, grad_2_B1, dev_var_in, B1Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);

    calc_ko_deriv42_x_bflag(tid, grad_0_B2, dev_var_in, B2Int, hx, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_y_bflag(tid, grad_1_B2, dev_var_in, B2Int, hy, host_sz_x, host_sz_y, host_sz_z, bflag);
    calc_ko_deriv42_z_bflag(tid, grad_2_B2, dev_var_in, B2Int, hz, host_sz_x, host_sz_y, host_sz_z, bflag);


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

    if (bflag!=0){
        calc_ko_derivs_bflag <<< number_of_blocks, 64, 0, stream>>> (dev_var_out, dev_var_in, hx, hy, hz, host_sz_x, host_sz_y, host_sz_z, bflag,
            #include "list_of_offset_args.h"
            ,
            #include "list_of_args.h"
        );
    }else{
        calc_ko_derivs <<< number_of_blocks, 64, 0, stream>>> (dev_var_out, dev_var_in, hx, hy, hz, host_sz_x, host_sz_y, host_sz_z, bflag,
            #include "list_of_offset_args.h"
            ,
            #include "list_of_args.h"
        );
    }
    
    CHECK_ERROR(cudaGetLastError(), "ko deriv Kernel launch failed");
}
