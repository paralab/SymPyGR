/**
 * Created on: Feb 12, 2018
 * 		Author: Eminda, Akila, Eranga, Ruwan
 **/

#ifndef RHS_CUDA_H_
#define RHS_CUDA_H_
#include <cmath>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "derivs_cuda.h"
#include "profile_param.h"

#define deriv_x cuda_deriv42_x
#define deriv_y cuda_deriv42_y
#define deriv_z cuda_deriv42_z

#define deriv_xx cuda_deriv42_xx
#define deriv_yy cuda_deriv42_yy
#define deriv_zz cuda_deriv42_zz

#define adv_deriv_x cuda_deriv42_adv_x
#define adv_deriv_y cuda_deriv42_adv_y
#define adv_deriv_z cuda_deriv42_adv_z

#define ko_deriv_x cuda_ko_deriv42_x
#define ko_deriv_y cuda_ko_deriv42_y
#define ko_deriv_z cuda_ko_deriv42_z

void cuda_bssnrhs(double * dev_var_out, double * dev_var_in, 
const unsigned int unzip_dof, const unsigned int& offset, 
const double *pmin,const double *pmax, const unsigned int *sz, 
const unsigned int& bflag);

void bssn_bcs(double * output, double * dev_var_in, int* dev_u_offset,
    double *dxf, double *dyf, double *dzf,
    double *pmin, double *pmax, const double f_falloff, const double f_asymptotic,
    const unsigned int *host_sz, int* dev_bflag, int* dev_sz);


void get_output (double* output,

        int *dev_alphaInt,
        int *dev_chiInt,
        int *dev_KInt,
        int *dev_gt0Int,
        int *dev_gt1Int,
        int *dev_gt2Int,
        int *dev_gt3Int,
        int *dev_gt4Int,
        int *dev_gt5Int,
        int *dev_beta0Int,
        int *dev_beta1Int,
        int *dev_beta2Int,
        int *dev_At0Int,
        int *dev_At1Int,
        int *dev_At2Int,
        int *dev_At3Int,
        int *dev_At4Int,
        int *dev_At5Int,
        int *dev_Gt0Int,
        int *dev_Gt1Int,
        int *dev_Gt2Int,
        int *dev_B0Int,
        int *dev_B1Int,
        int *dev_B2Int,

        double *grad_0_alpha,
        double *grad_1_alpha,
        double *grad_2_alpha,
        double *grad_0_beta0,
        double *grad_1_beta0,
        double *grad_2_beta0,
        double *grad_0_beta1,
        double *grad_1_beta1,
        double *grad_2_beta1,
        double *grad_0_beta2,
        double *grad_1_beta2,
        double *grad_2_beta2,

        double *grad_0_gt0,
        double *grad_1_gt0,
        double *grad_2_gt0,
        double *grad_0_gt1,
        double *grad_1_gt1,
        double *grad_2_gt1,
        double *grad_0_gt2,
        double *grad_1_gt2,
        double *grad_2_gt2,
        double *grad_0_gt3,
        double *grad_1_gt3,
        double *grad_2_gt3,
        double *grad_0_gt4,
        double *grad_1_gt4,
        double *grad_2_gt4,
        double *grad_0_gt5,
        double *grad_1_gt5,
        double *grad_2_gt5,  

        double *grad_0_chi,
        double *grad_1_chi,
        double *grad_2_chi,

        double *grad_0_At0,
        double *grad_1_At0,
        double *grad_2_At0,
        double *grad_0_At1,
        double *grad_1_At1,
        double *grad_2_At1,
        double *grad_0_At2,
        double *grad_1_At2,
        double *grad_2_At2,
        double *grad_0_At3,
        double *grad_1_At3,
        double *grad_2_At3,
        double *grad_0_At4,
        double *grad_1_At4,
        double *grad_2_At4,
        double *grad_0_At5,
        double *grad_1_At5,
        double *grad_2_At5,

        double *grad_0_K,
        double *grad_1_K,
        double *grad_2_K,

        double *grad_0_Gt0,
        double *grad_1_Gt0,
        double *grad_2_Gt0,
        double *grad_0_Gt1,
        double *grad_1_Gt1,
        double *grad_2_Gt1,
        double *grad_0_Gt2,
        double *grad_1_Gt2,
        double *grad_2_Gt2,

        double *grad_0_B0,
        double *grad_1_B0,
        double *grad_2_B0,
        double *grad_0_B1,
        double *grad_1_B1,
        double *grad_2_B1,
        double *grad_0_B2,
        double *grad_1_B2,
        double *grad_2_B2,
        
        int* dev_sz,
        const unsigned int* host_sz);

#endif
