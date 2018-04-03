#include "rhs_cuda.h"
#include "bssneqn_solve.h"

#include <iostream>
#include <stdio.h>

enum VAR_CU {U_ALPHA=0,U_CHI,U_K,U_GT0,U_GT1,U_GT2,U_BETA0,U_BETA1,U_BETA2,U_B0,U_B1,U_B2,U_SYMGT0,U_SYMGT1,U_SYMGT2,U_SYMGT3,U_SYMGT4,U_SYMGT5,U_SYMAT0,U_SYMAT1,U_SYMAT2,U_SYMAT3,U_SYMAT4,U_SYMAT5};

double ETA_CONST_CUDA=0.1;
double ETA_R0_CUDA=0.1;
double ETA_DAMPING_EXP_CUDA=0.1;
double KO_DISS_SIGMA_CUDA=1e-4;
unsigned int BSSN_LAMBDA_CUDA[4]={1,2,3,4};
double BSSN_LAMBDA_F_CUDA[2]={0.8,0.9};

void cuda_bssnrhs(double * dev_var_out, double * dev_var_in, const unsigned int unzip_dof, 
const unsigned int& offset, const double *pmin, const double *pmax, const unsigned int *sz, 
const unsigned int& bflag)
{
    cudaError_t cudaStatus;

    int alphaInt = (VAR_CU::U_ALPHA) * unzip_dof + offset;
    int chiInt = (VAR_CU::U_CHI) * unzip_dof + offset;
    int KInt = (VAR_CU::U_K) * unzip_dof + offset;
    int gt0Int = (VAR_CU::U_SYMGT0) * unzip_dof + offset;
    int gt1Int = (VAR_CU::U_SYMGT1) * unzip_dof + offset;
    int gt2Int =  (VAR_CU::U_SYMGT2) * unzip_dof + offset;
    int gt3Int =(VAR_CU::U_SYMGT3) * unzip_dof + offset;
    int gt4Int = (VAR_CU::U_SYMGT4) * unzip_dof + offset;
    int gt5Int = (VAR_CU::U_SYMGT5) * unzip_dof + offset;
    int beta0Int = (VAR_CU::U_BETA0) * unzip_dof + offset;
    int beta1Int = (VAR_CU::U_BETA1) * unzip_dof + offset;
    int beta2Int =(VAR_CU::U_BETA2) * unzip_dof + offset;
    int At0Int = (VAR_CU::U_SYMAT0) * unzip_dof + offset;
    int At1Int = (VAR_CU::U_SYMAT1) * unzip_dof + offset;
    int At2Int = (VAR_CU::U_SYMAT2) * unzip_dof + offset;
    int At3Int = (VAR_CU::U_SYMAT3) * unzip_dof + offset;
    int At4Int = (VAR_CU::U_SYMAT4) * unzip_dof + offset;
    int At5Int = (VAR_CU::U_SYMAT5) * unzip_dof + offset;
    int Gt0Int = (VAR_CU::U_GT0) * unzip_dof + offset;
    int Gt1Int = (VAR_CU::U_GT1) * unzip_dof + offset;
    int Gt2Int = (VAR_CU::U_GT2) * unzip_dof + offset;
    int B0Int = (VAR_CU::U_B0) * unzip_dof + offset;
    int B1Int = (VAR_CU::U_B1) * unzip_dof + offset;
    int B2Int = (VAR_CU::U_B2) * unzip_dof + offset;

    double hx = (pmax[0] - pmin[0]) / (sz[0] - 1);
    double hy = (pmax[1] - pmin[1]) / (sz[1] - 1);
    double hz = (pmax[2] - pmin[2]) / (sz[2] - 1);

    // Send above values to GPU memory
    #include "bssnrhs_cuda_offset_malloc.h"

    double * dev_dy_hx; //similar to hx in cpu code
    cudaStatus = cudaMalloc((void **) &dev_dy_hx, sizeof(double));
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "hx cudaMalloc failed!\n"); return;}
    cudaStatus = cudaMemcpy(dev_dy_hx, &hx, sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "hx cudaMemcpy failed!\n"); return;}

    double * dev_dy_hy;
    cudaStatus = cudaMalloc((void **) &dev_dy_hy, sizeof(double));
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "hy cudaMalloc failed!\n"); return;}
    cudaStatus = cudaMemcpy(dev_dy_hy, &hy, sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "hy cudaMemcpy failed!\n"); return;}

    double * dev_dy_hz;
    cudaStatus = cudaMalloc((void **) &dev_dy_hz, sizeof(double));
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "hz cudaMalloc failed!\n"); return;}
    cudaStatus = cudaMemcpy(dev_dy_hz, &hz, sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "hz cudaMemcpy failed!\n"); return;}

    int * dev_sz;
    cudaStatus = cudaMalloc((void **) &dev_sz, 3*sizeof(int));
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "sz cudaMalloc failed!\n"); return;}
    cudaStatus = cudaMemcpy(dev_sz, sz, 3*sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "sz cudaMemcpy failed!\n"); return;}

    int * dev_zero;
    cudaStatus = cudaMalloc((void **) &dev_zero, sizeof(int));
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "0 cudaMalloc failed!\n"); return;}

    double *dev_pmin;
    cudaStatus = cudaMalloc((void **) &dev_pmin, sizeof(pmin)*sizeof(double));
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "pmin cudaMalloc failed!\n"); return;}
    cudaStatus = cudaMemcpy(dev_pmin, pmin, sizeof(pmin)*sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "pmin cudaMemcpy failed!\n"); return;}

    double *dev_pmax;
    cudaStatus = cudaMalloc((void **) &dev_pmax, sizeof(pmax)*sizeof(double));
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "pmax cudaMalloc failed!\n"); return;}
    cudaStatus = cudaMemcpy(dev_pmax, pmax, sizeof(pmax)*sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "pmax cudaMemcpy failed!\n"); return;}

    const unsigned int lambda[4] = {BSSN_LAMBDA_CUDA[0], BSSN_LAMBDA_CUDA[1],
        BSSN_LAMBDA_CUDA[2], BSSN_LAMBDA_CUDA[3]};
    const double lambda_f[2] = {BSSN_LAMBDA_F_CUDA[0], BSSN_LAMBDA_F_CUDA[1]};

    unsigned int *dev_lambda;
    cudaStatus = cudaMalloc((void **) &dev_lambda, sizeof(lambda)*sizeof(int));
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "lambda cudaMalloc failed!\n"); return;}
    cudaStatus = cudaMemcpy(dev_lambda, lambda, sizeof(lambda)*sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "lambda cudaMemcpy failed!\n"); return;}

    double *dev_lambda_f;
    cudaStatus = cudaMalloc((void **) &dev_lambda_f, sizeof(lambda_f)*sizeof(double));
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "lambda_f cudaMalloc failed!\n"); return;}
    cudaStatus = cudaMemcpy(dev_lambda_f, lambda_f, sizeof(lambda_f)*sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "lambda_f cudaMemcpy failed!\n"); return;}

    // Allocate memory to store the output of derivs
    unsigned int n = sz[0]*sz[1]*sz[2];
    int size = n * sizeof(double);

    #include "bssnrhs_cuda_malloc.h"
    #include "bssnrhs_cuda_malloc_adv.h"

    bssn::timer::t_deriv.start();

    // Deriv calls are follows
    //#include "bssnrhs_cuda_derivs.h"
    #include "bssnrhs_cuda_derivs_adv.h"

    bssn::timer::t_deriv.stop();

//     int sizeArray=(sz[2]-3)*(sz[1]-3)*(sz[0]-3);
//     callculateBSSN_EQ(
//   grad_0_alpha,
//   grad_1_alpha,
//   grad_2_alpha,
//   grad_0_beta0,
//   grad_1_beta0,
//   grad_2_beta0,
//   grad_0_beta1,
//   grad_1_beta1,
//   grad_2_beta1,
//   grad_0_beta2,
//   grad_1_beta2,
//   grad_2_beta2,
//   grad_0_B0,
//   grad_1_B0,
//   grad_2_B0,
//   grad_0_B1,
//   grad_1_B1,
//   grad_2_B1,
//   grad_0_B2,
//   grad_1_B2,
//   grad_2_B2,
//   grad_0_chi,
//   grad_1_chi,
//   grad_2_chi,
//   grad_0_Gt0,
//   grad_1_Gt0,
//   grad_2_Gt0,
//   grad_0_Gt1,
//   grad_1_Gt1,
//   grad_2_Gt1,
//   grad_0_Gt2,
//   grad_1_Gt2,
//   grad_2_Gt2,
//   grad_0_K,
//   grad_1_K,
//   grad_2_K,
//   grad_0_gt0,
//   grad_1_gt0,
//   grad_2_gt0,
//   grad_0_gt1,
//   grad_1_gt1,
//   grad_2_gt1,
//   grad_0_gt2,
//   grad_1_gt2,
//   grad_2_gt2,
//   grad_0_gt3,
//   grad_1_gt3,
//   grad_2_gt3,
//   grad_0_gt4,
//   grad_1_gt4,
//   grad_2_gt4,
//   grad_0_gt5,
//   grad_1_gt5,
//   grad_2_gt5,
//   grad_0_At0,
//   grad_1_At0,
//   grad_2_At0,
//   grad_0_At1,
//   grad_1_At1,
//   grad_2_At1,
//   grad_0_At2,
//   grad_1_At2,
//   grad_2_At2,
//   grad_0_At3,
//   grad_1_At3,
//   grad_2_At3,
//   grad_0_At4,
//   grad_1_At4,
//   grad_2_At4,
//   grad_0_At5,
//   grad_1_At5,
//   grad_2_At5,
//   grad2_0_0_gt0,
//   grad2_0_1_gt0,
//   grad2_0_2_gt0,
//   grad2_1_1_gt0,
//   grad2_1_2_gt0,
//   grad2_2_2_gt0,
//   grad2_0_0_gt1,
//   grad2_0_1_gt1,
//   grad2_0_2_gt1,
//   grad2_1_1_gt1,
//   grad2_1_2_gt1,
//   grad2_2_2_gt1,
//   grad2_0_0_gt2,
//   grad2_0_1_gt2,
//   grad2_0_2_gt2,
//   grad2_1_1_gt2,
//   grad2_1_2_gt2,
//   grad2_2_2_gt2,
//   grad2_0_0_gt3,
//   grad2_0_1_gt3,
//   grad2_0_2_gt3,
//   grad2_1_1_gt3,
//   grad2_1_2_gt3,
//   grad2_2_2_gt3,
//   grad2_0_0_gt4,
//   grad2_0_1_gt4,
//   grad2_0_2_gt4,
//   grad2_1_1_gt4,
//   grad2_1_2_gt4,
//   grad2_2_2_gt4,
//   grad2_0_0_gt5,
//   grad2_0_1_gt5,
//   grad2_0_2_gt5,
//   grad2_1_1_gt5,
//   grad2_1_2_gt5,
//   grad2_2_2_gt5,
//   grad2_0_0_chi,
//   grad2_0_1_chi,
//   grad2_0_2_chi,
//   grad2_1_1_chi,
//   grad2_1_2_chi,
//   grad2_2_2_chi,
//   grad2_0_0_alpha,
//   grad2_0_1_alpha,
//   grad2_0_2_alpha,
//   grad2_1_1_alpha,
//   grad2_1_2_alpha,
//   grad2_2_2_alpha,
//   grad2_0_0_beta0,
//   grad2_0_1_beta0,
//   grad2_0_2_beta0,
//   grad2_1_1_beta0,
//   grad2_1_2_beta0,
//   grad2_2_2_beta0,
//   grad2_0_0_beta1,
//   grad2_0_1_beta1,
//   grad2_0_2_beta1,
//   grad2_1_1_beta1,
//   grad2_1_2_beta1,
//   grad2_2_2_beta1,
//   grad2_0_0_beta2,
//   grad2_0_1_beta2,
//   grad2_0_2_beta2,
//   grad2_1_1_beta2,
//   grad2_1_2_beta2,
//   grad2_2_2_beta2,
//   agrad_0_gt0,
//   agrad_1_gt0,
//   agrad_2_gt0,
//   agrad_0_gt1,
//   agrad_1_gt1,
//   agrad_2_gt1,
//   agrad_0_gt2,
//   agrad_1_gt2,
//   agrad_2_gt2,
//   agrad_0_gt3,
//   agrad_1_gt3,
//   agrad_2_gt3,
//   agrad_0_gt4,
//   agrad_1_gt4,
//   agrad_2_gt4,
//   agrad_0_gt5,
//   agrad_1_gt5,
//   agrad_2_gt5,
//   agrad_0_At0,
//   agrad_1_At0,
//   agrad_2_At0,
//   agrad_0_At1,
//   agrad_1_At1,
//   agrad_2_At1,
//   agrad_0_At2,
//   agrad_1_At2,
//   agrad_2_At2,
//   agrad_0_At3,
//   agrad_1_At3,
//   agrad_2_At3,
//   agrad_0_At4,
//   agrad_1_At4,
//   agrad_2_At4,
//   agrad_0_At5,
//   agrad_1_At5,
//   agrad_2_At5,
//   agrad_0_alpha,
//   agrad_1_alpha,
//   agrad_2_alpha,
//   agrad_0_beta0,
//   agrad_1_beta0,
//   agrad_2_beta0,
//   agrad_0_beta1,
//   agrad_1_beta1,
//   agrad_2_beta1,
//   agrad_0_beta2,
//   agrad_1_beta2,
//   agrad_2_beta2,
//   agrad_0_chi,
//   agrad_1_chi,
//   agrad_2_chi,
//   agrad_0_Gt0,
//   agrad_1_Gt0,
//   agrad_2_Gt0,
//   agrad_0_Gt1,
//   agrad_1_Gt1,
//   agrad_2_Gt1,
//   agrad_0_Gt2,
//   agrad_1_Gt2,
//   agrad_2_Gt2,
//   agrad_0_K,
//   agrad_1_K,
//   agrad_2_K,
//   agrad_0_B0,
//   agrad_1_B0,
//   agrad_2_B0,
//   agrad_0_B1,
//   agrad_1_B1,
//   agrad_2_B1,
//   agrad_0_B2,
//   agrad_1_B2,
//   agrad_2_B2,
//   dev_alphaInt,
//   dev_chiInt,
//   dev_KInt,
//   dev_gt0Int,
//   dev_gt1Int,
//   dev_gt2Int,
//   dev_gt3Int,
//   dev_gt4Int,
//   dev_gt5Int,
//   dev_beta0Int,
//   dev_beta1Int,
//   dev_beta2Int,
//   dev_At0Int,
//   dev_At1Int,
//   dev_At2Int,
//   dev_At3Int,
//   dev_At4Int,
//   dev_At5Int,
//   dev_Gt0Int,
//   dev_Gt1Int,
//   dev_Gt2Int,
//   dev_B0Int,
//   dev_B1Int,
//   dev_B2Int,
//   dev_lambda,
//   dev_lambda_f,
//   pmin,
//   dev_sz,
//   dev_dy_hx,
//   dev_dy_hy,
//   dev_dy_hz,
//   dev_var_in,
//   dev_var_out,
//   &sizeArray
// );
    // Free up GPU memory
    #include "bssnrhs_cuda_offset_demalloc.h"
    #include "bssnrhs_cuda_mdealloc.h"
    #include "bssnrhs_cuda_mdealloc_adv.h"
    cudaFree(dev_dy_hx);
    cudaFree(dev_dy_hy);
    cudaFree(dev_dy_hz);
    cudaFree(dev_sz);
    cudaFree(dev_zero);
    cudaFree(dev_pmin);
    cudaFree(dev_lambda);
    cudaFree(dev_lambda_f);
}

__global__ void cacl_bssn_bcs_x(double * output, double * dev_var_in, int* dev_u_offset,
    double *dxf, double *dyf, double *dzf,
    const double *pmin, const double f_falloff, const double f_asymptotic,
    const unsigned int *dev_sz, int* dev_bflag) {

        int x;
        int y;
        int z;

        if (*dev_bflag & (1u<<OCT_DIR_LEFT)) {
            x = pmin[0] + 3*((pmax[0] - pmin[0]) / (nx - 1));
            for (unsigned int k = kb; k < ke; k++) {
               z = pmin[2] + k*hz;
              for (unsigned int j = jb; j < je; j++) {
                 y = pmin[1] + j*hy;
                 pp = IDX(ib,j,k);
                 inv_r = 1.0 / sqrt(x*x + y*y + z*z);
        
                f_rhs[pp] = -  inv_r * (
                                 x * dxf[pp]
                               + y * dyf[pp]
                               + z * dzf[pp]
                               + f_falloff * (   f[pp] - f_asymptotic ) );
        
              }
            }
          }
        
          if (*dev_bflag & (1u<<OCT_DIR_RIGHT)) {
             x = pmin[0] + ie*hx;
            for (unsigned int k = kb; k < ke; k++) {
               z = pmin[2] + k*hz;
              for (unsigned int j = jb; j < je; j++) {
                 y = pmin[1] + j*hy;
                 pp = IDX(ie,j,k);
                 inv_r = 1.0 / sqrt(x*x + y*y + z*z);
        
                f_rhs[pp] = -  inv_r * (
                                 x * dxf[pp]
                               + y * dyf[pp]
                               + z * dzf[pp]
                               + f_falloff * (   f[pp] - f_asymptotic ) );
        
              }
            }
          }

}

__global__ void cacl_bssn_bcs_y(double * output, double * dev_var_in, int* dev_u_offset,
    double *dxf, double *dyf, double *dzf,
    const double *pmin, const double f_falloff, const double f_asymptotic,
    const unsigned int *host_sz, int* dev_bflag) {



}

__global__ void cacl_bssn_bcs_z(double * output, double * dev_var_in, int* dev_u_offset,
    double *dxf, double *dyf, double *dzf,
    const double *pmin, const double f_falloff, const double f_asymptotic,
    const unsigned int *host_sz, int* dev_bflag) {



}
void bssn_bcs(double * output, double * dev_var_in, int* dev_u_offset,
    double *dxf, double *dyf, double *dzf,
    const double *pmin,const double *pmax, const double f_falloff, const double f_asymptotic,
    const unsigned int *host_sz, int* dev_bflag, int* dev_sz) {

        const unsigned int nx = host_sz[0];
        const unsigned int ny = host_sz[1];
        const unsigned int nz = host_sz[2];

        double hx = (pmax[0] - pmin[0]) / (nx - 1);
        double hy = (pmax[1] - pmin[1]) / (ny - 1);
        double hz = (pmax[2] - pmin[2]) / (nz - 1);

        const int ie = nx - 3;//x direction
        const int je = ny - 3;//y direction
        const int ke = nz - 3;//z direction

        int temp_max = (ie>je)? ie : je;
        int maximumIterations = (temp_max>ke) ? temp_max: ke;
        
        int requiredBlocks = maximumIterations / 10;
        if (ie % 10 != 0 || je % 10 != 0 || ke % 10 != 0) {
            requiredBlocks++;
        }
        
        int threads_x = ie / requiredBlocks;
        int threads_y = je / requiredBlocks;
        int threads_z = ke / requiredBlocks;
        
        calc_co_deriv42_x <<< dim2(threads_x,threads_y), dim2(threads_x,threads_y) >>> (output, dev_var_in,
            dev_dx, dev_bflag, dev_sz, dev_u_offset);


    }