#include "rhs_cuda.h"

#include <iostream>
#include <stdio.h>

enum VAR_CU {U_ALPHA=0,U_CHI,U_K,U_GT0,U_GT1,U_GT2,U_BETA0,U_BETA1,U_BETA2,U_B0,U_B1,U_B2,U_SYMGT0,U_SYMGT1,U_SYMGT2,U_SYMGT3,U_SYMGT4,U_SYMGT5,U_SYMAT0,U_SYMAT1,U_SYMAT2,U_SYMAT3,U_SYMAT4,U_SYMAT5};

__global__ void example_kernal(double * val){
    // Eminda you can use this if it is required
    //test GPU mem values
}

void cuda_bssnrhs(double * dev_var_out, double * dev_var_in, const unsigned int unzip_dof, 
const unsigned int& offset, const double *pmin, const double *pmax, const unsigned int *sz, 
const unsigned int& bflag)
{
    cudaError_t cudaStatus;

    // Choose which GPU to run on, change this on a multi-GPU system.
     cudaStatus = cudaSetDevice(0);
     if (cudaStatus != cudaSuccess) {
         fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?\n");
         return;
     }

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

    double * dev_dy_hx;
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

    // Allocate memory to store the output of derivs
    unsigned int n = sz[0]*sz[1]*sz[2];
    #include "bssnrhs_cuda_malloc.h"

    bssn::timer::t_deriv.start();

    // Deriv calls are follows
    #include "bssnrhs_cuda_derivs.h"

    bssn::timer::t_deriv.stop();


    // Free up GPU memory
    // #include "bssnrhs_cuda_offset_demalloc.h"
    // #include "bssnrhs_cuda_mdealloc.h"
    // cudaFree(&dev_dy_hx);
    // cudaFree(&dev_dy_hy);
    // cudaFree(&dev_dy_hz);
    // cudaFree(&dev_sz);
}