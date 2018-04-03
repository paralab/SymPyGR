#include "rhs_cuda.h"
#include "bssneqn_solve.h"

#include <iostream>
#include <stdio.h>


enum VAR_CU {U_ALPHA=0,U_CHI,U_K,U_GT0,U_GT1,U_GT2,U_BETA0,U_BETA1,U_BETA2,U_B0,U_B1,U_B2,U_SYMGT0,U_SYMGT1,U_SYMGT2,U_SYMGT3,U_SYMGT4,U_SYMGT5,U_SYMAT0,U_SYMAT1,U_SYMAT2,U_SYMAT3,U_SYMAT4,U_SYMAT5};

void cuda_bssnrhs(double * dev_var_out, double * dev_var_in, const unsigned int unzip_dof, 
const unsigned int& offset, const double *pmin, const double *pmax, const unsigned int *sz, 
const unsigned int& bflag)
{
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
    cudaError_t cudaStatus;
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

    double * dev_pmin;
    cudaStatus = cudaMalloc((void **) &dev_pmin, 3*sizeof(double));
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "pmin cudaMalloc failed!\n"); return;}
    cudaStatus = cudaMemcpy(dev_pmin, pmin, 3*sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "pmin cudaMemcpy failed!\n"); return;}

    double *dev_pmax;
    cudaStatus = cudaMalloc((void **) &dev_pmax, sizeof(pmax)*sizeof(double));
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "pmax cudaMalloc failed!\n"); return;}
    cudaStatus = cudaMemcpy(dev_pmax, pmax, sizeof(pmax)*sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "pmax cudaMemcpy failed!\n"); return;}

    // Allocate memory to store the output of derivs
    unsigned int n = sz[0]*sz[1]*sz[2];
    int size = n * sizeof(double);

    bssn::timer::t_deriv.start();

    #include "bssnrhs_cuda_malloc.h"
    #include "bssnrhs_cuda_malloc_adv.h"
   

    // Deriv calls are follows
    #include "bssnrhs_cuda_derivs.h"
    #include "bssnrhs_cuda_derivs_adv.h"

    bssn::timer::t_deriv.stop();


    bssn::timer::t_rhs.start();
    calc_bssn_eqns(sz, dev_sz, dev_pmin, dev_dy_hz, dev_dy_hy, dev_dy_hx, dev_var_in, dev_var_out,
        #include "list_of_args.h"
    );
    bssn::timer::t_rhs.stop();

    #if test
    // Copying specified array to CPU for testing purpose
    double * host_array_cpu = (double *) malloc(size);
    cudaStatus = cudaMemcpy(host_array_cpu, grad_1_alpha, size, cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "TEST: host_array_cpu cudaMemcpy from GPU to CPU failed!\n"); return;}
    test_file_write::writeToFile("output_cuda.txt", host_array_cpu, n);
    free(host_array_cpu);

    // double * host_array_cpu = (double *) malloc(unzip_dof*24);
    // cudaStatus = cudaMemcpy(host_array_cpu, dev_var_out, unzip_dof*24, cudaMemcpyDeviceToHost);
    // if (cudaStatus != cudaSuccess) {fprintf(stderr, "TEST: host_array_cpu cudaMemcpy from GPU to CPU failed!\n"); return;}
    // test_file_write::writeToFile("output_cuda.txt", host_array_cpu, unzip_dof*24);
    // free(host_array_cpu);
    #endif

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
