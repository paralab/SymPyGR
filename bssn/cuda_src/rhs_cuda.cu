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

    bssn::timer::t_deriv_gpu.start();

    #include "bssnrhs_cuda_malloc.h"
    #include "bssnrhs_cuda_malloc_adv.h"
   

    // Deriv calls are follows
    #include "bssnrhs_cuda_derivs.h"
    #include "bssnrhs_cuda_derivs_adv.h"

    bssn::timer::t_deriv_gpu.stop();


    bssn::timer::t_rhs_gpu.start();
    calc_bssn_eqns(sz, dev_sz, dev_pmin, dev_dy_hz, dev_dy_hy, dev_dy_hx, dev_var_in, dev_var_out,
        #include "list_of_args.h"
    );
    bssn::timer::t_rhs_gpu.stop();

        
    if (bflag != 0) {
        bssn::timer::t_bdyc_gpu.start();

        bssn_bcs(dev_var_out, dev_var_in, dev_alphaInt, grad_0_alpha, grad_1_alpha, grad_2_alpha,
            dev_pmin, dev_pmax, 1.0, 1.0, sz, dev_bflag, dev_sz);
        bssn_bcs(dev_var_out, dev_var_in, dev_chiInt, grad_0_chi, grad_1_chi, grad_2_chi,
            dev_pmin, dev_pmax, 1.0, 1.0, sz, dev_bflag, dev_sz);
        bssn_bcs(dev_var_out, dev_var_in, dev_KInt, grad_0_K, grad_1_K, grad_2_K,
            dev_pmin, dev_pmax, 1.0, 0.0, sz, dev_bflag, dev_sz);

        bssn_bcs(dev_var_out, dev_var_in, dev_beta0Int, grad_0_beta0, grad_1_beta0, grad_2_beta0,
            dev_pmin, dev_pmax, 1.0, 0.0, sz, dev_bflag, dev_sz);
        bssn_bcs(dev_var_out, dev_var_in, dev_beta1Int, grad_0_beta1, grad_1_beta1, grad_2_beta1,
            dev_pmin, dev_pmax, 1.0, 0.0, sz, dev_bflag, dev_sz);
        bssn_bcs(dev_var_out, dev_var_in, dev_beta2Int, grad_0_beta2, grad_1_beta2, grad_2_beta2,
            dev_pmin, dev_pmax, 1.0, 0.0, sz, dev_bflag, dev_sz);

        bssn_bcs(dev_var_out, dev_var_in, dev_Gt0Int, grad_0_Gt0, grad_1_Gt0, grad_2_Gt0,
            dev_pmin, dev_pmax, 2.0, 0.0, sz, dev_bflag, dev_sz);
        bssn_bcs(dev_var_out, dev_var_in, dev_Gt1Int, grad_0_Gt1, grad_1_Gt1, grad_2_Gt1,
            dev_pmin, dev_pmax, 2.0, 0.0, sz, dev_bflag, dev_sz);
        bssn_bcs(dev_var_out, dev_var_in, dev_Gt2Int, grad_0_Gt2, grad_1_Gt2, grad_2_Gt2,
            dev_pmin, dev_pmax, 2.0, 0.0, sz, dev_bflag, dev_sz);

        bssn_bcs(dev_var_out, dev_var_in, dev_B0Int, grad_0_B0, grad_1_B0, grad_2_B0,
            dev_pmin, dev_pmax, 1.0, 0.0, sz, dev_bflag, dev_sz);
        bssn_bcs(dev_var_out, dev_var_in, dev_B1Int, grad_0_B1, grad_1_B1, grad_2_B1,
            dev_pmin, dev_pmax, 1.0, 0.0, sz, dev_bflag, dev_sz);
        bssn_bcs(dev_var_out, dev_var_in, dev_B2Int, grad_0_B2, grad_1_B2, grad_2_B2,
            dev_pmin, dev_pmax, 1.0, 0.0, sz, dev_bflag, dev_sz);

        bssn_bcs(dev_var_out, dev_var_in, dev_At0Int, grad_0_At0, grad_1_At0, grad_2_At0,
            dev_pmin, dev_pmax, 2.0, 0.0, sz, dev_bflag, dev_sz);
        bssn_bcs(dev_var_out, dev_var_in, dev_At1Int, grad_0_At1, grad_1_At1, grad_2_At1,
            dev_pmin, dev_pmax, 2.0, 0.0, sz, dev_bflag, dev_sz);
        bssn_bcs(dev_var_out, dev_var_in, dev_At2Int, grad_0_At2, grad_1_At2, grad_2_At2,
            dev_pmin, dev_pmax, 2.0, 0.0, sz, dev_bflag, dev_sz);
        bssn_bcs(dev_var_out, dev_var_in, dev_At3Int, grad_0_At3, grad_1_At3, grad_2_At3,
            dev_pmin, dev_pmax, 2.0, 0.0, sz, dev_bflag, dev_sz);
        bssn_bcs(dev_var_out, dev_var_in, dev_At4Int, grad_0_At4, grad_1_At4, grad_2_At4,
            dev_pmin, dev_pmax, 2.0, 0.0, sz, dev_bflag, dev_sz);
        bssn_bcs(dev_var_out, dev_var_in, dev_At5Int, grad_0_At5, grad_1_At5, grad_2_At5,
            dev_pmin, dev_pmax, 2.0, 0.0, sz, dev_bflag, dev_sz); 

        bssn_bcs(dev_var_out, dev_var_in, dev_gt0Int, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                dev_pmin, dev_pmax, 1.0, 1.0, sz, dev_bflag, dev_sz);
        bssn_bcs(dev_var_out, dev_var_in, dev_gt1Int, grad_0_gt1, grad_1_gt1, grad_2_gt1,
            dev_pmin, dev_pmax, 1.0, 0.0, sz, dev_bflag, dev_sz);
        bssn_bcs(dev_var_out, dev_var_in, dev_gt2Int, grad_0_gt2, grad_1_gt2, grad_2_gt2,
            dev_pmin, dev_pmax, 1.0, 0.0, sz, dev_bflag, dev_sz);
        bssn_bcs(dev_var_out, dev_var_in, dev_gt3Int, grad_0_gt3, grad_1_gt3, grad_2_gt3,
            dev_pmin, dev_pmax, 1.0, 1.0, sz, dev_bflag, dev_sz);
        bssn_bcs(dev_var_out, dev_var_in, dev_gt4Int, grad_0_gt4, grad_1_gt4, grad_2_gt4,
            dev_pmin, dev_pmax, 1.0, 0.0, sz, dev_bflag, dev_sz);
        bssn_bcs(dev_var_out, dev_var_in, dev_gt5Int, grad_0_gt5, grad_1_gt5, grad_2_gt5,
            dev_pmin, dev_pmax, 1.0, 1.0, sz, dev_bflag, dev_sz);
          

        bssn::timer::t_bdyc_gpu.stop();
    }

    bssn::timer::t_deriv_gpu.start();
    #include "bssnrhs_cuda_ko_derivs.h"
    bssn::timer::t_deriv_gpu.stop();

    #if test && 0
    // Copying specified array to CPU for testing purpose
    double * host_array_cpu = (double *) malloc(size);
    #include "test_GPU_derivs.h" // only one of both at a time
    #include "test_GPU_adv_derivs.h"
    free(host_array_cpu); 
    #endif

    bssn::timer::t_rhs_gpu.start();
    get_output(dev_var_out, dev_sz, sz,
        #include "list_of_args.h"
    );
    bssn::timer::t_rhs_gpu.stop();

    bssn::timer::t_deriv_gpu.start();
    #include "bssnrhs_cuda_offset_demalloc.h"
    #include "bssnrhs_cuda_mdealloc.h"
    #include "bssnrhs_cuda_mdealloc_adv.h"
    cudaFree(dev_dy_hx);
    cudaFree(dev_dy_hy);
    cudaFree(dev_dy_hz);
    cudaFree(dev_sz);
    cudaFree(dev_zero);
    cudaFree(dev_pmin);
    bssn::timer::t_deriv_gpu.stop();

}

__global__ void cacl_bssn_bcs_x(double * output, double * dev_var_in, int* dev_u_offset,
    double *dxf, double *dyf, double *dzf,
    double *pmin, double *pmax, const double f_falloff, const double f_asymptotic,
    int *dev_sz, int* dev_bflag) {

        int j = 3 + threadIdx.x + blockIdx.x * blockDim.x;
        int k = 3 + threadIdx.y + blockIdx.y * blockDim.y;
        int nx = dev_sz[0];
        int ny = dev_sz[1];
        int nz = dev_sz[2];

        if(j >= ny-3 || k >= nz-3) return;

        double inv_r;
        double hx = (pmax[0] - pmin[0]) / (nx - 1);
        double hy = (pmax[1] - pmin[1]) / (ny - 1);
        double hz = (pmax[2] - pmin[2]) / (nz - 1);
        double x, y, z;
        int pp;

        if (*dev_bflag & (1u<<OCT_DIR_LEFT)) {
            
            x = pmin[0] + 3*hx;
            z = pmin[2] + k*hz;
            y = pmin[1] + j*hy;
            pp = IDX(3,j,k);
            inv_r = 1.0 / sqrt(x*x + y*y + z*z);
   
            output[*dev_u_offset + pp] = -  inv_r * (
                            x * dxf[pp]
                          + y * dyf[pp]
                          + z * dzf[pp]
                          + f_falloff * (   dev_var_in[*dev_u_offset + pp] - f_asymptotic ) );
          }
        
          if (*dev_bflag & (1u<<OCT_DIR_RIGHT)) {
             x = pmin[0] + (nx - 3)*hx;
             z = pmin[2] + k*hz;
             y = pmin[1] + j*hy;
             pp = IDX((nx - 3),j,k);
             inv_r = 1.0 / sqrt(x*x + y*y + z*z);
    
             output[*dev_u_offset + pp] = -  inv_r * (
                             x * dxf[pp]
                           + y * dyf[pp]
                           + z * dzf[pp]
                           + f_falloff * (   dev_var_in[*dev_u_offset + pp] - f_asymptotic ) );
          }

}

__global__ void cacl_bssn_bcs_y(double * output, double * dev_var_in, int* dev_u_offset,
    double *dxf, double *dyf, double *dzf,
    double *pmin, double *pmax, const double f_falloff, const double f_asymptotic,
    int *dev_sz, int* dev_bflag) {

        int i = 3 + threadIdx.x + blockIdx.x * blockDim.x;
        int k = 3 + threadIdx.y + blockIdx.y * blockDim.y;
        int nx = dev_sz[0];
        int ny = dev_sz[1];
        int nz = dev_sz[2];

        if(i >= nx-3 || k >= nz-3) return;

        double inv_r;
        double hx = (pmax[0] - pmin[0]) / (nx - 1);
        double hy = (pmax[1] - pmin[1]) / (ny - 1);
        double hz = (pmax[2] - pmin[2]) / (nz - 1);
        double x, y, z;
        int pp;

        if (*dev_bflag & (1u<<OCT_DIR_DOWN)) {
            
            y = pmin[1] + 3*hy;
            z = pmin[2] + k*hz;
            x = pmin[0] + i*hx;
            pp = IDX(i,3,k);
            inv_r = 1.0 / sqrt(x*x + y*y + z*z);
   
            output[*dev_u_offset + pp] = -  inv_r * (
                            x * dxf[pp]
                          + y * dyf[pp]
                          + z * dzf[pp]
                          + f_falloff * (   dev_var_in[*dev_u_offset + pp] - f_asymptotic ) );
            
          }
        
          if (*dev_bflag & (1u<<OCT_DIR_UP)) {
             x = pmin[0] + i*hx;
             z = pmin[2] + k*hz;
             y = pmin[1] + (ny-3)*hy;
             pp = IDX(i,(ny - 3),k);
             inv_r = 1.0 / sqrt(x*x + y*y + z*z);
    
             output[*dev_u_offset + pp] = -  inv_r * (
                             x * dxf[pp]
                           + y * dyf[pp]
                           + z * dzf[pp]
                           + f_falloff * (   dev_var_in[*dev_u_offset + pp] - f_asymptotic ) );
               
          }
}

__global__ void cacl_bssn_bcs_z(double * output, double * dev_var_in, int* dev_u_offset,
    double *dxf, double *dyf, double *dzf,
    double *pmin, double *pmax, const double f_falloff, const double f_asymptotic,
    int *dev_sz, int* dev_bflag) {
        
            int i = 3 + threadIdx.x + blockIdx.x * blockDim.x;
            int j = 3 + threadIdx.y + blockIdx.y * blockDim.y;
            int nx = dev_sz[0];
            int ny = dev_sz[1];
            int nz = dev_sz[2];

            if(i >= nx-3 || j >= ny-3) return;

            double inv_r;
            double hx = (pmax[0] - pmin[0]) / (nx - 1);
            double hy = (pmax[1] - pmin[1]) / (ny - 1);
            double hz = (pmax[2] - pmin[2]) / (nz - 1);
            double x, y, z;
            int pp;

            if (*dev_bflag & (1u<<OCT_DIR_BACK)) {
            
            y = pmin[1] + j*hy;
            z = pmin[2] + 3*hz;
            x = pmin[0] + i*hx;
            pp = IDX(i,j,3);
            inv_r = 1.0 / sqrt(x*x + y*y + z*z);
   
            output[*dev_u_offset + pp] = -  inv_r * (
                            x * dxf[pp]
                          + y * dyf[pp]
                          + z * dzf[pp]
                          + f_falloff * (   dev_var_in[*dev_u_offset + pp] - f_asymptotic ) );
            
          }
        
          if (*dev_bflag & (1u<<OCT_DIR_FRONT)) {
            x = pmin[0] + i*hx;
            z = pmin[2] + (nz-3)*hz;
            y = pmin[1] + j*hy;
            pp = IDX(i,j,3);
            inv_r = 1.0 / sqrt(x*x + y*y + z*z);
    
            output[*dev_u_offset + pp] = -  inv_r * (
                             x * dxf[pp]
                           + y * dyf[pp]
                           + z * dzf[pp]
                           + f_falloff * (   dev_var_in[*dev_u_offset + pp] - f_asymptotic ) );
               
          }
}

void bssn_bcs(double * output, double * dev_var_in, int* dev_u_offset,
    double *dxf, double *dyf, double *dzf,
    double *pmin, double *pmax, const double f_falloff, const double f_asymptotic,
    const unsigned int *host_sz, int* dev_bflag, int* dev_sz) {
        
        cudaError_t cudaStatus;
        const unsigned int nx = host_sz[0];
        const unsigned int ny = host_sz[1];
        const unsigned int nz = host_sz[2];

        const int ie = nx - 3;//x direction
        const int je = ny - 3;//y direction
        const int ke = nz - 3;//z direction

        int maximumIterations = (je>ke) ? je: ke;
        
        int requiredBlocks = (9 + maximumIterations) / 10;
        
        int threads_y = (requiredBlocks-1+je) / requiredBlocks;
        int threads_z = (requiredBlocks-1+ke) / requiredBlocks;
        
        cacl_bssn_bcs_x <<< dim3(threads_y,threads_z), dim3(threads_y,threads_z) >>> (output, dev_var_in,
           dev_u_offset, dxf, dyf, dzf, pmin, pmax, f_falloff, f_asymptotic, dev_sz, dev_bflag );

        cudaStatus = cudaDeviceSynchronize();
           if (cudaStatus != cudaSuccess) {
                   fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching bssn_bcs_x kernal!\n", cudaStatus);
                   return;
           }
        
        maximumIterations = (ke>ie) ? ke : ie ;
        requiredBlocks = (9 + maximumIterations)/10;
        int threads_x = (requiredBlocks-1+ie) / requiredBlocks;
        threads_z = (requiredBlocks-1+ke) / requiredBlocks;
        cacl_bssn_bcs_y <<< dim3(threads_x,threads_z), dim3(threads_x,threads_z) >>> (output, dev_var_in,
            dev_u_offset, dxf, dyf, dzf, pmin, pmax, f_falloff, f_asymptotic, dev_sz, dev_bflag );
 
        cudaStatus = cudaDeviceSynchronize();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching bssn_bcs_y kernal!\n", cudaStatus);
            return;
        }

        maximumIterations = (je>ie) ? je : ie ;
        requiredBlocks = (9 + maximumIterations)/10;
        threads_x = (requiredBlocks-1+ie) / requiredBlocks;
        threads_y = (requiredBlocks-1+je) / requiredBlocks;
        cacl_bssn_bcs_z <<< dim3(threads_x,threads_y), dim3(threads_x,threads_y) >>> (output, dev_var_in,
            dev_u_offset, dxf, dyf, dzf, pmin, pmax, f_falloff, f_asymptotic, dev_sz, dev_bflag );
 
        cudaStatus = cudaDeviceSynchronize();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching bssn_bcs_z kernal!\n", cudaStatus);
            return;
        }
    }

    __global__ void kernal_get_output (double * output, int * dev_sz, 
        #include "list_of_para.h"
    ) 
    {

        int i = 3 + threadIdx.x + blockIdx.x * blockDim.x;
        int j = 3 + threadIdx.y + blockIdx.y * blockDim.y;
        int k = 3 + threadIdx.z + blockIdx.z * blockDim.z;

        int nx = dev_sz[0];
        int ny = dev_sz[1];

        if(i >= nx-3 || j >= ny-3 || k >= dev_sz[2]-3) return;

        const  double sigma = 1e-4;
        int pp = i + nx*(j + ny*k);

        output[*dev_alphaInt + pp] += sigma * (grad_0_alpha[pp] + grad_1_alpha[pp] + grad_2_alpha[pp]);
        output[*dev_beta0Int + pp] += sigma * (grad_0_beta0[pp] + grad_1_beta0[pp] + grad_2_beta0[pp]);
        output[*dev_beta1Int + pp] += sigma * (grad_0_beta1[pp] + grad_1_beta1[pp] + grad_2_beta1[pp]);
        output[*dev_beta2Int + pp] += sigma * (grad_0_beta2[pp] + grad_1_beta2[pp] + grad_2_beta2[pp]);

        output[*dev_gt0Int + pp] += sigma * (grad_0_gt0[pp] + grad_1_gt0[pp] + grad_2_gt0[pp]);
        output[*dev_gt1Int + pp] += sigma * (grad_0_gt1[pp] + grad_1_gt1[pp] + grad_2_gt1[pp]);
        output[*dev_gt2Int + pp] += sigma * (grad_0_gt2[pp] + grad_1_gt2[pp] + grad_2_gt2[pp]);
        output[*dev_gt3Int + pp] += sigma * (grad_0_gt3[pp] + grad_1_gt3[pp] + grad_2_gt3[pp]);
        output[*dev_gt4Int + pp] += sigma * (grad_0_gt4[pp] + grad_1_gt4[pp] + grad_2_gt4[pp]);
        output[*dev_gt5Int + pp] += sigma * (grad_0_gt5[pp] + grad_1_gt5[pp] + grad_2_gt5[pp]);

        output[*dev_chiInt + pp]  += sigma * (grad_0_chi[pp] + grad_1_chi[pp] + grad_2_chi[pp]);

        output[*dev_At0Int + pp] += sigma * (grad_0_At0[pp] + grad_1_At0[pp] + grad_2_At0[pp]);
        output[*dev_At1Int + pp] += sigma * (grad_0_At1[pp] + grad_1_At1[pp] + grad_2_At1[pp]);
        output[*dev_At2Int + pp] += sigma * (grad_0_At2[pp] + grad_1_At2[pp] + grad_2_At2[pp]);
        output[*dev_At3Int + pp] += sigma * (grad_0_At3[pp] + grad_1_At3[pp] + grad_2_At3[pp]);
        output[*dev_At4Int + pp] += sigma * (grad_0_At4[pp] + grad_1_At4[pp] + grad_2_At4[pp]);
        output[*dev_At5Int + pp] += sigma * (grad_0_At5[pp] + grad_1_At5[pp] + grad_2_At5[pp]);

        output[*dev_KInt + pp] += sigma * (grad_0_K[pp] + grad_1_K[pp] + grad_2_K[pp]);
        
        output[*dev_Gt0Int + pp] += sigma * (grad_0_Gt0[pp] + grad_1_Gt0[pp] + grad_2_Gt0[pp]);
        output[*dev_Gt1Int + pp] += sigma * (grad_0_Gt1[pp] + grad_1_Gt1[pp] + grad_2_Gt1[pp]);
        output[*dev_Gt2Int + pp] += sigma * (grad_0_Gt2[pp] + grad_1_Gt2[pp] + grad_2_Gt2[pp]);

        output[*dev_B0Int + pp] += sigma * (grad_0_B0[pp] + grad_1_B0[pp] + grad_2_B0[pp]);
        output[*dev_B1Int + pp] += sigma * (grad_0_B1[pp] + grad_1_B1[pp] + grad_2_B1[pp]);
        output[*dev_B2Int + pp] += sigma * (grad_0_B2[pp] + grad_1_B2[pp] + grad_2_B2[pp]);

    }

    void get_output (double* output, int* dev_sz, const unsigned int* host_sz, 
        #include "list_of_para.h"
    ) 
    {
            const int ie = host_sz[0] - 3;//x direction
            const int je = host_sz[1] - 3;//y direction
            const int ke = host_sz[2] - 3;//z direction
  
            int temp_max = (ie>je)? ie : je;
            int maximumIterations = (temp_max>ke) ? temp_max: ke;
            
            int requiredBlocks = (9+maximumIterations) / 10;

            kernal_get_output <<< dim3(requiredBlocks, requiredBlocks, requiredBlocks),
                      dim3((ie + requiredBlocks -1)/requiredBlocks,
                      (je + requiredBlocks -1)/requiredBlocks, 
                      (ke + requiredBlocks -1)/requiredBlocks) >>> (output, dev_sz, 
                        #include "list_of_args.h"
                      );
    }
