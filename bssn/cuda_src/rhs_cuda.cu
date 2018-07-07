#include "rhs_cuda.h"
#include "bssneqn_solve.h"

#include <iostream>
#include <stdio.h>


enum VAR_CU {U_ALPHA=0,U_CHI,U_K,U_GT0,U_GT1,U_GT2,U_BETA0,U_BETA1,U_BETA2,U_B0,U_B1,U_B2,U_SYMGT0,U_SYMGT1,U_SYMGT2,U_SYMGT3,U_SYMGT4,U_SYMGT5,U_SYMAT0,U_SYMAT1,U_SYMAT2,U_SYMAT3,U_SYMAT4,U_SYMAT5};

void cuda_bssnrhs(double * dev_var_out, double * dev_var_in, const unsigned int unzip_dof, 
const double * pmin, const double * pmax, const unsigned int * sz, 
const unsigned int& bflag, cudaStream_t stream,
#include "list_of_para.h"
)
{ 
    int alphaInt = (VAR_CU::U_ALPHA) * unzip_dof;
    int chiInt = (VAR_CU::U_CHI) * unzip_dof;
    int KInt = (VAR_CU::U_K) * unzip_dof;
    int gt0Int = (VAR_CU::U_SYMGT0) * unzip_dof;
    int gt1Int = (VAR_CU::U_SYMGT1) * unzip_dof;
    int gt2Int =  (VAR_CU::U_SYMGT2) * unzip_dof;
    int gt3Int = (VAR_CU::U_SYMGT3) * unzip_dof;
    int gt4Int = (VAR_CU::U_SYMGT4) * unzip_dof;
    int gt5Int = (VAR_CU::U_SYMGT5) * unzip_dof;
    int beta0Int = (VAR_CU::U_BETA0) * unzip_dof;
    int beta1Int = (VAR_CU::U_BETA1) * unzip_dof;
    int beta2Int = (VAR_CU::U_BETA2) * unzip_dof;
    int At0Int = (VAR_CU::U_SYMAT0) * unzip_dof;
    int At1Int = (VAR_CU::U_SYMAT1) * unzip_dof;
    int At2Int = (VAR_CU::U_SYMAT2) * unzip_dof;
    int At3Int = (VAR_CU::U_SYMAT3) * unzip_dof;
    int At4Int = (VAR_CU::U_SYMAT4) * unzip_dof;
    int At5Int = (VAR_CU::U_SYMAT5) * unzip_dof;
    int Gt0Int = (VAR_CU::U_GT0) * unzip_dof;
    int Gt1Int = (VAR_CU::U_GT1) * unzip_dof;
    int Gt2Int = (VAR_CU::U_GT2) * unzip_dof;
    int B0Int = (VAR_CU::U_B0) * unzip_dof;
    int B1Int = (VAR_CU::U_B1) * unzip_dof;
    int B2Int = (VAR_CU::U_B2) * unzip_dof;

    double hx = (pmax[0] - pmin[0]) / (sz[0] - 1);
    double hy = (pmax[1] - pmin[1]) / (sz[1] - 1);
    double hz = (pmax[2] - pmin[2]) / (sz[2] - 1);

    // Deriv calls are follows
    cuda_calc_all(dev_var_in, hx, hy, hz, sz[0], sz[1], sz[2], bflag, stream, 
        #include "list_of_args.h"
        ,
        #include "list_of_offset_args.h"
            );
        
    cuda_deriv_calc_all_adv(dev_var_in, hx, hy, hz, sz[0], sz[1], sz[2], bflag, stream,
        #include "list_of_args.h"
        ,
        #include "list_of_offset_args.h"
            );
    
    calc_bssn_eqns(dev_var_in, dev_var_out, sz, pmin, hz, hy, hx, stream,
    #include "list_of_offset_args.h"
    ,
    #include "list_of_args.h"
    );

    if (bflag != 0) {

        bssn_bcs(dev_var_out, dev_var_in, alphaInt, grad_0_alpha, grad_1_alpha, grad_2_alpha,
            pmin, pmax, 1.0, 1.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, chiInt, grad_0_chi, grad_1_chi, grad_2_chi,
            pmin, pmax, 1.0, 1.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, KInt, grad_0_K, grad_1_K, grad_2_K,
            pmin, pmax, 1.0, 0.0, sz, bflag, stream);

        bssn_bcs(dev_var_out, dev_var_in, beta0Int, grad_0_beta0, grad_1_beta0, grad_2_beta0,
            pmin, pmax, 1.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, beta1Int, grad_0_beta1, grad_1_beta1, grad_2_beta1,
            pmin, pmax, 1.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, beta2Int, grad_0_beta2, grad_1_beta2, grad_2_beta2,
            pmin, pmax, 1.0, 0.0, sz, bflag, stream);

        bssn_bcs(dev_var_out, dev_var_in, Gt0Int, grad_0_Gt0, grad_1_Gt0, grad_2_Gt0,
            pmin, pmax, 2.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, Gt1Int, grad_0_Gt1, grad_1_Gt1, grad_2_Gt1,
            pmin, pmax, 2.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, Gt2Int, grad_0_Gt2, grad_1_Gt2, grad_2_Gt2,
            pmin, pmax, 2.0, 0.0, sz, bflag, stream);

        bssn_bcs(dev_var_out, dev_var_in, B0Int, grad_0_B0, grad_1_B0, grad_2_B0,
            pmin, pmax, 1.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, B1Int, grad_0_B1, grad_1_B1, grad_2_B1,
            pmin, pmax, 1.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, B2Int, grad_0_B2, grad_1_B2, grad_2_B2,
            pmin, pmax, 1.0, 0.0, sz, bflag, stream);

        bssn_bcs(dev_var_out, dev_var_in, At0Int, grad_0_At0, grad_1_At0, grad_2_At0,
            pmin, pmax, 2.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, At1Int, grad_0_At1, grad_1_At1, grad_2_At1,
            pmin, pmax, 2.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, At2Int, grad_0_At2, grad_1_At2, grad_2_At2,
            pmin, pmax, 2.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, At3Int, grad_0_At3, grad_1_At3, grad_2_At3,
            pmin, pmax, 2.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, At4Int, grad_0_At4, grad_1_At4, grad_2_At4,
            pmin, pmax, 2.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, At5Int, grad_0_At5, grad_1_At5, grad_2_At5,
            pmin, pmax, 2.0, 0.0, sz, bflag, stream); 

        bssn_bcs(dev_var_out, dev_var_in, gt0Int, grad_0_gt0, grad_1_gt0, grad_2_gt0,
            pmin, pmax, 1.0, 1.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, gt1Int, grad_0_gt1, grad_1_gt1, grad_2_gt1,
            pmin, pmax, 1.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, gt2Int, grad_0_gt2, grad_1_gt2, grad_2_gt2,
            pmin, pmax, 1.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, gt3Int, grad_0_gt3, grad_1_gt3, grad_2_gt3,
            pmin, pmax, 1.0, 1.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, gt4Int, grad_0_gt4, grad_1_gt4, grad_2_gt4,
            pmin, pmax, 1.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, gt5Int, grad_0_gt5, grad_1_gt5, grad_2_gt5,
            pmin, pmax, 1.0, 1.0, sz, bflag, stream);
        
    }

    calc_ko_deriv_all(dev_var_in, hx, hy, hz, sz[0], sz[1], sz[2], bflag, stream,
        #include "list_of_args.h"
        ,
        #include "list_of_offset_args.h"
        );

    get_output(dev_var_out, sz, stream,
        #include "list_of_offset_args.h"
        ,
        #include "list_of_args.h"
    );

}

__global__ void cacl_bssn_bcs_x(double * dev_var_out, double * dev_var_in, 
    int u_offset,
    double * dxf, double * dyf, double * dzf,
    double pmin_x, double pmin_y, double pmin_z, double pmax_x, double pmax_y, double pmax_z,
    const double f_falloff, const double f_asymptotic,
    const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z,
    int bflag) {

        int j = 3 + threadIdx.x + blockIdx.x * blockDim.x;
        int k = 3 + threadIdx.y + blockIdx.y * blockDim.y;
        int nx = host_sz_x;
        int ny = host_sz_y;
        int nz = host_sz_z;

        if(j >= ny-3 || k >= nz-3) return;

        double inv_r;
        double hx = (pmax_x - pmin_x) / (nx - 1);
        double hy = (pmax_y - pmin_y) / (ny - 1);
        double hz = (pmax_z - pmin_z) / (nz - 1);
        double x, y, z;
        int pp;

        if (bflag & (1u<<OCT_DIR_LEFT)) {
            
            x = pmin_x + 3*hx;
            z = pmin_z + k*hz;
            y = pmin_y + j*hy;
            pp = IDX(3,j,k);
            inv_r = 1.0 / sqrt(x*x + y*y + z*z);
   
            dev_var_out[u_offset + pp] = -  inv_r * (
                            x * dxf[pp]
                          + y * dyf[pp]
                          + z * dzf[pp]
                          + f_falloff * (   dev_var_in[u_offset + pp] - f_asymptotic ) );
          }
        
          if (bflag & (1u<<OCT_DIR_RIGHT)) {
             x = pmin_x + (nx - 3)*hx;
             z = pmin_z + k*hz;
             y = pmin_y + j*hy;
             pp = IDX((nx - 3),j,k);
             inv_r = 1.0 / sqrt(x*x + y*y + z*z);
    
             dev_var_out[u_offset + pp] = -  inv_r * (
                             x * dxf[pp]
                           + y * dyf[pp]
                           + z * dzf[pp]
                           + f_falloff * (   dev_var_in[u_offset + pp] - f_asymptotic ) );
          }

}

__global__ void cacl_bssn_bcs_y(double * dev_var_out, double * dev_var_in, 
    int u_offset,
    double * dxf, double * dyf, double * dzf,
    double pmin_x, double pmin_y, double pmin_z, double pmax_x, double pmax_y, double pmax_z,
    const double f_falloff, const double f_asymptotic,
    const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z,
    int bflag) {

        int i = 3 + threadIdx.x + blockIdx.x * blockDim.x;
        int k = 3 + threadIdx.y + blockIdx.y * blockDim.y;
        int nx = host_sz_x;
        int ny = host_sz_y;
        int nz = host_sz_z;

        if(i >= nx-3 || k >= nz-3) return;

        double inv_r;
        double hx = (pmax_x - pmin_x) / (nx - 1);
        double hy = (pmax_y - pmin_y) / (ny - 1);
        double hz = (pmax_z - pmin_z) / (nz - 1);
        double x, y, z;
        int pp;

        if (bflag & (1u<<OCT_DIR_DOWN)) {
            
            y = pmin_y + 3*hy;
            z = pmin_z + k*hz;
            x = pmin_x + i*hx;
            pp = IDX(i,3,k);
            inv_r = 1.0 / sqrt(x*x + y*y + z*z);
   
            dev_var_out[u_offset + pp] = -  inv_r * (
                            x * dxf[pp]
                          + y * dyf[pp]
                          + z * dzf[pp]
                          + f_falloff * (   dev_var_in[u_offset + pp] - f_asymptotic ) );
            
          }
        
          if (bflag & (1u<<OCT_DIR_UP)) {
             x = pmin_x + i*hx;
             z = pmin_z + k*hz;
             y = pmin_y + (ny-3)*hy;
             pp = IDX(i,(ny - 3),k);
             inv_r = 1.0 / sqrt(x*x + y*y + z*z);
    
             dev_var_out[u_offset + pp] = -  inv_r * (
                             x * dxf[pp]
                           + y * dyf[pp]
                           + z * dzf[pp]
                           + f_falloff * (   dev_var_in[u_offset + pp] - f_asymptotic ) );
               
          }
}

__global__ void cacl_bssn_bcs_z(double * dev_var_out, double * dev_var_in, 
    int u_offset,
    double * dxf, double * dyf, double * dzf,
    double pmin_x, double pmin_y, double pmin_z, double pmax_x, double pmax_y, double pmax_z,
    const double f_falloff, const double f_asymptotic,
    const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z,
    int bflag) {
        
            int i = 3 + threadIdx.x + blockIdx.x * blockDim.x;
            int j = 3 + threadIdx.y + blockIdx.y * blockDim.y;
            int nx = host_sz_x;
            int ny = host_sz_y;
            int nz = host_sz_z;

            if(i >= nx-3 || j >= ny-3) return;

            double inv_r;
            double hx = (pmax_x - pmin_x) / (nx - 1);
            double hy = (pmax_y - pmin_y) / (ny - 1);
            double hz = (pmax_z - pmin_z) / (nz - 1);
            double x, y, z;
            int pp;

            if (bflag & (1u<<OCT_DIR_BACK)) {
            
            y = pmin_y + j*hy;
            z = pmin_z + 3*hz;
            x = pmin_x + i*hx;
            pp = IDX(i,j,3);
            inv_r = 1.0 / sqrt(x*x + y*y + z*z);
   
            dev_var_out[u_offset + pp] = -  inv_r * (
                            x * dxf[pp]
                          + y * dyf[pp]
                          + z * dzf[pp]
                          + f_falloff * (   dev_var_in[u_offset + pp] - f_asymptotic ) );
            
          }
        
          if (bflag & (1u<<OCT_DIR_FRONT)) {
            x = pmin_x + i*hx;
            z = pmin_z + (nz-3)*hz;
            y = pmin_y + j*hy;
            pp = IDX(i,j,3);
            inv_r = 1.0 / sqrt(x*x + y*y + z*z);
    
            dev_var_out[u_offset + pp] = -  inv_r * (
                             x * dxf[pp]
                           + y * dyf[pp]
                           + z * dzf[pp]
                           + f_falloff * (   dev_var_in[u_offset + pp] - f_asymptotic ) );
               
          }
}


// dev_var_out, dev_var_in, dev_At5Int, grad_0_At5, grad_1_At5, grad_2_At5,
//             dev_pmin, dev_pmax, 2.0, 0.0, sz, dev_bflag, dev_sz, stream


void bssn_bcs(double * dev_var_out, double * dev_var_in, 
    int u_offset, double * dxf, double * dyf, double * dzf,
    const double * pmin, const double * pmax, const double f_falloff, const double f_asymptotic,
    const unsigned int * host_sz, int bflag, cudaStream_t stream) {
        
        const unsigned int nx = host_sz[0];
        const unsigned int ny = host_sz[1];
        const unsigned int nz = host_sz[2];

        const int ie = nx - 3;//x direction
        const int je = ny - 3;//y direction
        const int ke = nz - 3;//z direction

        double pmin_x = pmin[0];
        double pmin_y = pmin[1];
        double pmin_z = pmin[2];

        double pmax_x = pmax[0];
        double pmax_y = pmax[1];
        double pmax_z = pmax[2];

        const unsigned int host_sz_x = host_sz[0];
        const unsigned int host_sz_y = host_sz[1];
        const unsigned int host_sz_z = host_sz[2];

        int maximumIterations = (je>ke) ? je: ke;
        
        int requiredBlocks = (9 + maximumIterations) / 10;
        
        int threads_y = (requiredBlocks-1+je) / requiredBlocks;
        int threads_z = (requiredBlocks-1+ke) / requiredBlocks;
        
        cacl_bssn_bcs_x <<< dim3(threads_y,threads_z), dim3(threads_y,threads_z), 0, stream >>> (
            dev_var_out, dev_var_in, 
            u_offset, dxf, dyf, dzf, 
            pmin_x, pmin_y, pmin_z, pmax_x, pmax_y, pmax_z, 
            f_falloff, f_asymptotic, 
            host_sz_x, host_sz_y, host_sz_z, 
            bflag );
        
        CHECK_ERROR(cudaGetLastError(), "cacl_bssn_bcs_x Kernel launch failed");
           
        maximumIterations = (ke>ie) ? ke : ie ;
        requiredBlocks = (9 + maximumIterations)/10;
        int threads_x = (requiredBlocks-1+ie) / requiredBlocks;
        threads_z = (requiredBlocks-1+ke) / requiredBlocks;
        cacl_bssn_bcs_y <<< dim3(threads_x,threads_z), dim3(threads_x,threads_z), 0, stream >>> (
            dev_var_out, dev_var_in, 
            u_offset, dxf, dyf, dzf, 
            pmin_x, pmin_y, pmin_z, pmax_x, pmax_y, pmax_z, 
            f_falloff, f_asymptotic, 
            host_sz_x, host_sz_y, host_sz_z, 
            bflag );
 
        CHECK_ERROR(cudaGetLastError(), "cacl_bssn_bcs_y Kernel launch failed");

        maximumIterations = (je>ie) ? je : ie ;
        requiredBlocks = (9 + maximumIterations)/10;
        threads_x = (requiredBlocks-1+ie) / requiredBlocks;
        threads_y = (requiredBlocks-1+je) / requiredBlocks;
        cacl_bssn_bcs_z <<< dim3(threads_x,threads_y), dim3(threads_x,threads_y), 0, stream >>> (
            dev_var_out, dev_var_in, 
            u_offset, dxf, dyf, dzf, 
            pmin_x, pmin_y, pmin_z, pmax_x, pmax_y, pmax_z, 
            f_falloff, f_asymptotic, 
            host_sz_x, host_sz_y, host_sz_z, 
            bflag );

        CHECK_ERROR(cudaGetLastError(), "cacl_bssn_bcs_z Kernel launch failed");
    }

__global__ void kernal_get_output (double * dev_var_out, 
    const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z,
    #include "list_of_offset_para.h"
    ,
    #include "list_of_para.h"
) 
{
    int thread_id = blockIdx.x*threads_per_block + threadIdx.x;

    for (int id = thread_id*thread_load_output; id<(thread_id+1)*thread_load_output; id++){
                

        int i = id%(host_sz_x-6) + 3;
        int j = ((id/(host_sz_x-6))%(host_sz_y-6)) + 3;
        int k = (id/(host_sz_y-6)/(host_sz_x-6)) + 3;

        int nx = host_sz_x;
        int ny = host_sz_y;

        if(i >= nx-3 || j >= ny-3 || k >= host_sz_z-3) return;

        const  double sigma = 1e-4;
        int pp = i + nx*(j + ny*k);

        dev_var_out[alphaInt + pp] += sigma * (grad_0_alpha[pp] + grad_1_alpha[pp] + grad_2_alpha[pp]);
        dev_var_out[beta0Int + pp] += sigma * (grad_0_beta0[pp] + grad_1_beta0[pp] + grad_2_beta0[pp]);
        dev_var_out[beta1Int + pp] += sigma * (grad_0_beta1[pp] + grad_1_beta1[pp] + grad_2_beta1[pp]);
        dev_var_out[beta2Int + pp] += sigma * (grad_0_beta2[pp] + grad_1_beta2[pp] + grad_2_beta2[pp]);

        dev_var_out[gt0Int + pp] += sigma * (grad_0_gt0[pp] + grad_1_gt0[pp] + grad_2_gt0[pp]);
        dev_var_out[gt1Int + pp] += sigma * (grad_0_gt1[pp] + grad_1_gt1[pp] + grad_2_gt1[pp]);
        dev_var_out[gt2Int + pp] += sigma * (grad_0_gt2[pp] + grad_1_gt2[pp] + grad_2_gt2[pp]);
        dev_var_out[gt3Int + pp] += sigma * (grad_0_gt3[pp] + grad_1_gt3[pp] + grad_2_gt3[pp]);
        dev_var_out[gt4Int + pp] += sigma * (grad_0_gt4[pp] + grad_1_gt4[pp] + grad_2_gt4[pp]);
        dev_var_out[gt5Int + pp] += sigma * (grad_0_gt5[pp] + grad_1_gt5[pp] + grad_2_gt5[pp]);

        dev_var_out[chiInt + pp]  += sigma * (grad_0_chi[pp] + grad_1_chi[pp] + grad_2_chi[pp]);

        dev_var_out[At0Int + pp] += sigma * (grad_0_At0[pp] + grad_1_At0[pp] + grad_2_At0[pp]);
        dev_var_out[At1Int + pp] += sigma * (grad_0_At1[pp] + grad_1_At1[pp] + grad_2_At1[pp]);
        dev_var_out[At2Int + pp] += sigma * (grad_0_At2[pp] + grad_1_At2[pp] + grad_2_At2[pp]);
        dev_var_out[At3Int + pp] += sigma * (grad_0_At3[pp] + grad_1_At3[pp] + grad_2_At3[pp]);
        dev_var_out[At4Int + pp] += sigma * (grad_0_At4[pp] + grad_1_At4[pp] + grad_2_At4[pp]);
        dev_var_out[At5Int + pp] += sigma * (grad_0_At5[pp] + grad_1_At5[pp] + grad_2_At5[pp]);

        dev_var_out[KInt + pp] += sigma * (grad_0_K[pp] + grad_1_K[pp] + grad_2_K[pp]);
        
        dev_var_out[Gt0Int + pp] += sigma * (grad_0_Gt0[pp] + grad_1_Gt0[pp] + grad_2_Gt0[pp]);
        dev_var_out[Gt1Int + pp] += sigma * (grad_0_Gt1[pp] + grad_1_Gt1[pp] + grad_2_Gt1[pp]);
        dev_var_out[Gt2Int + pp] += sigma * (grad_0_Gt2[pp] + grad_1_Gt2[pp] + grad_2_Gt2[pp]);

        dev_var_out[B0Int + pp] += sigma * (grad_0_B0[pp] + grad_1_B0[pp] + grad_2_B0[pp]);
        dev_var_out[B1Int + pp] += sigma * (grad_0_B1[pp] + grad_1_B1[pp] + grad_2_B1[pp]);
        dev_var_out[B2Int + pp] += sigma * (grad_0_B2[pp] + grad_1_B2[pp] + grad_2_B2[pp]);
    }
}

void get_output (double * dev_var_out, const unsigned int * host_sz, cudaStream_t stream,
    #include "list_of_offset_para.h"
    ,
    #include "list_of_para.h"
) 
{
        const int ie = host_sz[0] - 3;//x direction
        const int je = host_sz[1] - 3;//y direction
        const int ke = host_sz[2] - 3;//z direction

        const unsigned int host_sz_x = host_sz[0];
        const unsigned int host_sz_y = host_sz[1];
        const unsigned int host_sz_z = host_sz[2];

        int total_points = ceil(1.0*ie*je*ke/thread_load_output);
        int blocks = ceil(1.0*total_points/threads_per_block);

        kernal_get_output <<< blocks, threads_per_block, 0, stream >>> (dev_var_out, 
                    host_sz_x, host_sz_y, host_sz_z,
                    #include "list_of_offset_args.h"
                    ,
                    #include "list_of_args.h"
                    );
        
        CHECK_ERROR(cudaGetLastError(), "kernal_get_output Kernel launch failed");
}
