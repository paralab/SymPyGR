/**
 * Created on: Sep 21, 2018
 * 		Author: Akila, Eranga, Eminda, Ruwan
 **/
 
#include "getOutput.cuh"

__global__ void calc_get_output (double * dev_var_out, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z,
    #include "para_derivs_offsets.h"
) 
{
    int thread_id = blockIdx.x*1024 + threadIdx.x;

    int i = thread_id%(host_sz_x-6) + 3;
    int j = ((thread_id/(host_sz_x-6))%(host_sz_y-6)) + 3;
    int k = (thread_id/(host_sz_y-6)/(host_sz_x-6)) + 3;

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

void get_output_kernel_wrapper(double * dev_var_out, const unsigned int * host_sz, cudaStream_t stream,
    #include "para_derivs_offsets.h"
) 
{
    const int ie = host_sz[0] - 3;//x direction
    const int je = host_sz[1] - 3;//y direction
    const int ke = host_sz[2] - 3;//z direction

    const unsigned int host_sz_x = host_sz[0];
    const unsigned int host_sz_y = host_sz[1];
    const unsigned int host_sz_z = host_sz[2];

    int total_points = ceil(1.0*ie*je*ke);
    int blocks = ceil(1.0*total_points/1024);

    calc_get_output <<< blocks, 1024, 0, stream >>> (dev_var_out, 
                host_sz_x, host_sz_y, host_sz_z,
                #include "args_derivs_offsets.h"
                );
    
    CHECK_ERROR(cudaGetLastError(), "kernal_get_output Kernel launch failed");
}