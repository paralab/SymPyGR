#include "bssneqn_solve.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <iostream>

using namespace std;

int threads_per_block_cpu=250;
int blocks_cpu=50;

__constant__ int threads_per_block=250;
__constant__ int blocks=50;

__constant__ double ETA_CONST=0.1;
__constant__ double ETA_R0=0.1;
__constant__ double ETA_DAMPING_EXP=0.1;
__constant__ unsigned int lambda[4]={1,2,3,4};
__constant__ double lambda_f[2]={0.8,0.9};

__global__ void cuda_bssn_eqns_points(int * dev_offset, int * dev_sz, double * dev_pmin, double * dev_dy_hz, double * dev_dy_hy, double * dev_dy_hx, 
    double * dev_var_in, double * dev_var_out,
    #include "list_of_para.h"
    )
{
    int id = *dev_offset + blockIdx.x*threads_per_block + threadIdx.x;

    int i = id%(dev_sz[0]-6) + 3;
    int j = ((id/(dev_sz[0]-6))%(dev_sz[1]-6)) + 3;
    int k = (id/(dev_sz[2]-6)/(dev_sz[1]-6)) + 3;

    if (k>=dev_sz[2]-3) return;

    double z = dev_pmin[2] + *dev_dy_hz*k;
    double y = dev_pmin[1] + *dev_dy_hy*j;
    double x = dev_pmin[0] + *dev_dy_hx*i;

    int pp = i + (dev_sz[0])*(j + (dev_sz[1])*k);
    double r_coord = sqrt(x*x + y*y + z*z);
    double eta = ETA_CONST;
    if (r_coord >= ETA_R0) {
        eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
    }
    #include "cuda_bssneqs.h"
}

void calc_bssn_eqns(const unsigned int * sz, int * dev_sz, double * dev_pmin, double * dev_dy_hz, double * dev_dy_hy, double * dev_dy_hx, 
double * dev_var_in, double * dev_var_out, 
#include "list_of_para.h"
, cudaStream_t stream
)
{
    int total_points = (sz[2]-6)*(sz[1]-6)*(sz[0]-6);

    int points_at_once = threads_per_block_cpu*blocks_cpu;
    int loops = ceil(1.0*total_points/points_at_once);

    cudaError_t cudaStatus;
    // full
    for(int i=0; i<loops; i++){
        int offset = i*points_at_once;

        int * dev_offset;
        cudaStatus = cudaMalloc((void **) &dev_offset, sizeof(int));
        if (cudaStatus != cudaSuccess) {fprintf(stderr, "dev_offset cudaMalloc failed!\n"); return;}

        cudaStatus = cudaMemcpyAsync(dev_offset, &offset, sizeof(int), cudaMemcpyHostToDevice, stream);
        if (cudaStatus != cudaSuccess) {fprintf(stderr, "dev_offset cudaMemcpy failed!\n"); return;}

        cuda_bssn_eqns_points<<< blocks_cpu, threads_per_block_cpu, 0, stream >>>(dev_offset, dev_sz, dev_pmin, dev_dy_hz, dev_dy_hy, dev_dy_hx, dev_var_in, dev_var_out,
            #include "list_of_args.h"
        );
        
        // Check for any errors launching the kernel
        cudaError_t cudaStatus;
        cudaStatus = cudaGetLastError();
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cuda_bssn_eqns_points Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
            return;
        }
        cudaStatus = cudaFree(dev_offset);
        if (cudaStatus != cudaSuccess) { fprintf(stderr, "dev_offset cudaFree failed: %s\n", cudaGetErrorString(cudaStatus));};
    } 
    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching cuda_bssn_eqns_points kernals!\n", cudaStatus);
        return;
    }
}
