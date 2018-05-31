#include "bssneqn_solve.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <iostream>

using namespace std;

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
    int i = 3 + threadIdx.x + blockIdx.x * blockDim.x;
    int j = 3 + threadIdx.y + blockIdx.y * blockDim.y;
    int k = 3 + threadIdx.z + blockIdx.z * blockDim.z;

    int nx = dev_sz[0]; 
    int ny = dev_sz[1]; 
    if(i >= nx-3 || j >= ny-3 || k >= dev_sz[2]-3) return;

    double z = dev_pmin[2] + *dev_dy_hz*k;
    double y = dev_pmin[1] + *dev_dy_hy*j;
    double x = dev_pmin[0] + *dev_dy_hx*i;

    int pp = i + nx*(j + ny*k);
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
)
{
    const int ie = sz[0] - 3;//x direction
    const int je = sz[1] - 3;//y direction
    const int ke = sz[2] - 3;//z direction
  
    int temp_max = (ie>je)? ie : je;
    int maximumIterations = (temp_max>ke) ? temp_max: ke;
     
    int requiredBlocks = (5+maximumIterations) / 6;

    cuda_bssn_eqns_points <<< dim3(requiredBlocks, requiredBlocks, requiredBlocks),
                      dim3((ie + requiredBlocks -1)/requiredBlocks,
                      (je + requiredBlocks -1)/requiredBlocks, 
                      (ke + requiredBlocks -1)/requiredBlocks) >>> (0, dev_sz, dev_pmin, dev_dy_hz, dev_dy_hy, dev_dy_hx, dev_var_in, dev_var_out,
                        #include "list_of_args.h"
                    );
    CHECK_ERROR(cudaGetLastError(), "cuda_bssn_eqns_points Kernel launch failed");
    
}
