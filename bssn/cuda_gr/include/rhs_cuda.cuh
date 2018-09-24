// * Created by milinda on 8/10/18
/**
 * @author Milinda Fernando.
 * School of Computing, University of Utah.
 * @brief Computation of the rhs using cuda
 * Created by milinda on 8/10/18
 * */

#ifndef SFCSORTBENCH_CUDARHS_H
#define SFCSORTBENCH_CUDARHS_H

#include "cuda_runtime.h"
#include <device_launch_parameters.h>
#include "block_cu.h"
#include "bssn_rhs_deriv_mem_cuda.h"
#include "cudaUtils.h"
#include "params_cu.h"
#include "rhs_bssn.cuh"
#include "rhs_bssn.cuh"
#include "profile_gpu.h"
#include "cudaUtils.cuh"





namespace cuda
{

    /***
     * @brief performs kernel pre-launch tasks and launch the bssnrhs kernel
     *
     **/
     void computeRHS(double **unzipVarsRHS, 
     const double **uZipVars,const cuda::_Block* blkList,unsigned int numBlocks, cudaStream_t stream, cuda::_Block* blockListReference,
     double** tmp2D, double** referenceToInput,cuda::MemoryDerivs derivWorkSpace, cuda::MemoryDerivs* derivPointer, cudaDeviceProp* cudaDeviceProperties,
     cuda::BSSNComputeParams* bSSNComputeParams, double** outputReference);


}

#endif //SFCSORTBENCH_CUDARHS_H
