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
    void computeRHS(double **unzipVarsRHS, const double **uZipVars,const cuda::_Block* blkList,unsigned int numBlocks,const cuda::BSSNComputeParams* bssnPars,std::vector<unsigned int >& blockMap,dim3 gridDim,dim3 blockDim);


    /**
     * @brief Asynchronous version of the computeRHS() function call.
     * */
     void computeRHSAsync(double **OUTPUT_REFERENCE, double **INPUT_REFERENCE, 
        cuda::_Block* DENDRO_BLOCK_LIST, unsigned int numBlocks, cuda::BSSNComputeParams* bssnPars,
        std::vector<unsigned int >& blockMap,dim3 gridDim,dim3 blockDim,unsigned int numStreams, cudaStream_t stream);




}

#endif //SFCSORTBENCH_CUDARHS_H
