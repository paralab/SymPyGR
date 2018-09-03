//
// Created by milinda on 8/10/18.
//
#include "rhs_cuda.cuh"
#include "../include/bssn_rhs_deriv_mem_cuda.h"

namespace cuda
{



    void computeRHS(double **unzipVarsRHS, const double **uZipVars,const cuda::_Block* blkList,unsigned int numBlocks)
    {

        cuda::profile::t_overall.start();


        cuda::profile::t_H2D_Comm.start();

            //get GPU information.
            // assumes the if there are multiple gpus per node all have the same specification.
            cuda::__CUDA_DEVICE_PROPERTIES=getGPUDeviceInfo(0);

            // device properties for the host
            cudaDeviceProp deviceProp;
            cudaGetDeviceProperties(&deviceProp,0);

            const double GPU_BLOCK_SHARED_MEM_UTIL=0.8;
            const unsigned int BSSN_NUM_VARS=24;
            const unsigned int BSSN_CONSTRAINT_NUM_VARS=6;

            //send blocks to the gpu
            cuda::__DENDRO_BLOCK_LIST=cuda::copyArrayToDevice(blkList,numBlocks);
            cuda::__DENDRO_NUM_BLOCKS=cuda::copyValueToDevice(&numBlocks);

            cuda::__BSSN_NUM_VARS=cuda::copyValueToDevice(&BSSN_NUM_VARS);
            cuda::__BSSN_CONSTRAINT_NUM_VARS=cuda::copyValueToDevice(&BSSN_CONSTRAINT_NUM_VARS);

            cuda::__GPU_BLOCK_SHARED_MEM_UTIL=cuda::copyValueToDevice(&GPU_BLOCK_SHARED_MEM_UTIL);

        cuda::profile::t_H2D_Comm.stop();

        unsigned int maxBlkSz_1d=0;
        for(unsigned int blk=0;blk<numBlocks;blk++)
        {
            const unsigned int* sz=blkList[blk].getSz();
            if(maxBlkSz_1d<(sz[0]*sz[1]*sz[2]))
                maxBlkSz_1d=sz[0]*sz[1]*sz[2];
        }

        const unsigned int derivSz=(maxBlkSz_1d*maxBlkSz_1d*maxBlkSz_1d);
        cuda::__DENDRO_BLK_MAX_SZ=cuda::copyValueToDevice(&derivSz);
        const size_t deriv_mem_sz= derivSz*(deviceProp.multiProcessorCount);

        cuda::profile::t_cudaMalloc_derivs.start();

            cuda::MemoryDerivs derivWorkSpace;
            derivWorkSpace.allocateDerivMemory(derivSz);

        cuda::profile::t_cudaMalloc_derivs.stop();



        dim3 blockGrid(numBlocks,1);
        dim3 threadBlock(8,8,8);


        cuda::profile::t_derivs.start();


        cuda::__computeDerivPass1 <<<blockGrid,threadBlock>>> (uZipVars,derivWorkSpace,cuda::__DENDRO_BLOCK_LIST,cuda::__CUDA_DEVICE_PROPERTIES);
        CUDA_CHECK_ERROR();
        cuda::__computeDerivPass2 <<<blockGrid,threadBlock>>> (uZipVars,derivWorkSpace,cuda::__DENDRO_BLOCK_LIST,cuda::__CUDA_DEVICE_PROPERTIES);
        CUDA_CHECK_ERROR();
        cuda::__computeDerivPass3 <<<blockGrid,threadBlock>>> (uZipVars,derivWorkSpace,cuda::__DENDRO_BLOCK_LIST,cuda::__CUDA_DEVICE_PROPERTIES);
        CUDA_CHECK_ERROR();
        cuda::__computeDerivPass4 <<<blockGrid,threadBlock>>> (uZipVars,derivWorkSpace,cuda::__DENDRO_BLOCK_LIST,cuda::__CUDA_DEVICE_PROPERTIES);
        CUDA_CHECK_ERROR();
        cuda::__computeDerivPass5 <<<blockGrid,threadBlock>>> (uZipVars,derivWorkSpace,cuda::__DENDRO_BLOCK_LIST,cuda::__CUDA_DEVICE_PROPERTIES);
        CUDA_CHECK_ERROR();

        cudaDeviceSynchronize();

        cuda::profile::t_derivs.stop();


        cuda::profile::t_cudaMalloc_derivs.start();

            derivWorkSpace.deallocateDerivMemory();

        cuda::profile::t_cudaMalloc_derivs.stop();

        cudaFree(cuda::__CUDA_DEVICE_PROPERTIES);
        cudaFree(cuda::__DENDRO_BLOCK_LIST);
        cudaFree(cuda::__DENDRO_NUM_BLOCKS);
        cudaFree(cuda::__BSSN_NUM_VARS);
        cudaFree(cuda::__BSSN_CONSTRAINT_NUM_VARS);
        cudaFree(cuda::__GPU_BLOCK_SHARED_MEM_UTIL);




        cuda::profile::t_overall.stop();


    }

}
