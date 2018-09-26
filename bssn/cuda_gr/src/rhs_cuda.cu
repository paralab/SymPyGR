//
// Created by milinda on 8/10/18.
//
#include "rhs_cuda.cuh"


namespace cuda
{



    void computeRHS(double **unzipVarsRHS, const double **uZipVars,const cuda::_Block* blkList,
            unsigned int numBlocks, cudaStream_t stream, cuda::_Block* blockListReference, 
            double** referenceToInput, cuda::MemoryDerivs* derivPointer, 
            cudaDeviceProp* cudaDeviceProperties, cuda::BSSNComputeParams* bSSNComputeParams,
            double** outputReference)
    {
        cuda::profile::t_overall.start();

        cuda::profile::t_H2D_Comm.start();

            //get GPU information.
            cuda::__CUDA_DEVICE_PROPERTIES = cudaDeviceProperties;
            cuda::__UNZIP_OUTPUT = outputReference;
            cuda::__BSSN_COMPUTE_PARMS = bSSNComputeParams;
            cuda::__DENDRO_BLOCK_LIST = blockListReference;
            cuda::__UNZIP_INPUT = referenceToInput;
            cuda::__BSSN_DERIV_WORKSPACE = derivPointer;

        cuda::profile::t_H2D_Comm.stop();

        unsigned int maxBlkSz=0;
        for(unsigned int blk=0;blk<numBlocks;blk++)
        {
            const unsigned int* sz=blkList[blk].getSz();
            if(maxBlkSz<(sz[0]*sz[1]*sz[2]))
                maxBlkSz=sz[0]*sz[1]*sz[2];
        }

        dim3 blockGrid(numBlocks,1);
        dim3 threadBlock(32,4,1);

        cuda::profile::t_derivs.start();

        cuda::__RSWS_computeDerivs <<<blockGrid,threadBlock>>> ((const double**)cuda::__UNZIP_INPUT,cuda::__BSSN_DERIV_WORKSPACE,cuda::__DENDRO_BLOCK_LIST,cuda::__CUDA_DEVICE_PROPERTIES);
        CUDA_CHECK_ERROR();

        // cudaDeviceSynchronize();
        // CUDA_CHECK_ERROR();

        cuda::profile::t_derivs.stop();

        threadBlock=dim3(6,6,6);
        cuda::profile::t_rhs.start();

        cuda::__compute_a_rhs<<<blockGrid,threadBlock>>>(cuda::__UNZIP_OUTPUT,(const double**)cuda::__UNZIP_INPUT,cuda::__BSSN_DERIV_WORKSPACE,cuda::__DENDRO_BLOCK_LIST,cuda::__BSSN_COMPUTE_PARMS,cuda::__CUDA_DEVICE_PROPERTIES);
        CUDA_CHECK_ERROR();

        cuda::__compute_b_rhs<<<blockGrid,threadBlock, 0, stream>>>(cuda::__UNZIP_OUTPUT,(const double**)cuda::__UNZIP_INPUT,cuda::__BSSN_DERIV_WORKSPACE,cuda::__DENDRO_BLOCK_LIST,cuda::__BSSN_COMPUTE_PARMS,cuda::__CUDA_DEVICE_PROPERTIES);
        CUDA_CHECK_ERROR();

        cuda::__compute_gt_rhs<<<blockGrid,threadBlock, 0, stream>>>(cuda::__UNZIP_OUTPUT,(const double**)cuda::__UNZIP_INPUT,cuda::__BSSN_DERIV_WORKSPACE,cuda::__DENDRO_BLOCK_LIST,cuda::__BSSN_COMPUTE_PARMS,cuda::__CUDA_DEVICE_PROPERTIES);
        CUDA_CHECK_ERROR();

        cuda::__compute_chi_rhs<<<blockGrid,threadBlock, 0, stream>>>(cuda::__UNZIP_OUTPUT,(const double**)cuda::__UNZIP_INPUT,cuda::__BSSN_DERIV_WORKSPACE,cuda::__DENDRO_BLOCK_LIST,cuda::__BSSN_COMPUTE_PARMS,cuda::__CUDA_DEVICE_PROPERTIES);
        CUDA_CHECK_ERROR();

        cuda::__compute_At_rhs<<<blockGrid,threadBlock, 0, stream>>>(cuda::__UNZIP_OUTPUT,(const double**)cuda::__UNZIP_INPUT,cuda::__BSSN_DERIV_WORKSPACE,cuda::__DENDRO_BLOCK_LIST,cuda::__BSSN_COMPUTE_PARMS,cuda::__CUDA_DEVICE_PROPERTIES);
        CUDA_CHECK_ERROR();

        cuda::__compute_K_rhs<<<blockGrid,threadBlock, 0, stream>>>(cuda::__UNZIP_OUTPUT,(const double**)cuda::__UNZIP_INPUT,cuda::__BSSN_DERIV_WORKSPACE,cuda::__DENDRO_BLOCK_LIST,cuda::__BSSN_COMPUTE_PARMS,cuda::__CUDA_DEVICE_PROPERTIES);
        CUDA_CHECK_ERROR();

        cuda::__compute_Gt_rhs<<<blockGrid,threadBlock, 0, stream>>>(cuda::__UNZIP_OUTPUT,(const double**)cuda::__UNZIP_INPUT,cuda::__BSSN_DERIV_WORKSPACE,cuda::__DENDRO_BLOCK_LIST,cuda::__BSSN_COMPUTE_PARMS,cuda::__CUDA_DEVICE_PROPERTIES);
        CUDA_CHECK_ERROR();

        cuda::__compute_B_rhs<<<blockGrid, threadBlock, 0, stream>>>(cuda::__UNZIP_OUTPUT,(const double**)cuda::__UNZIP_INPUT,cuda::__BSSN_DERIV_WORKSPACE,cuda::__DENDRO_BLOCK_LIST,cuda::__BSSN_COMPUTE_PARMS,cuda::__CUDA_DEVICE_PROPERTIES);
        CUDA_CHECK_ERROR();

        cudaDeviceSynchronize();
        CUDA_CHECK_ERROR();

        cuda::profile::t_rhs.stop();

        // cuda::profile::t_cudaMalloc_derivs.start();
        //     derivWorkSpace.deallocateDerivMemory();
        //     CUDA_CHECK_ERROR();
        // cuda::profile::t_cudaMalloc_derivs.stop();

        // cudaFree(cuda::__CUDA_DEVICE_PROPERTIES);
        // cudaFree(cuda::__DENDRO_BLOCK_LIST);

        // cuda::copyBackMemory(cuda::__UNZIP_INPUT,BSSN_NUM_VARS, stream);

        cuda::profile::t_overall.stop();
    }

}
