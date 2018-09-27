//
// Created by milinda on 8/10/18.
//
#include "rhs_cuda.cuh"


namespace cuda
{



    void computeRHS(double **unzipVarsRHS, const double **uZipVars,const cuda::_Block* blkList,unsigned int numBlocks,const cuda::BSSNComputeParams* bssnPars)
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

            const unsigned int UNZIP_DOF_SZ=blkList[numBlocks-1].getOffset()+(blkList[numBlocks-1].getSz()[0]*blkList[numBlocks-1].getSz()[1]*blkList[numBlocks-1].getSz()[2]);


            //send blocks to the gpu
            cuda::__DENDRO_BLOCK_LIST=cuda::copyArrayToDevice(blkList,numBlocks);
            cuda::__DENDRO_NUM_BLOCKS=cuda::copyValueToDevice(&numBlocks);

            cuda::__BSSN_NUM_VARS=cuda::copyValueToDevice(&BSSN_NUM_VARS);
            cuda::__BSSN_CONSTRAINT_NUM_VARS=cuda::copyValueToDevice(&BSSN_CONSTRAINT_NUM_VARS);

            cuda::__GPU_BLOCK_SHARED_MEM_UTIL=cuda::copyValueToDevice(&GPU_BLOCK_SHARED_MEM_UTIL);

            //allocate memory for unzip vectors
            cuda::__UNZIP_INPUT=cuda::alloc2DCudaArray<double>(uZipVars,BSSN_NUM_VARS,UNZIP_DOF_SZ);
            cuda::__UNZIP_OUTPUT=cuda::alloc2DCudaArray<double>(BSSN_NUM_VARS,UNZIP_DOF_SZ);

            cuda::__BSSN_COMPUTE_PARMS=cuda::copyValueToDevice(&(*bssnPars));


        cuda::profile::t_H2D_Comm.stop();

        unsigned int maxBlkSz=0;
        for(unsigned int blk=0;blk<numBlocks;blk++)
        {
            const unsigned int* sz=blkList[blk].getSz();
            if(maxBlkSz<(sz[0]*sz[1]*sz[2]))
                maxBlkSz=sz[0]*sz[1]*sz[2];
        }

        const unsigned int derivSz=(maxBlkSz);
        cuda::__DENDRO_BLK_MAX_SZ=cuda::copyValueToDevice(&derivSz);
        //const size_t deriv_mem_sz= derivSz*(deviceProp.multiProcessorCount);
        const unsigned int numSM=deviceProp.multiProcessorCount;

        //std::cout<<"deriv alloc begin"<<std::endl;

        cuda::profile::t_cudaMalloc_derivs.start();

            cuda::MemoryDerivs derivWorkSpace;
            derivWorkSpace.allocateDerivMemory(maxBlkSz,numSM);
            CUDA_CHECK_ERROR();

            cuda::__BSSN_DERIV_WORKSPACE=cuda::copyValueToDevice(&derivWorkSpace);
            CUDA_CHECK_ERROR();

        cuda::profile::t_cudaMalloc_derivs.stop();



        dim3 blockGrid(numBlocks,1);
        dim3 threadBlock(2,2,2);


        cuda::profile::t_derivs.start();


        cuda::__RSWS_computeDerivs <<<blockGrid,threadBlock>>> ((const double**)cuda::__UNZIP_INPUT,cuda::__BSSN_DERIV_WORKSPACE,cuda::__DENDRO_BLOCK_LIST,cuda::__CUDA_DEVICE_PROPERTIES);
        CUDA_CHECK_ERROR();

        cudaDeviceSynchronize();
        CUDA_CHECK_ERROR();

        cuda::profile::t_derivs.stop();

        threadBlock=dim3(6,6,6);
        cuda::profile::t_rhs.start();

        /*cuda::__compute_a_rhs<<<blockGrid,threadBlock>>>(cuda::__UNZIP_OUTPUT,(const double**)cuda::__UNZIP_INPUT,cuda::__BSSN_DERIV_WORKSPACE,cuda::__DENDRO_BLOCK_LIST,cuda::__BSSN_COMPUTE_PARMS,cuda::__CUDA_DEVICE_PROPERTIES);
        CUDA_CHECK_ERROR();

        cuda::__compute_b_rhs<<<blockGrid,threadBlock>>>(cuda::__UNZIP_OUTPUT,(const double**)cuda::__UNZIP_INPUT,cuda::__BSSN_DERIV_WORKSPACE,cuda::__DENDRO_BLOCK_LIST,cuda::__BSSN_COMPUTE_PARMS,cuda::__CUDA_DEVICE_PROPERTIES);
        CUDA_CHECK_ERROR();

        cuda::__compute_gt_rhs<<<blockGrid,threadBlock>>>(cuda::__UNZIP_OUTPUT,(const double**)cuda::__UNZIP_INPUT,cuda::__BSSN_DERIV_WORKSPACE,cuda::__DENDRO_BLOCK_LIST,cuda::__BSSN_COMPUTE_PARMS,cuda::__CUDA_DEVICE_PROPERTIES);
        CUDA_CHECK_ERROR();

        cuda::__compute_chi_rhs<<<blockGrid,threadBlock>>>(cuda::__UNZIP_OUTPUT,(const double**)cuda::__UNZIP_INPUT,cuda::__BSSN_DERIV_WORKSPACE,cuda::__DENDRO_BLOCK_LIST,cuda::__BSSN_COMPUTE_PARMS,cuda::__CUDA_DEVICE_PROPERTIES);
        CUDA_CHECK_ERROR();

        cuda::__compute_At_rhs<<<blockGrid,threadBlock>>>(cuda::__UNZIP_OUTPUT,(const double**)cuda::__UNZIP_INPUT,cuda::__BSSN_DERIV_WORKSPACE,cuda::__DENDRO_BLOCK_LIST,cuda::__BSSN_COMPUTE_PARMS,cuda::__CUDA_DEVICE_PROPERTIES);
        CUDA_CHECK_ERROR();

        cuda::__compute_K_rhs<<<blockGrid,threadBlock>>>(cuda::__UNZIP_OUTPUT,(const double**)cuda::__UNZIP_INPUT,cuda::__BSSN_DERIV_WORKSPACE,cuda::__DENDRO_BLOCK_LIST,cuda::__BSSN_COMPUTE_PARMS,cuda::__CUDA_DEVICE_PROPERTIES);
        CUDA_CHECK_ERROR();

        cuda::__compute_Gt_rhs<<<blockGrid,threadBlock>>>(cuda::__UNZIP_OUTPUT,(const double**)cuda::__UNZIP_INPUT,cuda::__BSSN_DERIV_WORKSPACE,cuda::__DENDRO_BLOCK_LIST,cuda::__BSSN_COMPUTE_PARMS,cuda::__CUDA_DEVICE_PROPERTIES);
        CUDA_CHECK_ERROR();

        cuda::__compute_B_rhs<<<blockGrid,threadBlock>>>(cuda::__UNZIP_OUTPUT,(const double**)cuda::__UNZIP_INPUT,cuda::__BSSN_DERIV_WORKSPACE,cuda::__DENDRO_BLOCK_LIST,cuda::__BSSN_COMPUTE_PARMS,cuda::__CUDA_DEVICE_PROPERTIES);
        CUDA_CHECK_ERROR();

        cudaDeviceSynchronize();
        CUDA_CHECK_ERROR();*/

        cuda::profile::t_rhs.stop();

        cuda::profile::t_cudaMalloc_derivs.start();
            derivWorkSpace.deallocateDerivMemory();
            CUDA_CHECK_ERROR();
        cuda::profile::t_cudaMalloc_derivs.stop();

        cudaFree(cuda::__CUDA_DEVICE_PROPERTIES);
        cudaFree(cuda::__DENDRO_BLOCK_LIST);
        cudaFree(cuda::__DENDRO_NUM_BLOCKS);
        cudaFree(cuda::__BSSN_NUM_VARS);
        cudaFree(cuda::__BSSN_CONSTRAINT_NUM_VARS);
        cudaFree(cuda::__GPU_BLOCK_SHARED_MEM_UTIL);

        cuda::dealloc2DCudaArray(cuda::__UNZIP_INPUT,BSSN_NUM_VARS);
        cuda::dealloc2DCudaArray(cuda::__UNZIP_OUTPUT,BSSN_NUM_VARS);



        cuda::profile::t_overall.stop();


    }

}
