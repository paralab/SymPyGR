//
// Created by milinda on 8/9/18.
//
/**
 * @author Milinda Fernando
 * School of Computing, University of Utah
 * @brief Contains utility function for the host related to GPUs
 * */
#include "cudaUtils.h"

namespace cuda
{


    cudaDeviceProp* getGPUDeviceInfo(unsigned int device, cudaStream_t stream)
    {

        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp,device);

        cudaDeviceProp* __deviceProp;
        cudaMalloc(&__deviceProp,sizeof(cudaDeviceProp));
        CUDA_CHECK_ERROR();

        cudaMemcpyAsync(__deviceProp,&deviceProp,sizeof(cudaDeviceProp),cudaMemcpyHostToDevice, stream);
        CUDA_CHECK_ERROR();

        return __deviceProp;

    }







}