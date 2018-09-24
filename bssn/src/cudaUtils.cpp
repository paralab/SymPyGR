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


    cudaDeviceProp* getGPUDeviceInfo(unsigned int device)
    {

        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp,device);

        cudaDeviceProp* __deviceProp;
        cudaMalloc(&__deviceProp,sizeof(cudaDeviceProp));
        CUDA_CHECK_ERROR();

        cudaMemcpy(__deviceProp,&deviceProp,sizeof(cudaDeviceProp),cudaMemcpyHostToDevice);
        CUDA_CHECK_ERROR();

        return __deviceProp;

    }







}