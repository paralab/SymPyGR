//
// Created by milinda on 8/9/18.
//
/**
 * @author Milinda Fernando
 * School of Computing, University of Utah
 * @brief Contains utility function for the host related to GPUs
 * */

#include "cuda_runtime.h"
#include <iostream>
#include "block.h"



#ifndef SFCSORTBENCH_CUDAUTILS_H
#define SFCSORTBENCH_CUDAUTILS_H


//Macro for checking cuda errors following a cuda launch or api call
#define CUDA_CHECK_ERROR() {                                          \
 cudaError_t e=cudaGetLastError();                                 \
 if(e!=cudaSuccess) {                                              \
   printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e));           \
   exit(0); \
 }                                                                 \
}



namespace cuda
{

    /**
     * @brief send device information to the gpu
     * @param[in] device : gpu device ID
     * @return cudaDeviceProp allocated on the device
     * */
     cudaDeviceProp* getGPUDeviceInfo(unsigned int device);

     /**
      * @breif send mesh blocks to the gpu
      * @param[in] in : input array
      * @param[in] out: device pointer where the data is copied to.
      * */
      template<typename T>
       T * copyArrayToDevice(const T* in, unsigned int numElems);


    /**
     * @breif copy value to device
     * @param[in] in : input value
     * @param[in] out: device pointer where the data is copied to.
     * */
      template<typename T>
      inline T * copyValueToDevice(const T* in);


}




// templated functions

namespace cuda
{

    template<typename T>
    T * copyArrayToDevice(const T* in, unsigned int numElems)
    {

        T* __devicePtr;
        cudaMalloc(&__devicePtr,sizeof(T)*numElems);
        CUDA_CHECK_ERROR();

        cudaMemcpy(__devicePtr,in,sizeof(T)*numElems,cudaMemcpyHostToDevice);
        CUDA_CHECK_ERROR();

        return __devicePtr;

    }


    template<typename T>
    inline T * copyValueToDevice(const T* in)
    {

        T* __devicePtr;
        cudaMalloc(&__devicePtr,sizeof(T));
        CUDA_CHECK_ERROR();

        cudaMemcpy(__devicePtr,in,sizeof(T),cudaMemcpyHostToDevice);
        CUDA_CHECK_ERROR();

        return __devicePtr;

    }


}


#endif //SFCSORTBENCH_CUDAUTILS_H
