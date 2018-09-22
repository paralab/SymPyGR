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
     cudaDeviceProp* getGPUDeviceInfo(unsigned int device, cudaStream_t stream);

     /**
      * @breif send mesh blocks to the gpu
      * @param[in] in : input array
      * @param[in] out: device pointer where the data is copied to.
      * */
      template<typename T>
       T * copyArrayToDevice(const T* in, unsigned int numElems, cudaStream_t stream);


    /**
     * @breif copy value to device
     * @param[in] in : input value
     * @param[in] out: device pointer where the data is copied to.
     * */
      template<typename T>
      inline T * copyValueToDevice(const T* in, cudaStream_t stream);

      /**
       * @brief allocates a 2D cuda array on the device.
       * @param[in] sz1: dim 1 size
       * @param[in] sz2: dim 2 size
       * @returns the double pointer to the 2D array.
       * */
      template <typename T>
      T** alloc2DCudaArray(unsigned int sz1,  unsigned int sz2, cudaStream_t stream);

      /**
       * @brief allocates a 2D cuda array on the device and copy data.
       * @param[in] sz1: dim 1 size
       * @param[in] sz2: dim 2 size
       * @returns the double pointer to the 2D array.
       * */
      template <typename T>
      T** alloc2DCudaArray(const T** in,unsigned int sz1,  unsigned int sz2, cudaStream_t stream);


      /**
       * @brief deallocates the 2D cuda array.
       * @param[in] sz1: dim 1 size
       * */
      template <typename T>
      void dealloc2DCudaArray(T ** & __array2D,  unsigned int sz1, cudaStream_t stream);
}




// templated functions

namespace cuda
{

    template<typename T>
    T * copyArrayToDevice(const T* in, unsigned int numElems, cudaStream_t stream)
    {

        T* __devicePtr;
        cudaMalloc(&__devicePtr,sizeof(T)*numElems);
        CUDA_CHECK_ERROR();

        cudaMemcpyAsync(__devicePtr,in,sizeof(T)*numElems,cudaMemcpyHostToDevice, stream);
        CUDA_CHECK_ERROR();

        return __devicePtr;

    }


    template<typename T>
    inline T * copyValueToDevice(const T* in, cudaStream_t stream)
    {

        T* __devicePtr;
        cudaMalloc(&__devicePtr,sizeof(T));
        CUDA_CHECK_ERROR();

        cudaMemcpyAsync(__devicePtr,in,sizeof(T),cudaMemcpyHostToDevice, stream);
        CUDA_CHECK_ERROR();

        return __devicePtr;

    }

    template <typename T>
    T** alloc2DCudaArray(unsigned int sz1, unsigned int sz2, cudaStream_t stream)
    {

        T** __tmp2d;
        cudaMalloc(&__tmp2d,sizeof(T*)*sz1);
        CUDA_CHECK_ERROR();

        T** tmp2D=new T*[sz1];

        for(unsigned int i=0;i<sz1;i++)
        {
            cudaMalloc(&tmp2D[i],sizeof(T)*sz2);
            CUDA_CHECK_ERROR();
        }

        cudaMemcpyAsync(__tmp2d,tmp2D,sizeof(T*)*sz1,cudaMemcpyHostToDevice, stream);
        CUDA_CHECK_ERROR();
        delete [] tmp2D;

        return __tmp2d;

    }

    template <typename T>
    T** alloc2DCudaArray(const T** in,unsigned int sz1,  unsigned int sz2, cudaStream_t stream)
    {
        T** __tmp2d;
        cudaMalloc(&__tmp2d,sizeof(T*)*sz1);
        CUDA_CHECK_ERROR();

        T** tmp2D=new T*[sz1];

        for(unsigned int i=0;i<sz1;i++)
        {
            cudaMalloc(&tmp2D[i],sizeof(T)*sz2);
            CUDA_CHECK_ERROR();
            cudaMemcpyAsync(tmp2D[i],in[i], sizeof(T)*sz2 ,cudaMemcpyHostToDevice, stream);
            CUDA_CHECK_ERROR();
        }

        cudaMemcpyAsync(__tmp2d,tmp2D,sizeof(T*)*sz1,cudaMemcpyHostToDevice, stream);
        CUDA_CHECK_ERROR();
        delete [] tmp2D;

        return __tmp2d;
    }


    template <typename T>
    void dealloc2DCudaArray(T ** & __array2D, unsigned int sz1, cudaStream_t stream)
    {
        T** tmp2D=new T*[sz1];

        cudaMemcpyAsync(tmp2D,__array2D,sizeof(T*)*sz1,cudaMemcpyDeviceToHost, stream);
        CUDA_CHECK_ERROR();

        for(unsigned int i=0;i<sz1;i++)
        {
            cudaFree(tmp2D[i]);
            CUDA_CHECK_ERROR();
        }

        delete [] tmp2D;

        cudaFree(__array2D);
        CUDA_CHECK_ERROR();
    }




}


#endif //SFCSORTBENCH_CUDAUTILS_H
