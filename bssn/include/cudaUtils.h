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
       T * allocateDeviceMemory(const T* in, unsigned int numElems);

     /**
      * @breif send mesh blocks to the gpu
      * @param[in] in : input array
      * @param[in] out: device pointer where the data is copied to.
      * */
      template<typename T>
       T * copyArrayToDevice(T* __devicePtr, const T* in, unsigned int numElems, cudaStream_t stream);


    /**
     * @breif copy value to device
     * @param[in] in : input value
     * @param[in] out: device pointer where the data is copied to.
     * */
      template<typename T>
      inline T * copyValueToDevice(const T* in);

      /**
       * @brief allocates a 2D cuda array on the device.
       * @param[in] sz1: dim 1 size
       * @param[in] sz2: dim 2 size
       * @returns the double pointer to the 2D array.
       * */
    template <typename T>
    T** alloc2DCudaArray(unsigned int sz1,  unsigned int sz2);

      /**
       * @brief allocates a 2D cuda array on the device and copy data.
       * @param[in] sz1: dim 1 size
       * @param[in] sz2: dim 2 size
       * @returns the double pointer to the 2D array.
       * */
    template <typename T>
    T** alloc2DCudaArray(const T** in,unsigned int sz1,  unsigned int sz2);


    template <typename T>
    T** copy2DCudaArray(const T** in,  T** tmp2D, unsigned int sz1, unsigned int sz2,
             T** temp2d, cudaStream_t stream);

    template <typename T>
    T** getReferenceTo2DArray(unsigned int sz1);

      /**
       * @brief deallocates the 2D cuda array.
       * @param[in] sz1: dim 1 size
       * */
      template <typename T>
      void dealloc2DCudaArray(T ** & __array2D,  unsigned int sz1, cudaStream_t stream);

      template<typename T>
        inline T * copyDeriveValues(T* __devicePtr, const T* in, cudaStream_t stream);

    template<typename T>
    inline T * allocateDerivativeMemory(const T* in);

    template <typename T>
    void copyBackMemory(T ** & __array2D, unsigned int sz1, cudaStream_t stream);

    template <typename T>
      void copy2DArrayToHost(T** in, T** out, unsigned int sz1,  unsigned int sz2, cudaStream_t stream, unsigned int offset);
}




// templated functions

namespace cuda
{
    template<typename T>
       T * allocateDeviceMemory(const T* in, unsigned int numElems)
    {
         T* __devicePtr;
        cudaMalloc(&__devicePtr,sizeof(T)*numElems);
        CUDA_CHECK_ERROR();
        return __devicePtr;
    }

    template<typename T>
    T * copyArrayToDevice( T* __devicePtr, const T* in, unsigned int numElems, cudaStream_t stream)
    {

        cudaMemcpyAsync(__devicePtr,in,sizeof(T)*numElems,cudaMemcpyHostToDevice, stream);
        CUDA_CHECK_ERROR();

        return __devicePtr;

    }

    template<typename T>
    inline T * allocateDerivativeMemory(const T* in)
    {

        T* __devicePtr;
        cudaMalloc(&__devicePtr,sizeof(T));
        CUDA_CHECK_ERROR();

        return __devicePtr;

    }

    template<typename T>
    inline T * copyDeriveValues(T* __devicePtr, const T* in, cudaStream_t stream)
    {

        cudaMemcpyAsync(__devicePtr,in,sizeof(T),cudaMemcpyHostToDevice, stream);
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

    template <typename T>
    T** alloc2DCudaArray(unsigned int sz1, unsigned int sz2)
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

        cudaMemcpy(__tmp2d,tmp2D,sizeof(T*)*sz1,cudaMemcpyHostToDevice);
        CUDA_CHECK_ERROR();
        delete [] tmp2D;

        return __tmp2d;

    }

    template <typename T>
    T** alloc2DCudaArray(T** in,unsigned int sz1,  unsigned int sz2)
    {
        // sz1 == 24
        T** tmp2D=new T*[sz1];

        for(unsigned int i=0;i<sz1;i++)
        {
            cudaMalloc(&tmp2D[i],sizeof(T)*sz2);
            CUDA_CHECK_ERROR();

        }

        return tmp2D;
    }

    template <typename T>
    T** getReferenceTo2DArray(unsigned int sz1) {
         T** __tmp2d;
        cudaMalloc(&__tmp2d,sizeof(T*)*sz1);
        CUDA_CHECK_ERROR();

        return __tmp2d;
    }

    template <typename T>
    T** copy2DCudaArray(T** in,  T** tmp2D, unsigned int sz1,
         unsigned int sz2, T** temp2d, cudaStream_t stream)
    {

        for(unsigned int i=0;i<sz1;i++)
        {
            cudaMemcpyAsync(tmp2D[i], in[i], sizeof(T)*sz2,cudaMemcpyHostToDevice, stream);
            CUDA_CHECK_ERROR();
        }
        cudaMemcpyAsync(temp2d,tmp2D,sizeof(T*)*sz1,cudaMemcpyHostToDevice, stream);
        CUDA_CHECK_ERROR();
        return temp2d;
    }

    template <typename T>
    void copyBackMemory(T ** & __array2D, unsigned int sz1, cudaStream_t stream)
    {
        T** tmp2D=new T*[sz1];

        cudaMemcpyAsync(tmp2D,__array2D,sizeof(T*)*sz1,cudaMemcpyDeviceToHost, stream);
        CUDA_CHECK_ERROR();

        // for(unsigned int i=0;i<sz1;i++)
        // {
        //     cudaFree(tmp2D[i]);
        //     CUDA_CHECK_ERROR();
        // }

        // delete [] tmp2D;

        // cudaFree(__array2D);
        // CUDA_CHECK_ERROR();
    }

    template <typename T>
    void dealloc2DCudaArray(T ** & __array2D, unsigned int sz1, cudaStream_t stream)
    {
        T** tmp2D=new T*[sz1];

        cudaMemcpyAsync(tmp2D, __array2D, sizeof(T*)*sz1, cudaMemcpyDeviceToHost, stream);
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

    template <typename T>
    void copy2DArrayToHost(T** in, T** out, unsigned int sz1,  unsigned int sz2, cudaStream_t stream, unsigned int offset)
    {
        T** tmp2D=new T*[sz1];
        cudaMemcpyAsync(tmp2D, in, sizeof(T)*sz1, cudaMemcpyDeviceToHost, stream);
        CUDA_CHECK_ERROR();

        for(unsigned int i=0; i<sz1; i++)
        {
            cudaMemcpyAsync(&out[i][offset], tmp2D[i], sizeof(T)*sz2, cudaMemcpyDeviceToHost, stream);
            CUDA_CHECK_ERROR();
        }
    }
}


#endif //SFCSORTBENCH_CUDAUTILS_H
