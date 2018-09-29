//
// Created by milinda on 8/9/18.
//
/**
 * @author Milinda Fernando
 * School of Computing, University of Utah
 * @brief Contains utility function for the host related to GPUs
 * */




#ifndef SFCSORTBENCH_CUDAUTILS_H
#define SFCSORTBENCH_CUDAUTILS_H

#include "cuda_runtime.h"
#include <iostream>
#include "block.h"
#include <vector>


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
       T * copyArrayToDevice(T* __devicePtr, const T* in, unsigned int numElems);

    /**
      * @breif allocate mesh blocks in gpu
      * @param[in] in : input array
      * @param[in] out: device pointer where the data is copied to.
      * */
    template<typename T>
    T * allocateMemoryForArray(unsigned int numElems);

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



    /**
     * @breif send mesh blocks to the gpu (async)
     * @param[in] in : input array
     * @param[in] out: device pointer where the data is copied to.
     * */
    template<typename T>
    T * copyArrayToDeviceAsync(const T* in, unsigned int numElems,const cudaStream_t* stream,unsigned int * sendCounts,unsigned int * sendOffset,unsigned int numStreams);

    /**
       * @brief allocates a 2D cuda array on the device and copy data.
       * @param[in] sz1: dim 1 size
       * @param[in] sz2: dim 2 size
       * @returns the double pointer to the 2D array.
       * */
    template <typename T>
    T** alloc2DCudaArrayAsync(const T** in,unsigned int sz1,  unsigned int sz2,const cudaStream_t * stream,unsigned int * sendCounts,unsigned int * sendOffset,unsigned int numStreams);



    /**
     * @brief deallocates the 2D cuda array.
     * @param[in] sz1: dim 1 size
     * */
    template <typename T>
    void dealloc2DCudaArray(T ** & __array2D,  unsigned int sz1);

    /***
     * computes the how dendro blocks (octree blocks result in from unzip) to the gpu/
     * @param[in] blkList: list of dendro data blocks
     * @param[in] numBlocks: number of data blocks,
     * @param[out] blockMap: (blockMap[2*blocDim.x] , blockMap[2*blocDim.x+1]) begin & end of data block that is going to be process by the gpu block
     * */
    void computeDendroBlockToGPUMap(const ot::Block* blkList, unsigned int numBlocks, std::vector<unsigned int >& blockMap,dim3 & gridDim);






}



    /**
       * @copy a 2D array from GPU to CPU.
       * @param[in] sz1: dim 1 size
       * */
    template <typename T>
    void copy2DArrayToHost(T** in, T** out, unsigned int sz1,  unsigned int sz2, unsigned int offset);
}

// templated functions

namespace cuda
{

    template<typename T>
    T * copyArrayToDevice(T* __devicePtr, const T* in, unsigned int numElems)
    {

        cudaMemcpy(__devicePtr,in,sizeof(T)*numElems,cudaMemcpyHostToDevice);
        CUDA_CHECK_ERROR();

        return __devicePtr;

    }

    template<typename T>
    T * allocateMemoryForArray(unsigned int numElems)
    {

        T* __devicePtr;
        cudaMalloc(&__devicePtr,sizeof(T)*numElems);
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

        T** tmp2D=new T*[sz1];

        for(unsigned int i=0;i<sz1;i++)
        {
            cudaMalloc(&tmp2D[i],sizeof(T)*sz2);
            CUDA_CHECK_ERROR();

        }

        return tmp2D;
    }

    template <typename T>
    T** copy2DCudaArray(const T** in,  T** tmp2D, unsigned int sz1, unsigned int sz2, T** inputReferenceGPU,
                        unsigned int offset)
    {

        for(unsigned int i=0;i<sz1;i++)
        {
            cudaMemcpy(tmp2D[i], &in[i][offset], sizeof(T)*sz2, cudaMemcpyHostToDevice);
            CUDA_CHECK_ERROR();
        }
        
        cudaMemcpy(inputReferenceGPU, tmp2D, sizeof(T*)*sz1, cudaMemcpyHostToDevice);
        CUDA_CHECK_ERROR();

        return inputReferenceGPU;
    }

    template <typename T>
    T** alloc2DGPUArray(unsigned int sz1,  unsigned int sz2)
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
    void dealloc2DCudaArray(T ** & __array2D, unsigned int sz1)
    {
        T** tmp2D=new T*[sz1];

        cudaMemcpy(tmp2D,__array2D,sizeof(T*)*sz1,cudaMemcpyDeviceToHost);
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
    T** getReferenceTo2DArray(unsigned int sz1) {
         T** __tmp2d;
        cudaMalloc(&__tmp2d,sizeof(T*)*sz1);
        CUDA_CHECK_ERROR();

        return __tmp2d;
    }

    template<typename T>
    T * copyArrayToDeviceAsync(const T* in, unsigned int numElems,const cudaStream_t* stream,unsigned int * sendCounts,unsigned int * sendOffset,unsigned int numStreams)
    {


        T* __devicePtr;
        cudaMalloc(&__devicePtr,sizeof(T)*numElems);
        CUDA_CHECK_ERROR();

        sendOffset[0]=0;
        for(unsigned int i=0;i<numStreams;i++)
            sendCounts[i]=((((i+1)*numElems)/numStreams)-(((i)*numElems)/numStreams));

        for(unsigned int i=1;i<numStreams;i++)
            sendOffset[i]=sendOffset[i-1]+sendCounts[i-1];

        for(unsigned int i=0;i<numStreams;i++)
        {
            cudaMemcpyAsync(__devicePtr+sendOffset[i],in+sendOffset[i],sizeof(T)*sendCounts[i],cudaMemcpyHostToDevice,stream[i]);
            CUDA_CHECK_ERROR();
        }

        return __devicePtr;

    }


    template <typename T>
    T** alloc2DCudaArrayAsync(const T** in,unsigned int sz1,  unsigned int sz2,const cudaStream_t * stream,unsigned int * sendCounts,unsigned int * sendOffset,unsigned int numStreams)
    {
        T** __tmp2d;
        cudaMalloc(&__tmp2d,sizeof(T*)*sz1);
        CUDA_CHECK_ERROR();

        T** tmp2D=new T*[sz1];


        sendOffset[0]=0;
        for(unsigned int i=0;i<numStreams;i++)
            sendCounts[i]=((((i+1)*sz2)/numStreams)-(((i)*sz2)/numStreams));

        for(unsigned int i=1;i<numStreams;i++)
            sendOffset[i]=sendOffset[i-1]+sendCounts[i-1];


        for(unsigned int i=0;i<sz1;i++)
        {
            cudaMalloc(&tmp2D[i],sizeof(T)*sz2);
            CUDA_CHECK_ERROR();

            for(unsigned int j=0;j<numStreams;j++)
            {
                cudaMemcpyAsync(tmp2D[i]+sendOffset[j],in[i]+sendOffset[j], sizeof(T)*(sendCounts[j]) ,cudaMemcpyHostToDevice,stream[j]);
                CUDA_CHECK_ERROR();
            }

        }

        cudaMemcpyAsync(__tmp2d,tmp2D,sizeof(T*)*sz1,cudaMemcpyHostToDevice,stream[0]);
        CUDA_CHECK_ERROR();
        delete [] tmp2D;

        return __tmp2d;
    }

        for(unsigned int i=0; i<sz1; i++)
        {
            cudaMemcpy(&out[i][offset], tmp2D[i], sizeof(T)*sz2, cudaMemcpyDeviceToHost);
            CUDA_CHECK_ERROR();
        }
    }

}

#endif //SFCSORTBENCH_CUDAUTILS_H
