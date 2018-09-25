//
// Created by milinda on 8/9/18.
//
/**
* @author : Milinda Fernando
* @brief : contains utility functions for bssn cuda computations. 
*/
#ifndef BSSN_CUDAUTILS_H
#define BSSN_CUDAUTILS_H

#include "cuda_runtime.h"
#include <device_launch_parameters.h>
#include "derivs.cuh"




namespace cuda
{


    enum VAR {U_ALPHA=0,U_CHI,U_K,U_GT0,U_GT1,U_GT2,U_BETA0,U_BETA1,U_BETA2,U_B0,U_B1,U_B2,U_SYMGT0,U_SYMGT1,U_SYMGT2,U_SYMGT3,U_SYMGT4,U_SYMGT5,U_SYMAT0,U_SYMAT1,U_SYMAT2,U_SYMAT3,U_SYMAT4,U_SYMAT5};
    //enum VAR_PSI4 {C_PSI4_REAL, C_PSI4_IMG};
    enum VAR_CONSTRAINT {C_HAM=0, C_MOM0, C_MOM1, C_MOM2, C_PSI4_REAL, C_PSI4_IMG};



    /**
    * @brief load (1D) data from global memory to shared memory
    * @param[in] __globalIn : source global address to copy data from
    * @param[out] sharedOut : shared memory location
    * @param[in] ijk_lm : min max of index i,j,k of the dendro block included in the tile 3x2
    * @param[in] sz: dendro block size
    * @param[in] tile_sz: tile sz
    */
    template<typename T>
    __device__ void __loadGlobalToShared1D(const T* __globalIn,  T* sharedOut, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz);


    /**
    * @brief load (3D) data from global memory to shared memory
    * @param[in] __globalIn : source global address to copy data from
    * @param[out] sharedOut : shared memory location
    * @param[in] ijk_lm : min max of index i,j,k of the dendro block included in the tile 3x2
    * @param[in] sz: dendro block size
    * @param[in] tile_sz: tile sz
    */
    template<typename T>
    __device__ void __loadGlobalToShared3D(const T* __globalIn,  T* sharedOut, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz);


    /**
   * @brief extract the sign of the shift vectors
   * @param[in] __globalIn : source global address to copy data from
   * @param[out] sharedOut : shared memory location (bool value for of the sign)
   * @param[in] ijk_lm : min max of index i,j,k of the dendro block included in the tile 3x2
   * @param[in] sz: dendro block size
   * @param[in] tile_sz: tile sz
   */
    template<typename T>
    __device__ void __extractSign3D(const T* sharedIn,  bool * sharedOut, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz);



    /**
    * @brief store (1D) data from shared memory to global memory
    * @param[in] __globalout : dest. global address to copy data to
    * @param[out] sharedIn : shared memory location
    * @param[in] ijk_lm : min max of index i,j,k of the dendro block included in the tile
    * @param[in] sz: dendro block size
    * @param[in] tile_sz: tile sz
    */
    template <typename T>
    __device__ void __storeSharedToGlobal1D(const T* sharedIn, T* __globalOut,const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz);


    /**
    * @brief store (3D) data from shared memory to global memory
    * @param[in] __globalout : dest. global address to copy data to
    * @param[out] sharedIn : shared memory location
    * @param[in] ijk_lm : min max of index i,j,k of the dendro block included in the tile
    * @param[in] sz: dendro block size
    * @param[in] tile_sz: tile sz
    */
    template <typename T>
    __device__ void __storeSharedToGlobal3D(const T* sharedIn, T* __globalOut,const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz);





}






namespace cuda
{


    template<typename T>
    __device__ void __loadGlobalToShared1D(const T* __globalIn, T* sharedOut, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz)
    {

        const unsigned int i_b=ijk_lm[2*0+0];
        const unsigned int i_e=ijk_lm[2*0+1];

        for(unsigned int i=i_b+threadIdx.x;i<i_e;i+=blockDim.x)
            sharedOut[(i-i_b )] = __globalIn[i];

        return ;

    }



    template<typename T>
    __device__ void __loadGlobalToShared3D(const T* __globalIn, T* sharedOut, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz)
    {

        const unsigned int i_b=ijk_lm[2*0+0];
        const unsigned int i_e=ijk_lm[2*0+1];

        const unsigned int j_b=ijk_lm[2*1+0];
        const unsigned int j_e=ijk_lm[2*1+1];

        const unsigned int k_b=ijk_lm[2*2+0];
        const unsigned int k_e=ijk_lm[2*2+1];


        for(unsigned int k=k_b+threadIdx.z;k<k_e;k+=blockDim.z)
            for(unsigned int j=j_b+threadIdx.y;j<j_e;j+=blockDim.y)
                for(unsigned int i=i_b+threadIdx.x;i<i_e;i+=blockDim.x)
                    sharedOut[(k-k_b) * (tile_sz[0]*tile_sz[1]) + (j-j_b) * (tile_sz[0])+ (i-i_b )] = __globalIn[k*(sz[1]*sz[0])+j*(sz[0])+i];



        return ;

    }


    template<typename T>
    __device__ void __storeSharedToGlobal1D(const T* sharedIn, T* __globalOut,const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz)
    {

        const unsigned int i_b=ijk_lm[2*0+0];
        const unsigned int i_e=ijk_lm[2*0+1];

        for(unsigned int i=i_b+threadIdx.x;i<i_e;i+=blockDim.x)
           __globalOut[i]=sharedIn[(i-i_b )];

        return ;


    }




    template<typename T>
    __device__ void __storeSharedToGlobal3D(const T* sharedIn, T* __globalOut,const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz)
    {

        const unsigned int i_b=ijk_lm[2*0+0];
        const unsigned int i_e=ijk_lm[2*0+1];

        const unsigned int j_b=ijk_lm[2*1+0];
        const unsigned int j_e=ijk_lm[2*1+1];

        const unsigned int k_b=ijk_lm[2*2+0];
        const unsigned int k_e=ijk_lm[2*2+1];


        for(unsigned int k=k_b+threadIdx.z;k<k_e;k+=blockDim.z)
            for(unsigned int j=j_b+threadIdx.y;j<j_e;j+=blockDim.y)
                for(unsigned int i=i_b+threadIdx.x;i<i_e;i+=blockDim.x)
                    __globalOut[k*(sz[1]*sz[0])+j*(sz[0])+i]=sharedIn[(k-k_b) * (tile_sz[0]*tile_sz[1]) + (j-j_b) * (tile_sz[0])+ (i-i_b )];



        return ;


    }

    template<typename T>
    __device__ void __extractSign3D(const T* sharedIn,  bool * sharedOut, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz)
    {

        const unsigned int i_b=ijk_lm[2*0+0];
        const unsigned int i_e=ijk_lm[2*0+1];

        const unsigned int j_b=ijk_lm[2*1+0];
        const unsigned int j_e=ijk_lm[2*1+1];

        const unsigned int k_b=ijk_lm[2*2+0];
        const unsigned int k_e=ijk_lm[2*2+1];


        for(unsigned int k=k_b+threadIdx.z;k<k_e;k+=blockDim.z)
            for(unsigned int j=j_b+threadIdx.y;j<j_e;j+=blockDim.y)
                for(unsigned int i=i_b+threadIdx.x;i<i_e;i+=blockDim.x)
                {
                    (sharedIn[(k-k_b) * (tile_sz[0]*tile_sz[1]) + (j-j_b) * (tile_sz[0])+ (i-i_b )]>0) ?  sharedOut[(k-k_b) * (tile_sz[0]*tile_sz[1]) + (j-j_b) * (tile_sz[0])+ (i-i_b )]=true : sharedOut[(k-k_b) * (tile_sz[0]*tile_sz[1]) + (j-j_b) * (tile_sz[0])+ (i-i_b )]=false;
                }

    }

}


#endif //BSSN_CUDAUTILS_H





