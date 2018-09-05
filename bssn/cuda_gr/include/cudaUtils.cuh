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
    * @brief load data from global memory to shared memory
    * @param[in] __globalIn : source global address to copy data from
    * @param[out] sharedOut : shared memory location
    * @param[in] ijk_lm : min max of index i,j,k of the dendro block included in the tile 3x2
    * @param[in] sz: dendro block size
    * @param[in] tile_sz: tile sz
    */
    template<typename T>
    __device__ void __loadGlobalToShared(const T* __globalIn,  T* sharedOut, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz);
    

    /**
    * @brief store data from shared memory to global memory
    * @param[in] __globalout : dest. global address to copy data to
    * @param[out] sharedIn : shared memory location
    * @param[in] ijk_lm : min max of index i,j,k of the dendro block included in the tile
    * @param[in] sz: dendro block size
    * @param[in] tile_sz: tile sz
    */
    template <typename T>
    __device__ void __storeSharedToGlobal(const T* sharedIn, T* __globalOut,const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz);
}

#endif //SFCSORTBENCH_CUDAUTILS_H




namespace cuda
{


    template<typename T>
    __device__ void __loadGlobalToShared(const T* __globalIn, T* sharedOut, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz)
    {

        const unsigned int l_x=ijk_lm[2*0+1]-ijk_lm[2*0+0];
        const unsigned int l_y=ijk_lm[2*1+1]-ijk_lm[2*1+0];
        const unsigned int l_z=ijk_lm[2*2+1]-ijk_lm[2*2+0];

        const unsigned int i_b=ijk_lm[2*0+0];
        const unsigned int i_e=ijk_lm[2*0+1];

        const unsigned int j_b=ijk_lm[2*1+0];
        const unsigned int j_e=ijk_lm[2*1+1];

        const unsigned int k_b=ijk_lm[2*2+0];
        const unsigned int k_e=ijk_lm[2*2+1];

        if(threadIdx.x>=l_x || threadIdx.y >= l_y || threadIdx.z>=l_z) return;

        const unsigned int ix_b=i_b + (threadIdx.x * l_x)/blockDim.x;
        const unsigned int ix_e=i_b + ((threadIdx.x+1) * l_x)/blockDim.x;

        const unsigned int jy_b=j_b + (threadIdx.y * l_y)/blockDim.y;
        const unsigned int jy_e=j_b + ((threadIdx.y+1) * l_y)/blockDim.y;

        const unsigned int kz_b=k_b + (threadIdx.z * (l_z))/blockDim.z;
        const unsigned int kz_e=k_b + ((threadIdx.z+1) * (l_z))/blockDim.z;

        for(unsigned int k=kz_b;k<kz_e;k++)
            for(unsigned int j=jy_b;j<jy_e;j++)
                for(unsigned int i=ix_b;i<ix_e;i++)
                    sharedOut[(k-k_b) * (tile_sz[0]*tile_sz[1]) + (j-j_b) * (tile_sz[0])+ (i-i_b )] = __globalIn[k*(sz[1]*sz[0])+j*(sz[0])+i];

        return ;

    }




    template<typename T>
    __device__ void __storeSharedToGlobal(const T* sharedIn, T* __globalOut,const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz)
    {

        const unsigned int l_x=ijk_lm[2*0+1]-ijk_lm[2*0+0];
        const unsigned int l_y=ijk_lm[2*1+1]-ijk_lm[2*1+0];
        const unsigned int l_z=ijk_lm[2*2+1]-ijk_lm[2*2+0];

        const unsigned int i_b=ijk_lm[2*0+0];
        const unsigned int i_e=ijk_lm[2*0+1];

        const unsigned int j_b=ijk_lm[2*1+0];
        const unsigned int j_e=ijk_lm[2*1+1];

        const unsigned int k_b=ijk_lm[2*2+0];
        const unsigned int k_e=ijk_lm[2*2+1];

        if(threadIdx.x>=l_x || threadIdx.y >= l_y || threadIdx.z>=l_z) return;

        const unsigned int ix_b=i_b + (threadIdx.x * l_x)/blockDim.x;
        const unsigned int ix_e=i_b + ((threadIdx.x+1) * l_x)/blockDim.x;

        const unsigned int jy_b=j_b + (threadIdx.y * l_y)/blockDim.y;
        const unsigned int jy_e=j_b + ((threadIdx.y+1) * l_y)/blockDim.y;

        const unsigned int kz_b=k_b + (threadIdx.z * (l_z))/blockDim.z;
        const unsigned int kz_e=k_b + ((threadIdx.z+1) * (l_z))/blockDim.z;

        for(unsigned int k=kz_b;k<kz_e;k++)
            for(unsigned int j=jy_b;j<jy_e;j++)
                for(unsigned int i=ix_b;i<ix_e;i++)
                    __globalOut[k*(sz[1]*sz[0])+j*(sz[0])+i]=sharedIn[(k-k_b) * (tile_sz[0]*tile_sz[1]) + (j-j_b) * (tile_sz[0])+ (i-i_b )];


        return ;


    }

}






