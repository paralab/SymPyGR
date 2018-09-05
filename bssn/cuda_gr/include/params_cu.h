//
// Created by milinda on 8/9/18.
//

/***
 * @brief contains GPU related parameters for block sheduling data transfers etc.
 *
 *
 */

#ifndef SFCSORTBENCH_PARAMS_H
#define SFCSORTBENCH_PARAMS_H

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "cuda_runtime.h"
#include "block.h"
#include "block_cu.h"
#include <iostream>
#include "bssn_rhs_deriv_mem_cuda.h"


#define __Rx (cuda::_bssn::__BSSN_COMPD_MAX[0]-cuda::_bssn::__BSSN_COMPD_MIN[0])
#define __Ry (cuda::_bssn::__BSSN_COMPD_MAX[1]-cuda::_bssn::__BSSN_COMPD_MIN[1])
#define __Rz (cuda::_bssn::__BSSN_COMPD_MAX[2]-cuda::_bssn::__BSSN_COMPD_MIN[2])

#define __RgX (cuda::_bssn::__BSSN_OCTREE_MAX[0]-cuda::_bssn::__BSSN_OCTREE_MIN[0])
#define __RgY (cuda::_bssn::__BSSN_OCTREE_MAX[1]-cuda::_bssn::__BSSN_OCTREE_MIN[1])
#define __RgZ (cuda::_bssn::__BSSN_OCTREE_MAX[2]-cuda::_bssn::__BSSN_OCTREE_MIN[2])

#define __GRIDX_TO_X(xg) (((__Rx/__RgX)*(xg-cuda::_bssn::__BSSN_OCTREE_MIN[0]))+cuda::_bssn::__BSSN_COMPD_MIN[0])
#define __GRIDY_TO_Y(yg) (((__Ry/__RgY)*(yg-cuda::_bssn::__BSSN_OCTREE_MIN[1]))+cuda::_bssn::__BSSN_COMPD_MIN[1])
#define __GRIDZ_TO_Z(zg) (((__Rz/__RgZ)*(zg-cuda::_bssn::__BSSN_OCTREE_MIN[2]))+cuda::_bssn::__BSSN_COMPD_MIN[2])

#define __X_TO_GRIDX(xc) (((__RgX/__Rx)*(xc-cuda::_bssn::__BSSN_COMPD_MIN[0]))+cuda::_bssn::__BSSN_OCTREE_MIN[0])
#define __Y_TO_GRIDY(yc) (((__RgY/__Ry)*(yc-cuda::_bssn::__BSSN_COMPD_MIN[1]))+cuda::_bssn::__BSSN_OCTREE_MIN[1])
#define __Z_TO_GRIDZ(zc) (((__RgZ/__Rz)*(zc-cuda::_bssn::__BSSN_COMPD_MIN[2]))+cuda::_bssn::__BSSN_OCTREE_MIN[2])




namespace cuda
{


        /**stores the device properties*/
        extern  cudaDeviceProp* __CUDA_DEVICE_PROPERTIES;

        /**pointer to the block list*/
        extern _Block* __DENDRO_BLOCK_LIST;

        /** number of blocks*/
        extern  unsigned int*  __DENDRO_NUM_BLOCKS;

        /** number of evol vars */
        extern  unsigned int* __BSSN_NUM_VARS;

        /** number of constraint vars */
        extern  unsigned int* __BSSN_CONSTRAINT_NUM_VARS;


        /** x% of block shared memory utilised for the bssn computations*/
        extern double * __GPU_BLOCK_SHARED_MEM_UTIL;

        /**max block size*/
        extern unsigned int * __DENDRO_BLK_MAX_SZ;

        /**memory for deriv computations. */
        extern MemoryDerivs __BSSN_DERIV_WORKSPACE;

        /**unzip input */
        extern double ** __UNZIP_INPUT;

        /**unzip output*/
        extern double ** __UNZIP_OUTPUT;






}// end of namespace cuda


#endif //SFCSORTBENCH_PARAMS_H
