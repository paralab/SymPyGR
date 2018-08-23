//
// Created by milinda on 8/9/18.
//
/**
 * @brief Contians cuda derivs for bssn computation
 *
 * */
#ifndef SFCSORTBENCH_DERIVS_H
#define SFCSORTBENCH_DERIVS_H



#include "cuda_runtime.h"
#include <stdio.h>
#include <iostream>


#define IDX(i,j,k) ( (i) + nx * ( (j) + ny * (k) ) )

namespace cuda
{

        // first derivatives
        __device__ void deriv42_z(double * const  Dzu, const double * const  u, const double dz, const unsigned int** ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw ,unsigned bflag);
        __device__ void deriv42_y(double * const  Dyu, const double * const  u, const double dy, const unsigned int** ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw ,unsigned bflag);
        __device__ void deriv42_x(double * const  Dxu, const double * const  u, const double dx, const unsigned int** ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw ,unsigned bflag);

        // advective derivatives
        __device__ void deriv42adv_z(double * const  Dzu, const double * const  u,const double dz, const unsigned int** ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, const double * const betaz, unsigned int pw,unsigned bflag);
        __device__ void deriv42adv_y(double * const  Dyu, const double * const  u,const double dy, const unsigned int** ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, const double * const betay, unsigned int pw,unsigned bflag);
        __device__ void deriv42adv_x(double * const  Dxu, const double * const  u,const double dx, const unsigned int** ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, const double * const betax, unsigned int pw,unsigned bflag);

        // second order derivatives
        __device__ void deriv42_zz(double * const  Du, const double * const  u, const double dz, const unsigned int** ijk_lm, const unsigned int * sz, const unsigned int* tile_sz,unsigned int pw ,unsigned bflag);
        __device__ void deriv42_yy(double * const  Du, const double * const  u, const double dy, const unsigned int** ijk_lm, const unsigned int * sz, const unsigned int* tile_sz,unsigned int pw , unsigned bflag);
        __device__ void deriv42_xx(double * const  DxDxu, const double * const  u,const double dx, const unsigned int** ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw, unsigned bflag);


        // Kriess-Oliger derivatives
        __device__ void ko_deriv42_z(double * const Du, const double * const u, const double dz, const unsigned int** ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw,unsigned bflag);
        __device__ void ko_deriv42_y(double * const  Du, const double * const  u, const double dy, const unsigned int** ijk_lm, const unsigned int * sz, const unsigned int* tile_sz,unsigned int pw, unsigned bflag);
        __device__ void ko_deriv42_x(double * const  Du, const double * const  u, const double dx, const unsigned int** ijk_lm, const unsigned int * sz, const unsigned int* tile_sz,unsigned int pw, unsigned bflag);


}// end of namespace cuda
#endif //SFCSORTBENCH_DERIVS_H


