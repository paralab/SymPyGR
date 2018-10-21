//
// Created by milinda on 8/20/18.
//

#ifndef DENDRO_5_0_GPUBSSNEXAMPLE_H
#define DENDRO_5_0_GPUBSSNEXAMPLE_H

/**
 * @author Milinda Fernando
 * @brief simple gr mesh generation and launch cuda kernels to compute BSSN rhs
 *
 * */


#include <iostream>
#include <random>
#include "block.h"
#include <vector>
#include "profile_param.h"
#include <grUtils.h>
#include "rhs.h"
#include "block_cu.h"
#include "profile_gpu.h"
#include "parameters.h"
#include "point.h"
#include "grDef.h"
#include <chrono>
#ifdef BSSN_ENABLE_CUDA
    #include "rhs_cuda.cuh"
#endif




typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::milliseconds ms;
typedef std::chrono::duration<float> fsec;


#endif //DENDRO_5_0_GPUBSSNEXAMPLE_H
