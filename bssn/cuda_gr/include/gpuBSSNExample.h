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
#ifdef BSSN_ENABLE_CUDA
    #include "rhs_cuda.cuh"
#endif





#endif //DENDRO_5_0_GPUBSSNEXAMPLE_H
