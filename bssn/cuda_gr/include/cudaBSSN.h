#ifndef CUDA_BSSN_H
#define CUDA_BSSN_H

#include "dataGeneration.h"
#include "rhsMethods.h"

int steamCountToLevel[5] = {3, 3, 3, 2, 2};
int hybrid_divider = 3;

// Select RHS Method
#define parallelized 1
#define async 0
#define hybrid 0

#endif