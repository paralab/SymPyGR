// generated by Dendro-GR SymPyGR code gernation framework
//date: 2018-09-28 14:08:54
#ifndef DERIVS_BSSN_CUH 
#define DERIVS_BSSN_CUH 
#include<iostream>
#include"cuda_runtime.h"
#include<device_launch_parameters.h>
#include "block_cu.h"
#include "params_cu.h"
#include "bssn_rhs_deriv_mem_cuda.h"
#include "cudaUtils.cuh"
#include "derivs.cuh"
namespace cuda {

/**compute derivative kernel __RSWS_computeDerivs*/
__global__ void __RSWS_computeDerivs(const double**__unzipInVar,MemoryDerivs* __derivWorkspace, const cuda::_Block* __dendroBlkList, const unsigned int * __gpuBlockMap, const cudaDeviceProp* __deviceProperties);

}// end of namespace cuda


#endif
