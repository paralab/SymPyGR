//
// Created by milinda on 1/15/18.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief
*/
//

#ifndef BSSN_PROFILE_PARAM_H
#define BSSN_PROFILE_PARAM_H

#include "bssn_profiler.h"
#include "test_param.h"

namespace bssn
{
    namespace timer
    {

        extern bssn_profiler_t 	total_runtime;
        extern bssn_profiler_t t_cpu_runtime;
        extern bssn_profiler_t	t_deriv;
        extern bssn_profiler_t	t_rhs;
        extern bssn_profiler_t  t_bdyc;

        extern bssn_profiler_t t_gpu_runtime;
        extern bssn_profiler_t t_sorting;
        extern bssn_profiler_t t_malloc_free;
        extern bssn_profiler_t t_mem_handling_gpu;
        extern bssn_profiler_t t_memcopy_kernel;
        extern bssn_profiler_t t_memcopy;

        extern bssn_profiler_t flop_count;

        void profileInfo();
    }
}



#endif //BSSN_PROFILE_PARAM_H
