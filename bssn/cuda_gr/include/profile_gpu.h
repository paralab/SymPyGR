//
// Created by milinda on 8/30/18.
//

/**
 * @author Milinda Fernando
 * School Of Computing, University of Utah
 * @brief contains profile parameters for the derv computations.
 *
 * */

#ifndef DENDRO_5_0_PROFILE_GPU_H
#define DENDRO_5_0_PROFILE_GPU_H

#include "bssn_profiler.h"
#include "block.h"
#include <iostream>
#include <vector>



namespace cuda
{
    namespace profile
    {

        /** overall gpu exuxution time*/
        extern bssn_profiler_t t_overall;

        /**host to device communication*/
        extern bssn_profiler_t t_H2D_Comm;
        /**device to host communication*/
        extern bssn_profiler_t t_D2H_Comm;

        /**memory allocation deallocation time*/
        extern bssn_profiler_t t_cudaMalloc_derivs;

        /**deriv computation time*/
        extern bssn_profiler_t t_derivs;



        /**initialize the profile counters*/
        void initialize();

        /**output the timings*/
        void printOutput(const std::vector<ot::Block>& localBlkList);



    }// end of namespace profile

}// end of namespace cuda

#endif //DENDRO_5_0_PROFILE_GPU_H
