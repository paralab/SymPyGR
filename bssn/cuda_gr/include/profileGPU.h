//
// Created by milinda on 8/30/18.
//

/**
 * @author Milinda Fernando
 * School Of Computing, University of Utah
 * @brief contains profile parameters for the derv computations.
 *
 * */

/**
 * Created on: Sep 21, 2018
 * 		Adapted by: Akila, Eranga, Eminda, Ruwan
 **/

#ifndef PROFILE_GPU_H
#define PROFILE_GPU_H

#include "profilerClass.h"
#include <iostream>

namespace cuda
{
    namespace profile
    {
        extern bssn_profiler_t t_cpu;

        /** overall gpu exuxution time*/
        extern bssn_profiler_t t_overall;

        /** memory allocation */
        extern bssn_profiler_t t_malloc;

        /** releasing memory */
        extern bssn_profiler_t t_free;

        /** kernel + memcopy time */
        extern bssn_profiler_t t_kernel_memcopy;

        /**output the timings*/
        void printOutput();



    }// end of namespace profile

}// end of namespace cuda

#endif //DENDRO_5_0_PROFILE_GPU_H
