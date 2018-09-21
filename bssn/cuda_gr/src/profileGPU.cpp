//
// Created by milinda on 8/30/18.
//

/**
 * Created on: Sep 21, 2018
 * 		Adapted by: Akila, Eranga, Eminda, Ruwan
 **/

#include "profileGPU.h"

namespace cuda
{
    namespace profile
    {
        bssn_profiler_t t_cpu;

        /** overall gpu exuxution time*/
        bssn_profiler_t t_overall;

        /** memory allocation */
        bssn_profiler_t t_malloc;

        /** releasing memory */
        bssn_profiler_t t_free;

        /** kernel + memcopy time */
        bssn_profiler_t t_kernel_memcopy;

    } // end of namespace profile
} // end of namespace cuda

void cuda::profile::printOutput() {

    const char separator    = ' ';
    const int nameWidth     = 30;
    const int numWidth      = 10;
    double t_stat;

    #ifdef ENABLE_CUDA_TEST
    t_stat=t_cpu.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"+cpu_overall(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat<<std::endl;
    #endif

    t_stat=t_overall.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"+cuda_overall(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat<<std::endl;

    t_stat=t_malloc.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  -mem_alloc(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat<<std::endl;

    t_stat=t_free.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  -cuda_free(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat<<std::endl;

    t_stat=t_kernel_memcopy.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  -kernel_memcopy(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat<<std::endl;

    return;
}


