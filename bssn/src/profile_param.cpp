//
// Created by milinda on 1/15/18.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief
*/
//

#include "profile_param.h"

#define max_dp_flop_rate 147.8
#define max_bandwidth 12.05

namespace bssn
{
    namespace timer
    {

        bssn_profiler_t total_runtime;
        bssn_profiler_t t_cpu_runtime;
        bssn_profiler_t	t_deriv;
        bssn_profiler_t	t_rhs;
        bssn_profiler_t t_bdyc;

        bssn_profiler_t t_gpu_runtime;
        bssn_profiler_t t_sorting;
        bssn_profiler_t t_malloc_free;
        bssn_profiler_t t_mem_handling_gpu;
        bssn_profiler_t t_memcopy_kernel;
        bssn_profiler_t t_memcopy;

        bssn_profiler_t flop_count;
        bssn_profiler_t total_points;

    }
}


void bssn::timer::profileInfo()
{

    const char separator    = ' ';
    const int nameWidth     = 30;
    const int numWidth      = 10;

    double t_stat;

    t_stat=total_runtime.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"runtime(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat<<std::endl;

    std::cout<<"\nCPU Timings"<<std::endl;

    t_stat=t_cpu_runtime.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"--cpu_runtime(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat<<std::endl;


    std::cout<<"GPU Timings"<<std::endl;

    t_stat=t_gpu_runtime.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"--gpu_runtime(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat<<std::endl;

    t_stat=t_sorting.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"--block_sorting(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) << t_stat <<std::endl;

    t_stat=t_malloc_free.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"--malloc_free(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) << t_stat <<std::endl;

    t_stat=t_memcopy_kernel.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"--kernel_memcopy(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) << t_stat <<std::endl;

    t_stat=t_memcopy.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"--memcopy(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) << t_stat <<std::endl;

    std::cout<<"Speedup"<<std::endl;

    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  Overall speedup";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<< t_cpu_runtime.seconds/t_gpu_runtime.seconds <<std::endl;

    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  Approx. FlopRate (GFlop/s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<< flop_count.seconds/(t_memcopy_kernel.seconds) << "Percentage " << flop_count.seconds/(t_memcopy_kernel.seconds)*100/max_dp_flop_rate << "%" << std::endl;

    // std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  Appr. Bandwidth (GB/s)";
    // std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<< total_points.seconds*8/1024/1024/1024/t_memcopy_kernel.seconds << "Percentage " << (total_points.seconds*8/1024/1024/1024/t_memcopy_kernel.seconds*100)/max_bandwidth << "%" << std::endl;
}

