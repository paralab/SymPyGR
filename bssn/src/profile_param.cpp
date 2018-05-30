//
// Created by milinda on 1/15/18.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief
*/
//

#include "profile_param.h"

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
        bssn_profiler_t t_deriv_gpu;
        bssn_profiler_t t_rhs_gpu;
        bssn_profiler_t t_bdyc_gpu;
        bssn_profiler_t t_mem_handling_gpu;

    }
}


void bssn::timer::profileInfo()
{

    const char separator    = ' ';
    const int nameWidth     = 30;
    const int numWidth      = 10;

    double t_stat;

    t_stat=total_runtime.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"+runtime(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat<<std::endl;

    std::cout<<"CPU Timings"<<std::endl;
    t_stat=t_deriv.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"--deriv(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat<<std::endl;

    t_stat=t_rhs.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"--rhs(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat<<std::endl;

    t_stat=t_bdyc.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"--bdy cond.(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat<<std::endl;

    t_stat=t_cpu_runtime.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"--cpu_runtime(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat<<std::endl;


    std::cout<<"GPU Timings"<<std::endl;
    t_stat=t_mem_handling_gpu.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"--mem_handling(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat<<std::endl;

    t_stat=t_deriv_gpu.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"--deriv(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat<<std::endl;

    t_stat=t_rhs_gpu.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"--rhs(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat<<std::endl;

    t_stat=t_bdyc_gpu.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"--bdy cond.(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat<<std::endl;

    t_stat=t_gpu_runtime.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"--gpu_runtime(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat<<std::endl;

    #if isGPU && isCPU
    std::cout<<"Speedup"<<std::endl;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  deriv speedup";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<< t_deriv.seconds/t_deriv_gpu.seconds <<std::endl;

    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  rhs speedup";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<< t_rhs.seconds/t_rhs_gpu.seconds <<std::endl;

    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  Overall speedup";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<< t_cpu_runtime.seconds/t_gpu_runtime.seconds <<std::endl;
    #endif

    #if test
    double times[4] = {(double) total_runtime.seconds, (double) t_deriv.seconds, (double) t_rhs.seconds, (double) t_bdyc.seconds};
    test_file_write::appendToFile("performance.txt", times, 4);
    #endif

}
