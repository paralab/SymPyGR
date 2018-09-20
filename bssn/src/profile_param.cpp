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
        bssn_profiler_t	t_deriv;
        bssn_profiler_t	t_rhs;
        bssn_profiler_t t_bdyc;

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
}