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
#include <iostream>
#include "block.h"

namespace bssn
{
    namespace timer
    {

        extern bssn_profiler_t 	total_runtime;
        extern bssn_profiler_t	t_deriv;
        extern bssn_profiler_t	t_rhs;
        extern bssn_profiler_t  t_bdyc;



        void initialize();

        void profileInfo();
    }
}



#endif //BSSN_PROFILE_PARAM_H
