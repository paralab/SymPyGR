//
// Created by milinda on 10/20/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief simple profiler based on Hari's sort_profiler for bssn application.
*/
//

#ifndef SFCSORTBENCH_BSSN_PROFILER_H
#define SFCSORTBENCH_BSSN_PROFILER_H

#include <iostream>
#include <iomanip>
#include <fstream>

#ifdef HAVE_PAPI
#include <papi.h>
#endif

#include <omp.h>

namespace bssn
{

    namespace timer
    {


        class bssn_profiler_t
        {
        public:
            bssn_profiler_t ();
            virtual ~bssn_profiler_t ();

            void start();
            void stop();
            void clear();
            void snapshotclear();
            void setTime(float time);

        public:
            long double	seconds;  // openmp wall time
            long long p_flpops; // papi floating point operations
            long double snap; // snap shot of the cumilative time.

        private:
            void  flops_papi();

        protected:
            long double	  _pri_seconds;  // openmp wall time
            long long _pri_p_flpops; // papi floating point operations
        };



    } // end of namespace

} //end of namespace






#endif //SFCSORTBENCH_BSSN_PROFILER_H
