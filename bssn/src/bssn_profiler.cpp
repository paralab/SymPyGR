//
// Created by milinda on 10/20/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief simple profiler based on Hari's sort_profiler for bssn application.
*/
//


#include "bssn_profiler.h"

namespace bssn
{
    namespace timer
    {

        bssn_profiler_t::bssn_profiler_t () {
            seconds  = 0.0;   // openmp wall time
            p_flpops =   0;   // papi floating point operations

            _pri_seconds  = 0.0;
            _pri_p_flpops =   0;
        }

        bssn_profiler_t::~bssn_profiler_t () {

        }

        void bssn_profiler_t::start() {
            _pri_seconds = omp_get_wtime();
            flops_papi();
        }

        void bssn_profiler_t::stop() {
            seconds -= _pri_seconds;
            p_flpops -= _pri_p_flpops;
            snap-=_pri_seconds;

            _pri_seconds = omp_get_wtime();
            flops_papi();

            seconds  += _pri_seconds;
            p_flpops += _pri_p_flpops;
            snap     += _pri_seconds;
        }

        void bssn_profiler_t::snapshotclear()
        {
            snap=0.0;
        }

        void bssn_profiler_t::clear() {
            seconds  = 0.0;
            p_flpops =   0;

            _pri_seconds  = 0.0;
            _pri_p_flpops =   0;
        }


        void   bssn_profiler_t::flops_papi() {
#ifdef HAVE_PAPI
            int 		retval;
	float rtime, ptime, mflops;
	retval  = PAPI_flops(&rtime, &ptime, &_pri_p_flpops, &mflops);
	// assert (retval == PAPI_OK);
#else
            _pri_p_flpops =   0;
#endif
        }

        void bssn_profiler_t::setTime(float time) {
            seconds += time;
        }



    }
}

