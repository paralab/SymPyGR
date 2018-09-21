//
// Created by milinda on 10/20/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief simple profiler based on Hari's sort_profiler for bssn application.
*/
//

/**
 * Created on: Sep 21, 2018
 * 		Adapted by: Akila, Eranga, Eminda, Ruwan
 **/

#ifndef PROFILER_H
#define PROFILER_H

#include <iostream>
#include <iomanip>
#include <fstream>

#ifdef HAVE_PAPI
#include <papi.h>
#endif

#include <omp.h>



class bssn_profiler_t
{
public:
    bssn_profiler_t ();
    virtual ~bssn_profiler_t ();

    void start();
    void stop();
    void clear();
    void snapshotclear();

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


#endif 
