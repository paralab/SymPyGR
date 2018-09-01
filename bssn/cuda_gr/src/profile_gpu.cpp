//
// Created by milinda on 8/30/18.
//

#include "profile_gpu.h"

namespace cuda
{
    namespace profile
    {

        /** overall gpu exuxution time*/
        bssn_profiler_t t_overall;

        /**host to device communication*/
        bssn_profiler_t t_H2D_Comm;
        /**device to host communication*/
        bssn_profiler_t t_D2H_Comm;

        /**memory allocation deallocation time*/
        bssn_profiler_t t_cudaMalloc_derivs;

        /**deriv computation time*/
        bssn_profiler_t t_derivs;




    } // end of namespace profile
} // end of namespace cuda


void cuda::profile::initialize()
{
    t_overall.start();
    t_H2D_Comm.start();
    t_D2H_Comm.start();
    t_cudaMalloc_derivs.start();
    t_derivs.start();

}




void cuda::profile::printOutput(const std::vector<ot::Block>& localBlkList) {


    const unsigned int NUM_DERIV_OPS=3396;
    unsigned int unzipInternal=0;

    const char separator    = ' ';
    const int nameWidth     = 30;
    const int numWidth      = 10;

    for(unsigned int blk=0;blk<localBlkList.size();blk++)
        unzipInternal+=(localBlkList[blk].get1DArraySize()-2*localBlkList[blk].get1DPadWidth())*(localBlkList[blk].get1DArraySize()-2*localBlkList[blk].get1DPadWidth())*(localBlkList[blk].get1DArraySize()-2*localBlkList[blk].get1DPadWidth());

    double t_stat;
    t_stat=t_overall.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"+cuda_overall(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat<<std::endl;;



    t_stat=t_H2D_Comm.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  -H2D(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat<<std::endl;

    t_stat=t_D2H_Comm.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  -D2H(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat<<std::endl;


    t_stat=t_cudaMalloc_derivs.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  -deriv(malloc)(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat<<std::endl;

    t_stat=t_derivs.seconds;
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  -deriv(compute)(s)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat<<std::endl;

    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  -deriv(compute)(flops)";
    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<((unzipInternal*NUM_DERIV_OPS)/t_stat)<<std::endl;


    return;


}


