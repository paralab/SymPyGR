//
// Created by milinda on 1/15/18.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief
*/
//

#include <iostream>
#include "computeBSSN.h"


int main (int argc, char** argv)
{


    /**
     *
     * parameters:
     * blk_lb: block element 1d lower bound (int)
     * blk_up: block element 1d upper bound (int) (blk_up>=blk_lb)
     * numblks : number of blocks needed for each block sizes. (total_blks= (blk_up-blk_lb+1)*numblks)
     *
     * */


    if(argc<2)
    {
        std::cout<<"Usage: "<<argv[0]<<"blk_lb blk_up numblks"<<std::endl;
        exit(0);
    }


    // initialize profile counters.
    bssn::timer::total_runtime.start();
    bssn::timer::t_deriv.start();
    bssn::timer::t_rhs.start();
    bssn::timer::t_bdyc.start();



    bssn::timer::total_runtime.start();

    unsigned int num_blks=100;
    unsigned int blk_lb=0;
    unsigned int blk_up=5;

    blk_lb=atoi(argv[1]);
    blk_up=atoi(argv[2]);
    num_blks=atoi(argv[3]);

    const unsigned int total_blks=num_blks*(blk_up-blk_lb+1);

    //1. setup the blk offsets.
    Block * blkList=new Block[total_blks];
    unsigned long unzipSz=0;
    for(unsigned int lev=blk_lb; lev<=blk_up; lev++)
       for(unsigned int i=0;i<num_blks;i++)
       {
           blkList[(lev-blk_lb)*num_blks+i]=Block((1u<<lev));
           blkList[(lev-blk_lb)*num_blks+i].offset=unzipSz;
           unzipSz+=((blkList[(lev-blk_lb)*num_blks+i].node1D_x)*(blkList[(lev-blk_lb)*num_blks+i].node1D_y)*(blkList[(lev-blk_lb)*num_blks+i].node1D_z));

       }

    const unsigned long unzip_dof=unzipSz;

    //2. allocate memory for bssn computation.
    double ** var_in=new double*[BSSN_NUM_VARS];
    double ** var_out=new double*[BSSN_NUM_VARS];

    double ** pre_computed = new double*[20];

    for(unsigned int i=0;i<BSSN_NUM_VARS;i++)
    {
        var_in[i]=new double[unzip_dof];
        var_out[i]=new double[unzip_dof];
        pre_computed[i] = new double[unzip_dof];

        for(unsigned int j=0;j<unzip_dof;j++)
        {
            // some random initialization.
            var_in[i][j]=sqrt(2)*0.001*j;
            var_out[i][j]=0.0;
        }


    }

    unsigned int offset;
    double ptmin[3], ptmax[3];
    unsigned int sz[3];
    unsigned int bflag;
    double dx,dy,dz;


     //----- timer begin:
    //3. perform computation.
    for(unsigned int blk=0;blk<total_blks;blk++)
    {
        offset=blkList[blk].offset;
        sz[0]=blkList[blk].node1D_x;
        sz[1]=blkList[blk].node1D_y;
        sz[2]=blkList[blk].node1D_z;

        bflag=0; // indicates if the block is bdy block.

        dx=0.1;
        dy=0.1;
        dz=0.1;

        ptmin[0]=0.0;
        ptmin[1]=0.0;
        ptmin[2]=0.0;

        ptmax[0]=1.0;
        ptmax[1]=1.0;
        ptmax[2]=1.0;

        bssnrhs(var_out, (const double **)var_in, pre_computed, offset, ptmin, ptmax, sz, bflag);

    }
    //-- timer end
    // (time this part of the code. )


    for(unsigned int i=0;i<BSSN_NUM_VARS;i++)
    {
        delete [] var_in[i];
        delete [] var_out[i];
    }

    delete [] blkList;
    delete [] var_in;
    delete [] var_out;

    bssn::timer::total_runtime.stop();

    bssn::timer::profileInfo();

    return 0;



}