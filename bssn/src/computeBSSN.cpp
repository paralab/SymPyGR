//
// Created by milinda on 1/15/18.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief
*/
//
#include "computeBSSN.h" 


int main (int argc, char** argv)
{
    std::cout << "Let's start debuging" <<std::endl;


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
    unsigned int index=0;
    for(unsigned int lev=blk_lb; lev<=blk_up; lev++)
       for(unsigned int i=0;i<num_blks;i++)
       {
           index = (lev-blk_lb)*num_blks+i;

           blkList[index]=Block((1u<<lev)); // lu<<lev = 2**lev

           Block & blk=blkList[index];
           blk.offset=unzipSz;
           unzipSz+=(blk.node1D_x*blk.node1D_y*blk.node1D_z);
       }
    const unsigned long unzip_dof=unzipSz;

    // //2. allocate memory for bssn computation.
    // double ** var_in=new double*[BSSN_NUM_VARS];
    // double ** var_out=new double*[BSSN_NUM_VARS];
    // Allocate memory on GPU for bssn computation
    double * dev_var_in;
    double * dev_var_out;
    cudaMalloc((void**)&dev_var_in, unzip_dof*BSSN_NUM_VARS*sizeof(double));
    cudaMalloc((void**)&dev_var_out, unzip_dof*BSSN_NUM_VARS*sizeof(double));

    // for(unsigned int i=0;i<BSSN_NUM_VARS;i++)
    // {
    //     var_in[i]=new double[unzip_dof];
    //     var_out[i]=new double[unzip_dof];

    //     for(unsigned int j=0;j<unzip_dof;j++)
    //     {
    //         // some random initialization.
    //         var_in[i][j]=sqrt(2)*0.001*j;
    //         var_out[i][j]=0.0;
    //     }
    // }

    // GPU usage requirement
    double * host_var_in = new double[BSSN_NUM_VARS*unzip_dof];
    double * host_var_out = new double[BSSN_NUM_VARS*unzip_dof];

    unsigned int j = 0;
    for(unsigned int i=0;i<BSSN_NUM_VARS*unzip_dof;i++)
    {
        if (j==unzip_dof){
            j=0;
        }
        // some random initialization.
        host_var_in[i]=sqrt(2)*0.001*j;
        host_var_out[i]=0.0;     
        j++;  
    }
    
    

    // double *dev_var_in;
    // double *dev_var_out;
    // cudaMalloc((void**)&dev_var_in, BSSN_NUM_VARS*unzip_dof*sizeof(double));
    // cudaMalloc((void**)&dev_var_out, BSSN_NUM_VARS*unzip_dof*sizeof(double));
    cudaMemcpy(dev_var_in, host_var_in, BSSN_NUM_VARS*unzip_dof*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_var_out, host_var_out, BSSN_NUM_VARS*unzip_dof*sizeof(double), cudaMemcpyHostToDevice);
    // std::cout << cudaStatus <<std::endl;


    // std::cout << dev_var_out <<std::endl;

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

        // bssnrhs(dev_var_out, dev_var_in, unzip_dof , offset, ptmin, ptmax, sz, bflag); // required to send the output array also

        cuda_bssnrhs(dev_var_out, dev_var_in, unzip_dof , offset, ptmin, ptmax, sz, bflag);

    }
    //-- timer end
    // (time this part of the code. )


    // for(unsigned int i=0;i<BSSN_NUM_VARS;i++)
    // {
    //     delete [] var_in[i];
    //     delete [] var_out[i];
    // }

    delete [] blkList;
    // delete [] var_in;
    // delete [] var_out;

    bssn::timer::total_runtime.stop();

    bssn::timer::profileInfo();

    // Free up memory
    cudaFree(&dev_var_in);
    cudaFree(&dev_var_out);

    return 0;

}
