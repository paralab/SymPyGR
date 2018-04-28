//
// Created by milinda on 1/15/18.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief
*/
//
#include "computeBSSN.h" 
#include "test_param.h"


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

    // 2. a. allocate memory for bssn computation on CPU.
    #if isCPU
    double ** var_in=new double*[BSSN_NUM_VARS];
    double ** var_out=new double*[BSSN_NUM_VARS];

    for(unsigned int i=0;i<BSSN_NUM_VARS;i++)
    {
        var_in[i]=new double[unzip_dof];
        var_out[i]=new double[unzip_dof];

        for(unsigned int j=0;j<unzip_dof;j++)
        {
            // some random initialization.
            var_in[i][j]=sin((j%360)*PI/180) + sin(((j+180)%360)*PI/180) + cos(((j+60)%360)*PI/180);
            var_out[i][j]=0.0;
        }
    }
    #endif
    
    // 2. b. Allocate memory on GPU for bssn computation
    #if isGPU
    cudaError_t cudaStatus;
    // Choose which GPU to run on, change this on a multi-GPU system.
     cudaStatus = cudaSetDevice(0);
     if (cudaStatus != cudaSuccess) {
         fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?\n");
         return 0;
     }

    double * dev_var_in;
    double * dev_var_out;
    cudaStatus = cudaMalloc((void**)&dev_var_in, unzip_dof*BSSN_NUM_VARS*sizeof(double));
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "var_in cudaMalloc failed!\n"); return 0;}

    cudaStatus = cudaMalloc((void**)&dev_var_out, unzip_dof*BSSN_NUM_VARS*sizeof(double));
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "var_out cudaMalloc failed!\n"); return 0;}

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
        host_var_in[i]=sin((j%360)*PI/180) + sin(((j+180)%360)*PI/180) + cos(((j+60)%360)*PI/180);
        host_var_out[i]=0.0;     
        j++;  
    }
    
    cudaStatus = cudaMemcpy(dev_var_in, host_var_in, BSSN_NUM_VARS*unzip_dof*sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "var_in cudaMemcpy failed!\n"); return 0;}
    #endif

    unsigned int offset;
    double ptmin[3], ptmax[3];
    unsigned int sz[3];
    unsigned int bflag;
    double dx,dy,dz;

    std::cout << "Est. total RAM req (var_in, var_out): " << 2*24*unzip_dof*8/1024/1024 << " mb" << std::endl;

    //----- timer begin:
    //3. perform computation.
    for(unsigned int blk=0;blk<total_blks;blk++)
    {
        offset=blkList[blk].offset;
        sz[0]=blkList[blk].node1D_x; 
        sz[1]=blkList[blk].node1D_y;
        sz[2]=blkList[blk].node1D_z;
        int total_points = sz[0]*sz[1]*sz[2];

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

        // CPU bssnrhs call
        #if isCPU
        #include "rhs.h"
        bssnrhs(var_out, (const double **)var_in, offset, ptmin, ptmax, sz, bflag);
        #endif

          // CUDA_bssnrhs call
        #if isGPU
        #include "rhs_cuda.h"
        printf("blockNo: %d \t| TotalBlockPoints: %d \t| Est. GPU memory req: %d mb\n", blk, sz[0]*sz[1]*sz[2],  (sz[0]*sz[1]*sz[2]*210*8 + 24*4)/1024/1024);
        cuda_bssnrhs(dev_var_out, dev_var_in, unzip_dof , offset, ptmin, ptmax, sz, bflag);
        #endif




    }
    printf("-------------------------\n");
    //-- timer end
    // (time this part of the code. )

    // CPU code
    #if isCPU
    for(unsigned int i=0;i<BSSN_NUM_VARS;i++)
    {
        delete [] var_in[i];
        delete [] var_out[i];
    }
    delete [] var_in;
    delete [] var_out;
    #endif

    // GPU code
    #if isGPU
    delete [] host_var_in;
    delete [] host_var_out;
    // Free up GPU memory
    cudaFree(dev_var_in);
    cudaFree(dev_var_out);
    #endif

    delete [] blkList;

    bssn::timer::total_runtime.stop();

    bssn::timer::profileInfo();

    return 0;
}
