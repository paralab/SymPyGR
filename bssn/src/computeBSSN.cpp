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

    unsigned int num_blks=100;
    unsigned int blk_lb=0;
    unsigned int blk_up=5;

    blk_lb=atoi(argv[1]);
    blk_up=atoi(argv[2]);
    num_blks=atoi(argv[3]);

    const unsigned int total_blks= num_blks*(blk_up-blk_lb+1);

    //1. setup the blk offsets.
    Block * blkList=new Block[total_blks];

    unsigned long unzipSz=0;
    unsigned int index=0;
    const unsigned int maxDepth=12; //blk_up+1;

    for(unsigned int lev=blk_lb; lev<=blk_up; lev++)
       for(unsigned int i=0;i<num_blks;i++)
       {
           index = (lev-blk_lb)*num_blks+i;

           blkList[index]=Block(0, 0, 0, 2*lev, lev, maxDepth);

           Block & blk=blkList[index];
           blk.offset=unzipSz;
           unzipSz+=(blk.node1D_x*blk.node1D_y*blk.node1D_z);
       }
    const unsigned long unzip_dof=unzipSz;

    double coord[3];
    double u[BSSN_NUM_VARS];
    double x,y,z,hx,hy,hz;
    unsigned int offset;
    unsigned int size_x,size_y,size_z;
    Block tmpBlock;

    // 2. a. allocate memory for bssn computation on CPU.
    #if isCPU || test
    bssn::timer::t_cpu_runtime.start();

    double ** var_in=new double*[BSSN_NUM_VARS];
    double ** var_out=new double*[BSSN_NUM_VARS];


    for(unsigned int i=0;i<BSSN_NUM_VARS;i++)
    {
        var_in[i] = new double[unzip_dof];
        var_out[i] = new double[unzip_dof];
    }

    bssn::timer::t_cpu_runtime.stop();
    #endif

    // 2. b. Allocate memory on GPU for bssn computation
    #if isGPU || test
    bssn::timer::t_gpu_runtime.start();
    

    cudaError_t cudaStatus;
    // Choose which GPU to run on, change this on a multi-GPU system.
     cudaStatus = cudaSetDevice(0);
     if (cudaStatus != cudaSuccess) {
         fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?\n");
         return 0;
     }

    bssn::timer::t_mem_handling_gpu.start();
    bssn::timer::t_rhs_gpu.start();
    double * dev_var_in;
    double * dev_var_out;
    cudaStatus = cudaMalloc((void**)&dev_var_in, unzip_dof*BSSN_NUM_VARS*sizeof(double));
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "var_in cudaMalloc failed!\n"); return 0;}

    cudaStatus = cudaMalloc((void**)&dev_var_out, unzip_dof*BSSN_NUM_VARS*sizeof(double));
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "var_out cudaMalloc failed!\n"); return 0;}
    bssn::timer::t_rhs_gpu.stop();
    bssn::timer::t_mem_handling_gpu.stop();

    // GPU usage requirement
    double * host_var_in = new double[BSSN_NUM_VARS*unzip_dof];
    double * host_var_out = new double[BSSN_NUM_VARS*unzip_dof];

    bssn::timer::t_gpu_runtime.stop();
    #endif

    for(unsigned int blk=0;blk<total_blks;blk++)
    {
        tmpBlock=blkList[blk];
        x=(double)tmpBlock.x;
        y=(double)tmpBlock.y;
        z=(double)tmpBlock.z;

        hx=0.001;  //(1u<<(maxDepth-tmpBlock.regLevel))/(double)(ELE_ORDER);//(1u<<(tmpBlock.regLevel-tmpBlock.blkLevel))/(double)ELE_ORDER;
        hy=0.001;  //(1u<<(maxDepth-tmpBlock.regLevel))/(double)(ELE_ORDER);//(1u<<(tmpBlock.regLevel-tmpBlock.blkLevel))/(double)ELE_ORDER;
        hz=0.001;  //(1u<<(maxDepth-tmpBlock.regLevel))/(double)(ELE_ORDER);//(1u<<(tmpBlock.regLevel-tmpBlock.blkLevel))/(double)ELE_ORDER;

        offset=tmpBlock.offset;
        size_x=tmpBlock.node1D_x;
        size_y=tmpBlock.node1D_y;
        size_z=tmpBlock.node1D_z;

        for(unsigned int k=0;k<tmpBlock.node1D_z;k++)
            for(unsigned int j=0;j<tmpBlock.node1D_y;j++)
                for(unsigned int i=0;i<tmpBlock.node1D_x;i++)
                {
                    coord[0]=x+i*hx;
                    coord[1]=y+j*hy;
                    coord[2]=z+k*hz;

                    initial_data(u,coord);

                    for(unsigned int var=0;var<BSSN_NUM_VARS;var++)
                    {
                        #if isCPU || test
                        bssn::timer::t_cpu_runtime.start();

                        var_in[var][offset+k*size_y*size_x+j*size_y+i]=u[var];
                        var_out[var][offset+k*size_y*size_x+j*size_y+i]=0;

                        bssn::timer::t_cpu_runtime.stop();
                        #endif

                        #if isGPU || test
                        bssn::timer::t_gpu_runtime.start();

                        host_var_in[var*unzip_dof+offset+k*size_y*size_x+j*size_y+i]=u[var];
                        host_var_out[var*unzip_dof+offset+k*size_y*size_x+j*size_y+i]=0;

                        bssn::timer::t_gpu_runtime.stop();
                        #endif
                    }
                }
    }

    #if isGPU || test
    bssn::timer::t_gpu_runtime.start();
    bssn::timer::t_mem_handling_gpu.start();
    bssn::timer::t_rhs_gpu.start();

    cudaStatus = cudaMemcpy(dev_var_in, host_var_in, BSSN_NUM_VARS*unzip_dof*sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "var_in cudaMemcpy failed!\n"); return 0;}
    cudaStatus = cudaMemcpy(dev_var_out, host_var_out, BSSN_NUM_VARS*unzip_dof*sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "var_in cudaMemcpy failed!\n"); return 0;}

    bssn::timer::t_rhs_gpu.stop();
    bssn::timer::t_mem_handling_gpu.stop();
    bssn::timer::t_gpu_runtime.stop();
    #endif


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

        std::cout << "Block no: " << blk << " Total Points: " << total_points << std::endl;

        // CPU bssnrhs call
        #if isCPU || test
        bssn::timer::t_cpu_runtime.start();

        #include "rhs.h"
        bssnrhs(var_out, (const double **)var_in, offset, ptmin, ptmax, sz, bflag);

        bssn::timer::t_cpu_runtime.stop();
        #endif

          // CUDA_bssnrhs call
        #if isGPU || test
        bssn::timer::t_gpu_runtime.start();

        #include "rhs_cuda.h"
        cuda_bssnrhs(dev_var_out, dev_var_in, unzip_dof , offset, ptmin, ptmax, sz, bflag);

        bssn::timer::t_gpu_runtime.stop();
        #endif
    }
    printf("-------------------------\n");
    //-- timer end
    // (time this part of the code. )

    #if isGPU || test
    bssn::timer::t_gpu_runtime.start();
    bssn::timer::t_mem_handling_gpu.start();
    bssn::timer::t_rhs_gpu.start();

    // Copy back dev_var_out to host_var_out
    cudaStatus = cudaMemcpy(host_var_out, dev_var_out, BSSN_NUM_VARS*unzip_dof*sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "dev_var_out to host_var_out cudaMemcpy failed!\n"); return 0;}

    bssn::timer::t_rhs_gpu.stop();
    bssn::timer::t_mem_handling_gpu.stop();
    bssn::timer::t_gpu_runtime.stop();
    #endif

    #if test
    // this will verify the GPU output with CPU output
    unsigned int error_count = 0;
    for(unsigned int i=0;i<BSSN_NUM_VARS;i++){
        for(unsigned int j=0; j<unzip_dof; j++){
            unsigned int abs_index = i*unzip_dof + j;
            double diff = var_out[i][j] - host_var_out[abs_index];
            if (fabs(diff)>threshold){
                error_count++;
                const char separator    = ' ';
                const int nameWidth     = 6;
                const int numWidth      = NUM_DIGITS+10;

                std::cout << std::left << std::setw(nameWidth) << setfill(separator) << "GPU: ";
                std::cout <<std::setprecision(NUM_DIGITS)<< std::left << std::setw(numWidth) << setfill(separator)  << host_var_out[abs_index];

                std::cout << std::left << std::setw(nameWidth) << setfill(separator) << "CPU: ";
                std::cout <<std::setprecision(NUM_DIGITS)<< std::left << std::setw(numWidth) << setfill(separator)  <<var_out[i][j];

                std::cout << std::left << std::setw(nameWidth) << setfill(separator) << "DIFF: ";
                std::cout <<std::setprecision(NUM_DIGITS)<< std::left << std::setw(numWidth) << setfill(separator)  << diff;

                std::cout << std::left << std::setw(nameWidth) << setfill(separator) << "BSSN_VAR: ";
                std::cout << std::left << std::setw(4) << setfill(separator)  << i;

                std::cout << std::left << std::setw(nameWidth) << setfill(separator) << "PP(valid when total_blks=1): ";
                std::cout << std::left << std::setw(numWidth) << setfill(separator)  << j << std::endl;
            }
        }
    }
    std::cout << "Total errors: " << error_count << " | Total number of dof: " << BSSN_NUM_VARS*unzip_dof << std::endl;
    #endif

    // CPU code
    #if isCPU || test
    bssn::timer::t_cpu_runtime.start();

    for(unsigned int i=0;i<BSSN_NUM_VARS;i++)
    {
        delete [] var_in[i];
        delete [] var_out[i];
    }
    delete [] var_in;
    delete [] var_out;

    bssn::timer::t_cpu_runtime.stop();
    #endif

    // GPU code
    #if isGPU || test
    bssn::timer::t_gpu_runtime.start();
    
    delete [] host_var_in;
    delete [] host_var_out;
    
    // Free up GPU memory
    bssn::timer::t_mem_handling_gpu.start();
    cudaFree(dev_var_in);
    cudaFree(dev_var_out);
    bssn::timer::t_mem_handling_gpu.stop();

    bssn::timer::t_gpu_runtime.stop();
    #endif

    delete [] blkList;

    bssn::timer::total_runtime.stop();

    bssn::timer::profileInfo();

    return 0;
}
