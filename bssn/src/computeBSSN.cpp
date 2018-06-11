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
#include "rhs.h"
#include "rhs_cuda.h"

void data_generation_3D_GPU_Async_supported(const unsigned int blk_lb, const unsigned int blk_up, const unsigned int num_blks, double ** var_in_array, Block * blkList){

    // double ** var_in_array = new double*[num_blks*(blk_up-blk_lb+1)];
    // double ** var_out_array = new double*[num_blks*(blk_up-blk_lb+1)];
    // Block * blkList = new Block[num_blks*(blk_up-blk_lb+1)];

    #pragma omp parallel for num_threads(4)
    for (int index=0; index<num_blks*(blk_up-blk_lb+1); index++){

        int block_no = index%num_blks;
        int lev = ((index/(num_blks))%(blk_up-blk_lb+1))+blk_lb;

        const unsigned int maxDepth=12;

        Block & blk=blkList[index];
        blk=Block(0, 0, 0, 2*lev, lev, maxDepth);
        blk.offset=0;
        const unsigned long unzip_dof=(blk.node1D_x*blk.node1D_y*blk.node1D_z);

        double * var_in_per_block = new double[unzip_dof*BSSN_NUM_VARS];
        double * var_out_per_block = new double[unzip_dof*BSSN_NUM_VARS];
        var_in_array[index] = var_in_per_block;

        double coord[3];
        double u[BSSN_NUM_VARS];
        double x,y,z,hx,hy,hz;
        unsigned int offset;
        unsigned int size_x,size_y,size_z;

        x=(double)blk.x;
        y=(double)blk.y;
        z=(double)blk.z;

        hx=0.001; 
        hy=0.001; 
        hz=0.001; 

        offset=blk.offset;
        size_x=blk.node1D_x;
        size_y=blk.node1D_y;
        size_z=blk.node1D_z;

        for(unsigned int k=0;k<blk.node1D_z;k++){
            for(unsigned int j=0;j<blk.node1D_y;j++){
                for(unsigned int i=0;i<blk.node1D_x;i++){
                    coord[0]=x+i*hx;
                    coord[1]=y+j*hy;
                    coord[2]=z+k*hz;

                    initial_data(u,coord);

                    for(unsigned int var=0; var<BSSN_NUM_VARS; var++){
                        var_in_per_block[var*unzip_dof+offset+k*size_y*size_x+j*size_y+i]=u[var];
                    }
                }
            }
        }
    }
}

void data_generation_2D(const unsigned int blk_lb, const unsigned int blk_up, const unsigned int num_blks, double ** var_in, double ** var_out, Block * blkList){

    const unsigned int total_blks=num_blks*(blk_up-blk_lb+1);
    const unsigned int maxDepth=12;

    // var_in = new double*[BSSN_NUM_VARS];
    // var_out = new double*[BSSN_NUM_VARS];
    // blkList = new Block[total_blks];

    unsigned long unzipSz=0;
    unsigned int index=0;
    for(unsigned int lev=blk_lb; lev<=blk_up; lev++){
       for(unsigned int i=0;i<num_blks;i++){
           index = (lev-blk_lb)*num_blks+i;

           blkList[index]=Block(0, 0, 0, 2*lev, lev, maxDepth);

           Block & blk=blkList[index];
           blk.offset=unzipSz;
           unzipSz+=(blk.node1D_x*blk.node1D_y*blk.node1D_z);
       }
    }
    const unsigned long unzip_dof=unzipSz;

    for(unsigned int i=0;i<BSSN_NUM_VARS;i++){
        var_in[i] = new double[unzip_dof];
        var_out[i] = new double[unzip_dof];
    }

    #pragma omp parallel for
    for(unsigned int blk=0;blk<total_blks;blk++){
        
        double coord[3];
        double u[BSSN_NUM_VARS];

        Block tmpBlock=blkList[blk];
        double x=(double)tmpBlock.x;
        double y=(double)tmpBlock.y;
        double z=(double)tmpBlock.z;

        double hx=0.001;
        double hy=0.001;
        double hz=0.001;

        unsigned int offset=tmpBlock.offset;
        unsigned int size_x=tmpBlock.node1D_x;
        unsigned int size_y=tmpBlock.node1D_y;
        unsigned int size_z=tmpBlock.node1D_z;

        for(unsigned int k=0;k<tmpBlock.node1D_z;k++){
            for(unsigned int j=0;j<tmpBlock.node1D_y;j++){
                for(unsigned int i=0;i<tmpBlock.node1D_x;i++){
                    coord[0]=x+i*hx;
                    coord[1]=y+j*hy;
                    coord[2]=z+k*hz;

                    initial_data(u,coord);

                    for(unsigned int var=0;var<BSSN_NUM_VARS;var++){
                        var_in[var][offset+k*size_y*size_x+j*size_y+i]=u[var];
                        var_out[var][offset+k*size_y*size_x+j*size_y+i]=0;
                    }
                }
            }
        }
    }
}

double ** GPU_Async_Iteration_Wise(const unsigned int blk_lb, const unsigned int blk_up, const unsigned int num_blks, double ** var_in_array, Block * blkList){
    // This function accepts only 3D input format

    double ** dev_var_in_array = new double*[num_blks*(blk_up-blk_lb+1)];
    double ** dev_var_out_array = new double*[num_blks*(blk_up-blk_lb+1)];
    double ** host_var_out_array = new double*[num_blks*(blk_up-blk_lb+1)];

    Block blk;
    int current_block = 0;
    int init_block = 0;
    double current_usage, fixed_usage;
    int total_blks = num_blks*(blk_up-blk_lb+1);
    // Check number of blocks can handle at onece

    cudaStream_t stream;
    cudaStream_t streams[numberOfStreams];

    int max_unzip_dof;
    for (int index=init_block; index< total_blks; index++){
        blk = blkList[index];
        max_unzip_dof = blk.node1D_x*blk.node1D_y*blk.node1D_z;
        CHECK_ERROR(cudaMalloc((void**)&dev_var_in_array[index], max_unzip_dof*BSSN_NUM_VARS*sizeof(double)), "dev_var_in_array[index]");
        CHECK_ERROR(cudaMalloc((void**)&dev_var_out_array[index], max_unzip_dof*BSSN_NUM_VARS*sizeof(double)), "dev_var_out_array[index]");
        CHECK_ERROR(cudaMallocHost((void**)&host_var_out_array[index], max_unzip_dof*BSSN_NUM_VARS*sizeof(double)), "host_var_out_array[index]");
    }
    
    int size = numberOfStreams * max_unzip_dof * sizeof(double);

    #include "bssnrhs_cuda_variable_malloc.h"
    #include "bssnrhs_cuda_variable_malloc_adv.h"
    #include "bssnrhs_cuda_malloc.h"
    #include "bssnrhs_cuda_malloc_adv.h"
    
    
    for (int index=0; index<numberOfStreams; index++){
        CHECK_ERROR(cudaStreamCreate(&streams[index]), "cudaStream creation");
    }

    for(unsigned int index=0;index<total_blks;index++) {

        cudaStream_t stream;
        int streamIndex = index % numberOfStreams;
        stream = streams[streamIndex];

        int block_no = index%num_blks;
        int level = ((index/(num_blks))%(blk_up-blk_lb+1))+blk_lb;

        Block & blk=blkList[index];

        unsigned long unzip_dof=(blk.node1D_x*blk.node1D_y*blk.node1D_z);
        unsigned int offset=blk.offset;

        unsigned int sz[3];
        sz[0]=blk.node1D_x; 
        sz[1]=blk.node1D_y;
        sz[2]=blk.node1D_z;

        unsigned int  bflag=0;

        double dx=0.1;
        double dy=0.1;
        double dz=0.1;

        double ptmin[3], ptmax[3];
        ptmin[0]=0.0;
        ptmin[1]=0.0;
        ptmin[2]=0.0;

        ptmax[0]=1.0;
        ptmax[1]=1.0;
        ptmax[2]=1.0;

        std::cout << "GPU - Block no: " << std::setw(2) << index << "  Total Points: " << std::setw(7) << unzip_dof << "  level: " << level << std::endl;

        CHECK_ERROR(cudaMemcpyAsync(dev_var_in_array[index], var_in_array[index], BSSN_NUM_VARS*unzip_dof*sizeof(double), cudaMemcpyHostToDevice, stream), "dev_var_in_array[index] cudaMemcpyHostToDevice");

        cuda_bssnrhs(dev_var_out_array[index], dev_var_in_array[index], unzip_dof , offset, ptmin, ptmax, sz, bflag, stream, NULL,
            #include "list_of_args_per_blk.h"
        );

        CHECK_ERROR(cudaMemcpyAsync(host_var_out_array[index], dev_var_out_array[index], BSSN_NUM_VARS*unzip_dof*sizeof(double), cudaMemcpyDeviceToHost, stream), "dev_var_out_array[index] cudaMemcpyDeviceToHost");
    }

    #include "bssnrhs_cuda_mdealloc.h"
    #include "bssnrhs_cuda_mdealloc_adv.h"

    for (int index=0; index< total_blks; index++){
        delete [] var_in_array[index];
        CHECK_ERROR(cudaFree(dev_var_in_array[index]), "dev_var_in_array[index] cudaFree");
        CHECK_ERROR(cudaFree(dev_var_out_array[index]), "dev_var_out_array[index] cudaFree");
    }

    return host_var_out_array;
}

void CPU_sequence(const unsigned int blk_lb, const unsigned int blk_up, const unsigned int num_blks, double ** var_in, double ** var_out, Block * blkList){
    // This function accepts only 2D version of the input

    double ptmin[3], ptmax[3];
    unsigned int sz[3];
    unsigned int bflag;
    unsigned int offset;
    double dx, dy, dz;

    for(unsigned int blk=0;blk<num_blks*(blk_up-blk_lb+1);blk++){

        int block_no = blk%num_blks;
        int level = ((blk/(num_blks))%(blk_up-blk_lb+1))+blk_lb;

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

        std::cout << "CPU - Block no: " << std::setw(2) << blk << "  Total Points: " << std::setw(7) << total_points << "  level: " << level << std::endl;
        
        bssnrhs(var_out, (const double **)var_in, offset, ptmin, ptmax, sz, bflag);
    }
}

int main (int argc, char** argv){
    /**
     *
     * parameters:
     * blk_lb: block element 1d lower bound (int)
     * blk_up: block element 1d upper bound (int) (blk_up>=blk_lb)
     * numblks : number of blocks needed for each block sizes. (total_blks= (blk_up-blk_lb+1)*numblks)
     *
     * */

    if(argc<2){
        std::cout<<"Usage: "<<argv[0]<<"blk_lb blk_up numblks"<<std::endl;
        exit(0);
    }

    // asyncDataTransfetTest();
    unsigned int blk_lb=atoi(argv[1]);
    unsigned int blk_up=atoi(argv[2]);
    unsigned int num_blks=atoi(argv[3]);

    double ** var_in_array = new double*[num_blks*(blk_up-blk_lb+1)];
    double ** var_in = new double*[BSSN_NUM_VARS];
    double ** var_out = new double*[BSSN_NUM_VARS];
    Block * blkList = new Block[num_blks*(blk_up-blk_lb+1)];
    
    #if isGPU
    data_generation_3D_GPU_Async_supported(blk_lb, blk_up, num_blks, var_in_array, blkList);
    #include "rhs_cuda.h"

    bssn::timer::t_gpu_runtime.start();
    double ** var_out_array = GPU_Async_Iteration_Wise(blk_lb, blk_up, num_blks, var_in_array, blkList);
    bssn::timer::t_gpu_runtime.stop();
    #endif

    std::cout << "" << std::endl;

    #if isCPU
    data_generation_2D(blk_lb, blk_up, num_blks, var_in, var_out, blkList);
    #include "rhs.h"

    bssn::timer::t_cpu_runtime.start();
    CPU_sequence(blk_lb, blk_up, num_blks, var_in, var_out, blkList);
    bssn::timer::t_cpu_runtime.stop();
    #endif

    #if test
    // Verify outputs
    for (int blk=0; blk<num_blks*(blk_up-blk_lb+1); blk++){
        for(int bssn_var=0; bssn_var<BSSN_NUM_VARS; bssn_var++){
            int sizeofBlock = blkList[blk].node1D_x * blkList[blk].node1D_y * blkList[blk].node1D_z;
            for (int pointInd=0; pointInd<sizeofBlock; pointInd++){
                double diff = var_out_array[blk][bssn_var*sizeofBlock+pointInd] - var_out[bssn_var][blkList[blk].offset+pointInd];
                if (fabs(diff)>threshold){
                    const char separator    = ' ';
                    const int nameWidth     = 6;
                    const int numWidth      = NUM_DIGITS+10;

                    std::cout << std::left << std::setw(nameWidth) << setfill(separator) << "GPU: ";
                    std::cout <<std::setprecision(NUM_DIGITS)<< std::left << std::setw(numWidth) << setfill(separator)  << var_out_array[blk][bssn_var*sizeofBlock+pointInd];

                    std::cout << std::left << std::setw(nameWidth) << setfill(separator) << "CPU: ";
                    std::cout <<std::setprecision(NUM_DIGITS)<< std::left << std::setw(numWidth) << setfill(separator)  << var_out[bssn_var][blkList[blk].offset+pointInd];

                    std::cout << std::left << std::setw(nameWidth) << setfill(separator) << "DIFF: ";
                    std::cout <<std::setprecision(NUM_DIGITS)<< std::left << std::setw(numWidth) << setfill(separator)  << diff << std::endl;
                    exit(0);
                }
            }
        }
    }
    #endif

    #if isGPU
    //Free host memory
    for (int blk=0; blk<num_blks*(blk_up-blk_lb+1); blk++){
        CHECK_ERROR(cudaFreeHost(var_out_array[blk]), "free host memory");
    }
    #endif
    
    bssn::timer::profileInfo();
    return 0;
}