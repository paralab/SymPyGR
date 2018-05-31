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


void data_generation_3D(const unsigned int blk_lb, const unsigned int blk_up, const unsigned int num_blks, double ** var_in_array, double ** var_out_array, Block * blkList){

    // double ** var_in_array = new double*[num_blks*(blk_up-blk_lb+1)];
    // double ** var_out_array = new double*[num_blks*(blk_up-blk_lb+1)];
    // Block * blkList = new Block[num_blks*(blk_up-blk_lb+1)];

    #pragma omp parallel for
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
        var_out_array[index] = var_out_per_block;

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
                        var_out_per_block[var*unzip_dof+offset+k*size_y*size_x+j*size_y+i]=0;
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

void GPU_2D(const unsigned int blk_lb, const unsigned int blk_up, const unsigned int num_blks, double ** var_in, double ** var_out, Block * blkList){
    // This function accepts only 2D input format
    unsigned long unzip_dof = 0;
    for (int index=0; index<num_blks*(blk_up-blk_lb+1); index++){
        unzip_dof += blkList[index].node1D_x * blkList[index].node1D_y * blkList[index].node1D_z;
    }

    double * dev_var_in;
    double * dev_var_out;

    CHECK_ERROR(cudaSetDevice(0), "cudaSetDevice in computeBSSN");

    CHECK_ERROR(cudaMalloc((void**)&dev_var_in, unzip_dof*BSSN_NUM_VARS*sizeof(double)), "dev_var_in");
    CHECK_ERROR(cudaMalloc((void**)&dev_var_out, unzip_dof*BSSN_NUM_VARS*sizeof(double)), "dev_var_out");

    // #pragma omp parallel for DO NOT USE IT MAKES SLOW
    for (int bssn_var=0; bssn_var<BSSN_NUM_VARS; bssn_var++){
        CHECK_ERROR(cudaMemcpy(dev_var_in+bssn_var*unzip_dof, var_in[bssn_var], unzip_dof*sizeof(double), cudaMemcpyHostToDevice), "dev_var_in cudaMemcpyHostToDevice");
    }

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

        bflag=0;

        dx=0.1;
        dy=0.1;
        dz=0.1;

        ptmin[0]=0.0;
        ptmin[1]=0.0;
        ptmin[2]=0.0;

        ptmax[0]=1.0;
        ptmax[1]=1.0;
        ptmax[2]=1.0;

        std::cout << "GPU - Block no: " << std::setw(2) << blk << "  Total Points: " << std::setw(7) << total_points << "  level: " << level << std::endl;

        #include "bssnrhs_cuda_offset_variable_malloc.h"
        #include "bssnrhs_cuda_variable_malloc.h"
        #include "bssnrhs_cuda_variable_malloc_adv.h"
        double * dev_dy_hx;
        double * dev_dy_hy;
        double * dev_dy_hz;
        int * dev_sz;
        int * dev_zero;
        double * dev_pmin;
        double *dev_pmax;
        int * dev_bflag;

        int size = unzip_dof * sizeof(double);

        #include "bssnrhs_cuda_offset_malloc.h"
        CHECK_ERROR(cudaMalloc((void **) &dev_dy_hx, sizeof(double)), "dev_dy_hx");
        CHECK_ERROR(cudaMalloc((void **) &dev_dy_hy, sizeof(double)), "dev_dy_hy");
        CHECK_ERROR(cudaMalloc((void **) &dev_dy_hz, sizeof(double)), "dev_dy_hz");
        CHECK_ERROR(cudaMalloc((void **) &dev_sz, 3*sizeof(int)), "dev_sz");
        CHECK_ERROR(cudaMalloc((void **) &dev_zero, sizeof(int)), "dev_zero");
        CHECK_ERROR(cudaMalloc((void **) &dev_pmin, 3*sizeof(double)), "dev_pmin");
        CHECK_ERROR(cudaMalloc((void **) &dev_pmax, 3*sizeof(double)), "dev_pmax");
        CHECK_ERROR(cudaMalloc((void **) &dev_bflag, sizeof(int)), "dev_bflag");
        #include "bssnrhs_cuda_malloc.h"
        #include "bssnrhs_cuda_malloc_adv.h"

        cuda_bssnrhs(dev_var_out, dev_var_in, unzip_dof , offset, ptmin, ptmax, sz, bflag,
        #include "list_of_args.h"
        , dev_dy_hx, dev_dy_hy, dev_dy_hz, dev_sz, dev_zero, dev_pmin, dev_pmax, dev_bflag
        );

        #include "bssnrhs_cuda_mdealloc.h"
        #include "bssnrhs_cuda_mdealloc_adv.h"
        #include "bssnrhs_cuda_offset_demalloc.h"
        CHECK_ERROR(cudaFree(dev_dy_hx), "dev_dy_hx cudaFree");
        CHECK_ERROR(cudaFree(dev_dy_hy), "dev_dy_hy cudaFree");
        CHECK_ERROR(cudaFree(dev_dy_hz), "dev_dy_hz cudaFree");
        CHECK_ERROR(cudaFree(dev_sz), "dev_sz cudaFree");
        CHECK_ERROR(cudaFree(dev_zero), "dev_zero cudaFree");
        CHECK_ERROR(cudaFree(dev_pmin), "dev_pmin cudaFree");
        CHECK_ERROR(cudaFree(dev_pmax), "dev_pmax cudaFree");
    }

    for (int bssn_var=0; bssn_var<BSSN_NUM_VARS; bssn_var++){
        CHECK_ERROR(cudaMemcpy(var_out[bssn_var], dev_var_out+bssn_var*unzip_dof, unzip_dof*sizeof(double), cudaMemcpyDeviceToHost), "dev_var_out[bssn_var*unzip_dof] cudaMemcpyDeviceToHost");
    }

    CHECK_ERROR(cudaFree(dev_var_in), "dev_var_in cudaFree");
    CHECK_ERROR(cudaFree(dev_var_out), "dev_var_out cudaFree");

    CHECK_ERROR(cudaDeviceSynchronize(), "cudaDeviceSynchronize in computeBSSN");
    CHECK_ERROR(cudaDeviceReset(), "cudaDeviceReset in computeBSSN");
}

void CPU_2D(const unsigned int blk_lb, const unsigned int blk_up, const unsigned int num_blks, double ** var_in, double ** var_out, Block * blkList){
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

    unsigned int blk_lb=atoi(argv[1]);
    unsigned int blk_up=atoi(argv[2]);
    unsigned int num_blks=atoi(argv[3]);

    Block * blkList = new Block[num_blks*(blk_up-blk_lb+1)];
    double ** var_in = new double*[BSSN_NUM_VARS];
    double ** var_out_GPU = new double*[BSSN_NUM_VARS];
    double ** var_out_CPU = new double*[BSSN_NUM_VARS];

    data_generation_2D(blk_lb, blk_up, num_blks, var_in, var_out_GPU, blkList);
    GPU_2D(blk_lb, blk_up, num_blks, var_in, var_out_GPU, blkList);

    std::cout << std::endl;

    #if test || cpu
    data_generation_2D(blk_lb, blk_up, num_blks, var_in, var_out_CPU, blkList);
    CPU_2D(blk_lb, blk_up, num_blks, var_in, var_out_CPU, blkList);
    #endif

    #if test
    // Verify output
    unsigned int unzip_dof=0;
    for (int blk=0; blk<num_blks*(blk_up-blk_lb+1); blk++){
        unzip_dof += blkList[blk].node1D_x * blkList[blk].node1D_y * blkList[blk].node1D_z;
    }
    unsigned int error_count = 0;
    for(unsigned int i=0;i<BSSN_NUM_VARS;i++){
        for(unsigned int j=0; j<unzip_dof; j++){
            double diff = var_out_CPU[i][j] - var_out_GPU[i][j];
            if (fabs(diff)>threshold){
                error_count++;
                const char separator    = ' ';
                const int nameWidth     = 6;
                const int numWidth      = NUM_DIGITS+10;

                std::cout << std::left << std::setw(nameWidth) << setfill(separator) << "GPU: ";
                std::cout <<std::setprecision(NUM_DIGITS)<< std::left << std::setw(numWidth) << setfill(separator) << var_out_GPU[i][j];

                std::cout << std::left << std::setw(nameWidth) << setfill(separator) << "CPU: ";
                std::cout <<std::setprecision(NUM_DIGITS)<< std::left << std::setw(numWidth) << setfill(separator) << var_out_CPU[i][j];

                std::cout << std::left << std::setw(nameWidth) << setfill(separator) << "DIFF: ";
                std::cout <<std::setprecision(NUM_DIGITS)<< std::left << std::setw(numWidth) << setfill(separator) << diff << std::endl;
            }
        }
    }
    #endif
    return 0;
}
