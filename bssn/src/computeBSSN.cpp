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
#include "rhs_cuda.h"

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

    // Create pointers to hold arrays
    double ** var_in_array = new double*[num_blks*(blk_up-blk_lb+1)];
    double ** var_out_array = new double*[num_blks*(blk_up-blk_lb+1)];
    double ** dev_var_in_array = new double*[num_blks*(blk_up-blk_lb+1)];
    double ** dev_var_out_array = new double*[num_blks*(blk_up-blk_lb+1)];
    Block * blkList = new Block[num_blks*(blk_up-blk_lb+1)];

    // Generate sample data
    #pragma omp parallel for
    for (int index=0; index<num_blks*(blk_up-blk_lb+1); index++){
        int block_no = index%num_blks;
        int level = ((index/(num_blks))%(blk_up-blk_lb+1))+blk_lb;

        const unsigned int maxDepth=12;
        int lev = blk_lb+level;

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

        for(unsigned int k=0;k<blk.node1D_z;k++)
        {
            for(unsigned int j=0;j<blk.node1D_y;j++)
            {
                for(unsigned int i=0;i<blk.node1D_x;i++)
                {
                    coord[0]=x+i*hx;
                    coord[1]=y+j*hy;
                    coord[2]=z+k*hz;

                    initial_data(u,coord);

                    for(unsigned int var=0; var<BSSN_NUM_VARS; var++)
                    {
                        var_in_per_block[var*unzip_dof+offset+k*size_y*size_x+j*size_y+i]=u[var];
                        var_out_per_block[var*unzip_dof+offset+k*size_y*size_x+j*size_y+i]=0;
                    }
                }
            }
        }
    }

    // Start processing
    #pragma omp parallel for num_threads(2)
    for (int index=0; index<num_blks*(blk_up-blk_lb+1); index++){
        int block_no = index%num_blks;
        int level = ((index/(num_blks))%(blk_up-blk_lb+1))+blk_lb;

        Block & blk=blkList[index];
        const unsigned long unzip_dof=(blk.node1D_x*blk.node1D_y*blk.node1D_z);
        
        unsigned int offset;
        double ptmin[3], ptmax[3];
        unsigned int sz[3];
        unsigned int bflag;
        double dx,dy,dz;

        offset=blk.offset;
        sz[0]=blk.node1D_x; 
        sz[1]=blk.node1D_y;
        sz[2]=blk.node1D_z;
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

        std::cout << "Block no: " << index << " Total Points: " << total_points << " level: " << level << " index: " << index << std::endl;
        
        // Check for GPU
        cudaError_t cudaStatus;
        cudaStatus = cudaSetDevice(0);
        // if (cudaStatus != cudaSuccess) {
        //     fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?\n");
        //     return 0;
        // }

        // GPU Memory allocations
        cudaStatus = cudaMalloc((void**)&dev_var_in_array[index], unzip_dof*BSSN_NUM_VARS*sizeof(double));
        // if (cudaStatus != cudaSuccess) {fprintf(stderr, "dev_var_in_array cudaMalloc failed!\n"); return 0;}

        cudaStatus = cudaMalloc((void**)&dev_var_out_array[index], unzip_dof*BSSN_NUM_VARS*sizeof(double));
        // if (cudaStatus != cudaSuccess) {fprintf(stderr, "var_out cudaMalloc failed!\n"); return 0;}

        #include "bssnrhs_cuda_offset_malloc.h"

        double * dev_dy_hx; //similar to hx in cpu code
        cudaStatus = cudaMalloc((void **) &dev_dy_hx, sizeof(double));
        // if (cudaStatus != cudaSuccess) {fprintf(stderr, "hx cudaMalloc failed!\n"); return;}
        double * dev_dy_hy;
        cudaStatus = cudaMalloc((void **) &dev_dy_hy, sizeof(double));
        // if (cudaStatus != cudaSuccess) {fprintf(stderr, "hy cudaMalloc failed!\n"); return;}
        double * dev_dy_hz;
        cudaStatus = cudaMalloc((void **) &dev_dy_hz, sizeof(double));
        // if (cudaStatus != cudaSuccess) {fprintf(stderr, "hz cudaMalloc failed!\n"); return;}
        int * dev_sz;
        cudaStatus = cudaMalloc((void **) &dev_sz, 3*sizeof(int));
        // if (cudaStatus != cudaSuccess) {fprintf(stderr, "sz cudaMalloc failed!\n"); return;}
        int * dev_zero;
        cudaStatus = cudaMalloc((void **) &dev_zero, sizeof(int));
        // if (cudaStatus != cudaSuccess) {fprintf(stderr, "0 cudaMalloc failed!\n"); return;}
        double * dev_pmin;
        cudaStatus = cudaMalloc((void **) &dev_pmin, 3*sizeof(double));
        // if (cudaStatus != cudaSuccess) {fprintf(stderr, "pmin cudaMalloc failed!\n"); return;}
        double *dev_pmax;
        cudaStatus = cudaMalloc((void **) &dev_pmax, 3*sizeof(double));
        // if (cudaStatus != cudaSuccess) {fprintf(stderr, "pmax cudaMalloc failed!\n"); return;}
        int * dev_bflag;
        cudaStatus = cudaMalloc((void **) &dev_bflag, sizeof(int));
        // if (cudaStatus != cudaSuccess)fprintf(stderr, "bflag cudaMalloc failed!\n");

        int size = total_points * sizeof(double);

        #include "bssnrhs_cuda_malloc.h"
        #include "bssnrhs_cuda_malloc_adv.h"

        cudaStream_t stream;
        cudaError_t result;
        result = cudaStreamCreate(&stream);
        // if (result != cudaSuccess) {fprintf(stderr, "cudaStream creation failed!\n"); return 0;}

        cudaStatus = cudaMemcpyAsync(dev_var_in_array[index], var_in_array[index], BSSN_NUM_VARS*unzip_dof*sizeof(double), cudaMemcpyHostToDevice, stream);
        // if (cudaStatus != cudaSuccess) {fprintf(stderr, "var_in_per_block asyncCudaMemcpy failed!\n"); return 0;}

    
        cuda_bssnrhs(dev_var_out_array[index], dev_var_in_array[index], unzip_dof , offset, ptmin, ptmax, sz, bflag, stream,
        #include "list_of_args.h"
        , dev_dy_hx, dev_dy_hy, dev_dy_hz, dev_sz, dev_zero, dev_pmin, dev_pmax, dev_bflag
        );

        cudaStatus = cudaMemcpyAsync(var_out_array[index], dev_var_out_array[index], BSSN_NUM_VARS*unzip_dof*sizeof(double), cudaMemcpyDeviceToHost, stream);

        cudaStreamSynchronize(stream);

        cudaStreamDestroy(stream);

        // #include "bssnrhs_cuda_mdealloc.h"
        // #include "bssnrhs_cuda_mdealloc_adv.h"

        // #include "bssnrhs_cuda_offset_demalloc.h"
        // cudaFree(dev_dy_hx);
        // cudaFree(dev_dy_hy);
        // cudaFree(dev_dy_hz);
        // cudaFree(dev_sz);
        // cudaFree(dev_zero);
        // cudaFree(dev_pmin);
        // cudaFree(dev_pmax);

        // delete [] var_in_array[index];
        // cudaFree(dev_var_in_array[index]);
        // cudaFree(dev_var_out_array[index]);
    }

    cudaError_t cudaStatus;
    cudaStatus = cudaDeviceSynchronize();

    // delete [] blkList;

    printf("Done!!!\n");
    return 0;
}