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
    double ** dev_var_out_array = new double*[num_blks*(blk_up-blk_lb+1)];;

    // Generate sample data
    for (int index=0; index<num_blks*(blk_up-blk_lb+1); index++){
        int block_no = index%num_blks;
        int level = ((index/(num_blks))%(blk_up-blk_lb+1))+blk_lb;

        Block * blkList=new Block[num_blks*(blk_up-blk_lb+1)];
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










    #pragma omp parallel for
    for (int index=0; index<num_blks*(blk_up-blk_lb+1); index++){
            int block_no = index%num_blks;
            int level = ((index/(num_blks))%(blk_up-blk_lb+1))+blk_lb;

            Block * blkList=new Block[1];
            const unsigned int maxDepth=12;
            int lev = blk_lb+level;
            blkList[0]=Block(0, 0, 0, 2*lev, lev, maxDepth);
            Block & blk=blkList[0];
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
            Block tmpBlock;

            cudaError_t cudaStatus;
            cudaStatus = cudaSetDevice(0);
            // if (cudaStatus != cudaSuccess) {
            //     fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?\n");
            //     return 0;
            // }

            cudaStatus = cudaMalloc((void**)&dev_var_in_array[index], unzip_dof*BSSN_NUM_VARS*sizeof(double));
            // if (cudaStatus != cudaSuccess) {fprintf(stderr, "dev_var_in_array cudaMalloc failed!\n"); return 0;}

            cudaStatus = cudaMalloc((void**)&dev_var_out_array[index], unzip_dof*BSSN_NUM_VARS*sizeof(double));
            // if (cudaStatus != cudaSuccess) {fprintf(stderr, "var_out cudaMalloc failed!\n"); return 0;}

            tmpBlock=blkList[0];
            x=(double)tmpBlock.x;
            y=(double)tmpBlock.y;
            z=(double)tmpBlock.z;

            hx=0.001; 
            hy=0.001; 
            hz=0.001; 

            offset=tmpBlock.offset;
            size_x=tmpBlock.node1D_x;
            size_y=tmpBlock.node1D_y;
            size_z=tmpBlock.node1D_z;

            for(unsigned int k=0;k<tmpBlock.node1D_z;k++)
            {
                for(unsigned int j=0;j<tmpBlock.node1D_y;j++)
                {
                    for(unsigned int i=0;i<tmpBlock.node1D_x;i++)
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

            cudaStream_t stream;
            cudaError_t result;
            result = cudaStreamCreate(&stream);
            // if (result != cudaSuccess) {fprintf(stderr, "cudaStream creation failed!\n"); return 0;}

            cudaStatus = cudaMemcpyAsync(dev_var_in_array[index], var_in_per_block, BSSN_NUM_VARS*unzip_dof*sizeof(double), cudaMemcpyHostToDevice, stream);
            // if (cudaStatus != cudaSuccess) {fprintf(stderr, "var_in_per_block asyncCudaMemcpy failed!\n"); return 0;}

            double ptmin[3], ptmax[3];
            unsigned int sz[3];
            unsigned int bflag;
            double dx,dy,dz;

            offset=blkList[0].offset;
            sz[0]=blkList[0].node1D_x; 
            sz[1]=blkList[0].node1D_y;
            sz[2]=blkList[0].node1D_z;
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

            std::cout << "Block no: " << index << " Total Points: " << total_points << std::endl;
            
            printf("index=%d | level=%d | block=%d\n", index, level, block_no);
            cuda_bssnrhs(dev_var_out_array[index], dev_var_in_array[index], unzip_dof , offset, ptmin, ptmax, sz, bflag, stream);

            cudaStatus = cudaMemcpyAsync(var_out_per_block, dev_var_out_array[index], BSSN_NUM_VARS*unzip_dof*sizeof(double), cudaMemcpyDeviceToHost, stream);
    
            cudaFree(dev_var_in_array[index]);
            cudaFree(dev_var_out_array[index]);

    }
    cudaError_t cudaStatus;
    cudaStatus = cudaDeviceSynchronize();


    // delete [] host_var_in;
    // delete [] host_var_out;
    
    // cudaFree(dev_var_in);
    // cudaFree(dev_var_out);

    printf("Executed\n");
    return 0;
}