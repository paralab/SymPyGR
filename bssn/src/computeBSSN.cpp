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
#include "rhs.h"

double ** GPU_sequence(unsigned int blk_lb, unsigned int blk_up, unsigned int num_blks)
{
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
    #pragma omp parallel for
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
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?\n");
        }


        #include "bssnrhs_cuda_offset_variable_malloc.h"
        #include "bssnrhs_cuda_variable_malloc.h"
        #include "bssnrhs_cuda_variable_malloc_adv.h"
        double * dev_dy_hx; //similar to hx in cpu code
        double * dev_dy_hy;
        double * dev_dy_hz;
        int * dev_sz;
        int * dev_zero;
        double * dev_pmin;
        double *dev_pmax;
        int * dev_bflag;

        int size = total_points * sizeof(double);

        // #pragma omp critical (global_data_lock)
        // {
        // GPU Memory allocations
        cudaStatus = cudaMalloc((void**)&dev_var_in_array[index], unzip_dof*BSSN_NUM_VARS*sizeof(double));
        if (cudaStatus != cudaSuccess) {fprintf(stderr, "dev_var_in_array cudaMalloc failed!\n");}

        cudaStatus = cudaMalloc((void**)&dev_var_out_array[index], unzip_dof*BSSN_NUM_VARS*sizeof(double));
        if (cudaStatus != cudaSuccess) {fprintf(stderr, "dev_var_out_array cudaMalloc failed!\n");}

        #include "bssnrhs_cuda_offset_malloc.h"

        cudaStatus = cudaMalloc((void **) &dev_dy_hx, sizeof(double));
        if (cudaStatus != cudaSuccess) {fprintf(stderr, "hx cudaMalloc failed!\n");}

        cudaStatus = cudaMalloc((void **) &dev_dy_hy, sizeof(double));
        if (cudaStatus != cudaSuccess) {fprintf(stderr, "hy cudaMalloc failed!\n");}

        cudaStatus = cudaMalloc((void **) &dev_dy_hz, sizeof(double));
        if (cudaStatus != cudaSuccess) {fprintf(stderr, "hz cudaMalloc failed!\n");}

        cudaStatus = cudaMalloc((void **) &dev_sz, 3*sizeof(int));
        if (cudaStatus != cudaSuccess) {fprintf(stderr, "sz cudaMalloc failed!\n");}

        cudaStatus = cudaMalloc((void **) &dev_zero, sizeof(int));
        if (cudaStatus != cudaSuccess) {fprintf(stderr, "0 cudaMalloc failed!\n");}

        cudaStatus = cudaMalloc((void **) &dev_pmin, 3*sizeof(double));
        if (cudaStatus != cudaSuccess) {fprintf(stderr, "pmin cudaMalloc failed!\n");}

        cudaStatus = cudaMalloc((void **) &dev_pmax, 3*sizeof(double));
        if (cudaStatus != cudaSuccess) {fprintf(stderr, "pmax cudaMalloc failed!\n");}

        cudaStatus = cudaMalloc((void **) &dev_bflag, sizeof(int));
        if (cudaStatus != cudaSuccess)fprintf(stderr, "bflag cudaMalloc failed!\n");

        #include "bssnrhs_cuda_malloc.h"
        #include "bssnrhs_cuda_malloc_adv.h"
        // }

        cudaStream_t stream;
        cudaStatus = cudaStreamCreate(&stream);
        if (cudaStatus != cudaSuccess) {fprintf(stderr, "cudaStream creation failed!\n");}

        cudaStatus = cudaMemcpyAsync(dev_var_in_array[index], var_in_array[index], BSSN_NUM_VARS*unzip_dof*sizeof(double), cudaMemcpyHostToDevice, stream);
        if (cudaStatus != cudaSuccess) {fprintf(stderr, "var_in_array[index] cudaMemcpyAsync failed!\n");}

        cuda_bssnrhs(dev_var_out_array[index], dev_var_in_array[index], unzip_dof , offset, ptmin, ptmax, sz, bflag, stream,
        #include "list_of_args.h"
        , dev_dy_hx, dev_dy_hy, dev_dy_hz, dev_sz, dev_zero, dev_pmin, dev_pmax, dev_bflag
        );

        cudaStatus = cudaMemcpyAsync(var_out_array[index], dev_var_out_array[index], BSSN_NUM_VARS*unzip_dof*sizeof(double), cudaMemcpyDeviceToHost, stream);
        if (cudaStatus != cudaSuccess) {fprintf(stderr, "var_out_array[index] cudaMemcpyAsync failed!\n");}

        cudaStatus = cudaStreamSynchronize(stream);
        if (cudaStatus != cudaSuccess) {fprintf(stderr, "3 stream cudaStreamSynchronize failed!\n");}

        cudaStatus = cudaStreamDestroy(stream);
        if (cudaStatus != cudaSuccess) {fprintf(stderr, "stream cudaStreamDestroy failed!\n");}

        // // #pragma omp critical (global_data_lock)
        // // {

        //     #include "bssnrhs_cuda_mdealloc.h"
        //     #include "bssnrhs_cuda_mdealloc_adv.h"

        //     #include "bssnrhs_cuda_offset_demalloc.h"

        //     cudaStatus = cudaFree(dev_dy_hx);
        //     if (cudaStatus != cudaSuccess) {fprintf(stderr, "dev_dy_hx cudaFree failed!\n");}

        //     cudaStatus = cudaFree(dev_dy_hy);
        //     if (cudaStatus != cudaSuccess) {fprintf(stderr, "dev_dy_hy cudaFree failed!\n");}

        //     cudaStatus = cudaFree(dev_dy_hz);
        //     if (cudaStatus != cudaSuccess) {fprintf(stderr, "dev_dy_hz cudaFree failed!\n");}

        //     cudaStatus = cudaFree(dev_sz);
        //     if (cudaStatus != cudaSuccess) {fprintf(stderr, "dev_sz cudaFree failed!\n");}

        //     cudaStatus = cudaFree(dev_zero);
        //     if (cudaStatus != cudaSuccess) {fprintf(stderr, "dev_zero cudaFree failed!\n");}

        //     cudaStatus = cudaFree(dev_pmin);
        //     if (cudaStatus != cudaSuccess) {fprintf(stderr, "dev_pmin cudaFree failed!\n");}

        //     cudaStatus = cudaFree(dev_pmax);
        //     if (cudaStatus != cudaSuccess) {fprintf(stderr, "dev_pmax cudaFree failed!\n");}

        //     cudaStatus = cudaFree(dev_var_in_array[index]);
        //     if (cudaStatus != cudaSuccess) {fprintf(stderr, "dev_var_in_array[index] cudaFree failed!\n");}

        //     cudaStatus = cudaFree(dev_var_out_array[index]);
        //     if (cudaStatus != cudaSuccess) {fprintf(stderr, "dev_var_out_array[index] cudaFree failed!\n");}
        // // }
        delete [] var_in_array[index];
    }

    cudaError_t cudaStatus;
    cudaStatus = cudaDeviceSynchronize();

    // delete [] blkList;
    printf("GPU done\n\n");
    return var_out_array;
}

// Customized CPU code for async GPU result testing
int CPU_sequence(unsigned int blk_lb, unsigned int blk_up, unsigned int num_blks, double ** var_out, unsigned int block_index)
{
    const unsigned int total_blks= num_blks*(blk_up-blk_lb+1);
    Block * blkList=new Block[total_blks];
    unsigned long unzipSz=0;
    unsigned int index=0;
    const unsigned int maxDepth=12;

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

    double ** var_in=new double*[BSSN_NUM_VARS];

    for(unsigned int i=0;i<BSSN_NUM_VARS;i++)
    {
        var_in[i] = new double[unzip_dof];
        var_out[i] = new double[unzip_dof];
    }

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
                        var_in[var][offset+k*size_y*size_x+j*size_y+i]=u[var];
                        var_out[var][offset+k*size_y*size_x+j*size_y+i]=0;
                    }
                }
    }

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

        bssnrhs(var_out, (const double **)var_in, offset, ptmin, ptmax, sz, bflag);

    }
    for(unsigned int i=0;i<BSSN_NUM_VARS;i++)
    {
        delete [] var_in[i];
    }
    delete [] var_in;
    delete [] blkList;

    printf("CPU done block no=%d\n", block_index);
    return unzip_dof;
}


int main (int argc, char** argv)
{
    
    unsigned int blk_lb=0;
    unsigned int blk_up=3;
    unsigned int num_blks=2;

    blk_lb=atoi(argv[1]);
    blk_up=atoi(argv[2]);
    num_blks=atoi(argv[3]);

    // GPU Call
    #include "rhs_cuda.h"
    double ** var_out_array = GPU_sequence(blk_lb, blk_up, num_blks);

    //CPU Call for test GPU result
    #include "rhs.h"
    
    for (int index=0; index<num_blks*(blk_up-blk_lb+1); index++){
        double ** var_out=new double*[BSSN_NUM_VARS];

        int block_no = index%num_blks;
        int level = ((index/(num_blks))%(blk_up-blk_lb+1))+blk_lb;

        int unzip_dof = CPU_sequence(level, level, 1, var_out, index);

        unsigned int error_count = 0;
        for(unsigned int i=0;i<BSSN_NUM_VARS;i++){
            for(unsigned int j=0; j<unzip_dof; j++){
                unsigned int abs_index = i*unzip_dof + j;
                double diff = var_out[i][j] - var_out_array[index][abs_index];
                if (fabs(diff)>threshold){
                    error_count++;
                    const char separator    = ' ';
                    const int nameWidth     = 6;
                    const int numWidth      = NUM_DIGITS+10;

                    std::cout << std::left << std::setw(nameWidth) << setfill(separator) << "GPU: ";
                    std::cout <<std::setprecision(NUM_DIGITS)<< std::left << std::setw(numWidth) << setfill(separator)  << var_out_array[index][abs_index];

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
        for(unsigned int i=0;i<BSSN_NUM_VARS;i++)
        {
            delete [] var_out[i];
        }
        delete [] var_out; 
    }
    return 0;
}
