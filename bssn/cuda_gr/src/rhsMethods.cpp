/**
 * Created on: Sep 21, 2018
 * 		Author: Akila, Eranga, Eminda, Ruwan
 **/

#include "rhsMethods.h"

void GPU_parallelized(unsigned int numberOfBlocks, Block * blkList, unsigned int lower_bound, unsigned int upper_bound, double ** var_in_array, double ** var_out_array){ 

    int steamCountToLevel[5] = {3, 3, 3, 2, 2};

    CHECK_ERROR(cudaSetDevice(0), "cudaSetDevice in computeBSSN"); // Set the GPU that we are going to deal with

    // Creating cuda streams for the process
    cudaStream_t stream;
    cudaStream_t streams[steamCountToLevel[0]]; // usually steamCountToLevel[0] should be the max number of streams verify it before execute.
    for (int index=0; index<steamCountToLevel[0]; index++){
        CHECK_ERROR(cudaStreamCreate(&streams[index]), "cudaStream creation");
    }

    // Check for available GPU memory
    size_t free_bytes, total_bytes;
    CHECK_ERROR(cudaMemGetInfo(&free_bytes, &total_bytes), "Available GPU memory checking failed");
    double GPUCapacity = 1.0*free_bytes/1024/1024 - 10; // in MB
    std::cout << "Available GPU with buffer of " << 10 << ": " << GPUCapacity << " | Total GPU memory: " << total_bytes/1024/1024 << std::endl << std::endl;

    // Sort the data block list
    mergeSort(blkList, 0, numberOfBlocks-1); // O(nlog(n))

    // block level seperation
    int * offsets = new int[upper_bound-lower_bound];
    int counter = 0;
    int init_level = lower_bound;
    offsets[counter] = 0;
    for (int i=0; i<numberOfBlocks; i++){
        if(blkList[i].blkLevel!=init_level){
            init_level = blkList[i].blkLevel;
            counter++;
            offsets[counter] = i;
        }
    }
    counter++;
    offsets[counter] = numberOfBlocks;

    // Let's process level by level
    Block blk;
    double data_usage = 0;
    double fixed_usage = 0;
    int current_index = 0;
    int num_streams = 2;
    for (int level=lower_bound; level<=upper_bound; level++){

        num_streams = steamCountToLevel[level];
        if (steamCountToLevel[level]==2) level=upper_bound;
        
        // Check GPU memory requirement
        blk = Block(0, 0, 0, 2*level, level, 12);
        fixed_usage = num_streams*(138+72)*blk.blkSize*sizeof(double)/1024/1024 + 210*4 + num_streams*4;
        data_usage = num_streams*(blk.blkSize*BSSN_NUM_VARS*sizeof(double)*2)/1024/1024;

        if (data_usage+fixed_usage>(GPUCapacity)){
            std::cout << "Required GPU memory = " << data_usage+fixed_usage << " Failed to allocate enough memory. Program terminated..." << std::endl;
            exit(0);
        }
        std::cout << "Required GPU memory = " << data_usage+fixed_usage << std::endl;

        // Allocating device memory to hold input and output
        double ** dev_var_in_array = new double*[num_streams];
        double ** dev_var_out_array = new double*[num_streams];

        for (int i=0; i<num_streams; i++){
            CHECK_ERROR(cudaMalloc((void**)&dev_var_in_array[i], blk.blkSize*BSSN_NUM_VARS*sizeof(double)), "dev_var_in_array[i]");
            CHECK_ERROR(cudaMalloc((void**)&dev_var_out_array[i], blk.blkSize*BSSN_NUM_VARS*sizeof(double)), "dev_var_out_array[i]");
        }

        int largest_intermediate_array = num_streams*blk.blkSize;
            
        // Allocation intermediate arrays
        int size = largest_intermediate_array * sizeof(double);
        #include "bssnrhs_cuda_variable_malloc.h"
        #include "bssnrhs_cuda_variable_malloc_adv.h"
        #include "bssnrhs_cuda_malloc.h"
        #include "bssnrhs_cuda_malloc_adv.h"

        
        std::cout << current_index << " \t " << offsets[level+1-lower_bound] << std::endl;
        
        // Start block processing
        int streamIndex = 0;
        int unzip_dof=blk.blkSize;
        bool isAsyncStarted = false;

        unsigned int sz[3];
        double ptmin[3], ptmax[3];
        unsigned int bflag = 0;
        double dx, dy, dz;

        bflag=1;

        dx=0.1;
        dy=0.1;
        dz=0.1;

        ptmin[0]=0.0;
        ptmin[1]=0.0;
        ptmin[2]=0.0;

        ptmax[0]=1.0;
        ptmax[1]=1.0;
        ptmax[2]=1.0;
        for (int i=current_index; i<offsets[level+1-lower_bound]; i++ ){
            blk=blkList[i];

            sz[0]=blk.node1D_x; 
            sz[1]=blk.node1D_y;
            sz[2]=blk.node1D_z;

            streamIndex = i%num_streams;

            CHECK_ERROR(cudaMemcpyAsync(dev_var_in_array[i%num_streams], var_in_array[blk.block_no], BSSN_NUM_VARS*blk.blkSize*sizeof(double), cudaMemcpyHostToDevice, streams[i%num_streams]), "dev_var_in_array[index] cudaMemcpyHostToDevice");

            cuda_bssnrhs(dev_var_out_array[i%num_streams], dev_var_in_array[i%num_streams], blk.blkSize, ptmin, ptmax, sz, bflag, streams[i%num_streams],
            #include "list_of_args_per_blk.h"
            );

            CHECK_ERROR(cudaMemcpyAsync(var_out_array[blk.block_no], dev_var_out_array[i%num_streams], BSSN_NUM_VARS*blk.blkSize*sizeof(double), cudaMemcpyDeviceToHost, streams[i%num_streams]), "dev_var_out_array[index] cudaMemcpyDeviceToHost");
        }
        CHECK_ERROR(cudaDeviceSynchronize(), "device sync in computeBSSN");

        current_index = offsets[level+1-lower_bound];

        // Release GPU memory
        #include "bssnrhs_cuda_mdealloc.h"
        #include "bssnrhs_cuda_mdealloc_adv.h"
        for (int i=0; i<num_streams; i++){
            CHECK_ERROR(cudaFree(dev_var_in_array[i]), "dev_var_in_array[0] cudaFree");
            CHECK_ERROR(cudaFree(dev_var_out_array[i]), "dev_var_out_array[0] cudaFree");
        }
    }

    // destroying cuda streams
    for (int i=0; i<steamCountToLevel[0]; i++){
        CHECK_ERROR(cudaStreamDestroy(streams[i]), "cudaStream destruction");
    }    
    return;
}

void GPU_async(unsigned int numberOfBlocks, Block * blkList, unsigned int lower_bound, unsigned int upper_bound, double ** var_in_array, double ** var_out_array){ 
    
    int steamCountToLevel[5] = {3, 3, 3, 2, 2};
    int hybrid_divider = 3;
    
    CHECK_ERROR(cudaSetDevice(0), "cudaSetDevice in computeBSSN"); // Set the GPU that we are going to deal with

    // creating cuda streams for the process
    cudaStream_t streams[3];
    CHECK_ERROR(cudaStreamCreate(&streams[0]), "cudaStream 1 creation");
    CHECK_ERROR(cudaStreamCreate(&streams[1]), "cudaStream 2 creation");
    CHECK_ERROR(cudaStreamCreate(&streams[2]), "cudaStream 3 creation");

    // Check for available GPU memory
    size_t free_bytes, total_bytes;
    CHECK_ERROR(cudaMemGetInfo(&free_bytes, &total_bytes), "Available GPU memory checking failed");
    double GPUCapacity = 1.0*free_bytes/1024/1024 - 10; // In MB
    std::cout << "Available GPU with buffer of " << 10 << "MB : " << GPUCapacity << " | Total GPU memory: " << total_bytes/1024/1024 << std::endl << std::endl;

    // Sort the data block list
    mergeSort(blkList, 0, numberOfBlocks-1); // O(nlog(n))

    // GPU memory requirement
    Block blk = blkList[numberOfBlocks-1];
    double fixed_usage, data_usage = 0;
    fixed_usage = (138+72)*blk.blkSize*sizeof(double)/1024/1024 + 210*4 + 2*4; // 4*210 means allocation overhead
    data_usage = blk.blkSize*BSSN_NUM_VARS*sizeof(double)*2/1024/1024;

    if (data_usage+fixed_usage>(GPUCapacity)){
        std::cout << "Required GPU memory = " << data_usage+fixed_usage << " Failed to allocate enough memory. Program terminated..." << std::endl;
        exit(0);
    }
    std::cout << "Required GPU memory = " << data_usage+fixed_usage << std::endl;

    // Allocating device memory to hold input and output
    double ** dev_var_in_array = new double*[2];
    double ** dev_var_out_array = new double*[2];

    CHECK_ERROR(cudaMalloc((void**)&dev_var_in_array[0], blk.blkSize*BSSN_NUM_VARS*sizeof(double)), "dev_var_in_array[0]");
    CHECK_ERROR(cudaMalloc((void**)&dev_var_in_array[1], blk.blkSize*BSSN_NUM_VARS*sizeof(double)), "dev_var_in_array[1]");

    CHECK_ERROR(cudaMalloc((void**)&dev_var_out_array[0], blk.blkSize*BSSN_NUM_VARS*sizeof(double)), "dev_var_out_array[0]");
    CHECK_ERROR(cudaMalloc((void**)&dev_var_out_array[1], blk.blkSize*BSSN_NUM_VARS*sizeof(double)), "dev_var_out_array[1]");

    int largest_intermediate_array = blk.blkSize;
        
    // Allocation intermediate arrays
    int size = largest_intermediate_array * sizeof(double);
    #include "bssnrhs_cuda_variable_malloc.h"
    #include "bssnrhs_cuda_variable_malloc_adv.h"
    #include "bssnrhs_cuda_malloc.h"
    #include "bssnrhs_cuda_malloc_adv.h"
        
    // Start block processing
    int streamIndex = 0;
    int unzip_dof=blk.blkSize;
    bool isAsyncStarted = false;

    unsigned int sz[3];
    double ptmin[3], ptmax[3];
    unsigned int bflag = 0;
    double dx, dy, dz;

    bflag=1;

    dx=0.1;
    dy=0.1;
    dz=0.1;

    ptmin[0]=0.0;
    ptmin[1]=0.0;
    ptmin[2]=0.0;

    ptmax[0]=1.0;
    ptmax[1]=1.0;
    ptmax[2]=1.0;

    for(int index=0; index<=numberOfBlocks-1; index++) {
        blk=blkList[index];

        sz[0]=blk.node1D_x; 
        sz[1]=blk.node1D_y;
        sz[2]=blk.node1D_z;

        // Starting async process
        CHECK_ERROR(cudaMemcpyAsync(dev_var_in_array[index%2], var_in_array[blk.block_no], BSSN_NUM_VARS*blk.blkSize*sizeof(double), cudaMemcpyHostToDevice, streams[index%3]), "dev_var_in_array[index] cudaMemcpyHostToDevice");

        if (isAsyncStarted){
            cudaStreamSynchronize(streams[(index+2)%3]);
        }

        cuda_bssnrhs(dev_var_out_array[index%2], dev_var_in_array[index%2], blk.blkSize, ptmin, ptmax, sz, bflag, streams[index%3],
            #include "list_of_args_per_blk.h"
        );
     
        if (isAsyncStarted){
            CHECK_ERROR(cudaMemcpyAsync(var_out_array[blkList[index-1].block_no], dev_var_out_array[(index+1)%2], BSSN_NUM_VARS*blkList[index-1].blkSize*sizeof(double), cudaMemcpyDeviceToHost, streams[(index+2)%3]), "dev_var_out_array[index] cudaMemcpyDeviceToHost");
        }

        //handle last DtoH memcopy
        if(index==numberOfBlocks-1){
            CHECK_ERROR(cudaMemcpyAsync(var_out_array[blk.block_no], dev_var_out_array[index%2], BSSN_NUM_VARS*blk.blkSize*sizeof(double), cudaMemcpyDeviceToHost, streams[index%3]), "dev_var_out_array[index] cudaMemcpyDeviceToHost");
        }

        if (!isAsyncStarted) isAsyncStarted=true;
    }

    CHECK_ERROR(cudaDeviceSynchronize(), "device sync in computeBSSN");

    // Release GPU memory
    #include "bssnrhs_cuda_mdealloc.h"
    #include "bssnrhs_cuda_mdealloc_adv.h"
    CHECK_ERROR(cudaFree(dev_var_in_array[0]), "dev_var_in_array[0] cudaFree");
    CHECK_ERROR(cudaFree(dev_var_in_array[1]), "dev_var_in_array[1] cudaFree");

    CHECK_ERROR(cudaFree(dev_var_out_array[0]), "dev_var_out_array[0] cudaFree");
    CHECK_ERROR(cudaFree(dev_var_out_array[1]), "dev_var_out_array[1] cudaFree");

    // destroying cuda streams
    CHECK_ERROR(cudaStreamDestroy(streams[0]), "cudaStream 1 destruction");
    CHECK_ERROR(cudaStreamDestroy(streams[1]), "cudaStream 2 destruction");
    CHECK_ERROR(cudaStreamDestroy(streams[2]), "cudaStream 3 destruction");

    return;
}

void GPU_hybrid(unsigned int numberOfBlocks, Block * blkList, unsigned int lower_bound, unsigned int upper_bound, double ** var_in_array, double ** var_out_array){ 

    int steamCountToLevel[5] = {3, 3, 3, 2, 2};
    int hybrid_divider = 3;

    CHECK_ERROR(cudaSetDevice(0), "cudaSetDevice in computeBSSN"); // Set the GPU that we are going to deal with

    // Creating cuda streams for the process
    cudaStream_t stream;
    cudaStream_t streams[steamCountToLevel[0]]; // usually steamCountToLevel[0] should be the max number of streams verify it before execute.
    for (int index=0; index<steamCountToLevel[0]; index++){
        CHECK_ERROR(cudaStreamCreate(&streams[index]), "cudaStream creation");
    }

    // Check for available GPU memory
    size_t free_bytes, total_bytes;
    CHECK_ERROR(cudaMemGetInfo(&free_bytes, &total_bytes), "Available GPU memory checking failed");
    double GPUCapacity = 1.0*free_bytes/1024/1024 - 10; // in MB
    std::cout << "Available GPU with buffer of " << 10 << ": " << GPUCapacity << " | Total GPU memory: " << total_bytes/1024/1024 << std::endl << std::endl;

    // Sort the data block list
    mergeSort(blkList, 0, numberOfBlocks-1); // O(nlog(n))

    // block level seperation
    int * offsets = new int[upper_bound-lower_bound];
    int counter = 0;
    int init_level = lower_bound;
    offsets[counter] = 0;
    for (int i=0; i<numberOfBlocks; i++){
        if(blkList[i].blkLevel!=init_level){
            init_level = blkList[i].blkLevel;
            counter++;
            offsets[counter] = i;
        }
    }
    counter++;
    offsets[counter] = numberOfBlocks;

    // Let's process level by level
    Block blk;
    double data_usage = 0;
    double fixed_usage = 0;
    int current_index = 0;
    int num_streams = 2;
    for (int level=lower_bound; level<=upper_bound; level++){

        num_streams = steamCountToLevel[level];
        if (steamCountToLevel[level]==2) level=upper_bound;
        
        // Check GPU memory requirement
        blk = Block(0, 0, 0, 2*level, level, 12);
        fixed_usage = num_streams*(138+72)*blk.blkSize*sizeof(double)/1024/1024 + 210*4 + num_streams*4;
        data_usage = num_streams*(blk.blkSize*BSSN_NUM_VARS*sizeof(double)*2)/1024/1024;
        if (num_streams==2){
            fixed_usage = (138+72)*blk.blkSize*sizeof(double)/1024/1024 + 210*4 + 2*4; // 4*210 means allocation overhead
            data_usage = blk.blkSize*BSSN_NUM_VARS*sizeof(double)*2/1024/1024;
        }

        if (data_usage+fixed_usage>(GPUCapacity)){
            std::cout << "Required GPU memory = " << data_usage+fixed_usage << " Failed to allocate enough memory. Program terminated..." << std::endl;
            exit(0);
        }
        std::cout << "Required GPU memory = " << data_usage+fixed_usage << std::endl;

        // Allocating device memory to hold input and output
        double ** dev_var_in_array = new double*[num_streams];
        double ** dev_var_out_array = new double*[num_streams];

        for (int i=0; i<num_streams; i++){
            CHECK_ERROR(cudaMalloc((void**)&dev_var_in_array[i], blk.blkSize*BSSN_NUM_VARS*sizeof(double)), "dev_var_in_array[i]");
            CHECK_ERROR(cudaMalloc((void**)&dev_var_out_array[i], blk.blkSize*BSSN_NUM_VARS*sizeof(double)), "dev_var_out_array[i]");
        }

        int largest_intermediate_array = num_streams*blk.blkSize;
        if (num_streams==2){
            largest_intermediate_array = blk.blkSize;
        }
            
        // Allocation intermediate arrays
        int size = largest_intermediate_array * sizeof(double);
        #include "bssnrhs_cuda_variable_malloc.h"
        #include "bssnrhs_cuda_variable_malloc_adv.h"
        #include "bssnrhs_cuda_malloc.h"
        #include "bssnrhs_cuda_malloc_adv.h"

        std::cout << current_index << " \t " << offsets[level+1-lower_bound] << std::endl;
        
        // Start block processing
        int streamIndex = 0;
        int unzip_dof=blk.blkSize;
        bool isAsyncStarted = false;

        unsigned int sz[3];
        double ptmin[3], ptmax[3];
        unsigned int bflag = 0;
        double dx, dy, dz;

        bflag=1;

        dx=0.1;
        dy=0.1;
        dz=0.1;

        ptmin[0]=0.0;
        ptmin[1]=0.0;
        ptmin[2]=0.0;

        ptmax[0]=1.0;
        ptmax[1]=1.0;
        ptmax[2]=1.0;
        for (int i=current_index; i<offsets[level+1-lower_bound]; i++ ){
            blk=blkList[i];

            sz[0]=blk.node1D_x; 
            sz[1]=blk.node1D_y;
            sz[2]=blk.node1D_z;

            if (num_streams==2){
                // Starting async process
                CHECK_ERROR(cudaMemcpyAsync(dev_var_in_array[i%2], var_in_array[blk.block_no], BSSN_NUM_VARS*blk.blkSize*sizeof(double), cudaMemcpyHostToDevice, streams[i%3]), "dev_var_in_array[index] cudaMemcpyHostToDevice");

                if (isAsyncStarted){
                    cudaStreamSynchronize(streams[(i+2)%3]);
                }

                cuda_bssnrhs(dev_var_out_array[i%2], dev_var_in_array[i%2], blk.blkSize, ptmin, ptmax, sz, bflag, streams[i%3],
                    #include "list_of_args_per_blk.h"
                );
            
                if (isAsyncStarted){
                    CHECK_ERROR(cudaMemcpyAsync(var_out_array[blkList[i-1].block_no], dev_var_out_array[(i+1)%2], BSSN_NUM_VARS*blkList[i-1].blkSize*sizeof(double), cudaMemcpyDeviceToHost, streams[(i+2)%3]), "dev_var_out_array[index] cudaMemcpyDeviceToHost");
                }

                //handle last DtoH memcopy
                if(i==offsets[level+1-lower_bound]-1){
                    CHECK_ERROR(cudaMemcpyAsync(var_out_array[blk.block_no], dev_var_out_array[i%2], BSSN_NUM_VARS*blk.blkSize*sizeof(double), cudaMemcpyDeviceToHost, streams[i%3]), "dev_var_out_array[index] cudaMemcpyDeviceToHost");
                }

                if (!isAsyncStarted) isAsyncStarted=true;
            }else{
                streamIndex = i%num_streams;

                CHECK_ERROR(cudaMemcpyAsync(dev_var_in_array[i%num_streams], var_in_array[blk.block_no], BSSN_NUM_VARS*blk.blkSize*sizeof(double), cudaMemcpyHostToDevice, streams[i%num_streams]), "dev_var_in_array[index] cudaMemcpyHostToDevice");

                cuda_bssnrhs(dev_var_out_array[i%num_streams], dev_var_in_array[i%num_streams], blk.blkSize, ptmin, ptmax, sz, bflag, streams[i%num_streams],
                #include "list_of_args_per_blk.h"
                );

                CHECK_ERROR(cudaMemcpyAsync(var_out_array[blk.block_no], dev_var_out_array[i%num_streams], BSSN_NUM_VARS*blk.blkSize*sizeof(double), cudaMemcpyDeviceToHost, streams[i%num_streams]), "dev_var_out_array[index] cudaMemcpyDeviceToHost");
            }
        }
        CHECK_ERROR(cudaDeviceSynchronize(), "device sync in computeBSSN");

        current_index = offsets[level+1-lower_bound];

        // Release GPU memory
        #include "bssnrhs_cuda_mdealloc.h"
        #include "bssnrhs_cuda_mdealloc_adv.h"
        for (int i=0; i<num_streams; i++){
            CHECK_ERROR(cudaFree(dev_var_in_array[i]), "dev_var_in_array[0] cudaFree");
            CHECK_ERROR(cudaFree(dev_var_out_array[i]), "dev_var_out_array[0] cudaFree");
        }
    }

    // destroying cuda streams
    for (int i=0; i<steamCountToLevel[0]; i++){
        CHECK_ERROR(cudaStreamDestroy(streams[i]), "cudaStream destruction");
    }    
    return;
}

#include "rhs.h"
void CPU_sequence(unsigned int numberOfLevels, Block * blkList, double ** var_in, double ** var_out){

    double ptmin[3], ptmax[3];
    unsigned int sz[3];
    unsigned int bflag;
    unsigned int offset;
    double dx, dy, dz;

    for(unsigned int blk=0; blk<numberOfLevels; blk++){

        offset=blkList[blk].offset;
        sz[0]=blkList[blk].node1D_x; 
        sz[1]=blkList[blk].node1D_y;
        sz[2]=blkList[blk].node1D_z;

        bflag=1;

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
}
