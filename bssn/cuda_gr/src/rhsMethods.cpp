#include "rhsMethods.h"

void GPU_parallelized(unsigned int numberOfLevels, Block * blkList, unsigned int lower_bound, unsigned int upper_bound, double ** var_in_array, double ** var_out_array){ 

    const int steamCountToLevel[5] = {3, 3, 3, 2, 2};
    const int hybrid_divider = 3;

    CHECK_ERROR(cudaSetDevice(0), "cudaSetDevice for the process"); // Set the GPU that we are going to deal with

    // Creating cuda streams for the process
    cudaStream_t stream;
    cudaStream_t streams[steamCountToLevel[0]]; // usually steamCountToLevel[0] should be the max number of streams verify it before execute.
    for (int index=0; index<steamCountToLevel[0]; index++){
        CHECK_ERROR(cudaStreamCreate(&streams[index]), "cudaStream creation");
    }

    // Check for available GPU memory
    size_t free_bytes, total_bytes;
    CHECK_ERROR(cudaMemGetInfo(&free_bytes, &total_bytes), "Available GPU memory checking failed");
    double GPU_capacity_buffer = 10;
    double GPUCapacity = 1.0*free_bytes/1024/1024 - GPU_capacity_buffer;
    std::cout << "Available GPU with buffer of " << GPU_capacity_buffer << ": " << GPUCapacity << " | Total GPU memory: " << total_bytes/1024/1024 << std::endl << std::endl;

    // Sort the data block list
    mergeSort(blkList, 0, numberOfLevels-1); // O(nlog(n))

    // Calculate how much data block can be executed at once in GPU
    Block blk;
    int init_block = 0;
    int current_block = 0;
    int total_points = 0;
    int numberOfStreams = 0;
    double current_usage = 0;
    double fixed_usage = 0;
    double prev_usage = 0;
    double actual_usage = 0;
    double malloc_overhead = 3;

    while (current_block<numberOfLevels){
        current_usage=0;
        fixed_usage=0;
        init_block=current_block;

        while ((current_usage+fixed_usage<(GPUCapacity)) && (current_block<numberOfLevels)){
            prev_usage = current_usage+fixed_usage; // usage of previous iteration
            blk = blkList[current_block];
            total_points = blk.blkSize;
            if (blk.blkLevel<5) {
                numberOfStreams = steamCountToLevel[blk.blkLevel];
            }else{
                numberOfStreams = 2;
            }
            if (fixed_usage<numberOfStreams*(138+72)*total_points*sizeof(double)/1024/1024 + 210*malloc_overhead + (current_block-init_block)*malloc_overhead){
                fixed_usage = numberOfStreams*(138+72)*total_points*sizeof(double)/1024/1024 + 210*malloc_overhead + (current_block-init_block)*malloc_overhead; // 5*210 means allocation overhead
            }
            current_usage += (total_points*BSSN_NUM_VARS*sizeof(double)*2)/1024/1024;
            current_block++;
        }
        actual_usage = current_usage+fixed_usage;
        if (current_usage+fixed_usage>(GPUCapacity)){
            current_block--;
            if (init_block>current_block-1) {
                std::cout << "Required GPU memory = " << actual_usage << " Failed to allocate enough memory. Program terminated..." << std::endl;
                exit(0);
            }
            actual_usage = prev_usage;
        }

        // Display the set of blocks selected to process with their GPU usage
        std::cout << "start: " << init_block << " end: " << current_block-1 << "| usage: " << actual_usage << std::endl;

        // Allocating device memory to hold input and output
        double ** dev_var_in_array = new double*[numberOfLevels];
        double ** dev_var_out_array = new double*[numberOfLevels];

        int largest_intermediate_array = 0;
        for (int index=init_block; index<=current_block-1; index++){
            blk = blkList[index];

            if (blk.blkLevel<5) {
                numberOfStreams = steamCountToLevel[blk.blkLevel];
            }else{
                numberOfStreams = 2;
            }

            if (largest_intermediate_array<blk.blkSize*numberOfStreams){
                largest_intermediate_array = blk.blkSize*numberOfStreams;
            }
            CHECK_ERROR(cudaMalloc((void**)&dev_var_in_array[index], blk.blkSize*BSSN_NUM_VARS*sizeof(double)), "dev_var_in_array[index]");
            CHECK_ERROR(cudaMalloc((void**)&dev_var_out_array[index], blk.blkSize*BSSN_NUM_VARS*sizeof(double)), "dev_var_out_array[index]");
        }
        
        // Allocation intermediate arrays
        int size = largest_intermediate_array * sizeof(double);
        #include "bssnrhs_cuda_variable_malloc.h"
        #include "bssnrhs_cuda_variable_malloc_adv.h"
        #include "bssnrhs_cuda_malloc.h"
        #include "bssnrhs_cuda_malloc_adv.h"

        // Start block processing
        int streamIndex;
        int unzip_dof;
        unsigned int sz[3];
        double ptmin[3], ptmax[3];
        unsigned int bflag;
        double dx, dy, dz;
        for(int index=init_block; index<=current_block-1; index++) {
            blk=blkList[index];

            // Call cudasync in the case of block level change occured
            if (index!=init_block && blkList[index-1].blkLevel!=blk.blkLevel) CHECK_ERROR(cudaDeviceSynchronize(), "device sync in computeBSSN level change");

            // Identify stream to schedule
            if (blk.blkLevel<5) {
                numberOfStreams = steamCountToLevel[blk.blkLevel];
            }else{
                numberOfStreams = 2;
            }
            streamIndex = index % numberOfStreams;
            stream = streams[streamIndex];

            // Configure the block specific values
            unzip_dof=blk.blkSize;

            sz[0]=blk.node1D_x; 
            sz[1]=blk.node1D_y;
            sz[2]=blk.node1D_z;

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
            
            // std::cout << "GPU - Count: " << std::setw(3) << index << " - Block no: " << std::setw(3) << blk.block_no << " - Bock level: " << std::setw(1) << blk.blkLevel << " - Block size: " << blk.blkSize << std::endl;
            
            CHECK_ERROR(cudaMemcpyAsync(dev_var_in_array[index], var_in_array[blk.block_no], BSSN_NUM_VARS*unzip_dof*sizeof(double), cudaMemcpyHostToDevice, stream), "dev_var_in_array[index] cudaMemcpyHostToDevice");
            // CHECK_ERROR(cudaDeviceSynchronize(), "device sync in computeBSSN");
            
            
            cuda_bssnrhs(dev_var_out_array[index], dev_var_in_array[index], unzip_dof, ptmin, ptmax, sz, bflag, stream,
            #include "list_of_args_per_blk.h"
            );
            // CHECK_ERROR(cudaDeviceSynchronize(), "device sync in computeBSSN");
            

            CHECK_ERROR(cudaMemcpyAsync(var_out_array[blk.block_no], dev_var_out_array[index], BSSN_NUM_VARS*unzip_dof*sizeof(double), cudaMemcpyDeviceToHost, stream), "dev_var_out_array[index] cudaMemcpyDeviceToHost");
            // CHECK_ERROR(cudaDeviceSynchronize(), "device sync in computeBSSN");
        }

        CHECK_ERROR(cudaDeviceSynchronize(), "device sync in computeBSSN");

       

        // Release GPU memory
        #include "bssnrhs_cuda_mdealloc.h"
        #include "bssnrhs_cuda_mdealloc_adv.h"

        for (int index=init_block; index<=current_block-1; index++){
            CHECK_ERROR(cudaFree(dev_var_in_array[index]), "dev_var_in_array[index] cudaFree");
            CHECK_ERROR(cudaFree(dev_var_out_array[index]), "dev_var_out_array[index] cudaFree");
        }
    }
    for (int index=0; index<steamCountToLevel[0]; index++){
        CHECK_ERROR(cudaStreamDestroy(streams[index]), "cudaStream destruction");
    }
    return;
}

// void GPU_parallelized_async_hybrid(unsigned int numberOfLevels, Block * blkList, unsigned int lower_bound, unsigned int upper_bound, double ** var_in_array, double ** var_out_array){ 
//     CHECK_ERROR(cudaSetDevice(0), "cudaSetDevice in computeBSSN"); // Set the GPU that we are going to deal with

//     // Creating cuda streams for the process
//     int num_streams = 2;
//     if (num_streams<steamCountToLevel[0]) num_streams=steamCountToLevel[0];
//     cudaStream_t stream;
//     cudaStream_t streams[num_streams]; // usually steamCountToLevel[0] should be the max number of streams verify it before execute.
//     for (int index=0; index<num_streams; index++){
//         CHECK_ERROR(cudaStreamCreate(&streams[index]), "cudaStream creation");
//     }

//     // Check for available GPU memory
//     size_t free_bytes, total_bytes;
//     CHECK_ERROR(cudaMemGetInfo(&free_bytes, &total_bytes), "Available GPU memory checking failed");
//     double GPU_capacity_buffer = 10;
//     double GPUCapacity = 1.0*free_bytes/1024/1024 - GPU_capacity_buffer;
//     std::cout << "Available GPU with buffer of " << GPU_capacity_buffer << ": " << GPUCapacity << " | Total GPU memory: " << total_bytes/1024/1024 << std::endl << std::endl;

//     // Sort the data block list
//     mergeSort(blkList, 0, numberOfLevels-1); // O(nlog(n))

//     // Calculate how much data block can be executed at once in GPU
//     Block blk;
//     int init_block = 0;
//     int current_block = 0;
//     int total_points = 0;
//     int numberOfStreams = 0;
//     double current_usage = 0;
//     double fixed_usage = 0;
//     double prev_usage = 0;
//     double actual_usage = 0;
//     double malloc_overhead = 3;

//     while (current_block<numberOfLevels){
//         current_usage=0;
//         fixed_usage=0;
//         init_block=current_block;

//         while ((current_usage+fixed_usage<(GPUCapacity)) && (current_block<numberOfLevels)){
//             prev_usage = current_usage+fixed_usage; // usage of previous iteration
//             blk = blkList[current_block];
//             total_points = blk.blkSize;
//             if (blk.blkLevel<hybrid_divider) {
//                 numberOfStreams = steamCountToLevel[blk.blkLevel];
//             }else{
//                 numberOfStreams = 1;
//             }
//             if (fixed_usage<numberOfStreams*(138+72)*total_points*sizeof(double)/1024/1024 + 210*malloc_overhead + (current_block-init_block)*malloc_overhead){
//                 fixed_usage = numberOfStreams*(138+72)*total_points*sizeof(double)/1024/1024 + 210*malloc_overhead + (current_block-init_block)*malloc_overhead; // 5*210 means allocation overhead
//             }
//             current_usage += (total_points*BSSN_NUM_VARS*sizeof(double)*2)/1024/1024;
//             current_block++;
//         }
//         actual_usage = current_usage+fixed_usage;
//         if (current_usage+fixed_usage>(GPUCapacity)){
//             current_block--;
//             if (init_block>current_block-1) {
//                 std::cout << "Required GPU memory = " << actual_usage << " Failed to allocate enough memory. Program terminated..." << std::endl;
//                 exit(0);
//             }
//             actual_usage = prev_usage;
//         }

//         // Display the set of blocks selected to process with their GPU usage
//         std::cout << "start: " << init_block << " end: " << current_block-1 << "| usage: " << actual_usage << std::endl;

//         // Allocating device memory to hold input and output
//         double ** dev_var_in_array = new double*[numberOfLevels];
//         double ** dev_var_out_array = new double*[numberOfLevels];

//         int largest_intermediate_array = 0;
//         for (int index=init_block; index<=current_block-1; index++){
//             blk = blkList[index];

//             if (blk.blkLevel<hybrid_divider) {
//                 numberOfStreams = steamCountToLevel[blk.blkLevel];
//             }else{
//                 numberOfStreams = 1;
//             }

//             if (largest_intermediate_array<blk.blkSize*numberOfStreams){
//                 largest_intermediate_array = blk.blkSize*numberOfStreams;
//             }
//             CHECK_ERROR(cudaMalloc((void**)&dev_var_in_array[index], blk.blkSize*BSSN_NUM_VARS*sizeof(double)), "dev_var_in_array[index]");
//             CHECK_ERROR(cudaMalloc((void**)&dev_var_out_array[index], blk.blkSize*BSSN_NUM_VARS*sizeof(double)), "dev_var_out_array[index]");
//         }
        
//         // Allocation intermediate arrays
//         int size = largest_intermediate_array * sizeof(double);
//         #include "bssnrhs_cuda_variable_malloc.h"
//         #include "bssnrhs_cuda_variable_malloc_adv.h"
//         #include "bssnrhs_cuda_malloc.h"
//         #include "bssnrhs_cuda_malloc_adv.h"

        
//         // Start block processing
//         int streamIndex;
//         int unzip_dof;
//         unsigned int sz[3];
//         double ptmin[3], ptmax[3];
//         unsigned int bflag;
//         double dx, dy, dz;
//         bool isAsyncStarted = false;
//         for(int index=init_block; index<=current_block-1; index++) {
//             blk=blkList[index];

//             // Configure the block specific values
//             unzip_dof=blk.blkSize;

//             sz[0]=blk.node1D_x; 
//             sz[1]=blk.node1D_y;
//             sz[2]=blk.node1D_z;

//             bflag=0;

//             dx=0.1;
//             dy=0.1;
//             dz=0.1;

//             ptmin[0]=0.0;
//             ptmin[1]=0.0;
//             ptmin[2]=0.0;

//             ptmax[0]=1.0;
//             ptmax[1]=1.0;
//             ptmax[2]=1.0;

//             // Block parallel processing
//             if (blk.blkLevel<hybrid_divider){
//                 // Call cudasync in the case of block level change occured
//                 if (index!=init_block && blkList[index-1].blkLevel!=blk.blkLevel) {
//                     CHECK_ERROR(cudaDeviceSynchronize(), "device sync in computeBSSN level change");
//                 }

//                 // Identify stream to schedule
//                 numberOfStreams = steamCountToLevel[blk.blkLevel];
//                 streamIndex = index % numberOfStreams;
//                 stream = streams[streamIndex];

//                 // std::cout << "GPU - Count: " << std::setw(3) << index << " - Block no: " << std::setw(3) << blk.block_no << " - Bock level: " << std::setw(1) << blk.blkLevel << " - Block size: " << blk.blkSize << std::endl;

//                 CHECK_ERROR(cudaMemcpyAsync(dev_var_in_array[index], var_in_array[blk.block_no], BSSN_NUM_VARS*unzip_dof*sizeof(double), cudaMemcpyHostToDevice, stream), "dev_var_in_array[index] cudaMemcpyHostToDevice");

//                 // cuda_bssnrhs(dev_var_out_array[index], dev_var_in_array[index], unzip_dof, ptmin, ptmax, sz, bflag, stream,
//                 //     #include "list_of_args_per_blk.h"
//                 // );

//                 CHECK_ERROR(cudaMemcpyAsync(var_out_array[blk.block_no], dev_var_out_array[index], BSSN_NUM_VARS*unzip_dof*sizeof(double), cudaMemcpyDeviceToHost, stream), "dev_var_out_array[index] cudaMemcpyDeviceToHost");

//             }else{
//                 // Starting async process

//                 // Sync any remaining process
//                 if (!isAsyncStarted) {
//                     streamIndex = 0;
//                     CHECK_ERROR(cudaDeviceSynchronize(), "device sync in computeBSSN");
//                 }

//                 // std::cout << "GPU - Count: " << std::setw(3) << index << " - Block no: " << std::setw(3) << blk.block_no << " - Bock level: " << std::setw(1) << blk.blkLevel << " - Block size: " << blk.blkSize << std::endl;

//                 CHECK_ERROR(cudaMemcpyAsync(dev_var_in_array[index], var_in_array[blk.block_no], BSSN_NUM_VARS*unzip_dof*sizeof(double), cudaMemcpyHostToDevice, streams[(index)%2]), "dev_var_in_array[index] cudaMemcpyHostToDevice");

//                 if (isAsyncStarted){
//                     cudaStreamSynchronize(streams[(index+1)%2]);
//                 }

//                 // cuda_bssnrhs(dev_var_out_array[index], dev_var_in_array[index], unzip_dof, ptmin, ptmax, sz, bflag, streams[(index)%2],
//                 //     #include "list_of_args_per_blk.h"
//                 // );

//                 if (isAsyncStarted){
//                     CHECK_ERROR(cudaMemcpyAsync(var_out_array[blkList[index-1].block_no], dev_var_out_array[index-1], BSSN_NUM_VARS*blkList[index-1].blkSize*sizeof(double), cudaMemcpyDeviceToHost, streams[(index+1)%2]), "dev_var_out_array[index] cudaMemcpyDeviceToHost");
//                 }

//                 //handle last DtoH memcopy
//                 if(index==current_block-1){
//                     CHECK_ERROR(cudaMemcpyAsync(var_out_array[blk.block_no], dev_var_out_array[index], BSSN_NUM_VARS*unzip_dof*sizeof(double), cudaMemcpyDeviceToHost, streams[(index)%2]), "dev_var_out_array[index] cudaMemcpyDeviceToHost");
//                 }

//                 if (!isAsyncStarted) isAsyncStarted=true;
//             }
            
//         }

//         CHECK_ERROR(cudaDeviceSynchronize(), "device sync in computeBSSN");

//         // Release GPU memory
//         #include "bssnrhs_cuda_mdealloc.h"
//         #include "bssnrhs_cuda_mdealloc_adv.h"

//         for (int index=init_block; index<=current_block-1; index++){
//             CHECK_ERROR(cudaFree(dev_var_in_array[index]), "dev_var_in_array[index] cudaFree");
//             CHECK_ERROR(cudaFree(dev_var_out_array[index]), "dev_var_out_array[index] cudaFree");
//         }
//     }
//     for (int index=0; index<num_streams; index++){
//         CHECK_ERROR(cudaStreamDestroy(streams[index]), "cudaStream destruction");
//     }
//     return;
// }

// void GPU_pure_async_htod_dtoH_overlap(unsigned int numberOfLevels, Block * blkList, unsigned int lower_bound, unsigned int upper_bound, double ** var_in_array, double ** var_out_array){ 
//     CHECK_ERROR(cudaSetDevice(0), "cudaSetDevice in computeBSSN"); // Set the GPU that we are going to deal with

//     // Creating cuda streams for the process
//     int num_streams = 3;
//     cudaStream_t stream;
//     cudaStream_t streams[num_streams]; // usually steamCountToLevel[0] should be the max number of streams verify it before execute.
//     for (int index=0; index<num_streams; index++){
//         CHECK_ERROR(cudaStreamCreate(&streams[index]), "cudaStream creation");
//     }

//     // Check for available GPU memory
//     size_t free_bytes, total_bytes;
//     CHECK_ERROR(cudaMemGetInfo(&free_bytes, &total_bytes), "Available GPU memory checking failed");
//     double GPU_capacity_buffer = 10;
//     double GPUCapacity = 1.0*free_bytes/1024/1024 - GPU_capacity_buffer;
//     std::cout << "Available GPU with buffer of " << GPU_capacity_buffer << ": " << GPUCapacity << " | Total GPU memory: " << total_bytes/1024/1024 << std::endl << std::endl;

//     // Sort the data block list
//     mergeSort(blkList, 0, numberOfLevels-1); // O(nlog(n))

//     // Calculate how much data block can be executed at once in GPU
//     Block blk;
//     int init_block = 0;
//     int current_block = 0;
//     int total_points = 0;
//     int numberOfStreams = 0;
//     double current_usage = 0;
//     double fixed_usage = 0;
//     double prev_usage = 0;
//     double actual_usage = 0;
//     double malloc_overhead = 3;

//     while (current_block<numberOfLevels){
//         current_usage=0;
//         fixed_usage=0;
//         init_block=current_block;

//         while ((current_usage+fixed_usage<(GPUCapacity)) && (current_block<numberOfLevels)){
//             prev_usage = current_usage+fixed_usage; // usage of previous iteration
//             blk = blkList[current_block];
//             total_points = blk.blkSize;
//             numberOfStreams = 1;
//             if (fixed_usage<numberOfStreams*(138+72)*total_points*sizeof(double)/1024/1024 + 210*malloc_overhead + (current_block-init_block)*malloc_overhead){
//                 fixed_usage = numberOfStreams*(138+72)*total_points*sizeof(double)/1024/1024 + 210*malloc_overhead + (current_block-init_block)*malloc_overhead; // 5*210 means allocation overhead
//             }
//             current_usage += (total_points*BSSN_NUM_VARS*sizeof(double)*2)/1024/1024;
//             current_block++;
//         }
//         actual_usage = current_usage+fixed_usage;
//         if (current_usage+fixed_usage>(GPUCapacity)){
//             current_block--;
//             if (init_block>current_block-1) {
//                 std::cout << "Required GPU memory = " << actual_usage << " Failed to allocate enough memory. Program terminated..." << std::endl;
//                 exit(0);
//             }
//             actual_usage = prev_usage;
//         }

//         // Display the set of blocks selected to process with their GPU usage
//         std::cout << "start: " << init_block << " end: " << current_block-1 << "| usage: " << actual_usage << std::endl;

//         // Allocating device memory to hold input and output
//         double ** dev_var_in_array = new double*[numberOfLevels];
//         double ** dev_var_out_array = new double*[numberOfLevels];

//         int largest_intermediate_array = 0;
//         for (int index=init_block; index<=current_block-1; index++){
//             blk = blkList[index];
//             numberOfStreams = 1;

//             if (largest_intermediate_array<blk.blkSize*numberOfStreams){
//                 largest_intermediate_array = blk.blkSize*numberOfStreams;
//             }
//             CHECK_ERROR(cudaMalloc((void**)&dev_var_in_array[index], blk.blkSize*BSSN_NUM_VARS*sizeof(double)), "dev_var_in_array[index]");
//             CHECK_ERROR(cudaMalloc((void**)&dev_var_out_array[index], blk.blkSize*BSSN_NUM_VARS*sizeof(double)), "dev_var_out_array[index]");
//         }
        
//         // Allocation intermediate arrays
//         int size = largest_intermediate_array * sizeof(double);
//         #include "bssnrhs_cuda_variable_malloc.h"
//         #include "bssnrhs_cuda_variable_malloc_adv.h"
//         #include "bssnrhs_cuda_malloc.h"
//         #include "bssnrhs_cuda_malloc_adv.h"
        
//         // Start block processing
//         int streamIndex = 0;
//         int unzip_dof;
//         unsigned int sz[3];
//         double ptmin[3], ptmax[3];
//         unsigned int bflag;
//         double dx, dy, dz;
//         bool isAsyncStarted = false;
//         for(int index=init_block; index<=current_block-1; index++) {
//             blk=blkList[index];

//             // Configure the block specific values
//             unzip_dof=blk.blkSize;

//             sz[0]=blk.node1D_x; 
//             sz[1]=blk.node1D_y;
//             sz[2]=blk.node1D_z;

//             bflag=0;

//             dx=0.1;
//             dy=0.1;
//             dz=0.1;

//             ptmin[0]=0.0;
//             ptmin[1]=0.0;
//             ptmin[2]=0.0;

//             ptmax[0]=1.0;
//             ptmax[1]=1.0;
//             ptmax[2]=1.0;

//             // Starting async process

//             // Sync any remaining process
//             // std::cout << "GPU - Count: " << std::setw(3) << index << " - Block no: " << std::setw(3) << blk.block_no << " - Bock level: " << std::setw(1) << blk.blkLevel << " - Block size: " << blk.blkSize << std::endl;

//             CHECK_ERROR(cudaMemcpyAsync(dev_var_in_array[index], var_in_array[blk.block_no], BSSN_NUM_VARS*unzip_dof*sizeof(double), cudaMemcpyHostToDevice, streams[(index)%3]), "dev_var_in_array[index] cudaMemcpyHostToDevice");

//             if (isAsyncStarted){
//                 cudaStreamSynchronize(streams[(index+2)%3]);
//             }

//             // cuda_bssnrhs(dev_var_out_array[index], dev_var_in_array[index], unzip_dof, ptmin, ptmax, sz, bflag, streams[(index)%3],
//             //     #include "list_of_args_per_blk.h"
//             // );

//             if (isAsyncStarted){
//                 CHECK_ERROR(cudaMemcpyAsync(var_out_array[blkList[index-1].block_no], dev_var_out_array[index-1], BSSN_NUM_VARS*blkList[index-1].blkSize*sizeof(double), cudaMemcpyDeviceToHost, streams[(index+2)%3]), "dev_var_out_array[index] cudaMemcpyDeviceToHost");
//             }

//             //handle last DtoH memcopy
//             if(index==current_block-1){
//                 CHECK_ERROR(cudaMemcpyAsync(var_out_array[blk.block_no], dev_var_out_array[index], BSSN_NUM_VARS*unzip_dof*sizeof(double), cudaMemcpyDeviceToHost, streams[(index)%3]), "dev_var_out_array[index] cudaMemcpyDeviceToHost");
//             }

//             if (!isAsyncStarted) isAsyncStarted=true;
//         }

//         CHECK_ERROR(cudaDeviceSynchronize(), "device sync in computeBSSN");

//         // Release GPU memory
//         #include "bssnrhs_cuda_mdealloc.h"
//         #include "bssnrhs_cuda_mdealloc_adv.h"

//         for (int index=init_block; index<=current_block-1; index++){
//             CHECK_ERROR(cudaFree(dev_var_in_array[index]), "dev_var_in_array[index] cudaFree");
//             CHECK_ERROR(cudaFree(dev_var_out_array[index]), "dev_var_out_array[index] cudaFree");
//         }
//     }
//     for (int index=0; index<num_streams; index++){
//         CHECK_ERROR(cudaStreamDestroy(streams[index]), "cudaStream destruction");
//     }
//     return;
// }
