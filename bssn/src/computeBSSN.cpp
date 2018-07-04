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

void data_generation_blockwise_mixed(double mean, double std, unsigned int numberOfLevels, unsigned int lower_bound, unsigned int upper_bound, bool isRandom, 
    Block * blkList, double ** var_in_array, double ** var_out_array){

    const unsigned int maxDepth=12;

    // Calculate block sizes
    std::map<int, double> block_sizes;
    for(int i=lower_bound; i<=upper_bound; i++){
        Block blk = Block(0, 0, 0, 2*i, i, maxDepth);
        block_sizes[i] = 1.0*blk.blkSize*BSSN_NUM_VARS*sizeof(double)/1024/1024;
    }

    // Create distribution
    int seed = 100;
    if (isRandom) seed=time(0);
    std::default_random_engine generator(seed);
    std::normal_distribution<double> distribution(mean, std);

    // Generating levels and block structure
    int total_grid_points = 0;
    int p[upper_bound-lower_bound+1]={};
    int level;
    int block_no = 0;
    while (block_no<numberOfLevels){
        level = int(distribution(generator)); // generate a number

        Block & blk=blkList[block_no];
        blk=Block(0, 0, 0, 2*level, level, maxDepth);

        total_grid_points += blk.blkSize;
           
        // distribution representation requirement
        if ((level>=lower_bound)&&(level<=upper_bound)) {
            block_no++;
            p[level-lower_bound]++;
        }
    }

    // distribution representation requirement
    std::cout << "---Level distribution---" << std::endl;
    double ram_requirement = 0.0;
    for (int i=0; i<=upper_bound-lower_bound; i++) {
        ram_requirement += block_sizes[i+lower_bound]*p[i]; // this is calculated in MBs

        std::cout << i+lower_bound << ": ";
        std::cout << std::setw(3) << p[i]*100/numberOfLevels << "% " << std::string(p[i]*100/numberOfLevels, '*') << std::endl;
    }
    std::cout << "Total blocks: " << block_no << " | Total grid points: " << 1.0*total_grid_points/1000000 << "x10^6" << std::endl;
    std::cout << "Total RAM requiremnet for input: " << ram_requirement << "MB" << std::endl;
    std::cout << "Total RAM requiremnet for both input/output: " << ram_requirement*2 << "MB" << std::endl;

    // Checking for available ram
    struct sysinfo myinfo; 
    sysinfo(&myinfo); 
    double total_MB = 1.0*myinfo.totalram/1024/1024; 
    double available_MB = 1.0*myinfo.freeram/1024/1024; 
    std::cout << "Total RAM: " << total_MB << "MB Available RAM: " << available_MB << "MB" << std::endl;

    double remaining = available_MB - ram_requirement*2;
    std::cout << "\nRemaining RAM if the process is going to continue: " << remaining << "MB" << std::endl;
    if (remaining>400){
        std::cout << "Data block generating starts..." << std::endl;
    } else if(remaining>0){
        std::cout << "Are you sure that you want to continue(Y/n):";
        char sure;
        cin >> sure;
        if (tolower(sure)!='y') {
            std::cout << "Program terminated..." << std::endl;
            exit(0);
        }
        std::cout << "Data block generating starts..." << std::endl;
    }else{
        std::cout << "Not enough RAM to process. Program terminated..." << std::endl;
        exit(0);
    }

    // Data block generation
    #pragma omp parallel for num_threads(4)
    for (int index=0; index<numberOfLevels; index++){
        Block & blk=blkList[index];
        const unsigned long unzip_dof=blk.blkSize;

        // Allocating pinned memory in RAM
        double * var_in_per_block;
        CHECK_ERROR(cudaMallocHost((void**)&var_in_per_block, unzip_dof*BSSN_NUM_VARS*sizeof(double)), "var_in_per_block");

        double * var_out_per_block = new double[unzip_dof*BSSN_NUM_VARS];
        CHECK_ERROR(cudaMallocHost((void**)&var_out_per_block, unzip_dof*BSSN_NUM_VARS*sizeof(double)), "var_out_per_block");

        var_in_array[index] = var_in_per_block;
        var_out_array[index] = var_out_per_block;

        double coord[3];
        double u[BSSN_NUM_VARS];
        double x,y,z,hx,hy,hz;
        unsigned int size_x,size_y,size_z;

        x=(double)blk.x;
        y=(double)blk.y;
        z=(double)blk.z;

        hx=0.001; 
        hy=0.001; 
        hz=0.001; 

        size_x=blk.node1D_x;
        size_y=blk.node1D_y;
        size_z=blk.node1D_z;

        for(unsigned int k=0; k<blk.node1D_z; k++){
            for(unsigned int j=0; j<blk.node1D_y; j++){
                for(unsigned int i=0; i<blk.node1D_x; i++){
                    coord[0]=x+i*hx;
                    coord[1]=y+j*hy;
                    coord[2]=z+k*hz;

                    initial_data(u,coord);

                    for(unsigned int var=0; var<BSSN_NUM_VARS; var++){
                        var_in_per_block[var*unzip_dof+k*size_y*size_x+j*size_y+i]=u[var];
                    }
                }
            }
        }
    }
}

void data_generation_blockwise_and_bssn_var_wise_mixed(double mean, double std, unsigned int numberOfLevels, unsigned int lower_bound, unsigned int upper_bound, bool isRandom,
    Block * blkList, double ** var_in_array, double ** var_out_array, double ** var_in, double ** var_out){

    const unsigned int maxDepth=12;

    // Calculate block sizes
    std::map<int, double> block_sizes;
    for(int i=lower_bound; i<=upper_bound; i++){
        Block blk = Block(0, 0, 0, 2*i, i, maxDepth);
        block_sizes[i] = 1.0*(blk.node1D_x*blk.node1D_y*blk.node1D_z)*BSSN_NUM_VARS*sizeof(double)/1024/1024;
    }

    // Create distribution
    int seed = 100;
    if (isRandom) seed=time(0);
    std::default_random_engine generator(seed);
    std::normal_distribution<double> distribution(mean, std);

    // Generating levels and block structure
    int p[upper_bound-lower_bound+1]={};
    int level;
    int block_no = 0;

    // RAM ---
    unsigned long unzipSz=0;
    while (block_no<numberOfLevels){
        level = int(distribution(generator)); // generate a number

        Block & blk=blkList[block_no];
        blk=Block(0, 0, 0, 2*level, level, maxDepth);

        // RAM ---
           blk.offset=unzipSz;
           unzipSz+=(blk.node1D_x*blk.node1D_y*blk.node1D_z);
           
        // distribution representation requirement
        if ((level>=lower_bound)&&(level<=upper_bound)) {
            block_no++;
            p[level-lower_bound]++;
        }
       }
    // RAM ---
    const unsigned long unzip_dof_cpu=unzipSz;

    // distribution representation requirement
    std::cout << "---Level distribution---" << std::endl;
    double ram_requirement = 0.0;
    for (int i=0; i<=upper_bound-lower_bound; i++) {
        ram_requirement += block_sizes[i+lower_bound]*p[i]; // this is calculated in MBs

        std::cout << i+lower_bound << ": ";
        std::cout << std::setw(3) << p[i]*100/numberOfLevels << "% " << std::string(p[i]*100/numberOfLevels, '*') << std::endl;
    }
    std::cout << "Total blocks: " << block_no << std::endl;
    std::cout << "Total RAM requiremnet for input: " << ram_requirement << "MB" << std::endl;
    std::cout << "Total RAM requiremnet for both input/output: " << ram_requirement*2 << "MB" << std::endl;
    // RAM ---
    std::cout << "Total RAM requiremnet for both GPU and CPU version: " << ram_requirement*2*2 + block_sizes[upper_bound]*210/BSSN_NUM_VARS << "MB" << std::endl;

    // Checking for available ram
    struct sysinfo myinfo; 
    sysinfo(&myinfo); 
    double total_MB = 1.0*myinfo.totalram/1024/1024; 
    double available_MB = 1.0*myinfo.freeram/1024/1024; 
    std::cout << "Total RAM: " << total_MB << "MB Available RAM: " << available_MB << "MB" << std::endl;

    double remaining = available_MB - ram_requirement*2*2 - block_sizes[upper_bound]*210/BSSN_NUM_VARS;
    std::cout << "\nRemaining RAM if the process is going to continue: " << remaining << "MB" << std::endl;
    if (remaining>400){
        std::cout << "Data block generating starts..." << std::endl;
    } else if(remaining>0){
        std::cout << "Are you sure that you want to continue(Y/n):";
        char sure;
        cin >> sure;
        if (tolower(sure)!='y') {
            std::cout << "Program terminated..." << std::endl;
            exit(0);
        }
        std::cout << "Data block generating starts..." << std::endl;
    }else{
        std::cout << "Not enough RAM to process. Program terminated..." << std::endl;
        exit(0);
    }

    // Data block generation
    #pragma omp parallel for num_threads(4)
    for (int index=0; index<numberOfLevels; index++){
        Block & blk=blkList[index];
        const unsigned long unzip_dof=(blk.node1D_x*blk.node1D_y*blk.node1D_z);

        // Allocating pinned memory in RAM
        double * var_in_per_block;
        CHECK_ERROR(cudaMallocHost((void**)&var_in_per_block, unzip_dof*BSSN_NUM_VARS*sizeof(double)), "var_in_per_block");

        double * var_out_per_block = new double[unzip_dof*BSSN_NUM_VARS];
        CHECK_ERROR(cudaMallocHost((void**)&var_out_per_block, unzip_dof*BSSN_NUM_VARS*sizeof(double)), "var_out_per_block");

        var_in_array[index] = var_in_per_block;
        var_out_array[index] = var_out_per_block;

        double coord[3];
        double u[BSSN_NUM_VARS];
        double x,y,z,hx,hy,hz;
        unsigned int size_x,size_y,size_z;
        
        x=(double)blk.x;
        y=(double)blk.y;
        z=(double)blk.z;

        hx=0.001; 
        hy=0.001; 
        hz=0.001; 

        size_x=blk.node1D_x;
        size_y=blk.node1D_y;
        size_z=blk.node1D_z;

        for(unsigned int k=0; k<blk.node1D_z; k++){
            for(unsigned int j=0; j<blk.node1D_y; j++){
                for(unsigned int i=0; i<blk.node1D_x; i++){
                    coord[0]=x+i*hx;
                    coord[1]=y+j*hy;
                    coord[2]=z+k*hz;

                    initial_data(u,coord);

                    for(unsigned int var=0; var<BSSN_NUM_VARS; var++){
                        var_in_per_block[var*unzip_dof+k*size_y*size_x+j*size_y+i]=u[var];
    }
                }
            }
        }
    }

    // RAM ---
    // Data generation for CPU version
    for(int i=0;i<BSSN_NUM_VARS;i++){
        var_in[i] = new double[unzip_dof_cpu];
        var_out[i] = new double[unzip_dof_cpu];
    }
        
    #pragma omp parallel for num_threads(4)
    for(unsigned int blk=0; blk<numberOfLevels; blk++){
        
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

void GPU_Async_Iteration_Wise(unsigned int numberOfLevels, Block * blkList, double ** var_in_array, double ** var_out_array){

    CHECK_ERROR(cudaSetDevice(0), "cudaSetDevice in computeBSSN");

    cudaEvent_t start[2], end[2];
    for (int i=0; i<2; i++){
        cudaEventCreate(&start[i]);
        cudaEventCreate(&end[i]);
    }

    // Check for available GPU memory
    size_t free_bytes, total_bytes;
    CHECK_ERROR(cudaMemGetInfo(&free_bytes, &total_bytes), "Available GPU memory checking failed");
    double GPUCapacity = 1.0*free_bytes/1024/1024 - 600;
    std::cout << "Available GPU with buffer of 600MB: " << GPUCapacity << " Total GPU memory: " << total_bytes/1024/1024 << std::endl << std::endl;

    double ** dev_var_in_array = new double*[numberOfLevels];
    double ** dev_var_out_array = new double*[numberOfLevels];

    Block blk;
    int current_block = 0;
    int init_block = 0;
    double current_usage, fixed_usage;
    int total_points;
    cudaStream_t stream;
    cudaStream_t streams[numberOfStreams];

    for (int index=0; index<numberOfStreams; index++){
        CHECK_ERROR(cudaStreamCreate(&streams[index]), "cudaStream creation");
    }

    // Check number of blocks can handle at once
    while (current_block<numberOfLevels){
        current_usage=0;
        fixed_usage=0;
        init_block=current_block;

        while ((current_usage+fixed_usage<(GPUCapacity)) && (current_block<numberOfLevels)){
            blk = blkList[current_block];
            total_points = blk.blkSize;

            if (fixed_usage<numberOfStreams*(138+72)*total_points*sizeof(double)/1024/1024){
                fixed_usage = numberOfStreams*(138+72)*total_points*sizeof(double)/1024/1024;
            }

            current_usage += (total_points*BSSN_NUM_VARS*sizeof(double)*2)/1024/1024;

            // std::cout << "\tUsage: " << current_usage+fixed_usage << " BlockNumber: " << current_block << " Fixed Usage: " << fixed_usage << std::endl;
            current_block++;
        }
        if (current_usage+fixed_usage>(GPUCapacity)){
            current_block--;
            if (init_block>current_block-1) {
                std::cout << "Required GPU memory = " << current_usage+fixed_usage << " Failed to allocate enough memory. Program terminated..." << std::endl;
                exit(0);
            }
        }
        
        std::cout << "start: " << init_block << " end: " << current_block-1 << "| usage: " << current_usage+fixed_usage << std::endl;

        int dof;
        int max_unzip_dof = 0;
        for (int index=init_block; index<=current_block-1; index++){
            blk = blkList[index];
            dof = blk.blkSize;
            if (max_unzip_dof<dof){
                max_unzip_dof=dof;
            }
            CHECK_ERROR(cudaMalloc((void**)&dev_var_in_array[index], dof*BSSN_NUM_VARS*sizeof(double)), "dev_var_in_array[index]");
            CHECK_ERROR(cudaMalloc((void**)&dev_var_out_array[index], dof*BSSN_NUM_VARS*sizeof(double)), "dev_var_out_array[index]");
        }
        
        int size = numberOfStreams * max_unzip_dof * sizeof(double);

        #include "bssnrhs_cuda_variable_malloc.h"
        #include "bssnrhs_cuda_variable_malloc_adv.h"
        #include "bssnrhs_cuda_malloc.h"
        #include "bssnrhs_cuda_malloc_adv.h"
        
        bssn::timer::t_memcopy_kernel.start();
        for(int index=init_block; index<=current_block-1; index++) {

            cudaStream_t stream;
            int streamIndex = index % numberOfStreams;
            stream = streams[streamIndex];

            Block & blk=blkList[index];

            unsigned long unzip_dof=blk.blkSize;

            unsigned int sz[3];
            sz[0]=blk.node1D_x; 
            sz[1]=blk.node1D_y;
            sz[2]=blk.node1D_z;

            unsigned int bflag=0;

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

            std::cout << "GPU - Block no: " << std::setw(3) << index << " - Bock level: " << std::setw(1) << blk.blkLevel << " - Block size: " << blk.blkSize << std::endl;

            if (index==init_block) cudaEventRecord(start[0], stream);
            CHECK_ERROR(cudaMemcpyAsync(dev_var_in_array[index], var_in_array[index], BSSN_NUM_VARS*unzip_dof*sizeof(double), cudaMemcpyHostToDevice, stream), "dev_var_in_array[index] cudaMemcpyHostToDevice");
            if (index==init_block) cudaEventRecord(end[0], stream);

            cuda_bssnrhs(dev_var_out_array[index], dev_var_in_array[index], unzip_dof, ptmin, ptmax, sz, bflag, stream,
                #include "list_of_args_per_blk.h"
            );

            if (index==current_block-1) cudaEventRecord(start[1], stream);
            CHECK_ERROR(cudaMemcpyAsync(var_out_array[index], dev_var_out_array[index], BSSN_NUM_VARS*unzip_dof*sizeof(double), cudaMemcpyDeviceToHost, stream), "dev_var_out_array[index] cudaMemcpyDeviceToHost");
            if (index==current_block-1) cudaEventRecord(end[1], stream);
        }

        CHECK_ERROR(cudaDeviceSynchronize(), "device sync in computeBSSN");

        float HtoD_ms, DtoH_ms;
        cudaEventElapsedTime(&HtoD_ms, start[0], end[0]);
        cudaEventElapsedTime(&DtoH_ms, start[1], end[1]);
        
        bssn::timer::t_memcopy.setTime((HtoD_ms+DtoH_ms)/1000);
        bssn::timer::t_memcopy_kernel.stop();

        // Release GPU memory
        #include "bssnrhs_cuda_mdealloc.h"
        #include "bssnrhs_cuda_mdealloc_adv.h"

        for (int index=init_block; index<=current_block-1; index++){
            CHECK_ERROR(cudaFree(dev_var_in_array[index]), "dev_var_in_array[index] cudaFree");
            CHECK_ERROR(cudaFree(dev_var_out_array[index]), "dev_var_out_array[index] cudaFree");
        }
    }
    return;
}

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

        std::cout << "CPU - Block no: " << std::setw(3) << blk << " - Bock level: " << std::setw(1) << blkList[blk].blkLevel << " - Block size: " << blkList[blk].blkSize << std::endl;
        
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

    if(argc<8){
        std::cout<<"Usage: " << argv[0] << " mean(double) std(double) numberOfBlocks(int) lowerLevel(int) upperLevel(int) isRandom(1/0) isTest(1/0)" <<std::endl;
        exit(0);
    }

    double mean = atoi(argv[1]);
    double std = atoi(argv[2]);
    unsigned int numberOfLevels = atoi(argv[3]);
    unsigned int lower_bound = atoi(argv[4]);
    unsigned int upper_bound = atoi(argv[5]);
    bool isRandom = atoi(argv[6]);
    bool isTest = atoi(argv[7]);

    Block * blkList = new Block[numberOfLevels];
    double ** var_in_array = new double*[numberOfLevels];
    double ** var_out_array = new double*[numberOfLevels];

    double ** var_in;
    double ** var_out;
    if (isTest){
        var_in = new double*[BSSN_NUM_VARS];
        var_out = new double*[BSSN_NUM_VARS];
        data_generation_blockwise_and_bssn_var_wise_mixed(mean, std, numberOfLevels, lower_bound, upper_bound, isRandom, blkList, var_in_array, var_out_array, var_in, var_out);
    }else{
        data_generation_blockwise_mixed(mean, std, numberOfLevels, lower_bound, upper_bound, isRandom, blkList, var_in_array, var_out_array);
    }
    
    #include "rhs_cuda.h"
    bssn::timer::t_gpu_runtime.start();
    GPU_Async_Iteration_Wise(numberOfLevels, blkList, var_in_array, var_out_array);
    bssn::timer::t_gpu_runtime.stop();

    std::cout << "" << std::endl;

    if (isTest){
        #include "rhs.h"
        bssn::timer::t_cpu_runtime.start();
        CPU_sequence(numberOfLevels, blkList, var_in, var_out);
        bssn::timer::t_cpu_runtime.stop();

        // Verify outputs
            for (int blk=0; blk<numberOfLevels; blk++){
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

        for (int var=0; var<BSSN_NUM_VARS; var++){
            delete [] var_in[var];
            delete [] var_out[var];
        }
    }

    //Free host memory
    for (int blk=0; blk<numberOfLevels; blk++){
        CHECK_ERROR(cudaFreeHost(var_in_array[blk]), "free host memory");
        CHECK_ERROR(cudaFreeHost(var_out_array[blk]), "free host memory");
    }

    bssn::timer::profileInfo();
    return 0;
}