/**
 * Created on: Sep 21, 2018
 * 		Author: Akila, Eranga, Eminda, Ruwan
 **/

#include "dataGeneration.h"

void data_generation_blockwise_mixed(double mean, double std, unsigned int numberOfLevels, unsigned int lower_bound, unsigned int upper_bound, bool isRandom, Block * blkList, double ** var_in_array, double ** var_out_array)
{
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
    std::normal_distribution<double> distributionN(mean, std);
    std::uniform_int_distribution<int> distributionU(lower_bound, upper_bound);

    // Generating levels and block structure
    int total_grid_points = 0;
    int p[upper_bound-lower_bound+1];
    for (int i=0; i<upper_bound-lower_bound+1; i++) p[i] = 0;

    int level;
    int block_no = 0;
    while (block_no<numberOfLevels){
        level = int(distributionN(generator)); // generate a number

        // distribution representation requirement
        if ((level>=lower_bound)&&(level<=upper_bound)) {
            Block & blk=blkList[block_no];
            blk=Block(0, 0, 0, 2*level, level, maxDepth);
            blk.block_no = block_no;

            total_grid_points += blk.blkSize;

            block_no++;
            p[level-lower_bound]++;
        }
    }

    // distribution representation requirement
    std::cout << "---Level distribution---" << "\033[0;33m" << std::endl;
    double ram_requirement = 0.0;
    for (int i=0; i<=upper_bound-lower_bound; i++) {
        ram_requirement += block_sizes[i+lower_bound]*p[i]; // this is calculated in MBs

        std::cout << i+lower_bound << ": ";
        std::cout << std::setw(3) << p[i]*100/numberOfLevels << "% " << std::string(p[i]*100/numberOfLevels, '*') << std::endl;
    }
    std::cout << "\033[0m" << "Total blocks: " << block_no << " | Total grid points: " << 1.0*total_grid_points/1000000 << "x10^6" << std::endl;
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
        std::cin >> sure;
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
    #pragma omp parallel for num_threads(20)
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


void data_generation_blockwise_and_bssn_var_wise_mixed(double mean, double std, unsigned int numberOfLevels, unsigned int lower_bound, unsigned int upper_bound, bool isRandom, Block * blkList, double ** var_in_array, double ** var_out_array, double ** var_in, double ** var_out)
{
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
    std::normal_distribution<double> distributionN(mean, std);
    std::uniform_int_distribution<int> distributionU(lower_bound, upper_bound);

    // Generating levels and block structure
    int total_grid_points = 0;
    int p[upper_bound-lower_bound+1];
    for (int i=0; i<upper_bound-lower_bound+1; i++) p[i] = 0;

    int level;
    int block_no = 0;

    // RAM ---
    unsigned long unzipSz=0;
    while (block_no<numberOfLevels){
        level = int(distributionN(generator)); // generate a number

        // distribution representation requirement
        if ((level>=lower_bound)&&(level<=upper_bound)) {
            Block & blk=blkList[block_no];
            blk=Block(0, 0, 0, 2*level, level, maxDepth);
            blk.block_no = block_no;
            // RAM ---
            blk.offset=unzipSz;
            unzipSz+=blk.blkSize;

            total_grid_points += blk.blkSize;

            block_no++;
            p[level-lower_bound]++;
        }
    }

    // RAM ---
    const unsigned long unzip_dof_cpu=unzipSz;

    // distribution representation requirement
    std::cout << "---Level distribution---" << "\033[0;33m" << std::endl;
    double ram_requirement = 0.0;
    for (int i=0; i<=upper_bound-lower_bound; i++) {
        ram_requirement += block_sizes[i+lower_bound]*p[i]; // this is calculated in MBs

        std::cout << i+lower_bound << ": ";
        std::cout << std::setw(3) << p[i]*100/numberOfLevels << "% " << std::string(p[i]*100/numberOfLevels, '*') << std::endl;
    }
    std::cout << "\033[0m" << "Total blocks: " << block_no << " | Total grid points: " << 1.0*total_grid_points/1000000 << "x10^6" << std::endl;
    std::cout << "Total RAM requiremnet for input: " << ram_requirement << "MB" << std::endl;
    std::cout << "Total RAM requiremnet for both input/output: " << ram_requirement*2 << "MB" << std::endl;
    // RAM ---
    std::cout << "Total RAM requiremnet for both GPU and CPU version: " << ram_requirement*2*2 + block_sizes[upper_bound]*210/BSSN_NUM_VARS << "MB" << std::endl; // block_sizes[upper_bound]*210/BSSN_NUM_VARS this is for allocating intermediate arrays in CPU

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
        std::cin >> sure;
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
    #pragma omp parallel for num_threads(20)
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
        
    #pragma omp parallel for num_threads(20)
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
