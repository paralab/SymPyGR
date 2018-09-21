/**
 * Created on: Sep 21, 2018
 * 		Author: Akila, Eranga, Eminda, Ruwan
 **/

#include "cudaBSSN.h" 

int main (int argc, char** argv){
    /**
     *
     * parameters:
     * lowerLevel: block element 1d lower bound (int)
     * upperLevel: block element 1d upper bound (int) (blk_up>=blk_lb)
     * numberOfBlocks : total blocks going to process 
     * isTest(1/0; optional): default=0
     * isRandom(1/0; optional): default=0
     * mean(double optional) 
     * std(double optional) 
     *
     * */

    int lower_bound;
    int upper_bound;
    int numberOfBlocks;
    bool isTest = 0;
    bool isRandom = 0;
    double mean;
    double std = 2;

    if(argc>=4){
        lower_bound=atoi(argv[1]);
        upper_bound=atoi(argv[2]);
        numberOfBlocks=atoi(argv[3]);

        mean = lower_bound + (upper_bound-lower_bound)/2;

        if (argc>=5) isTest=atoi(argv[4]);
        if (argc>=6) isRandom=atoi(argv[5]);
        if (argc>=7) mean=atoi(argv[6]);
        if (argc>=8) std=atoi(argv[7]);

    }else{
        std::cerr << "Correct usage: " << argv[0] << " lowerLevel(int) upperLevel(int) numberOfBlocks(int) isTest(1/0; optional) isRandom(1/0; optional) mean(double optional) std(double optional)" << std::endl;
        exit(0);
    }

    Block * blkList = new Block[numberOfBlocks];
    double ** var_in_array = new double*[numberOfBlocks];
    double ** var_out_array = new double*[numberOfBlocks];

    double ** var_in;
    double ** var_out;

    if (isTest){
        var_in = new double*[BSSN_NUM_VARS];
        var_out = new double*[BSSN_NUM_VARS];
        data_generation_blockwise_and_bssn_var_wise_mixed(mean, std, numberOfBlocks, lower_bound, upper_bound, isRandom, blkList, var_in_array, var_out_array, var_in, var_out);
    }else{
        data_generation_blockwise_mixed(mean, std, numberOfBlocks, lower_bound, upper_bound, isRandom, blkList, var_in_array, var_out_array);
    }

    #if parallelized
        GPU_parallelized(numberOfBlocks, blkList, lower_bound, upper_bound, var_in_array, var_out_array);
    #endif
    #if async
        GPU_async(numberOfBlocks, blkList, lower_bound, upper_bound, var_in_array, var_out_array);
    #endif
    #if hybrid
        GPU_hybrid(numberOfBlocks, blkList, lower_bound, upper_bound, var_in_array, var_out_array);
    #endif


    if (isTest){
        CPU_sequence(numberOfBlocks, blkList, var_in, var_out);

        // Verify outputs
        for (int blk=0; blk<numberOfBlocks; blk++){
            for(int bssn_var=0; bssn_var<BSSN_NUM_VARS; bssn_var++){
                int sizeofBlock = blkList[blk].blkSize;
                for (int pointInd=0; pointInd<sizeofBlock; pointInd++){
                    double diff = var_out_array[blkList[blk].block_no][bssn_var*sizeofBlock+pointInd] - var_out[bssn_var][blkList[blk].offset+pointInd];
                    if (fabs(diff)>1e-5){
                        char separator    = ' ';
                        const int nameWidth     = 6;
                        const int numWidth      = 20;
                        const int NUM_DIGITS = 10;

                        std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) << "GPU: ";
                        std::cout <<std::setprecision(NUM_DIGITS)<< std::left << std::setw(numWidth) << std::setfill(separator)  << var_out_array[blkList[blk].block_no][bssn_var*sizeofBlock+pointInd];

                        std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) << "CPU: ";
                        std::cout <<std::setprecision(NUM_DIGITS)<< std::left << std::setw(numWidth) << std::setfill(separator)  << var_out[bssn_var][blkList[blk].offset+pointInd];

                        std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) << "DIFF: ";
                        std::cout <<std::setprecision(NUM_DIGITS)<< std::left << std::setw(numWidth) << std::setfill(separator)  << diff << std::endl;

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
    for (int blk=0; blk<numberOfBlocks; blk++){
        CHECK_ERROR(cudaFreeHost(var_in_array[blk]), "free host memory");
        CHECK_ERROR(cudaFreeHost(var_out_array[blk]), "free host memory");
    }
    return 0;
}