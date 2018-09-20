#ifndef DATA_GENERATION_H
#define DATA_GENERATION_H

#include <random>
#include <map>
#include <sys/sysinfo.h> // Getting available RAM capacity
#include <iomanip> // setw formatting
#include "cuda_runtime.h"

#include "utils.h" // Block initialize with fake data code
#include "block.h"

#define CHECK_ERROR( err, msg ) if( err != cudaSuccess ) { std::cerr << "ERROR:" << cudaGetErrorName ( err ) << "  |  " << "ERROR DES: " << cudaGetErrorString( err ) << "  |  " << "User msg: " << msg << std::endl; exit( 0 ); }

void data_generation_blockwise_mixed(bool isNormal, double mean, double std, unsigned int numberOfLevels, 
unsigned int lower_bound, unsigned int upper_bound, bool isRandom, Block * blkList, double ** var_in_array, double ** var_out_array);


void data_generation_blockwise_and_bssn_var_wise_mixed(bool isNormal, double mean, double std, unsigned int numberOfLevels, 
unsigned int lower_bound, unsigned int upper_bound, bool isRandom, Block * blkList, double ** var_in_array, double ** var_out_array, 
double ** var_in, double ** var_out);

#endif