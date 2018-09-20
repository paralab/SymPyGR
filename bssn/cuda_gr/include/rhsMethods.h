
#include <iostream>
#include <stdio.h>

#include <random>
#include <map>
#include <sys/sysinfo.h> // Getting available RAM capacity
#include <iomanip> // setw formatting
#include "cuda_runtime.h"

#include "block.h"
#include "mergeSort.h"

#define CHECK_ERROR( err, msg ) if( err != cudaSuccess ) { std::cerr << "ERROR:" << cudaGetErrorName ( err ) << "  |  " << "ERROR DES: " << cudaGetErrorString( err ) << "  |  " << "User msg: " << msg << std::endl; exit( 0 ); }

int steamCountToLevel[5] = {3, 3, 3, 2, 2};
int hybrid_divider = 3;

void GPU_parallelized(unsigned int numberOfLevels, Block * blkList, unsigned int lower_bound, unsigned int upper_bound, 
double ** var_in_array, double ** var_out_array);

void GPU_parallelized_async_hybrid(unsigned int numberOfLevels, Block * blkList, unsigned int lower_bound, unsigned int upper_bound, 
double ** var_in_array, double ** var_out_array);

void GPU_pure_async_htod_dtoH_overlap(unsigned int numberOfLevels, Block * blkList, unsigned int lower_bound, unsigned int upper_bound, 
double ** var_in_array, double ** var_out_array);
