/**
 * Created on: Sep 21, 2018
 * 		Author: Akila, Eranga, Eminda, Ruwan
 **/

#ifndef RHS_METHODS_H
#define RHS_METHODS_H

#include <iostream>
#include <stdio.h>

#include <random>
#include <map>
#include <sys/sysinfo.h> // Getting available RAM capacity
#include <iomanip> // setw formatting
#include "cuda_runtime.h"

#include "block.h"
#include "mergeSort.h"
#include "cudaRHS.cuh"
#include "profileGPU.h"

#define CHECK_ERROR( err, msg ) if( err != cudaSuccess ) { std::cerr << "\033[1;31mERROR:" << cudaGetErrorName ( err ) << "  |  " << "ERROR DES: " << cudaGetErrorString( err ) << "  |  " << "User msg: " << msg << "\033[0m" << std::endl; exit( 0 ); }

void GPU_parallelized(unsigned int numberOfLevels, Block * blkList, unsigned int lower_bound, unsigned int upper_bound, 
double ** var_in_array, double ** var_out_array);

void GPU_hybrid(unsigned int numberOfLevels, Block * blkList, unsigned int lower_bound, unsigned int upper_bound, 
double ** var_in_array, double ** var_out_array);

void GPU_async(unsigned int numberOfLevels, Block * blkList, unsigned int lower_bound, unsigned int upper_bound, 
double ** var_in_array, double ** var_out_array);

#ifdef ENABLE_CUDA_TEST
void CPU_sequential(unsigned int numberOfLevels, Block * blkList, double ** var_in, double ** var_out);
#endif

#endif