/*
 *
 *  Created on: March 12, 2018
 *      Author: Eminda, Akila
 */

#ifndef BSSNEQN_SOLVE_H
#define BSSNEQN_SOLVE_H

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define CHECK_ERROR( err, msg ) if( err != cudaSuccess ) { std::cerr << "ERROR:" << cudaGetErrorName ( err ) << "  |  " << "ERROR DES: " << cudaGetErrorString( err ) << "  |  " << "User msg: " << msg << std::endl; exit( 0 ); }

void calc_bssn_eqns(double * dev_var_in, double * dev_var_out, const unsigned int * sz, const double * pmin, double hz, double hy, double hx, cudaStream_t stream,
#include "list_of_offset_para.h"
, 
#include "list_of_para.h"
);

#endif /* BSSNEQN_SOLVE_H_ */
