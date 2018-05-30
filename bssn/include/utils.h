//
// Created by milinda on 5/3/18.
//

#ifndef BSSN_UTILS_H
#define BSSN_UTILS_H

#include <cmath>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <time.h>


enum VAR {U_ALPHA=0,U_CHI,U_K,U_GT0,U_GT1,U_GT2,U_BETA0,U_BETA1,U_BETA2,U_B0,U_B1,U_B2,U_SYMGT0,U_SYMGT1,U_SYMGT2,U_SYMGT3,U_SYMGT4,U_SYMGT5,U_SYMAT0,U_SYMAT1,U_SYMAT2,U_SYMAT3,U_SYMAT4,U_SYMAT5};



/**
 * @brief : intialize the data for bssn computation. (not correct initial data)
 * @param[in,out] u  : pointer to the begin of the unzip block.
 * @param[in] xi : block left most corner coordinates
 * @param[in] shp; size of the block dimentions.
 * */
void initial_data(double *u, const double *xi);




#endif //BSSN_UTILS_H
