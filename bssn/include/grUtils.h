//
// Created by milinda on 8/31/18.
//

#ifndef BSSN_GRUTILS_H
#define BSSN_GRUTILS_H

#include <iostream>
#include <math.h>

namespace bssn {
/**
 * @brief These variable indexes are based on the variables defined in rkBSSN.h
 * */
    enum VAR {
        U_ALPHA = 0,
        U_CHI,
        U_K,
        U_GT0,
        U_GT1,
        U_GT2,
        U_BETA0,
        U_BETA1,
        U_BETA2,
        U_B0,
        U_B1,
        U_B2,
        U_SYMGT0,
        U_SYMGT1,
        U_SYMGT2,
        U_SYMGT3,
        U_SYMGT4,
        U_SYMGT5,
        U_SYMAT0,
        U_SYMAT1,
        U_SYMAT2,
        U_SYMAT3,
        U_SYMAT4,
        U_SYMAT5
    };

//enum VAR_PSI4 {C_PSI4_REAL, C_PSI4_IMG};
    enum VAR_CONSTRAINT {
        C_HAM = 0, C_MOM0, C_MOM1, C_MOM2, C_PSI4_REAL, C_PSI4_IMG
    };


    static const unsigned int BSSN_LAMBDA[4]={0,1,2,1};
    static const double BSSN_ETA_POWER[2]={2.0,2.0};
    static const double BSSN_LAMBDA_F[2]={0.1,0.1};
    static const double ETA_CONST=2.0;
    static const double ETA_R0=20.0;
    static const double ETA_DAMPING_EXP=1.0;
    static const double KO_DISS_SIGMA=0.1;
    static const double ETA_DAMPING=1.0;

    void fake_initial_data(double x, double y, double z, double *u);

    const unsigned int BSSN_NUM_VARS=24;



}

#endif //BSSN_GRUTILS_H
