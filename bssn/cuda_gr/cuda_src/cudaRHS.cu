/**
 * Created on: Sep 21, 2018
 * 		Author: Akila, Eranga, Eminda, Ruwan
 **/
 
#include "cudaRHS.cuh"

enum VAR_CU {U_ALPHA=0,U_CHI,U_K,U_GT0,U_GT1,U_GT2,U_BETA0,U_BETA1,U_BETA2,U_B0,U_B1,U_B2,U_SYMGT0,U_SYMGT1,U_SYMGT2,U_SYMGT3,U_SYMGT4,U_SYMGT5,U_SYMAT0,U_SYMAT1,U_SYMAT2,U_SYMAT3,U_SYMAT4,U_SYMAT5};

void cuda_bssnrhs(double * dev_var_out, double * dev_var_in, const unsigned int unzip_dof, 
const double * pmin, const double * pmax, const unsigned int * sz, 
const unsigned int& bflag, cudaStream_t stream,
#include "list_of_para.h"
)
{ 
    CHECK_ERROR(cudaMemsetAsync(dev_var_out, 0, 24*unzip_dof*sizeof(double), stream), "output array cleaning call"); // Clean output array

    int alphaInt = (VAR_CU::U_ALPHA) * unzip_dof;
    int chiInt = (VAR_CU::U_CHI) * unzip_dof;
    int KInt = (VAR_CU::U_K) * unzip_dof;
    int gt0Int = (VAR_CU::U_SYMGT0) * unzip_dof;
    int gt1Int = (VAR_CU::U_SYMGT1) * unzip_dof;
    int gt2Int =  (VAR_CU::U_SYMGT2) * unzip_dof;
    int gt3Int = (VAR_CU::U_SYMGT3) * unzip_dof;
    int gt4Int = (VAR_CU::U_SYMGT4) * unzip_dof;
    int gt5Int = (VAR_CU::U_SYMGT5) * unzip_dof;
    int beta0Int = (VAR_CU::U_BETA0) * unzip_dof;
    int beta1Int = (VAR_CU::U_BETA1) * unzip_dof;
    int beta2Int = (VAR_CU::U_BETA2) * unzip_dof;
    int At0Int = (VAR_CU::U_SYMAT0) * unzip_dof;
    int At1Int = (VAR_CU::U_SYMAT1) * unzip_dof;
    int At2Int = (VAR_CU::U_SYMAT2) * unzip_dof;
    int At3Int = (VAR_CU::U_SYMAT3) * unzip_dof;
    int At4Int = (VAR_CU::U_SYMAT4) * unzip_dof;
    int At5Int = (VAR_CU::U_SYMAT5) * unzip_dof;
    int Gt0Int = (VAR_CU::U_GT0) * unzip_dof;
    int Gt1Int = (VAR_CU::U_GT1) * unzip_dof;
    int Gt2Int = (VAR_CU::U_GT2) * unzip_dof;
    int B0Int = (VAR_CU::U_B0) * unzip_dof;
    int B1Int = (VAR_CU::U_B1) * unzip_dof;
    int B2Int = (VAR_CU::U_B2) * unzip_dof;

    double hx = (pmax[0] - pmin[0]) / (sz[0] - 1);
    double hy = (pmax[1] - pmin[1]) / (sz[1] - 1);
    double hz = (pmax[2] - pmin[2]) / (sz[2] - 1);

    calc_deriv_wrapper(dev_var_out, dev_var_in, hx, hy, hz, sz, bflag, stream,
        #include "list_of_offset_args.h"
        ,
        #include "list_of_args.h"
    );

    calc_bssn_eqns(dev_var_in, dev_var_out, sz, pmin, hz, hy, hx, stream,
    #include "list_of_offset_args.h"
    ,
    #include "list_of_args.h"
    );

    if (bflag!=0) {
        bssn_bcs(dev_var_out, dev_var_in, alphaInt, grad_0_alpha, grad_1_alpha, grad_2_alpha,
            pmin, pmax, 1.0, 1.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, chiInt, grad_0_chi, grad_1_chi, grad_2_chi,
            pmin, pmax, 1.0, 1.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, KInt, grad_0_K, grad_1_K, grad_2_K,
            pmin, pmax, 1.0, 0.0, sz, bflag, stream);

        bssn_bcs(dev_var_out, dev_var_in, beta0Int, grad_0_beta0, grad_1_beta0, grad_2_beta0,
            pmin, pmax, 1.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, beta1Int, grad_0_beta1, grad_1_beta1, grad_2_beta1,
            pmin, pmax, 1.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, beta2Int, grad_0_beta2, grad_1_beta2, grad_2_beta2,
            pmin, pmax, 1.0, 0.0, sz, bflag, stream);

        bssn_bcs(dev_var_out, dev_var_in, Gt0Int, grad_0_Gt0, grad_1_Gt0, grad_2_Gt0,
            pmin, pmax, 2.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, Gt1Int, grad_0_Gt1, grad_1_Gt1, grad_2_Gt1,
            pmin, pmax, 2.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, Gt2Int, grad_0_Gt2, grad_1_Gt2, grad_2_Gt2,
            pmin, pmax, 2.0, 0.0, sz, bflag, stream);

        bssn_bcs(dev_var_out, dev_var_in, B0Int, grad_0_B0, grad_1_B0, grad_2_B0,
            pmin, pmax, 1.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, B1Int, grad_0_B1, grad_1_B1, grad_2_B1,
            pmin, pmax, 1.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, B2Int, grad_0_B2, grad_1_B2, grad_2_B2,
            pmin, pmax, 1.0, 0.0, sz, bflag, stream);

        bssn_bcs(dev_var_out, dev_var_in, At0Int, grad_0_At0, grad_1_At0, grad_2_At0,
            pmin, pmax, 2.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, At1Int, grad_0_At1, grad_1_At1, grad_2_At1,
            pmin, pmax, 2.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, At2Int, grad_0_At2, grad_1_At2, grad_2_At2,
            pmin, pmax, 2.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, At3Int, grad_0_At3, grad_1_At3, grad_2_At3,
            pmin, pmax, 2.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, At4Int, grad_0_At4, grad_1_At4, grad_2_At4,
            pmin, pmax, 2.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, At5Int, grad_0_At5, grad_1_At5, grad_2_At5,
            pmin, pmax, 2.0, 0.0, sz, bflag, stream); 

        bssn_bcs(dev_var_out, dev_var_in, gt0Int, grad_0_gt0, grad_1_gt0, grad_2_gt0,
            pmin, pmax, 1.0, 1.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, gt1Int, grad_0_gt1, grad_1_gt1, grad_2_gt1,
            pmin, pmax, 1.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, gt2Int, grad_0_gt2, grad_1_gt2, grad_2_gt2,
            pmin, pmax, 1.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, gt3Int, grad_0_gt3, grad_1_gt3, grad_2_gt3,
            pmin, pmax, 1.0, 1.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, gt4Int, grad_0_gt4, grad_1_gt4, grad_2_gt4,
            pmin, pmax, 1.0, 0.0, sz, bflag, stream);
        bssn_bcs(dev_var_out, dev_var_in, gt5Int, grad_0_gt5, grad_1_gt5, grad_2_gt5,
            pmin, pmax, 1.0, 1.0, sz, bflag, stream); 
    }

    calc_ko_deriv_wrapper(dev_var_out, dev_var_in, hx, hy, hz, sz, bflag, stream,
        #include "list_of_offset_args.h"
        ,
        #include "list_of_args.h"
    );

    get_output(dev_var_out, sz, stream,
        #include "list_of_offset_args.h"
        ,
        #include "list_of_args.h"
    );
    return;
}
