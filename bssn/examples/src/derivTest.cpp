/**
 * @brief Evaluate the performance of the deriv routines. 
 * 
 */

#include<iostream>
#include "derivs.h"
#include <immintrin.h>
#include <chrono>
#include "memory_pool.h"
#include "grDef.h"
#include "grUtils.h"
#include "parameters.h"


typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::milliseconds ms;
typedef std::chrono::duration<float> fsec;

int main(int argc, char **argv)
{

    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " blk_sz iter" << std::endl;
        return 0;
    }

    const unsigned int blk_sz = atoi(argv[1]);
    const unsigned int iter   = atoi(argv[2]);
    const unsigned int sz[3]={blk_sz,blk_sz,blk_sz};

    const unsigned int nx=sz[0];
    const unsigned int ny=sz[1];
    const unsigned int nz=sz[2];
    const unsigned int NN =nx*ny*nz;
    const unsigned int bflag=0;
    
    const double hx=0.01;
    const double hy=0.01;
    const double hz=0.01;
    double init_var[bssn::BSSN_NUM_VARS];

    mem::memory_pool<double>* __mem_pool = new mem::memory_pool<double>(0,16);
    const unsigned int n=NN;
    double *alpha = __mem_pool->allocate(n);
    double *chi   = __mem_pool->allocate(n);
    double *K     = __mem_pool->allocate(n);
    double *gt0   = __mem_pool->allocate(n);
    double *gt1   = __mem_pool->allocate(n);
    double *gt2   = __mem_pool->allocate(n);
    double *gt3   = __mem_pool->allocate(n);
    double *gt4   = __mem_pool->allocate(n);
    double *gt5   = __mem_pool->allocate(n);
    double *beta0 = __mem_pool->allocate(n);
    double *beta1 = __mem_pool->allocate(n);
    double *beta2 = __mem_pool->allocate(n);
    double *At0   = __mem_pool->allocate(n);
    double *At1   = __mem_pool->allocate(n);
    double *At2   = __mem_pool->allocate(n);
    double *At3   = __mem_pool->allocate(n);
    double *At4   = __mem_pool->allocate(n);
    double *At5   = __mem_pool->allocate(n);
    double *Gt0   = __mem_pool->allocate(n);
    double *Gt1   = __mem_pool->allocate(n);
    double *Gt2   = __mem_pool->allocate(n);
    double *B0    = __mem_pool->allocate(n);
    double *B1    = __mem_pool->allocate(n);
    double *B2    = __mem_pool->allocate(n);


    double *a_rhs    = __mem_pool->allocate(n);
    double *chi_rhs  = __mem_pool->allocate(n);
    double *K_rhs    = __mem_pool->allocate(n);
    double *gt_rhs00 = __mem_pool->allocate(n);
    double *gt_rhs01 = __mem_pool->allocate(n);
    double *gt_rhs02 = __mem_pool->allocate(n);
    double *gt_rhs11 = __mem_pool->allocate(n);
    double *gt_rhs12 = __mem_pool->allocate(n);
    double *gt_rhs22 = __mem_pool->allocate(n);
    double *b_rhs0   = __mem_pool->allocate(n);
    double *b_rhs1   = __mem_pool->allocate(n);
    double *b_rhs2   = __mem_pool->allocate(n);
    double *At_rhs00 = __mem_pool->allocate(n);
    double *At_rhs01 = __mem_pool->allocate(n);
    double *At_rhs02 = __mem_pool->allocate(n);
    double *At_rhs11 = __mem_pool->allocate(n);
    double *At_rhs12 = __mem_pool->allocate(n);
    double *At_rhs22 = __mem_pool->allocate(n);
    double *Gt_rhs0  = __mem_pool->allocate(n);
    double *Gt_rhs1  = __mem_pool->allocate(n);
    double *Gt_rhs2  = __mem_pool->allocate(n);
    double *B_rhs0   = __mem_pool->allocate(n);
    double *B_rhs1   = __mem_pool->allocate(n);
    double *B_rhs2   = __mem_pool->allocate(n);

    const unsigned int lambda[4] = {bssn::BSSN_LAMBDA[0], bssn::BSSN_LAMBDA[1],
                                    bssn::BSSN_LAMBDA[2], bssn::BSSN_LAMBDA[3]
                                   };
    const double lambda_f[2] = {bssn::BSSN_LAMBDA_F[0], bssn::BSSN_LAMBDA_F[1]};

    for(unsigned int k=0; k < nz; k++)
    {
        const double z  =  k * hz;
        for(unsigned int j=0; j < ny; j++)
        {
            const double y  =  j * hy;
            for(unsigned int i=0; i < nx; i++)
            {
                const double x  =  i * hx;
                const unsigned int pp = k*ny*nx + j*nx + i;
                bssn::fake_initial_data(x,y,z,init_var);
                alpha[pp]   = init_var[bssn::VAR::U_ALPHA];
                
                beta0[pp]   = init_var[bssn::VAR::U_BETA0];
                beta1[pp]   = init_var[bssn::VAR::U_BETA1];
                beta2[pp]   = init_var[bssn::VAR::U_BETA2];

                B0[pp]      = init_var[bssn::VAR::U_B0];
                B1[pp]      = init_var[bssn::VAR::U_B1];
                B2[pp]      = init_var[bssn::VAR::U_B2];
                
                Gt0[pp]     = init_var[bssn::VAR::U_GT0];
                Gt1[pp]     = init_var[bssn::VAR::U_GT1];
                Gt2[pp]     = init_var[bssn::VAR::U_GT2];
                
                chi[pp]     = init_var[bssn::VAR::U_CHI];
                K[pp]       = init_var[bssn::VAR::U_K];
                
                gt0[pp]     = init_var[bssn::VAR::U_SYMGT0];
                gt1[pp]     = init_var[bssn::VAR::U_SYMGT3];
                gt2[pp]     = init_var[bssn::VAR::U_SYMGT5];
                gt3[pp]     = init_var[bssn::VAR::U_SYMGT1];
                gt4[pp]     = init_var[bssn::VAR::U_SYMGT2];
                gt5[pp]     = init_var[bssn::VAR::U_SYMGT4];

                At0[pp]     = init_var[bssn::VAR::U_SYMAT0];
                At1[pp]     = init_var[bssn::VAR::U_SYMAT1];
                At2[pp]     = init_var[bssn::VAR::U_SYMAT2];
                At3[pp]     = init_var[bssn::VAR::U_SYMAT3];
                At4[pp]     = init_var[bssn::VAR::U_SYMAT4];
                At5[pp]     = init_var[bssn::VAR::U_SYMAT5];
                

            }
        }
    }

    // allocate deriv vars. 
    #include "bssnrhs_memalloc.h"
    #include "bssnrhs_memalloc_adv.h"
    
    // compute deriv vars. 
    //#include "bssnrhs_derivs.h"
    //#include "bssnrhs_derivs_adv.h"

    // auto t1=Time::now();
    // for(unsigned int rr=0; rr < iter; rr++)
    // {
    //     //#include
    //     //#include "bssnrhs_derivs.h"
    //     // deriv64_x(grad_0_alpha, alpha, hx, sz, bflag);
    //     // deriv64_x(grad_0_beta0, beta0, hx, sz, bflag);
    //     // deriv64_x(grad_0_beta1, beta1, hx, sz, bflag);
    //     // deriv64_x(grad_0_beta2, beta2, hx, sz, bflag);
    //     //deriv64_xx(grad2_0_0_alpha, alpha, hx, sz, bflag);
    //     // deriv64_y(grad_1_alpha, alpha, hy, sz, bflag);
    //     // deriv64_z(grad_2_alpha, alpha, hz, sz, bflag);
        
    // }

    // auto t2=Time::now();
    // fsec fs = t2 - t1;
    // std::cout<<"Time: "<<fs.count()<<std::endl;

    //std::cout<<"iter:"<<iter<<std::endl;

    auto t1=Time::now();
    for(unsigned int rr=0; rr < iter; rr++)
    {

        for(unsigned int k=3; k < nz-3; k++)
            for(unsigned int j=3; j < ny-3; j++)
                for(unsigned int i=3; i < nx-3; i++)
                {
                    const double x  =  i * hx;
                    const double y  =  j * hy;
                    const double z  =  k * hz;

                    const double r_coord =x*x + y*y + z*z;
                    const double eta=2.0;
                    const double sigma=0.4;

                    
                    const unsigned int pp = k*ny*nx + j*nx + i;
                    #include "../../src/bssneqs_eta_const_standard_gauge.cpp"

                    // KO dissipation. 
                    a_rhs[pp]    += sigma * (grad_0_alpha[pp] + grad_1_alpha[pp] + grad_2_alpha[pp]);
                    b_rhs0[pp]   += sigma * (grad_0_beta0[pp] + grad_1_beta0[pp] + grad_2_beta0[pp]);
                    b_rhs1[pp]   += sigma * (grad_0_beta1[pp] + grad_1_beta1[pp] + grad_2_beta1[pp]);
                    b_rhs2[pp]   += sigma * (grad_0_beta2[pp] + grad_1_beta2[pp] + grad_2_beta2[pp]);
                    gt_rhs00[pp] += sigma * (grad_0_gt0[pp] + grad_1_gt0[pp] + grad_2_gt0[pp]);
                    gt_rhs01[pp] += sigma * (grad_0_gt1[pp] + grad_1_gt1[pp] + grad_2_gt1[pp]);
                    gt_rhs02[pp] += sigma * (grad_0_gt2[pp] + grad_1_gt2[pp] + grad_2_gt2[pp]);
                    gt_rhs11[pp] += sigma * (grad_0_gt3[pp] + grad_1_gt3[pp] + grad_2_gt3[pp]);
                    gt_rhs12[pp] += sigma * (grad_0_gt4[pp] + grad_1_gt4[pp] + grad_2_gt4[pp]);
                    gt_rhs22[pp] += sigma * (grad_0_gt5[pp] + grad_1_gt5[pp] + grad_2_gt5[pp]);
                    chi_rhs[pp]  += sigma * (grad_0_chi[pp] + grad_1_chi[pp] + grad_2_chi[pp]);
                    At_rhs00[pp] += sigma * (grad_0_At0[pp] + grad_1_At0[pp] + grad_2_At0[pp]);
                    At_rhs01[pp] += sigma * (grad_0_At1[pp] + grad_1_At1[pp] + grad_2_At1[pp]);
                    At_rhs02[pp] += sigma * (grad_0_At2[pp] + grad_1_At2[pp] + grad_2_At2[pp]);
                    At_rhs11[pp] += sigma * (grad_0_At3[pp] + grad_1_At3[pp] + grad_2_At3[pp]);
                    At_rhs12[pp] += sigma * (grad_0_At4[pp] + grad_1_At4[pp] + grad_2_At4[pp]);
                    At_rhs22[pp] += sigma * (grad_0_At5[pp] + grad_1_At5[pp] + grad_2_At5[pp]);
                    K_rhs[pp]    += sigma * (grad_0_K[pp] + grad_1_K[pp] + grad_2_K[pp]);
                    Gt_rhs0[pp]  += sigma * (grad_0_Gt0[pp] + grad_1_Gt0[pp] + grad_2_Gt0[pp]);
                    Gt_rhs1[pp]  += sigma * (grad_0_Gt1[pp] + grad_1_Gt1[pp] + grad_2_Gt1[pp]);
                    Gt_rhs2[pp]  += sigma * (grad_0_Gt2[pp] + grad_1_Gt2[pp] + grad_2_Gt2[pp]);
                    B_rhs0[pp]   += sigma * (grad_0_B0[pp] + grad_1_B0[pp] + grad_2_B0[pp]);
                    B_rhs1[pp]   += sigma * (grad_0_B1[pp] + grad_1_B1[pp] + grad_2_B1[pp]);
                    B_rhs2[pp]   += sigma * (grad_0_B2[pp] + grad_1_B2[pp] + grad_2_B2[pp]);

                }
    }

    auto t2=Time::now();
    auto fs = t2 - t1;
    std::cout<<"Time changed order: "<<fs.count()<<std::endl;


    for(unsigned int k=3; k < nz-3; k+=10)
            for(unsigned int j=3; j < ny-3; j+=10)
                for(unsigned int i=3; i < nx-3; i+=10)
                {
                    const unsigned int pp = k*ny*nx + j*nx + i;
                    std::cout<<"rhs: "<<a_rhs[pp]<<std::endl;
                    std::cout<<"rhs: "<<b_rhs0[pp]<<std::endl;
                    std::cout<<"rhs: "<<b_rhs1[pp]<<std::endl;
                    std::cout<<"rhs: "<<b_rhs2[pp]<<std::endl;
                    std::cout<<"rhs: "<<gt_rhs00[pp]<<std::endl;
                    std::cout<<"rhs: "<<gt_rhs01[pp]<<std::endl;
                    std::cout<<"rhs: "<<gt_rhs02[pp]<<std::endl;
                    std::cout<<"rhs: "<<gt_rhs11[pp]<<std::endl;
                    std::cout<<"rhs: "<<gt_rhs12[pp]<<std::endl;
                    std::cout<<"rhs: "<<gt_rhs22[pp]<<std::endl;
                    std::cout<<"rhs: "<<chi_rhs[pp]<<std::endl;
                    std::cout<<"rhs: "<<At_rhs00[pp]<<std::endl;
                    std::cout<<"rhs: "<<At_rhs01[pp]<<std::endl;
                    std::cout<<"rhs: "<<At_rhs02[pp]<<std::endl;
                    std::cout<<"rhs: "<<At_rhs11[pp]<<std::endl;
                    std::cout<<"rhs: "<<At_rhs12[pp]<<std::endl;
                    std::cout<<"rhs: "<<K_rhs[pp]<<std::endl;
                    std::cout<<"rhs: "<<Gt_rhs0[pp]<<std::endl;
                    std::cout<<"rhs: "<<Gt_rhs1[pp]<<std::endl;
                    std::cout<<"rhs: "<<Gt_rhs2[pp]<<std::endl;
                    std::cout<<"rhs: "<<B_rhs0[pp]<<std::endl;
                    std::cout<<"rhs: "<<B_rhs1[pp]<<std::endl;
                    std::cout<<"rhs: "<<B_rhs2[pp]<<std::endl;
                    

                    
                }
                    

    // de-alloaction.
    {
        __mem_pool->purge();
    }



    

    return 0;
    
}