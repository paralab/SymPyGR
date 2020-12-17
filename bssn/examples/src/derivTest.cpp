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

    #include "bssnrhs_memalloc.h"

    auto t1=Time::now();
    for(unsigned int rr=0; rr < iter; rr++)
    {
        #include "bssnrhs_derivs.h"
        // deriv64_x(grad_0_alpha, alpha, hx, sz, bflag);
        // deriv64_x(grad_0_beta0, beta0, hx, sz, bflag);
        // deriv64_x(grad_0_beta1, beta1, hx, sz, bflag);
        // deriv64_x(grad_0_beta2, beta2, hx, sz, bflag);
        //deriv64_xx(grad2_0_0_alpha, alpha, hx, sz, bflag);
        // deriv64_y(grad_1_alpha, alpha, hy, sz, bflag);
        // deriv64_z(grad_2_alpha, alpha, hz, sz, bflag);
        
    }

    auto t2=Time::now();
    fsec fs = t2 - t1;
    std::cout<<"Time: "<<fs.count()<<std::endl;


    t1=Time::now();
    for(unsigned int rr=0; rr < iter; rr++)
    {
        #include "bssnrhs_derivs.h"
        // deriv64_x(grad_0_alpha, alpha, hx, sz, bflag);
        // deriv64_x(grad_0_beta0, beta0, hx, sz, bflag);
        // deriv64_x(grad_0_beta1, beta1, hx, sz, bflag);
        // deriv64_x(grad_0_beta2, beta2, hx, sz, bflag);
        //deriv64_xx(grad2_0_0_alpha, alpha, hx, sz, bflag);
        // deriv64_y(grad_1_alpha, alpha, hy, sz, bflag);
        // deriv64_z(grad_2_alpha, alpha, hz, sz, bflag);
        
    }

    t2=Time::now();
    fs = t2 - t1;
    std::cout<<"Time changed order: "<<fs.count()<<std::endl;

    // de-alloaction.
    {
        __mem_pool->purge();
    }



    // double* u  = new double[NN];
    // double* Dxu = new double[NN];
    // double* Dyu = new double[NN];
    // double* Dzu = new double[NN];
    // // double* u=nullptr;
    // // double* Du=nullptr;
    // // posix_memalign((void**)&u,128,sizeof(double)*NN);
    // // posix_memalign((void**)&Du,128,sizeof(double)*NN);

    // for(unsigned int i=0; i< NN; i++)
    //     u[i] = i*0.01;

    // auto t1=Time::now();
    // for(unsigned int rr=0; rr < iter; rr++)
    // {
    //     deriv64_x(Dxu,u,hx,sz,0);
    //     // deriv64_y(Dyu,u,hx,sz,0);
    //     // deriv64_z(Dzu,u,hx,sz,0);
    // }
    // auto t2=Time::now();
    // fsec fs = t2 - t1;
    // std::cout<<"Time: "<<fs.count()<<std::endl;

    // delete []  u;
    // delete [] Dxu;
    // delete [] Dyu;
    // delete [] Dzu;

    return 0;
    
}