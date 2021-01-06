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
#include <mpi.h>


typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::milliseconds ms;
typedef std::chrono::duration<float> fsec;

int main(int argc, char **argv)
{

    unsigned int block_size = 16;
    unsigned int iter       = 10;   

    bool vectorize = false;
    bool stages = false;

    int i = 0;
    while(i < argc)
    {
        std::string arg = argv[i];

        if (arg == "--help" || arg == "-h")
        {
            std::cout << "Usage       : " << argv[0] << " --block 128 --iter 20 --vectorize --stages" << std::endl;
            std::cout << "--block     : block size. Defaults to " << block_size << std::endl;
            std::cout << "--iter      : Number of times to repeat RHS compute. Defaults to " << iter << std::endl;
            std::cout << "--vectorize : If present, use single vectorized loop equivalent to triple nested loop for RHS compute" << std::endl;
            std::cout << "--stages    : If present" << std::endl;

            return 0;
        }

        if(arg == "--vectorize")
            vectorize = true;
        if(arg == "--stages")
            stages = true;

        if(i+1 < argc && arg == "--block")
        {
            i += 1;
            block_size = std::stoi(argv[i]);
        }

        if(i+1 < argc && arg == "--iter")
        {
            i += 1;
            iter = std::stoi(argv[i]);
        }
  
        i++;       
    }

    MPI_Init(&argc,&argv);
    int rank, npes;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm,&npes);
    MPI_Comm_rank(comm,&rank);

    const unsigned int sz[3]={block_size,block_size,block_size};

    const unsigned int nx=sz[0];
    const unsigned int ny=sz[1];
    const unsigned int nz=sz[2];
    const unsigned int NN =nx*ny*nz;
    const unsigned int bflag=0;
    
    const double hx=0.01;
    const double hy=0.01;
    const double hz=0.01;

    double init_var[bssn::BSSN_NUM_VARS];
    

    mem::memory_pool<double>* __mem_pool = new mem::memory_pool<double>(0,64);
    const unsigned int n=NN;

    double * unzipIn   = new double[NN*bssn::BSSN_NUM_VARS];
    double * unzipOut  = new double[NN*bssn::BSSN_NUM_VARS];

    for (unsigned int pp =0; pp < NN * bssn::BSSN_NUM_VARS ; pp++)
    {
        unzipIn[pp]=0.0;
        unzipOut[pp]=0.0;
    }

    double *alpha = unzipIn + (bssn::VAR::U_ALPHA*NN);    //__mem_pool->allocate(n);
    double *chi   = unzipIn + (bssn::VAR::U_CHI*NN);      //__mem_pool->allocate(n);
    double *K     = unzipIn + (bssn::VAR::U_K*NN);        //__mem_pool->allocate(n);
    double *gt0   = unzipIn + (bssn::VAR::U_SYMGT0*NN);   //__mem_pool->allocate(n);
    double *gt1   = unzipIn + (bssn::VAR::U_SYMGT1*NN);   //__mem_pool->allocate(n);
    double *gt2   = unzipIn + (bssn::VAR::U_SYMGT2*NN);   //__mem_pool->allocate(n);
    double *gt3   = unzipIn + (bssn::VAR::U_SYMGT3*NN);   //__mem_pool->allocate(n);
    double *gt4   = unzipIn + (bssn::VAR::U_SYMGT4*NN);   //__mem_pool->allocate(n);
    double *gt5   = unzipIn + (bssn::VAR::U_SYMGT5*NN);   //__mem_pool->allocate(n);
    double *beta0 = unzipIn + (bssn::VAR::U_BETA0*NN);    //__mem_pool->allocate(n);
    double *beta1 = unzipIn + (bssn::VAR::U_BETA1*NN);    //__mem_pool->allocate(n);
    double *beta2 = unzipIn + (bssn::VAR::U_BETA2*NN);    //__mem_pool->allocate(n);
    double *At0   = unzipIn + (bssn::VAR::U_SYMAT0*NN);   //__mem_pool->allocate(n);
    double *At1   = unzipIn + (bssn::VAR::U_SYMAT1*NN);   //__mem_pool->allocate(n);
    double *At2   = unzipIn + (bssn::VAR::U_SYMAT2*NN);   //__mem_pool->allocate(n);
    double *At3   = unzipIn + (bssn::VAR::U_SYMAT3*NN);   //__mem_pool->allocate(n);
    double *At4   = unzipIn + (bssn::VAR::U_SYMAT4*NN);   //__mem_pool->allocate(n);
    double *At5   = unzipIn + (bssn::VAR::U_SYMAT5*NN);   //__mem_pool->allocate(n);
    double *Gt0   = unzipIn + (bssn::VAR::U_GT0*NN);      //__mem_pool->allocate(n);
    double *Gt1   = unzipIn + (bssn::VAR::U_GT1*NN);      //__mem_pool->allocate(n);
    double *Gt2   = unzipIn + (bssn::VAR::U_GT2*NN);      //__mem_pool->allocate(n);
    double *B0    = unzipIn + (bssn::VAR::U_B0*NN);       //__mem_pool->allocate(n);
    double *B1    = unzipIn + (bssn::VAR::U_B1*NN);       //__mem_pool->allocate(n);
    double *B2    = unzipIn + (bssn::VAR::U_B2*NN);       //__mem_pool->allocate(n);
  
    const unsigned int lambda[4] = {bssn::BSSN_LAMBDA[0], bssn::BSSN_LAMBDA[1],
                                    bssn::BSSN_LAMBDA[2], bssn::BSSN_LAMBDA[3]
                                   };
    const double lambda_f[2] = {bssn::BSSN_LAMBDA_F[0], bssn::BSSN_LAMBDA_F[1]};
    
    const double LAMBDA_F0 = lambda_f[0];
    const double LAMBDA_F1 = lambda_f[1];
    
    const double LAMBDA0 = lambda[0];
    const double LAMBDA1 = lambda[1];
    const double LAMBDA2 = lambda[2];
    const double LAMBDA3 = lambda[3];

    for(unsigned int k=0; k < nz; k++)
    {
        for(unsigned int j=0; j < ny; j++)
        {
            for(unsigned int i=0; i < nx; i++)
            {
                const double x  =  i * hx;
                const double z  =  k * hz;
                const double y  =  j * hy;
                    
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
    #include "bssnrhs_derivs.h"
    #include "bssnrhs_derivs_adv.h"

    double *a_rhs    =  unzipOut + (bssn::VAR::U_ALPHA*NN); 
    double *chi_rhs  =  unzipOut + (bssn::VAR::U_CHI*NN);   
    double *K_rhs    =  unzipOut + (bssn::VAR::U_K*NN);     
    double *gt_rhs00 =  unzipOut + (bssn::VAR::U_SYMGT0*NN);
    double *gt_rhs01 =  unzipOut + (bssn::VAR::U_SYMGT1*NN);
    double *gt_rhs02 =  unzipOut + (bssn::VAR::U_SYMGT2*NN);
    double *gt_rhs11 =  unzipOut + (bssn::VAR::U_SYMGT3*NN);
    double *gt_rhs12 =  unzipOut + (bssn::VAR::U_SYMGT4*NN);
    double *gt_rhs22 =  unzipOut + (bssn::VAR::U_SYMGT5*NN);
    double *b_rhs0   =  unzipOut + (bssn::VAR::U_BETA0*NN); 
    double *b_rhs1   =  unzipOut + (bssn::VAR::U_BETA1*NN); 
    double *b_rhs2   =  unzipOut + (bssn::VAR::U_BETA2*NN); 
    double *At_rhs00 =  unzipOut + (bssn::VAR::U_SYMAT0*NN);
    double *At_rhs01 =  unzipOut + (bssn::VAR::U_SYMAT1*NN);
    double *At_rhs02 =  unzipOut + (bssn::VAR::U_SYMAT2*NN);
    double *At_rhs11 =  unzipOut + (bssn::VAR::U_SYMAT3*NN);
    double *At_rhs12 =  unzipOut + (bssn::VAR::U_SYMAT4*NN);
    double *At_rhs22 =  unzipOut + (bssn::VAR::U_SYMAT5*NN);
    double *Gt_rhs0  =  unzipOut + (bssn::VAR::U_GT0*NN);   
    double *Gt_rhs1  =  unzipOut + (bssn::VAR::U_GT1*NN);   
    double *Gt_rhs2  =  unzipOut + (bssn::VAR::U_GT2*NN);   
    double *B_rhs0   =  unzipOut + (bssn::VAR::U_B0*NN);    
    double *B_rhs1   =  unzipOut + (bssn::VAR::U_B1*NN);    
    double *B_rhs2   =  unzipOut + (bssn::VAR::U_B2*NN);

    auto time_start=Time::now();
    for (unsigned int r=0; r < iter; r++)
    {
        for(unsigned int k=3; k < nz-3; k++)
        {
            for(unsigned int j=3; j < ny-3; j++)
            {
                #pragma novector if(!vectorize)
                #pragma vector if(vectorize)
                #pragma ivdep if(vectorize)
                for(unsigned int i=3; i < nx-3; i++)
                {
                        const double x  =  i * hx;
                        const double y  =  j * hy;
                        const double z  =  k * hz;

                        const double r_coord =sqrt(x*x + y*y + z*z);
                        const double sigma=0.4;
                        double eta=bssn::ETA_CONST;
                        if (r_coord >= bssn::ETA_R0) {
                            eta *= pow( (bssn::ETA_R0/r_coord), bssn::ETA_DAMPING_EXP);
                        }

                        const unsigned int pp = k*ny*nx + j* nx + i;

                        #include "../../src/bssneqs_eta_const_standard_gauge.cpp"
                }
            }
        }
        //prevent compiler from optimizing out the iter loop
        #pragma novector
        for (unsigned int pp=0; pp < NN*bssn::BSSN_NUM_VARS; pp+=10)
        {
            unzipOut[pp] += 1e-6;
            unzipIn[pp] +=1e-6;
        }
    }  

    auto time_end=Time::now();
    fsec fs = time_end - time_start;
    std::cout<<"Time(s): "<<fs.count()<<std::endl;

    #include "bssnrhs_dealloc.h"
    #include "bssnrhs_dealloc_adv.h"

    MPI_Finalize();
    return 0;
    
}