/**
 * @brief Evaluate the performance of the deriv routines. 
 * 
 */

#include<iostream>
#include "derivs.h"
#include <immintrin.h>
#include <chrono>
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

    double* u  = new double[NN];
    double* Du = new double[NN];
    // double* u=nullptr;
    // double* Du=nullptr;
    // posix_memalign((void**)&u,128,sizeof(double)*NN);
    // posix_memalign((void**)&Du,128,sizeof(double)*NN);

    for(unsigned int i=0; i< NN; i++)
        u[i] = i*0.01;

    auto t1=Time::now();
    for(unsigned int rr=0; rr < iter; rr++)
    {
        const double dx = hx;
        const double idx = 1.0/dx;
        const double idx_by_12 = idx/12.0;
        const double idx_by_60 = idx/60.0;

        const int nx = sz[0];
        const int ny = sz[1];
        const int nz = sz[2];
        const int ib = 3;
        const int jb = 3;
        const int kb = 3;
        const int ie = sz[0]-3;
        const int je = sz[1]-3;
        const int ke = sz[2]-3;
        const int n=1;
        
        
        for (int k = kb; k < ke; k++) 
        for (int j = jb; j < je; j++)
        for (int i = ib; i < ie; i++) 
        {
            int pp =IDX(i,j,k);
            Du[pp] = ( - 1.0   * u[pp-3*n] 
                        + 9.0  * u[pp-2*n]
                        -45.0  * u[pp-1*n]
                        +45.0  * u[pp+1*n]
                        - 9.0  * u[pp+2*n]
                        + 1.0  * u[pp+3*n] ) * idx_by_60;
        }
            
        
    }
    auto t2=Time::now();
    fsec fs = t2 - t1;
    std::cout<<"Time non tiled: "<<fs.count()<<std::endl;

    t1=Time::now();
    for(unsigned int rr=0; rr < iter; rr++)
    {
        const double dx = hx;
        const double idx = 1.0/dx;
        const double idx_by_12 = idx/12.0;
        const double idx_by_60 = idx/60.0;

        const int nx = sz[0];
        const int ny = sz[1];
        const int nz = sz[2];
        const int ib = 3;
        const int jb = 3;
        const int kb = 3;
        const int ie = sz[0]-3;
        const int je = sz[1]-3;
        const int ke = sz[2]-3;
        const int n=nx;
        const int tsz[] = {ie,je,ke};

        //for(int tk=kb; tk < ke; tk+=tsz[2])
        //for(int tj=jb; tj < je; tj+=tsz[1])
        //for(int ti=ib; ti < ie; ti+=tsz[0])
        //for (int k = tk; k < std::min( ke , tk + tsz[2]); k++)
        //for (int j = tj; j < std::min( je , tj + tsz[1]); j++)
        //for (int i = ti; i < std::min( ie , ti + tsz[0]); i++)
        for (int k = kb; k < ke; k++) 
        for (int i = ib; i < ie; i++) 
        for (int j = jb; j < je; j++)
        {
                int pp = IDX(i,j,k);
                Du[pp] = ( - 1.0   * u[pp-3*n] 
                            + 9.0  * u[pp-2*n]
                            -45.0  * u[pp-1*n]
                            +45.0  * u[pp+1*n]
                            - 9.0  * u[pp+2*n]
                            + 1.0  * u[pp+3*n] ) * idx_by_60;
        }

    }
    t2=Time::now();
    fs = t2 - t1;
    std::cout<<"Time tiled: "<<fs.count()<<std::endl;

    delete []  u;
    delete [] Du;

    return 0;
    
}