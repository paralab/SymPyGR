/**
 * @brief : Load compute store benchmark. 
*/

#include <iostream>
//#include <mpi.h>
#include<chrono>

typedef double DSCALAR;
typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::milliseconds ms;
typedef std::chrono::duration<double> fsec;



int main(int argc, char **argv)
{

    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " N R" << std::endl;
        return 0;
    }

    const unsigned int N = atoi(argv[1]);
    const unsigned int R   = atoi(argv[2]);
    std::cout<<" array size : "<<N<< " iter : "<<R<<std::endl;

    #pragma omp parallel
    {
        DSCALAR* A = (DSCALAR *)aligned_alloc(64,N * sizeof(DSCALAR));
        DSCALAR* B = (DSCALAR *)aligned_alloc(64,N * sizeof(DSCALAR));
        DSCALAR* C = (DSCALAR *)aligned_alloc(64,N * sizeof(DSCALAR));


        for (int n=0; n<N; n++){
            A[n] = (DSCALAR)1e-4*n;
            B[n] = (DSCALAR)1e-4*n+0.1;
            C[n] = (DSCALAR)0;
        }

        auto t1=Time::now();
        
        for (int iter = 0; iter < R; iter ++){
            for (int i = 0; i < N; i++){
            C[i] += A[i] * B[i]; 
            }
            A[43] += 0.0001;
        }
        
        auto t2=Time::now();
        fsec fs = t2 - t1;
        std::cout<<"Time(s) max: "<<fs.count()<<std::endl;

        for (int n=0; n<N; n++){
            C[n] *= (DSCALAR)1.00143;
        }


    }




    #if 0
        MPI_Init(&argc,&argv);

        MPI_Comm comm = MPI_COMM_WORLD;
        int rank, npes;

        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);

        DSCALAR* A = (DSCALAR *)aligned_alloc(64,N * sizeof(DSCALAR));
        DSCALAR* B = (DSCALAR *)aligned_alloc(64,N * sizeof(DSCALAR));
        DSCALAR* C = (DSCALAR *)aligned_alloc(64,N * sizeof(DSCALAR));


        for (int n=0; n<N; n++){
            A[n] = (DSCALAR)1e-4*n;
            B[n] = (DSCALAR)1e-4*n+0.1;
            C[n] = (DSCALAR)0;
        }

        auto t1=Time::now();
        
        for (int iter = 0; iter < R; iter ++){
        for (int i = 0; i < N; i++){
        C[i] += A[i] * B[i]; 
        }
        }
        
        auto t2=Time::now();
        fsec fs = t2 - t1;
        double time_s  = fs.count();
        double time_sg = 0; 
        
        for (int n=0; n<N; n++){
            C[n] *= (DSCALAR)1.00143;
        }
        
        MPI_Reduce(&time_s,&time_sg,1,MPI_DOUBLE,MPI_MAX,0,comm);
        if(!rank)
            std::cout<<"Time(s) max: "<<time_sg<<std::endl;

        MPI_Finalize();

    #endif

    return 0; 

    
}