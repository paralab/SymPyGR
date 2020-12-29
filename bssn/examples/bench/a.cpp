#ifdef LIKWID_PERFMON
#include <likwid.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif


#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <iostream>


double (*A);
double (*B);
double (*C);



double rtclock(void) {
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday(&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d", stat);
    return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}


using namespace std;
int main() {

  A = (double *)aligned_alloc(64,N * sizeof(A[0]));
  B = (double *)aligned_alloc(64,N * sizeof(B[0]));
  C = (double *)aligned_alloc(64,N * sizeof(C[0]));


  for (int n=0; n<N; n++){
      A[n] = 1e-4*n;
      B[n] = 1e-4*n+0.1;
  }

  double t1 = rtclock();

    LIKWID_MARKER_INIT;

#pragma omp parallel
{
  LIKWID_MARKER_THREADINIT;
}
    LIKWID_MARKER_REGISTER("Compute");
    
    /// code execute ///
    //#pragma omp parallel
{
    LIKWID_MARKER_START("Compute");

    for (int iter = 0; iter < R; iter ++){
	    #pragma vector
		for (int i = 0; i < N; i++){
	    	C[i] += A[i] * B[i]; 
	    }
	}

    LIKWID_MARKER_STOP("Compute");
}

    
    LIKWID_MARKER_CLOSE;


    double t2 = rtclock();
   
    for (int i = 0; i < N; i += 10){
	    	std::cout << C[i];
	   }
	
	std::cout << "\ntime:" << (t2-t1);
    return 0;
}


