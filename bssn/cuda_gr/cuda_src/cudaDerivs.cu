/**
 * Created on: Sep 21, 2018
 * 		Author: Akila, Eranga, Eminda, Ruwan
 **/
 
#include "cudaDerivs.cuh"

__global__ void calc_derivs1(double * dev_var_out, double * dev_var_in, double hx, double hy, double hz, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag,
    #include "list_of_offset_para.h"
    ,
    #include "list_of_para.h"
){
    int tid = blockIdx.x*64 + threadIdx.x;

    #include "calc_deriv_calls_1.cuh"
}

__global__ void calc_derivs2(double * dev_var_out, double * dev_var_in, double hx, double hy, double hz, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag,
    #include "list_of_offset_para.h"
    ,
    #include "list_of_para.h"
){
    int tid = blockIdx.x*64 + threadIdx.x;

    #include "calc_deriv_calls_2.cuh"
}

__global__ void calc_derivs1_bflag(double * dev_var_out, double * dev_var_in, double hx, double hy, double hz, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag,
    #include "list_of_offset_para.h"
    ,
    #include "list_of_para.h"
){
    int tid = blockIdx.x*64 + threadIdx.x;

    #include "calc_deriv_calls_1_bflag.cuh"
}

__global__ void calc_derivs2_bflag(double * dev_var_out, double * dev_var_in, double hx, double hy, double hz, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag,
    #include "list_of_offset_para.h"
    ,
    #include "list_of_para.h"
){
    int tid = blockIdx.x*64 + threadIdx.x;

    #include "calc_deriv_calls_2_bflag.cuh"
}

void calc_deriv_kernel_wrapper(double * dev_var_out, double * dev_var_in, double hx, double hy, double hz, const unsigned int * host_sz, int bflag, cudaStream_t stream,
    #include "list_of_offset_para.h"
    ,
    #include "list_of_para.h"
){
    const int ib = 1;
    const int jb = 1;
    const int kb = 1;
    const int ie = host_sz[0] - 1;
    const int je = host_sz[1] - 1;
    const int ke = host_sz[2] - 1;
    const unsigned int host_sz_x = host_sz[0];
    const unsigned int host_sz_y = host_sz[1];
    const unsigned int host_sz_z = host_sz[2];

    int number_of_threads_required;
    int number_of_blocks;

    if (bflag!=0){
        number_of_threads_required=ceil((ie-ib)*(je-jb)*(ke-kb));
        number_of_blocks=ceil(1.0*number_of_threads_required/64);
        calc_derivs1_bflag <<< number_of_blocks, 64, 0, stream>>> (dev_var_out, dev_var_in, hx, hy, hz, host_sz_x, host_sz_y, host_sz_z, bflag,
            #include "list_of_offset_args.h"
            ,
            #include "list_of_args.h"
        );
    
        number_of_threads_required=ceil((ie-ib)*(je-jb)*(ke-kb));
        number_of_blocks=ceil(1.0*number_of_threads_required/64);
        calc_derivs2_bflag <<< number_of_blocks, 64, 0, stream>>> (dev_var_out, dev_var_in, hx, hy, hz, host_sz_x, host_sz_y, host_sz_z, bflag,
            #include "list_of_offset_args.h"
            ,
            #include "list_of_args.h"
        );
    }else{
        number_of_threads_required=ceil((ie-ib)*(je-jb)*(ke-kb));
        number_of_blocks=ceil(1.0*number_of_threads_required/64);
        calc_derivs1 <<< number_of_blocks, 64, 0, stream>>> (dev_var_out, dev_var_in, hx, hy, hz, host_sz_x, host_sz_y, host_sz_z, bflag,
            #include "list_of_offset_args.h"
            ,
            #include "list_of_args.h"
        );
    
        number_of_threads_required=ceil((ie-ib)*(je-jb)*(ke-kb));
        number_of_blocks=ceil(1.0*number_of_threads_required/64);
        calc_derivs2 <<< number_of_blocks, 64, 0, stream>>> (dev_var_out, dev_var_in, hx, hy, hz, host_sz_x, host_sz_y, host_sz_z, bflag,
            #include "list_of_offset_args.h"
            ,
            #include "list_of_args.h"
        );
    }
    CHECK_ERROR(cudaGetLastError(), "deriv Kernel launch failed");
}



__global__ void calc_ko_derivs(double * dev_var_out, double * dev_var_in, double hx, double hy, double hz, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag,
    #include "list_of_offset_para.h"
    ,
    #include "list_of_para.h"
){
    int tid = blockIdx.x*64 + threadIdx.x;

    #include "calc_ko_deriv_calls.cuh"
}


__global__ void calc_ko_derivs_bflag(double * dev_var_out, double * dev_var_in, double hx, double hy, double hz, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag,
    #include "list_of_offset_para.h"
    ,
    #include "list_of_para.h"
){
    int tid = blockIdx.x*64 + threadIdx.x;

    #include "calc_ko_deriv_calls_bflag.cuh"
}

void calc_ko_deriv_kernel_wrapper(double * dev_var_out, double * dev_var_in, double hx, double hy, double hz, const unsigned int * host_sz, int bflag, cudaStream_t stream,
    #include "list_of_offset_para.h"
    ,
    #include "list_of_para.h"
    )
{
    const int ib = 1;
    const int jb = 1;
    const int kb = 1;
    const int ie = host_sz[0] - 1;
    const int je = host_sz[1] - 1;
    const int ke = host_sz[2] - 1;
    const unsigned int host_sz_x = host_sz[0];
    const unsigned int host_sz_y = host_sz[1];
    const unsigned int host_sz_z = host_sz[2];

    int number_of_threads_required;
    int number_of_blocks;

    number_of_threads_required=ceil((ie-ib)*(je-jb)*(ke-kb));
    number_of_blocks=ceil(1.0*number_of_threads_required/64);

    if (bflag!=0){
        calc_ko_derivs_bflag <<< number_of_blocks, 64, 0, stream>>> (dev_var_out, dev_var_in, hx, hy, hz, host_sz_x, host_sz_y, host_sz_z, bflag,
            #include "list_of_offset_args.h"
            ,
            #include "list_of_args.h"
        );
    }else{
        calc_ko_derivs <<< number_of_blocks, 64, 0, stream>>> (dev_var_out, dev_var_in, hx, hy, hz, host_sz_x, host_sz_y, host_sz_z, bflag,
            #include "list_of_offset_args.h"
            ,
            #include "list_of_args.h"
        );
    }
    
    CHECK_ERROR(cudaGetLastError(), "ko deriv Kernel launch failed");
}
