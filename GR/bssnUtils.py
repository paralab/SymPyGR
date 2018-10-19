##########################################################################
# author: Milinda Fernando
# email:  milinda@cs.utah.edu, 
# date: 08/13/2018
#
# python module to generate bssn derivative calls and support function 
# calls to call the generated code by bssn.py 
# (python code for perl script written by David)
#
##########################################################################

# Note: gbx, gby, gbz are not needed for the RHS, but the derivatives
# are needed for the boundary conditions.  The allocation of derivatives
# and calls to derivative routines for the boundaries uses the functions
# required for the rhs, so I include them here.

from collections import namedtuple
from datetime import datetime
from time import strftime

import dendro as dendro
import math as math
#import bssn_stages as bssn_stages
import bssn as bssn
import sympy as sympy
import re as re
import os as os

import cudaSharedMemManager as SharedMemManager


## ==== BSSN GPU code generation paramerters

D = ["alpha", "beta0", "beta1", "beta2",
      "B0", "B1", "B2",
      "chi", "Gt0", "Gt1", "Gt2", "K",
      "gt0", "gt1", "gt2", "gt3", "gt4", "gt5",
      "At0", "At1", "At2", "At3", "At4", "At5" ]

# variable names, to access the 2D array.
VAR_ENUM=["cuda::VAR::U_ALPHA",
          "cuda::VAR::U_BETA0",
          "cuda::VAR::U_BETA1",
          "cuda::VAR::U_BETA2",
          "cuda::VAR::U_B0",
          "cuda::VAR::U_B1",
          "cuda::VAR::U_B2",
          "cuda::VAR::U_CHI",
          "cuda::VAR::U_GT0",
          "cuda::VAR::U_GT1",
          "cuda::VAR::U_GT2",
          "cuda::VAR::U_K",
          "cuda::VAR::U_SYMGT0",
          "cuda::VAR::U_SYMGT1",
          "cuda::VAR::U_SYMGT2",
          "cuda::VAR::U_SYMGT3",
          "cuda::VAR::U_SYMGT4",
          "cuda::VAR::U_SYMGT5",
          "cuda::VAR::U_SYMAT0",
          "cuda::VAR::U_SYMAT1",
          "cuda::VAR::U_SYMAT2",
          "cuda::VAR::U_SYMAT3",
          "cuda::VAR::U_SYMAT4",
          "cuda::VAR::U_SYMAT5"]

# enum to symbolic input vars dictionary
VAR_ENUM_TO_INPUT_SYM =  { "alpha" : "cuda::VAR::U_ALPHA",
                       "beta0" : "cuda::VAR::U_BETA0",
                       "beta1" : "cuda::VAR::U_BETA1",
                       "beta2" : "cuda::VAR::U_BETA2",
                       "B0"    : "cuda::VAR::U_B0",
                       "B1"    : "cuda::VAR::U_B1",
                       "B2"    : "cuda::VAR::U_B2",
                       "chi"   : "cuda::VAR::U_CHI",
                       "Gt0"   : "cuda::VAR::U_GT0",
                       "Gt1"   : "cuda::VAR::U_GT1",
                       "Gt2"   : "cuda::VAR::U_GT2",
                       "K"     : "cuda::VAR::U_K",
                       "gt0"   : "cuda::VAR::U_SYMGT0",
                       "gt1"   : "cuda::VAR::U_SYMGT1",
                       "gt2"   : "cuda::VAR::U_SYMGT2",
                       "gt3"   : "cuda::VAR::U_SYMGT3",
                       "gt4"   : "cuda::VAR::U_SYMGT4",
                       "gt5"   : "cuda::VAR::U_SYMGT5",
                       "At0"   :"cuda::VAR::U_SYMAT0",
                       "At1"   :"cuda::VAR::U_SYMAT1",
                       "At2"   :"cuda::VAR::U_SYMAT2",
                       "At3"   :"cuda::VAR::U_SYMAT3",
                       "At4"   :"cuda::VAR::U_SYMAT4",
                       "At5"   :"cuda::VAR::U_SYMAT5"
                      }

# enum to symbolic output vars dictionary
VAR_ENUM_TO_OUTPUT_SYM={   "a_rhs" : "cuda::VAR::U_ALPHA",
                           "b_rhs0" : "cuda::VAR::U_BETA0",
                           "b_rhs1" : "cuda::VAR::U_BETA1",
                           "b_rhs2" : "cuda::VAR::U_BETA2",
                           "B_rhs0"    : "cuda::VAR::U_B0",
                           "B_rhs1"    : "cuda::VAR::U_B1",
                           "B_rhs2"    : "cuda::VAR::U_B2",
                           "chi_rhs"   : "cuda::VAR::U_CHI",
                           "Gt_rhs0"   : "cuda::VAR::U_GT0",
                           "Gt_rhs1"   : "cuda::VAR::U_GT1",
                           "Gt_rhs2"   : "cuda::VAR::U_GT2",
                           "K_rhs"     : "cuda::VAR::U_K",
                           "gt_rhs00"   : "cuda::VAR::U_SYMGT0",
                           "gt_rhs01"   : "cuda::VAR::U_SYMGT1",
                           "gt_rhs02"   : "cuda::VAR::U_SYMGT2",
                           "gt_rhs11"   : "cuda::VAR::U_SYMGT3",
                           "gt_rhs12"   : "cuda::VAR::U_SYMGT4",
                           "gt_rhs22"   : "cuda::VAR::U_SYMGT5",
                           "At_rhs00"   :"cuda::VAR::U_SYMAT0",
                           "At_rhs01"   :"cuda::VAR::U_SYMAT1",
                           "At_rhs02"   :"cuda::VAR::U_SYMAT2",
                           "At_rhs11"   :"cuda::VAR::U_SYMAT3",
                           "At_rhs12"   :"cuda::VAR::U_SYMAT4",
                           "At_rhs22"   :"cuda::VAR::U_SYMAT5"
                       }

# custom functions for code generation in cse.
custom_functions = {'grad': 'grad', 'grad2': 'grad2', 'agrad': 'agrad', 'kograd': 'kograd'}

# second derivs required for RHS
DD = ["gt0", "gt1", "gt2", "gt3", "gt4", "gt5", "chi",
       "alpha", "beta0", "beta1", "beta2" ]

# advective derivatives
AD = ["gt0", "gt1", "gt2", "gt3", "gt4", "gt5",
       "At0", "At1", "At2", "At3", "At4", "At5",
       "alpha", "beta0", "beta1", "beta2", "chi", "Gt0", "Gt1", "Gt2", "K",
       "B0", "B1", "B2"] 

KO=AD

# first derivs required for constraints--no gauge variables
CONSTRAINT_D = [ "chi", "Gt0", "Gt1", "Gt2", "K",
           "gt0", "gt1", "gt2", "gt3", "gt4", "gt5",
           "At0", "At1", "At2", "At3", "At4", "At5" ]

# second derivs required for constraints--no gauge variables
CONSTRAINT_DD = ["gt0", "gt1", "gt2", "gt3", "gt4", "gt5", "chi"]


PREFIX_D   = ["grad_0_", "grad_1_", "grad_2_"]
PREFIX_AD  = ["agrad_0_", "agrad_1_", "agrad_2_"]
PREFIX_KOD = ["kograd_0_", "kograd_1_", "kograd_2_"]
PREFIX_DD  = ["grad2_0_0_", "grad2_0_1_", "grad2_0_2_", "grad2_1_1_", "grad2_1_2_", "grad2_2_2_"]

# first derivative in i direction
FUNC_D_I=[]
for f in D:
    for p in PREFIX_D:
        FUNC_D_I.append(p+f)

# second derivative in ij direction
FUNC_D_IJ=[]
for f in DD:
    for p in PREFIX_DD:
        FUNC_D_IJ.append(p+f)

#advective derivative in i direction
FUNC_AD_I=[]
for f in AD:
    for p in PREFIX_AD:
        FUNC_AD_I.append(p+f)


#Kriess-Oliger derivative in i direction
FUNC_KOD_I=[]
for f in D:
    for p in PREFIX_KOD:
        FUNC_KOD_I.append(p+f)


# cuda utility functions 
## Note all the device vars which is global starts with __
FUNC_LOAD_VAR="cuda::__loadGlobalToShared3D<double>"
FUNC_STORE_VAR="cuda::__storeSharedToGlobal3D<double>"
FUNC_SIGN_EXT="cuda::__extractSign3D<double>"

VAR_UNZIP_IN="__unzipInVar"
VAR_UNZIP_OUT="__unzipOutVar"


## shift vector block shared variables to compute advective derivs
VAR_BETA0="beta0"
VAR_BETA1="beta1"
VAR_BETA2="beta2"

VAR_BETA0_BOOL="beta0_bool"
VAR_BETA1_BOOL="beta1_bool"
VAR_BETA2_BOOL="beta2_bool"


# shared input variable name for derivative kernels
VAR_IN_SHARED="unzipVarInShared"
# shared output variable name for derivative kernels
VAR_OUT_SHARED_0="unzipVarOutShared0"
VAR_OUT_SHARED_1="unzipVarOutShared1"


# block ids
VAR_BLK_ID_X="blockIdx.x"
VAR_BLK_ID_Y="blockIdx.y"
VAR_BLK_ID_Z="blockIdx.z"

# thread ids
VAR_TRD_ID_X="threadIdx.x"
VAR_TRD_ID_Y="threadIdx.y"
VAR_TRD_ID_Z="threadIdx.z"

# block dim
VAR_BLK_DIM_X="blockDim.x"
VAR_BLK_DIM_Y="blockDim.y"
VAR_BLK_DIM_Z="blockDim.z"


# x,y,z bounds of the time i_lm[0] is the min and i_lm[1] is the max. 
VAR_TILE_SZ="tile_sz"
VAR_DENDRO_BLK_ALIGNED_SZ="alignedSz"
VAR_TILE_LIMITS="ijk_lm"
VAR_TILE_LIMITS_STORE="tile_lm"


TYPE_DERIV_STRUCT="MemoryDerivs"
TYPE_BLK_CU="cuda::_Block"
TYPE_BSSN_COMP_PARS="BSSNComputeParams"


##
# generate the code to allocate derivative memory variables (allocated size unzip_dof)
##
def cudaDerivAllocDeallocHeader(fname,headers=[]):
    
    func_i=FUNC_D_I
    func_ij=FUNC_D_IJ
    afunc_i=FUNC_AD_I
    kofunc_i=FUNC_KOD_I

    with open(fname, 'w') as ofile:
        ofile.write("// generated by Dendro-GR SymPyGR code gernation framework\n")
        ofile.write("//date: "+str(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+"\n")
        fileName,fileExt=os.path.splitext(os.path.basename(fname))
        ofile.write("#ifndef "+fileName.upper()+"_"+fileExt[1:].upper()+" \n")
        ofile.write("#define "+fileName.upper()+"_"+fileExt[1:].upper()+" \n")
        ofile.write("\n")
        ofile.write("#include <iostream>\n")
        ofile.write("#include \"cuda_runtime.h\"\n")
        for header in headers:
            ofile.write("#include \""+header+"\"\n")
        ofile.write("\n")
        ofile.write("namespace cuda {\n")

        ofile.write("\tstruct "+TYPE_DERIV_STRUCT+"{\n\n")
        ofile.write("/**@brief upper bound of the block size processed by the GPU*/\n")
        ofile.write("\t unsigned int __maxBlkSz;\n")

        ofile.write("/**@brief number of streams the kernel get executed*/\n")
        ofile.write("\t unsigned int __numStream;\n")

        ofile.write("/**@brief size per stream*/\n")
        ofile.write("\t unsigned int __szPerStream;\n")


        for deriv in func_i:
            ofile.write("\t double* __"+deriv+";\n")

        for deriv in func_ij:
            ofile.write("\t double* __"+deriv+";\n")

        for deriv in afunc_i:
            ofile.write("\t double* __"+deriv+";\n")

        for deriv in kofunc_i:
            ofile.write("\t double* __"+deriv+";\n")

        ofile.write("\n")
        ofile.write("\n")

        ofile.write("\t/**@brief memory allocation for deriv variables */\n")
        ofile.write("\tvoid allocateDerivMemory(unsigned int maxBlkSz, unsigned int numSM,unsigned int numStreams=1); \n")

        ofile.write("\n")
        ofile.write("\n")

        ofile.write("\t/**@brief memory deallocation for deriv variables */\n")
        ofile.write("\tvoid deallocateDerivMemory(); \n")

        ofile.write("\n")
        ofile.write("\n")


        ofile.write("\t};\n\n")
        ofile.write("}// end of namespace cuda\n")
        ofile.write("\n")

        ofile.write("#endif\n")
    
    ofile.close()

##
# generate the code to allocate derivative memory variables (allocated size unzip_dof)
##
def cudaDerivAllocDeallocSource(fname,headers=[]):
    
    func_i=FUNC_D_I
    func_ij=FUNC_D_IJ
    afunc_i=FUNC_AD_I
    kofunc_i=FUNC_KOD_I


    with open(fname, 'w') as ofile:
        ofile.write("// generated by Dendro-GR SymPyGR code gernation framework\n")
        ofile.write("//date: "+str(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+"\n")
        for header in headers:
            ofile.write("#include \""+header+"\"\n")
        ofile.write("\n")
        ofile.write("namespace cuda {\n")

        ofile.write("\n")
        ofile.write("\n")

        ofile.write("\t/**@brief memory allocation for deriv variables */\n")
        ofile.write("\tvoid cuda::"+TYPE_DERIV_STRUCT+"::allocateDerivMemory(unsigned int maxBlkSz, unsigned int numSM,unsigned int numStreams){ \n")

        ofile.write("\t\t __maxBlkSz=maxBlkSz;\n")
        ofile.write("\t\t __numStream=numStreams;\n")
        ofile.write("\t\t __szPerStream=numSM*maxBlkSz;\n")

        ofile.write("\t\t const size_t bytes=sizeof(double)*numSM*maxBlkSz*numStreams;\n")

        for deriv in func_i:
            ofile.write("\t\t cudaMalloc((void **)&__"+deriv+",bytes);\n")
            #ofile.write("\t\tCUDA_CHECK_ERROR();\n")

        for deriv in func_ij:
            ofile.write("\t\t cudaMalloc((void **)&__"+deriv+",bytes);\n")
            #ofile.write("\t\tCUDA_CHECK_ERROR();\n")

        for deriv in afunc_i:
            ofile.write("\t\t cudaMalloc((void **)&__"+deriv+",bytes);\n")
            #ofile.write("\t\tCUDA_CHECK_ERROR();\n")

        for deriv in kofunc_i:
            ofile.write("\t\t cudaMalloc((void **)&__"+deriv+",bytes);\n")
            #ofile.write("\t\tCUDA_CHECK_ERROR()\n")

        ofile.write("\n")
        ofile.write("\n")

        ofile.write("} \n")

        ofile.write("\n")
        ofile.write("\n")

        ofile.write("\t/**@brief memory deallocation for deriv variables */\n")
        ofile.write("\tvoid cuda::"+TYPE_DERIV_STRUCT+"::deallocateDerivMemory(){ \n")

        for deriv in func_i:
            ofile.write("\t\t cudaFree((void *)__"+deriv+");\n")
            #ofile.write("\t\tCUDA_CHECK_ERROR()\n")

        for deriv in func_ij:
            ofile.write("\t\t cudaFree((void *)__"+deriv+");\n")
            #ofile.write("\t\tCUDA_CHECK_ERROR()\n")

        for deriv in afunc_i:
            ofile.write("\t\t cudaFree((void *)__"+deriv+");\n")
            #ofile.write("\t\tCUDA_CHECK_ERROR()\n")

        for deriv in kofunc_i:
            ofile.write("\t\t cudaFree((void *)__"+deriv+");\n")
            #ofile.write("\t\tCUDA_CHECK_ERROR()\n")
        
        ofile.write("\n")
        ofile.write("\n")

        ofile.write("} \n")

        ofile.write("\n")
        ofile.write("\n")

        ofile.write("}// end of namespace cuda\n")

    ofile.close()


def computeTileStore(dir,out,padWidth=3):

    padWidth=str(padWidth)
    out.write("\n")
    if(dir=="x"):

        out.write("\t\t"+VAR_TILE_LIMITS_STORE+"[0]="+padWidth+";\n")
        out.write("\t\t"+VAR_TILE_LIMITS_STORE+"[1]=("+VAR_TILE_LIMITS+"[1]-"+VAR_TILE_LIMITS+"[0]);\n")
        out.write("\n")
        out.write("\t\t"+VAR_TILE_LIMITS_STORE+"[2]=(iter_y)? "+padWidth+": 0;\n")
        out.write("\t\t"+VAR_TILE_LIMITS_STORE+"[3]=(iter_y==(BLK_ITERATIONS_Y-1)) ? ("+VAR_TILE_LIMITS+"[3]-"+VAR_TILE_LIMITS+"[2]) : "+"("+VAR_TILE_LIMITS+"[3]-"+VAR_TILE_LIMITS+"[2]-"+padWidth+")" +";\n")
        out.write("\n")
        out.write("\t\t"+VAR_TILE_LIMITS_STORE+"[4]=(iter_z)? "+padWidth+": 0;\n")
        out.write("\t\t"+VAR_TILE_LIMITS_STORE+"[5]=(iter_z==(BLK_ITERATIONS_Z-1)) ? ("+VAR_TILE_LIMITS+"[5]-"+VAR_TILE_LIMITS+"[4]) : "+"("+VAR_TILE_LIMITS+"[5]-"+VAR_TILE_LIMITS+"[4]-"+padWidth+")" +";\n")
        out.write("\n")

    elif(dir=="y"):

        out.write("\t\t"+VAR_TILE_LIMITS_STORE+"[0]=(iter_x)? "+padWidth+": 0;\n")
        out.write("\t\t"+VAR_TILE_LIMITS_STORE+"[1]=(iter_x==(BLK_ITERATIONS_X-1)) ? ("+VAR_TILE_LIMITS+"[1]-"+VAR_TILE_LIMITS+"[0]) : "+"("+VAR_TILE_LIMITS+"[1]-"+VAR_TILE_LIMITS+"[0]-"+padWidth+")" +";\n")
        out.write("\n")
        out.write("\t\t"+VAR_TILE_LIMITS_STORE+"[2]="+padWidth+";\n")
        out.write("\t\t"+VAR_TILE_LIMITS_STORE+"[3]=("+VAR_TILE_LIMITS+"[3]-"+VAR_TILE_LIMITS+"[2]);\n")
        out.write("\n")
        out.write("\t\t"+VAR_TILE_LIMITS_STORE+"[4]=(iter_z)? "+padWidth+": 0;\n")
        out.write("\t\t"+VAR_TILE_LIMITS_STORE+"[5]=(iter_z==(BLK_ITERATIONS_Z-1)) ? ("+VAR_TILE_LIMITS+"[5]-"+VAR_TILE_LIMITS+"[4]) : "+"("+VAR_TILE_LIMITS+"[5]-"+VAR_TILE_LIMITS+"[4]-"+padWidth+")" +";\n")
        out.write("\n")

    elif(dir=="z"):
        out.write("\t\t"+VAR_TILE_LIMITS_STORE+"[0]=(iter_x)? "+padWidth+": 0;\n")
        out.write("\t\t"+VAR_TILE_LIMITS_STORE+"[1]=(iter_x==(BLK_ITERATIONS_X-1)) ? ("+VAR_TILE_LIMITS+"[1]-"+VAR_TILE_LIMITS+"[0]) : "+"("+VAR_TILE_LIMITS+"[1]-"+VAR_TILE_LIMITS+"[0]-"+padWidth+")" +";\n")
        out.write("\n")
        out.write("\t\t"+VAR_TILE_LIMITS_STORE+"[2]=(iter_y)? "+padWidth+": 0;\n")
        out.write("\t\t"+VAR_TILE_LIMITS_STORE+"[3]=(iter_y==(BLK_ITERATIONS_Y-1)) ? ("+VAR_TILE_LIMITS+"[3]-"+VAR_TILE_LIMITS+"[2]) : "+"("+VAR_TILE_LIMITS+"[3]-"+VAR_TILE_LIMITS+"[2]-"+padWidth+")" +";\n")
        out.write("\n")
        out.write("\t\t"+VAR_TILE_LIMITS_STORE+"[4]="+padWidth+";\n")
        out.write("\t\t"+VAR_TILE_LIMITS_STORE+"[5]=("+VAR_TILE_LIMITS+"[5]-"+VAR_TILE_LIMITS+"[4]);\n")
        out.write("\n")

    out.write("\n")




def cudaCompute(fname_cuh,fname_cu,derivs,outs,varnames,kernelName,headers=[]):
     # cuda device properties
    VAR_CUDA_DEVICE="__deviceProperties"
    # dendro block list parameters
    VAR_DENDRO_BLK_LIST="__dendroBlkList"
    VAR_NUM_BLOCKS="cuda::__DENDRO_NUM_BLOCKS"
    VAR_DERIV_WORKSPACE="__derivWorkspace"
    VAR_GPU_BLOCK_MAP="__gpuBlockMap"
    VAR_MAX_DENDRO_BLK_SZ=VAR_DERIV_WORKSPACE+"->__maxBlkSz"
    VAR_DW_SZ_PER_STREAM=VAR_DERIV_WORKSPACE+"->__szPerStream"

    VAR_BSSN_PARAMS="__bssnParams"
    TYPE_BSSN_PARAMS="cuda::BSSNComputeParams"

    FUNC_DERIV_COMP="__compute_derivatives"
    FUNC_COMP_RHS_PRE="__compute"
    FUNC_KO_DISS="__ko_dissipation"

    VAR_SHARED_MEM="__sm_base"
    TYPE_SHARED_MEM="double"

    VAR_ADV_COMPRESS_0=VAR_BETA0_BOOL
    VAR_ADV_COMPRESS_1=VAR_BETA1_BOOL
    VAR_ADV_COMPRESS_2=VAR_BETA2_BOOL
    TYPE_ADV_COMPRESS="bool"

    VAR_DBLOCK="dblock"
    VAR_STREAM_ID="stream_id"
    VAR_DERIV_WORKSPACE_OFFSET=VAR_STREAM_ID+"*("+VAR_DW_SZ_PER_STREAM+") + SM_ID*("+VAR_MAX_DENDRO_BLK_SZ+")"



    ######################################################
    # Writing the header
    ######################################################
    with open(fname_cuh, 'w') as ofile:
        ofile.write("// generated by Dendro-GR SymPyGR code gernation framework\n")
        ofile.write("//date: "+str(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+"\n")
        fileName,fileExt=os.path.splitext(os.path.basename(fname_cuh))
        ofile.write("#ifndef "+fileName.upper()+"_"+fileExt[1:].upper()+" \n")
        ofile.write("#define "+fileName.upper()+"_"+fileExt[1:].upper()+" \n")
        ofile.write("#include<iostream>\n")
        ofile.write("#include\"cuda_runtime.h\"\n")
        ofile.write("#include<device_launch_parameters.h>\n")
        for header in headers:
            ofile.write("#include \""+header+"\"\n")

        ofile.write("namespace cuda {\n")
        ofile.write("\n")

        ofile.write("/**@brief compute derivs \n")
        ofile.write(" @param[in] "+VAR_UNZIP_IN+": unzipped input array (global memory) \n")
        ofile.write(" @param[in] "+TYPE_DERIV_STRUCT+": allocated workspace for derivative computations \n")
        ofile.write(" @param[in] "+VAR_DENDRO_BLK_LIST+": dendro block list \n")
        ofile.write(" @param[in] "+VAR_GPU_BLOCK_MAP+": gpu block map  \n")
        ofile.write(" @param[in] "+VAR_CUDA_DEVICE+": cuda device properties  \n")
        #ofile.write(" @param[out] "+VAR_UNZIP_OUT+": unzipped output computed rhs  \n")
        ofile.write("*/ \n")
        ofile.write("__device__ void "+FUNC_DERIV_COMP+"(const double**"+VAR_UNZIP_IN+","+TYPE_DERIV_STRUCT+"* __derivWorkspace, const "+ TYPE_BLK_CU+ "* dblock, const unsigned int * "+VAR_GPU_BLOCK_MAP+",const "+TYPE_BSSN_PARAMS+" * "+VAR_BSSN_PARAMS+",const cudaDeviceProp* "+VAR_CUDA_DEVICE+", "+TYPE_SHARED_MEM+"* "+VAR_SHARED_MEM+", "+TYPE_ADV_COMPRESS+"* "+VAR_ADV_COMPRESS_0+", "+TYPE_ADV_COMPRESS+"* "+VAR_ADV_COMPRESS_1+", "+TYPE_ADV_COMPRESS+"* "+VAR_ADV_COMPRESS_2+",unsigned int "+VAR_STREAM_ID+");\n")
        ofile.write("\n")

        for var in varnames:
            ofile.write("/**@brief compute "+var+" \n")
            ofile.write(" @param[in] "+VAR_UNZIP_IN+": unzipped input array (global memory) \n")
            ofile.write(" @param[in] "+TYPE_DERIV_STRUCT+": allocated workspace for derivative computations \n")
            ofile.write(" @param[in] "+VAR_DENDRO_BLK_LIST+": dendro block list \n")
            ofile.write(" @param[in] "+VAR_GPU_BLOCK_MAP+": gpu block map  \n")
            ofile.write(" @param[in] "+VAR_CUDA_DEVICE+": cuda device properties  \n")
            ofile.write(" @param[out] "+VAR_UNZIP_OUT+": unzipped output computed rhs  \n")
            ofile.write("*/ \n")
            ofile.write("__device__ void "+FUNC_COMP_RHS_PRE+"_"+var+"(double **"+VAR_UNZIP_OUT+", const double**"+VAR_UNZIP_IN+","+TYPE_DERIV_STRUCT+"* __derivWorkspace, const "+ TYPE_BLK_CU+ "* dblock, const unsigned int * "+VAR_GPU_BLOCK_MAP+",const "+TYPE_BSSN_PARAMS+" * "+VAR_BSSN_PARAMS+",const cudaDeviceProp* "+VAR_CUDA_DEVICE+", "+TYPE_SHARED_MEM+"* "+VAR_SHARED_MEM+", unsigned int "+VAR_STREAM_ID+");\n")
            ofile.write("\n")

        ofile.write("/**@brief apply KO dissipation \n")
        ofile.write(" @param[in] "+VAR_UNZIP_IN+": unzipped input array (global memory) \n")
        ofile.write(" @param[in] "+TYPE_DERIV_STRUCT+": allocated workspace for derivative computations \n")
        ofile.write(" @param[in] "+VAR_DENDRO_BLK_LIST+": dendro block list \n")
        ofile.write(" @param[in] "+VAR_GPU_BLOCK_MAP+": gpu block map  \n")
        ofile.write(" @param[in] "+VAR_CUDA_DEVICE+": cuda device properties  \n")
        ofile.write(" @param[out] "+VAR_UNZIP_OUT+": unzipped output computed rhs  \n")
        ofile.write("*/ \n")
        ofile.write("__device__ void "+FUNC_KO_DISS+"(double **"+VAR_UNZIP_OUT+", const double**"+VAR_UNZIP_IN+","+TYPE_DERIV_STRUCT+"* __derivWorkspace, const "+ TYPE_BLK_CU+ "* dblock, const unsigned int * "+VAR_GPU_BLOCK_MAP+",const "+TYPE_BSSN_PARAMS+" * "+VAR_BSSN_PARAMS+",const cudaDeviceProp* "+VAR_CUDA_DEVICE+", "+TYPE_SHARED_MEM+"* "+VAR_SHARED_MEM+", unsigned int "+VAR_STREAM_ID+");\n")
        ofile.write("\n")

        ofile.write("/**@brief compute RHS \n")
        ofile.write(" @param[in] "+VAR_UNZIP_IN+": unzipped input array (global memory) \n")
        ofile.write(" @param[in] "+TYPE_DERIV_STRUCT+": allocated workspace for derivative computations \n")
        ofile.write(" @param[in] "+VAR_DENDRO_BLK_LIST+": dendro block list \n")
        ofile.write(" @param[in] "+VAR_GPU_BLOCK_MAP+": gpu block map  \n")
        ofile.write(" @param[in] "+VAR_CUDA_DEVICE+": cuda device properties  \n")
        ofile.write(" @param[out] "+VAR_UNZIP_OUT+": unzipped output computed rhs  \n")
        ofile.write("*/ \n")
        ofile.write("__global__ void "+kernelName+"(double **"+VAR_UNZIP_OUT+", const double**"+VAR_UNZIP_IN+","+TYPE_DERIV_STRUCT+"* __derivWorkspace, const "+ TYPE_BLK_CU+ "* "+VAR_DENDRO_BLK_LIST+", const unsigned int * "+VAR_GPU_BLOCK_MAP+",const "+TYPE_BSSN_PARAMS+" * "+VAR_BSSN_PARAMS+",const cudaDeviceProp* "+VAR_CUDA_DEVICE+", unsigned int "+VAR_STREAM_ID+");\n")
        ofile.write("\n")
        ofile.write("}// end of namespace cuda\n")

        ofile.write("\n")
        ofile.write("\n")
        ofile.write("#endif\n")


    ofile.close()

    
    ######################################################
    # Writing the source
    ######################################################

    with open(fname_cu, 'w') as ofile:
        ofile.write("// generated by Dendro-GR SymPyGR code gernation framework\n")
        ofile.write("//date: "+str(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+"\n")
        ofile.write("#include \""+os.path.basename(fname_cuh)+"\"\n")
        # namespace begin
        ofile.write("namespace cuda {\n\n")


        ofile.write("/**@brief compute RHS \n")
        ofile.write(" @param[in] "+VAR_UNZIP_IN+": unzipped input array (global memory) \n")
        ofile.write(" @param[in] "+TYPE_DERIV_STRUCT+": allocated workspace for derivative computations \n")
        ofile.write(" @param[in] "+VAR_DENDRO_BLK_LIST+": dendro block list \n")
        ofile.write(" @param[in] "+VAR_GPU_BLOCK_MAP+": gpu block map  \n")
        ofile.write(" @param[in] "+VAR_CUDA_DEVICE+": cuda device properties  \n")
        ofile.write(" @param[out] "+VAR_UNZIP_OUT+": unzipped output computed rhs  \n")
        ofile.write("*/ \n")
        # function begin
        ofile.write("__global__ void "+kernelName+"(double **"+VAR_UNZIP_OUT+", const double**"+VAR_UNZIP_IN+","+TYPE_DERIV_STRUCT+"* __derivWorkspace, const "+ TYPE_BLK_CU+ "* "+VAR_DENDRO_BLK_LIST+", const unsigned int * "+VAR_GPU_BLOCK_MAP+",const "+TYPE_BSSN_PARAMS+" * "+VAR_BSSN_PARAMS+" ,const cudaDeviceProp* "+VAR_CUDA_DEVICE+", unsigned int "+VAR_STREAM_ID+"){\n\n")

        ofile.write("// shared memory allocation for deriv and rhs computation\n")
        memManager=SharedMemManager.MemoryManager(maxMemSz=48*1024,memUsable=41*1024,cout=ofile,baseName=VAR_SHARED_MEM,varType=TYPE_SHARED_MEM)

        deriv_tile_sz_1d=0
        deriv_req_pad=0
        deriv_max_pad=3

        for deriv in derivs:
            if deriv_tile_sz_1d <deriv.DerivTile1D :
                deriv_tile_sz_1d=deriv.DerivTile1D
            if deriv_req_pad < deriv.padWidth :
                deriv_req_pad =deriv.padWidth

        deriv_tile_sz=deriv_tile_sz_1d**3

        ofile.write("\t__shared__ bool "+VAR_BETA0_BOOL+ "["+str(deriv_tile_sz)+"];\n")
        ofile.write("\t__shared__ bool "+VAR_BETA1_BOOL+ "["+str(deriv_tile_sz)+"];\n")
        ofile.write("\t__shared__ bool "+VAR_BETA2_BOOL+ "["+str(deriv_tile_sz)+"];\n\n")


        ofile.write("\tfor(unsigned int blk="+VAR_GPU_BLOCK_MAP+"[2*"+VAR_BLK_ID_X+"];blk<"+VAR_GPU_BLOCK_MAP+"[2*"+VAR_BLK_ID_X+"+1];++blk){\n\n\n")
        ofile.write("\t// blocks assigned to each gpu block \n")

        ofile.write("\tconst _Block * "+VAR_DBLOCK+"=&"+VAR_DENDRO_BLK_LIST+"[blk];\n")


        ofile.write("\t// compute the derivatives\n")
        ofile.write("\t"+FUNC_DERIV_COMP+"("+VAR_UNZIP_IN+","+VAR_DERIV_WORKSPACE+","+VAR_DBLOCK+","+VAR_GPU_BLOCK_MAP+","+VAR_BSSN_PARAMS+","+VAR_CUDA_DEVICE+","+VAR_SHARED_MEM+","+VAR_ADV_COMPRESS_0+","+VAR_ADV_COMPRESS_1+","+VAR_ADV_COMPRESS_2+","+VAR_STREAM_ID+");\n")
        ofile.write("\t__syncthreads();\n")

        ofile.write("\t// compute the RHS\n")
        for var in varnames:
            ofile.write("\t"+FUNC_COMP_RHS_PRE+"_"+var+"("+VAR_UNZIP_OUT+","+VAR_UNZIP_IN+","+VAR_DERIV_WORKSPACE+","+VAR_DBLOCK+","+VAR_GPU_BLOCK_MAP+","+VAR_BSSN_PARAMS+","+VAR_CUDA_DEVICE+","+VAR_SHARED_MEM+","+VAR_STREAM_ID+");\n")
            ofile.write("\t__syncthreads();\n")

        ofile.write("\t"+FUNC_KO_DISS+"("+VAR_UNZIP_OUT+","+VAR_UNZIP_IN+","+VAR_DERIV_WORKSPACE+","+VAR_DBLOCK+","+VAR_GPU_BLOCK_MAP+","+VAR_BSSN_PARAMS+","+VAR_CUDA_DEVICE+","+VAR_SHARED_MEM+","+VAR_STREAM_ID+");\n")
        ofile.write("\t__syncthreads();\n")

        ofile.write("\t}// end of the block loop\n")

        ofile.write("} // end of kernel \n\n")
        ofile.write("\n")



        ofile.write("/**@brief compute derivs \n")
        ofile.write(" @param[in] "+VAR_UNZIP_IN+": unzipped input array (global memory) \n")
        ofile.write(" @param[in] "+TYPE_DERIV_STRUCT+": allocated workspace for derivative computations \n")
        ofile.write(" @param[in] "+VAR_DENDRO_BLK_LIST+": dendro block list \n")
        ofile.write(" @param[in] "+VAR_GPU_BLOCK_MAP+": gpu block map  \n")
        ofile.write(" @param[in] "+VAR_CUDA_DEVICE+": cuda device properties  \n")
        #ofile.write(" @param[out] "+VAR_UNZIP_OUT+": unzipped output computed rhs  \n")
        ofile.write("*/ \n")
        ofile.write("__device__ void "+FUNC_DERIV_COMP+"(const double**"+VAR_UNZIP_IN+","+TYPE_DERIV_STRUCT+"* __derivWorkspace, const "+ TYPE_BLK_CU+ "* dblock, const unsigned int * "+VAR_GPU_BLOCK_MAP+",const "+TYPE_BSSN_PARAMS+" * "+VAR_BSSN_PARAMS+",const cudaDeviceProp* "+VAR_CUDA_DEVICE+", "+TYPE_SHARED_MEM+"* "+VAR_SHARED_MEM+", "+TYPE_ADV_COMPRESS+"* "+VAR_ADV_COMPRESS_0+", "+TYPE_ADV_COMPRESS+"* "+VAR_ADV_COMPRESS_1+", "+TYPE_ADV_COMPRESS+"* "+VAR_ADV_COMPRESS_2+", unsigned int "+VAR_STREAM_ID+"){\n")
        ofile.write("\n")


        ofile.write("\tconst unsigned int NUM_SM_UNITS="+VAR_CUDA_DEVICE+"->multiProcessorCount;\n")
        ofile.write("\tconst unsigned int SM_ID=get_smid();//"+VAR_BLK_ID_X+"%NUM_SM_UNITS;\n")
        ofile.write("\tconst unsigned int offset=dblock->getOffset();\n")
        ofile.write("\tconst unsigned int *sz=dblock->getSz();\n")
        ofile.write("\tconst unsigned int *"+VAR_DENDRO_BLK_ALIGNED_SZ+"=dblock->getAlignedSz();\n")
        ofile.write("\tconst double* hx=dblock->getDx();\n")

        ofile.write("\tconst double dx=hx[0];\n")
        ofile.write("\tconst double dy=hx[1];\n")
        ofile.write("\tconst double dz=hx[2];\n")

        ofile.write("\tconst double* ptmin=dblock->getPtMin();\n")
        ofile.write("\tconst double* ptmax=dblock->getPtMax();\n")
        ofile.write("\tconst unsigned int bflag=dblock->getBFlag();\n")
        ofile.write("\n")

        if(deriv_req_pad>deriv_max_pad):
            print("code generation error : maxPadwith for derivatives is larger than the dendro block pad width\n")
            os.sys.exit(0)

        ofile.write("\tconst unsigned int "+VAR_TILE_SZ+"[3]={"+str(deriv_tile_sz_1d)+","+str(deriv_tile_sz_1d)+","+str(deriv_tile_sz_1d)+"};\n")

        memManager.malloc(VAR_IN_SHARED,deriv_tile_sz,ofile,prefix="\t")
        memManager.malloc(VAR_OUT_SHARED_0,deriv_tile_sz,ofile,prefix="\t")
        memManager.malloc(VAR_OUT_SHARED_1,deriv_tile_sz,ofile,prefix="\t")

        ofile.write("\tconst unsigned int Lb = "+str(deriv_max_pad-deriv_req_pad)+";// load begin bound\n")
        ofile.write("\tconst unsigned int Le = sz[0]-"+str(deriv_max_pad-deriv_req_pad)+";// load end bound\n")


        # !! Note that we assume tile size are cubic.
        ofile.write("//!! Note that we assume tile size are cubic.\n")
        ofile.write("\tconst unsigned int BLK_ITERATIONS_X = ((Le-Lb)<"+VAR_TILE_SZ+"[0])? 1: ((int)ceil((double)(Le-Lb-"+VAR_TILE_SZ+"[0])/("+VAR_TILE_SZ+"[0]-2*" +str(deriv_req_pad)+")))+1;\n")
        ofile.write("\tconst unsigned int BLK_ITERATIONS_Y = BLK_ITERATIONS_X;\n")
        ofile.write("\tconst unsigned int BLK_ITERATIONS_Z = BLK_ITERATIONS_X;\n")
        ofile.write("\n")

        ofile.write("\tunsigned int "+VAR_TILE_LIMITS+"[3*2];\n")
        ofile.write("\tunsigned int "+VAR_TILE_LIMITS_STORE+"[3*2];\n")


        ofile.write("\tfor(unsigned int iter_z=0;iter_z<BLK_ITERATIONS_Z;iter_z++){\n\n")
        ofile.write("\t\t "+VAR_TILE_LIMITS+"[2*2+0]=max("+str(deriv_max_pad-deriv_req_pad)+",(int)("+str(deriv_max_pad-deriv_req_pad)+" + "+VAR_TILE_SZ+"[2]*iter_z -2*iter_z*"+str(deriv_req_pad)+"));\n")
        ofile.write("\t\t "+VAR_TILE_LIMITS+"[2*2+1]=min("+VAR_TILE_LIMITS+"[2*2+0]+"+VAR_TILE_SZ+"[2]-1,sz[2]-"+str(deriv_max_pad-deriv_req_pad)+ "-1);\n")
        ofile.write("\n")

        ofile.write("\n")
        ofile.write("\t\t if(("+VAR_TILE_LIMITS+"[5]-"+VAR_TILE_LIMITS+"[4]+1)<="+str(2*deriv_req_pad+3)+") \n\t\t  "+VAR_TILE_LIMITS+"[4]="+VAR_TILE_LIMITS+"[4]-("+str(2*deriv_req_pad+3)+"-("+VAR_TILE_LIMITS+"[5]-"+VAR_TILE_LIMITS+"[4]+1)) ; \n ")
        ofile.write("\n")

        ofile.write("\t  for(unsigned int iter_y=0;iter_y<BLK_ITERATIONS_Y;iter_y++){\n\n")
        ofile.write("\t\t "+VAR_TILE_LIMITS+"[2*1+0]=max("+str(deriv_max_pad-deriv_req_pad)+",(int)("+str(deriv_max_pad-deriv_req_pad)+" + "+VAR_TILE_SZ+"[1]*iter_y -2*iter_y*"+str(deriv_req_pad)+"));\n")
        ofile.write("\t\t "+VAR_TILE_LIMITS+"[2*1+1]=min("+VAR_TILE_LIMITS+"[2*1+0]+"+VAR_TILE_SZ+"[1]-1,sz[1]-"+str(deriv_max_pad-deriv_req_pad)+ "-1);\n")
        ofile.write("\n")

        ofile.write("\t\t if(("+VAR_TILE_LIMITS+"[3]-"+VAR_TILE_LIMITS+"[2]+1)<="+str(2*deriv_req_pad+3)+") \n\t\t  "+VAR_TILE_LIMITS+"[2]="+VAR_TILE_LIMITS+"[2]-("+str(2*deriv_req_pad+3)+"-("+VAR_TILE_LIMITS+"[3]-"+VAR_TILE_LIMITS+"[2]+1)) ; \n ")
        ofile.write("\n")

        ofile.write("\t    for(unsigned int iter_x=0;iter_x<BLK_ITERATIONS_X;iter_x++){\n")

        ofile.write("\t\t "+VAR_TILE_LIMITS+"[2*0+0]=max("+str(deriv_max_pad-deriv_req_pad)+",(int)("+str(deriv_max_pad-deriv_req_pad)+" + "+VAR_TILE_SZ+"[0]*iter_x -2*iter_x*"+str(deriv_req_pad)+"));\n")
        ofile.write("\t\t "+VAR_TILE_LIMITS+"[2*0+1]=min("+VAR_TILE_LIMITS+"[2*0+0]+"+VAR_TILE_SZ+"[0]-1,sz[0]-"+str(deriv_max_pad-deriv_req_pad)+ "-1);\n")

        ofile.write("\n")

        ofile.write("\t\t if(("+VAR_TILE_LIMITS+"[1]-"+VAR_TILE_LIMITS+"[0]+1)<="+str(2*deriv_req_pad+3)+") \n\t\t  "+VAR_TILE_LIMITS+"[0]="+VAR_TILE_LIMITS+"[0]-("+str(2*deriv_req_pad+3)+"-("+VAR_TILE_LIMITS+"[1]-"+VAR_TILE_LIMITS+"[0]+1)) ; \n ")
        ofile.write("\n")

        ofile.write("\n")
        ofile.write("\t\t //if(threadIdx.x ==0 && threadIdx.y==0 && threadIdx.z==0)\n")
        ofile.write("\t\t //printf(\" iter %d %d %d : threadid (%d,%d,%d) tile begin: (%d,%d,%d) tile end: (%d,%d,%d) \\n\",iter_x,iter_y,iter_z, threadIdx.x,threadIdx.y,threadIdx.z,ijk_lm[0],ijk_lm[2],ijk_lm[4],ijk_lm[1],ijk_lm[3],ijk_lm[5]);\n\n")
        ofile.write("\n")

        ofile.write("\t\t"+FUNC_LOAD_VAR+"(&"+VAR_UNZIP_IN+"[cuda::VAR::U_BETA0][offset],(double *) "+VAR_IN_SHARED+",(const unsigned int *) "+VAR_TILE_LIMITS+",(const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+",(const unsigned int *) "+VAR_TILE_SZ+");\n")
        ofile.write("\t\t"+FUNC_LOAD_VAR+"(&"+VAR_UNZIP_IN+"[cuda::VAR::U_BETA1][offset],(double *) "+VAR_OUT_SHARED_0+",(const unsigned int *) "+VAR_TILE_LIMITS+",(const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+",(const unsigned int *) "+VAR_TILE_SZ+");\n")
        ofile.write("\t\t"+FUNC_LOAD_VAR+"(&"+VAR_UNZIP_IN+"[cuda::VAR::U_BETA2][offset],(double *) "+VAR_OUT_SHARED_1+",(const unsigned int *) "+VAR_TILE_LIMITS+",(const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+",(const unsigned int *) "+VAR_TILE_SZ+");\n")
        ofile.write("\t\t__syncthreads();\n")
        ofile.write("\t\t"+FUNC_SIGN_EXT+"((double *)"+VAR_IN_SHARED+",(bool *) "+VAR_BETA0_BOOL+",(const unsigned int *) "+VAR_TILE_LIMITS+",(const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+",(const unsigned int *) "+VAR_TILE_SZ+");\n")
        ofile.write("\t\t"+FUNC_SIGN_EXT+"((double *)"+VAR_OUT_SHARED_0+",(bool *) "+VAR_BETA1_BOOL+",(const unsigned int *) "+VAR_TILE_LIMITS+",(const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+",(const unsigned int *) "+VAR_TILE_SZ+");\n")
        ofile.write("\t\t"+FUNC_SIGN_EXT+"((double *)"+VAR_OUT_SHARED_1+",(bool *) "+VAR_BETA2_BOOL+",(const unsigned int *) "+VAR_TILE_LIMITS+",(const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+",(const unsigned int *) "+VAR_TILE_SZ+");\n")
        ofile.write("\t\t__syncthreads();\n")

        for e in D:
            enumStr=VAR_ENUM[D.index(e)]
            ofile.write("\n")
            ofile.write("\t\t//load input data from global to shared memory\n")

            ofile.write("\t\t"+FUNC_LOAD_VAR+"(&"+VAR_UNZIP_IN+"["+enumStr+"][offset],(double *) "+VAR_IN_SHARED+",(const unsigned int *) "+VAR_TILE_LIMITS+",(const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+",(const unsigned int *) "+VAR_TILE_SZ+");\n")
            ofile.write("\t\t__syncthreads();\n")
            ofile.write("\t\t//sync to make sure all the data is loaded\n")

            for deriv in derivs:
                if((deriv.DerivType=="d") and (deriv.DerivDir=="x")):
                    ofile.write("\t\t// computing deriv "+deriv.DerivDir+" for variable "+e+" \n")
                    ofile.write("\t\t"+deriv.DerivFuncCall+"\n")
                    ofile.write("\t\t__syncthreads();\n") # not essential if each thread writes only the points it has computed in the block.
                    for deriv1 in derivs:
                        if((e in DD) and (deriv1.DerivType=="dd") and ((deriv1.DerivDir=="xy") or (deriv1.DerivDir=="xz"))):
                            ofile.write("\t\t// computing deriv "+deriv1.DerivDir+" for variable "+e+" \n")
                            ofile.write("\t\t"+deriv1.DerivFuncCall+"\n")
                            ofile.write("\t\t__syncthreads();\n") # not essential if each thread writes only the points it has computed in the block.

                            if(deriv1.DerivDir=="xy"):
                                computeTileStore("y",ofile,deriv_req_pad)
                            elif(deriv1.DerivDir=="xz"):
                                computeTileStore("z",ofile,deriv_req_pad)

                            #!!!! NOTE that for mixed derivs you need to store the padding region as well.
                            ofile.write("\t\t"+FUNC_STORE_VAR+"((double *) "+VAR_OUT_SHARED_1+",&("+VAR_DERIV_WORKSPACE+"->__"+deriv1.DerivOutput+"_"+ e +"[("+VAR_DERIV_WORKSPACE_OFFSET+")]),(const unsigned int *) "+VAR_TILE_LIMITS+",(const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+", (const unsigned int *) "+VAR_TILE_LIMITS_STORE+",(const unsigned int *) "+VAR_TILE_SZ+");\n")
                            if((deriv1.DerivDir=="xy")):
                                ofile.write("\t\t__syncthreads();\n")
                            ofile.write("\n")

                    computeTileStore("x",ofile,deriv_req_pad)
                    ofile.write("\t\t"+FUNC_STORE_VAR+"((double *) "+VAR_OUT_SHARED_0+",&("+VAR_DERIV_WORKSPACE+"->__"+deriv.DerivOutput+"_"+ e +"[("+VAR_DERIV_WORKSPACE_OFFSET+")]),(const unsigned int *) "+VAR_TILE_LIMITS+",(const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+", (const unsigned int *) "+VAR_TILE_LIMITS_STORE+",(const unsigned int *) "+VAR_TILE_SZ+");\n")
                    ofile.write("\t\t__syncthreads();\n")
                    ofile.write("\n")

                if((deriv.DerivType=="d") and (deriv.DerivDir=="y")):
                    ofile.write("\t\t// computing deriv "+deriv.DerivDir+" for variable "+e+" \n")
                    ofile.write("\t\t"+deriv.DerivFuncCall+"\n")
                    ofile.write("\t\t__syncthreads();\n") # not essential if each thread writes only the points it has computed in the block.
                    for deriv1 in derivs:
                        if((e in DD) and (deriv1.DerivType=="dd") and (deriv1.DerivDir=="yz")):
                            ofile.write("\t\t// computing deriv "+deriv1.DerivDir+" for variable "+e+" \n")
                            ofile.write("\t\t"+deriv1.DerivFuncCall+"\n")
                            ofile.write("\t\t__syncthreads();\n") # not essential if each thread writes only the points it has computed in the block.
                            computeTileStore("z",ofile,deriv_req_pad)
                            ofile.write("\t\t"+FUNC_STORE_VAR+"((double *) "+VAR_OUT_SHARED_1+",&("+VAR_DERIV_WORKSPACE+"->__"+deriv1.DerivOutput+"_"+ e +"[("+VAR_DERIV_WORKSPACE_OFFSET+")]),(const unsigned int *) "+VAR_TILE_LIMITS+",(const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+", (const unsigned int *) "+VAR_TILE_LIMITS_STORE+",(const unsigned int *) "+VAR_TILE_SZ+");\n")
                            #ofile.write("\t\t__syncthreads();\n")
                            ofile.write("\n")
                    #write the x, y,z derivs.

                    computeTileStore("y",ofile,deriv_req_pad)
                    ofile.write("\t\t"+FUNC_STORE_VAR+"((double *) "+VAR_OUT_SHARED_0+",&("+VAR_DERIV_WORKSPACE+"->__"+deriv.DerivOutput+"_"+ e +"[("+VAR_DERIV_WORKSPACE_OFFSET+")]),(const unsigned int *) "+VAR_TILE_LIMITS+",(const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+", (const unsigned int *) "+VAR_TILE_LIMITS_STORE+",(const unsigned int *) "+VAR_TILE_SZ+");\n")
                    ofile.write("\t\t__syncthreads();\n")
                    ofile.write("\n")

                if((deriv.DerivType=="d") and (deriv.DerivDir=="z")):
                    ofile.write("\t\t// computing deriv "+deriv.DerivDir+" for variable "+e+" \n")
                    ofile.write("\t\t"+deriv.DerivFuncCall+"\n")
                    ofile.write("\t\t__syncthreads();\n") # not essential if each thread writes only the points it has computed in the block.
                    computeTileStore("z",ofile,deriv_req_pad)
                    ofile.write("\t\t"+FUNC_STORE_VAR+"((double *) "+VAR_OUT_SHARED_0+",&("+VAR_DERIV_WORKSPACE+"->__"+deriv.DerivOutput+"_"+ e +"[("+VAR_DERIV_WORKSPACE_OFFSET+")]),(const unsigned int *) "+VAR_TILE_LIMITS+",(const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+", (const unsigned int *) "+VAR_TILE_LIMITS_STORE+",(const unsigned int *) "+VAR_TILE_SZ+");\n")
                    ofile.write("\t\t__syncthreads();\n")

            ofile.write("\n")

            derivCount=0

            for deriv in derivs:
                if((e in DD) and (deriv.DerivType=="dd")  and ((deriv.DerivDir=="xx") or (deriv.DerivDir=="yy") or (deriv.DerivDir=="zz"))):
                    ofile.write("\t\t"+deriv.DerivFuncCall+"\n")
                    ofile.write("\t\t__syncthreads();\n") # not essential if each thread writes only the points it has computed in the block.
                    if(deriv.DerivDir=="xx"):
                        computeTileStore("x",ofile,deriv_req_pad)
                    elif(deriv.DerivDir=="yy"):
                        computeTileStore("y",ofile,deriv_req_pad)
                    elif(deriv.DerivDir=="zz"):
                        computeTileStore("z",ofile,deriv_req_pad)
                    ofile.write("\t\t"+FUNC_STORE_VAR+"((double *) "+VAR_OUT_SHARED_0+",&("+VAR_DERIV_WORKSPACE+"->__"+deriv.DerivOutput+"_"+ e +"[("+VAR_DERIV_WORKSPACE_OFFSET+")]),(const unsigned int *) "+VAR_TILE_LIMITS+",(const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+", (const unsigned int *) "+VAR_TILE_LIMITS_STORE+",(const unsigned int *) "+VAR_TILE_SZ+");\n")
                    ofile.write("\t\t__syncthreads();\n")

            ofile.write("\n")

            for deriv in derivs:
                if( (e in KO) and (deriv.DerivType=="ko")):
                    ofile.write("\t\t"+deriv.DerivFuncCall+"\n")
                    ofile.write("\t\t__syncthreads();\n") # not essential if each thread writes only the points it has computed in the block.
                    if(deriv.DerivDir=="x"):
                        computeTileStore("x",ofile,deriv_req_pad)
                    elif(deriv.DerivDir=="y"):
                        computeTileStore("y",ofile,deriv_req_pad)
                    elif(deriv.DerivDir=="z"):
                        computeTileStore("z",ofile,deriv_req_pad)
                    ofile.write("\t\t"+FUNC_STORE_VAR+"((double *) "+VAR_OUT_SHARED_0+",&("+VAR_DERIV_WORKSPACE+"->__"+deriv.DerivOutput+"_"+ e +"[("+VAR_DERIV_WORKSPACE_OFFSET+")]),(const unsigned int *) "+VAR_TILE_LIMITS+",(const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+", (const unsigned int *) "+VAR_TILE_LIMITS_STORE+",(const unsigned int *) "+VAR_TILE_SZ+");\n")
                    ofile.write("\t\t__syncthreads();\n")

            for deriv in derivs:
                if( (e in AD) and (deriv.DerivType=="ad")):
                    ofile.write("\t\t"+deriv.DerivFuncCall+"\n")
                    ofile.write("\t\t__syncthreads();\n") # not essential if each thread writes only the points it has computed in the block.
                    if(deriv.DerivDir=="x"):
                        computeTileStore("x",ofile,deriv_req_pad)
                    elif(deriv.DerivDir=="y"):
                        computeTileStore("y",ofile,deriv_req_pad)
                    elif(deriv.DerivDir=="z"):
                        computeTileStore("z",ofile,deriv_req_pad)
                    ofile.write("\t\t"+FUNC_STORE_VAR+"((double *) "+VAR_OUT_SHARED_0+",&("+VAR_DERIV_WORKSPACE+"->__"+deriv.DerivOutput+"_"+ e +"[("+VAR_DERIV_WORKSPACE_OFFSET+")]),(const unsigned int *) "+VAR_TILE_LIMITS+",(const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+", (const unsigned int *) "+VAR_TILE_LIMITS_STORE+",(const unsigned int *) "+VAR_TILE_SZ+");\n")
                    ofile.write("\t\t__syncthreads();\n")



        ofile.write("\t\t  } // end of block tile loop x\n")
        ofile.write("\t\t } // end of block tile loop y\n")
        ofile.write("\t\t} // end of block tile loop z\n\n")

        ofile.write("} // end of function "+FUNC_DERIV_COMP+"\n\n")


        ##############################################################################
        ## RHS code generation
        ##############################################################################

        for var_id in range(0,len(varnames)):

            memManager.deallocAll()
            memManager.clearScopeVariables()

            ofile.write("/**@brief compute "+varnames[var_id]+" \n")
            ofile.write(" @param[in] "+VAR_UNZIP_IN+": unzipped input array (global memory) \n")
            ofile.write(" @param[in] "+TYPE_DERIV_STRUCT+": allocated workspace for derivative computations \n")
            ofile.write(" @param[in] "+VAR_DENDRO_BLK_LIST+": dendro block list \n")
            ofile.write(" @param[in] "+VAR_GPU_BLOCK_MAP+": gpu block map  \n")
            ofile.write(" @param[in] "+VAR_CUDA_DEVICE+": cuda device properties  \n")
            ofile.write(" @param[out] "+VAR_UNZIP_OUT+": unzipped output computed rhs  \n")
            ofile.write("*/ \n")
            ofile.write("__device__ void "+FUNC_COMP_RHS_PRE+"_"+varnames[var_id]+"(double **"+VAR_UNZIP_OUT+", const double**"+VAR_UNZIP_IN+","+TYPE_DERIV_STRUCT+"* __derivWorkspace, const "+ TYPE_BLK_CU+ "* dblock, const unsigned int * "+VAR_GPU_BLOCK_MAP+",const "+TYPE_BSSN_PARAMS+" * "+VAR_BSSN_PARAMS+",const cudaDeviceProp* "+VAR_CUDA_DEVICE+", "+TYPE_SHARED_MEM+"* "+VAR_SHARED_MEM+", unsigned int "+VAR_STREAM_ID+"){\n")
            ofile.write("\n\n")


            mi = [0, 1, 2, 4, 5, 8]
            midx = ['00', '01', '02', '11', '12', '22']

            ofile.write("\n")
            idx="[pp]"

            ofile.write("\t///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n")
            ofile.write("\t//                             generated code for "+varnames[var_id]+"              begin   \n")
            ofile.write("\t///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n")
            varOut=varnames[var_id]
            exp=outs[var_id]

            ofile.write("\t// bssn compute parameters \n")

            ofile.write("\tconst double lambda[4]={"+VAR_BSSN_PARAMS+"->BSSN_LAMBDA[0],"+VAR_BSSN_PARAMS+"->BSSN_LAMBDA[1],"+VAR_BSSN_PARAMS+"->BSSN_LAMBDA[2],"+VAR_BSSN_PARAMS+"->BSSN_LAMBDA[3]};\n")
            ofile.write("\tconst double lambda_f[2]={"+VAR_BSSN_PARAMS+"->BSSN_LAMBDA_F[0],"+VAR_BSSN_PARAMS+"->BSSN_LAMBDA_F[1]};\n")
            ofile.write("\tconst double kosigma="+VAR_BSSN_PARAMS+"->KO_DISS_SIGMA;\n")
            ofile.write("\tconst double ETA_R0="+VAR_BSSN_PARAMS+"->ETA_R0;\n")
            ofile.write("\tconst double R0="+VAR_BSSN_PARAMS+"->ETA_R0;\n")
            ofile.write("\tconst double ETA_DAMPING="+VAR_BSSN_PARAMS+"->ETA_DAMPING;\n")
            ofile.write("\tconst double ETA_DAMPING_EXP="+VAR_BSSN_PARAMS+"->ETA_DAMPING_EXP;\n")
            ofile.write("\tconst double ETA_CONST="+VAR_BSSN_PARAMS+"->ETA_CONST;\n")
            ofile.write("\tconst double eta_power[2]={"+VAR_BSSN_PARAMS+"->BSSN_ETA_POWER[0],"+VAR_BSSN_PARAMS+"->BSSN_ETA_POWER[1]};\n")

            print("code generation for : "+varOut)


            ofile.write("\tconst unsigned int NUM_SM_UNITS="+VAR_CUDA_DEVICE+"->multiProcessorCount;\n")
            ofile.write("\tconst unsigned int SM_ID=get_smid();//"+VAR_BLK_ID_X+"%NUM_SM_UNITS;\n")
            ofile.write("\tconst unsigned int offset=dblock->getOffset();\n")
            ofile.write("\tconst unsigned int *sz=dblock->getSz();\n")
            ofile.write("\tconst unsigned int *"+VAR_DENDRO_BLK_ALIGNED_SZ+"=dblock->getAlignedSz();\n")
            ofile.write("\tconst double* hx=dblock->getDx();\n")

            ofile.write("\tconst double dx=hx[0];\n")
            ofile.write("\tconst double dy=hx[1];\n")
            ofile.write("\tconst double dz=hx[2];\n")

            ofile.write("\tconst double* ptmin=dblock->getPtMin();\n")
            ofile.write("\tconst double* ptmax=dblock->getPtMax();\n")
            ofile.write("\tconst unsigned int bflag=dblock->getBFlag();\n")

            ofile.write("\n")

            num_e = 0
            lexp = []
            lname = []

            if type(exp) == list:
                num_e = num_e + len(exp)
                for j, ev in enumerate(exp):
                    lexp.append(ev)
                    lname.append(varOut+repr(j)+idx)
            elif type(exp) == sympy.Matrix:
                num_e = num_e + len(exp)
                for j, k in enumerate(mi):
                    lexp.append(exp[k])
                    lname.append(varOut+midx[j]+idx)
            else:
                num_e = num_e + 1
                lexp.append(exp)
                lname.append(varOut+idx)

            print("cse tree build begin")
            ee_name = 'DENDRO_' #''.join(random.choice(string.ascii_uppercase) for _ in range(5))
            ee_syms = sympy.utilities.numbered_symbols(prefix=ee_name)
            _v = sympy.cse(lexp, symbols=ee_syms, optimizations='basic')

            print("cse tree build completed")

            # bssn variables needed for rhs computation.
            bssnInputVars=[]

            # bssn variables output
            bssnOutputVars=[]

            # derivative variables needed for rhs computation
            derivVars=[]

            # staged bssn variables.
            bssnStagedVars=[]

            if type(exp) == list:
                #print("list \n")
                for j, ev in enumerate(exp):
                    regm=re.findall(re.compile(r"([A-Z,a-z,0-9,_]*\[pp\])"),dendro.change_deriv_names(str(ev)))
                    for varDep in regm:
                        if varDep[0:-4] in VAR_ENUM_TO_INPUT_SYM.keys():
                            bssnInputVars.append(varDep[0:-4])
                        elif varDep[0:-4] in VAR_ENUM_TO_OUTPUT_SYM.keys():
                            bssnOutputVars.append(varDep[0:-4])
                        else:
                            for key,value in custom_functions.items():
                                if value in varDep[0:-4]:
                                    derivVars.append(varDep[0:-4])
                                    break


            elif type(exp)==sympy.Matrix:
                #print(dendro.change_deriv_names(str(exp)))
                #print(exp.free_symbols)
                regm=re.findall(re.compile(r"([A-Z,a-z,0-9,_]*\[pp\])"),dendro.change_deriv_names(str(exp)))
                for varDep in regm:
                    if varDep[0:-4] in VAR_ENUM_TO_INPUT_SYM.keys():
                        bssnInputVars.append(varDep[0:-4])
                    elif varDep[0:-4] in VAR_ENUM_TO_OUTPUT_SYM.keys():
                        bssnOutputVars.append(varDep[0:-4])
                    else:
                        for key,value in custom_functions.items():
                            if value in varDep[0:-4]:
                                derivVars.append(varDep[0:-4])
                                break


            else:
                #print(dendro.change_deriv_names(str(exp)))
                regm=re.findall(re.compile(r"([A-Z,a-z,0-9,_]*\[pp\])"),dendro.change_deriv_names(str(exp)))
                for varDep in regm:
                    #print (varDep[0:-4])
                    if varDep[0:-4] in VAR_ENUM_TO_INPUT_SYM.keys():
                        bssnInputVars.append(varDep[0:-4])
                    elif varDep[0:-4] in VAR_ENUM_TO_OUTPUT_SYM.keys():
                        bssnOutputVars.append(varDep[0:-4])
                    else:
                        for key,value in custom_functions.items():
                            if value in varDep[0:-4]:
                                derivVars.append(varDep[0:-4])
                                break


            for lvar in lname:
                if lvar[0:-4] in VAR_ENUM_TO_OUTPUT_SYM.keys():
                    bssnOutputVars.append(lvar[0:-4])
                else:
                    bssnStagedVars.append(lvar[0:-4])

            bssnInputVars=list(set(bssnInputVars))
            bssnOutputVars=list(set(bssnOutputVars))
            bssnStagedVars=list(set(bssnStagedVars))
            derivVars=list(set(derivVars))


            total_dep=len(bssnInputVars)+len(bssnStagedVars)+len(derivVars)+len(bssnOutputVars)
            rhs_tile_size_1d=math.floor(((memManager.getMemUsable())/(total_dep*8))**(1.0/3.0))

            ofile.write("\tconst unsigned int "+VAR_TILE_SZ+"[3]={"+str(rhs_tile_size_1d)+","+str(rhs_tile_size_1d)+","+str(rhs_tile_size_1d)+"};\n")
            rhs_tile_size=rhs_tile_size_1d**3

            ofile.write("\t\n")
            # no padding region required for rhs computation
            rhs_req_pad=0

            ofile.write("\t //input vars begin\n")
            for var in bssnInputVars:
                memManager.malloc(var,rhs_tile_size,ofile,prefix="\t")
            ofile.write("\t //input vars end\n")

            ofile.write("\t // staged vars begin\n")
            for var in bssnStagedVars:
                memManager.malloc(var,rhs_tile_size,ofile,prefix="\t")

            ofile.write("\t // staged vars end\n")

            ofile.write("\t // deriv vars begin\n")
            for var in derivVars:
                memManager.malloc(var,rhs_tile_size,ofile,prefix="\t")
            ofile.write("\t // deriv vars end\n")

            ofile.write("\t // output vars begin\n")
            for var in bssnOutputVars:
                memManager.malloc(var,rhs_tile_size,ofile,prefix="\t")
            ofile.write("\t // output vars end\n")


            ofile.write("\tconst unsigned int Lb = "+str(deriv_max_pad-rhs_req_pad)+";// load begin bound\n")
            ofile.write("\tconst unsigned int Le = sz[0]-"+str(deriv_max_pad-rhs_req_pad)+";// load end bound\n")

            # !! Note that we assume tile size are cubic.
            ofile.write("//!! Note that we assume tile size are cubic.\n")
            ofile.write("\tconst unsigned int BLK_ITERATIONS_X = ((Le-Lb)<"+VAR_TILE_SZ+"[0])? 1: ((int)ceil((double)(Le-Lb-"+VAR_TILE_SZ+"[0])/("+VAR_TILE_SZ+"[0]-2*" +str(rhs_req_pad)+")))+1;\n")
            ofile.write("\tconst unsigned int BLK_ITERATIONS_Y = BLK_ITERATIONS_X;\n")
            ofile.write("\tconst unsigned int BLK_ITERATIONS_Z = BLK_ITERATIONS_X;\n")
            ofile.write("\n")

            ofile.write("\tunsigned int "+VAR_TILE_LIMITS+"[3*2];\n")
            ofile.write("\tunsigned int "+VAR_TILE_LIMITS_STORE+"[3*2];\n")

            ofile.write("\tfor(unsigned int iter_z=0;iter_z<BLK_ITERATIONS_Z;iter_z++){\n\n")
            ofile.write("\t\t "+VAR_TILE_LIMITS+"[2*2+0]=max("+str(deriv_max_pad-rhs_req_pad)+",(int)("+str(deriv_max_pad-rhs_req_pad)+" + "+VAR_TILE_SZ+"[2]*iter_z -2*iter_z*"+str(rhs_req_pad)+"));\n")
            ofile.write("\t\t "+VAR_TILE_LIMITS+"[2*2+1]=min("+VAR_TILE_LIMITS+"[2*2+0]+"+VAR_TILE_SZ+"[2]-1,sz[2]-"+str(deriv_max_pad-rhs_req_pad)+"-1);\n")
            ofile.write("\n")

            ofile.write("\t  for(unsigned int iter_y=0;iter_y<BLK_ITERATIONS_Y;iter_y++){\n\n")
            ofile.write("\t\t "+VAR_TILE_LIMITS+"[2*1+0]=max("+str(deriv_max_pad-rhs_req_pad)+",(int)("+str(deriv_max_pad-rhs_req_pad)+" + "+VAR_TILE_SZ+"[1]*iter_y -2*iter_y*"+str(rhs_req_pad)+"));\n")
            ofile.write("\t\t "+VAR_TILE_LIMITS+"[2*1+1]=min("+VAR_TILE_LIMITS+"[2*1+0]+"+VAR_TILE_SZ+"[1]-1,sz[1]-"+str(deriv_max_pad-rhs_req_pad)+"-1);\n")
            ofile.write("\n")

            ofile.write("\t    for(unsigned int iter_x=0;iter_x<BLK_ITERATIONS_X;iter_x++){\n")

            ofile.write("\t\t "+VAR_TILE_LIMITS+"[2*0+0]=max("+str(deriv_max_pad-rhs_req_pad)+",(int)("+str(deriv_max_pad-rhs_req_pad)+" + "+VAR_TILE_SZ+"[0]*iter_x -2*iter_x*"+str(rhs_req_pad)+"));\n")
            ofile.write("\t\t "+VAR_TILE_LIMITS+"[2*0+1]=min("+VAR_TILE_LIMITS+"[2*0+0]+"+VAR_TILE_SZ+"[0]-1,sz[0]-"+str(deriv_max_pad-rhs_req_pad)+"-1);\n")
            ofile.write("\n")
            ofile.write("\n")

            ofile.write("\t\t"+VAR_TILE_LIMITS_STORE+"[0]=0;\n")
            ofile.write("\t\t"+VAR_TILE_LIMITS_STORE+"[1]="+VAR_TILE_LIMITS+"[1] - "+VAR_TILE_LIMITS+"[0];\n")
            ofile.write("\n")
            ofile.write("\t\t"+VAR_TILE_LIMITS_STORE+"[2]=0;\n")
            ofile.write("\t\t"+VAR_TILE_LIMITS_STORE+"[3]="+VAR_TILE_LIMITS+"[3] - "+VAR_TILE_LIMITS+"[2];\n")
            ofile.write("\n")
            ofile.write("\t\t"+VAR_TILE_LIMITS_STORE+"[4]=0;\n")
            ofile.write("\t\t"+VAR_TILE_LIMITS_STORE+"[5]="+VAR_TILE_LIMITS+"[5] - "+VAR_TILE_LIMITS+"[4];\n")

            ofile.write("\t\t //if(threadIdx.x ==0 && threadIdx.y==0 && threadIdx.z==0)\n")
            ofile.write("\t\t //printf(\" iter %d %d %d : threadid (%d,%d,%d) tile begin: (%d,%d,%d) tile end: (%d,%d,%d) \\n\",iter_x,iter_y,iter_z, threadIdx.x,threadIdx.y,threadIdx.z,ijk_lm[0],ijk_lm[2],ijk_lm[4],ijk_lm[1],ijk_lm[3],ijk_lm[5]);\n\n")

            ofile.write("\n\n")


            ofile.write("\t\t //load data from global to shared memory\n")
            for var in bssnInputVars:
                ofile.write("\t\t "+FUNC_LOAD_VAR+"(&"+VAR_UNZIP_IN+"["+VAR_ENUM_TO_INPUT_SYM[var]+"][offset],(double *) "+var+",(const unsigned int *) "+VAR_TILE_LIMITS+",(const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+",(const unsigned int *) "+VAR_TILE_SZ+");\n")

            for var in derivVars:
                ofile.write("\t\t "+FUNC_LOAD_VAR+"(&("+VAR_DERIV_WORKSPACE+"->__"+var+"[("+VAR_DERIV_WORKSPACE_OFFSET+")]),(double *) "+var+",(const unsigned int *) "+VAR_TILE_LIMITS+",(const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+",(const unsigned int *) "+VAR_TILE_SZ+");\n")

            ofile.write("\t\t __syncthreads();\n\n")

            ofile.write("\n\n")
            #/*|| ("+VAR_TRD_ID_Z+">=("+VAR_TILE_LIMITS+"[5]-"+VAR_TILE_LIMITS+"[4]))*/
            ofile.write("\tif(!(("+VAR_TRD_ID_X+">("+VAR_TILE_LIMITS+"[1]-"+VAR_TILE_LIMITS+"[0])) || ("+VAR_TRD_ID_Y+">("+VAR_TILE_LIMITS+"[3]-"+VAR_TILE_LIMITS+"[2]))) ){ \n\n")

            ofile.write("\t\t double x,y,z,r_coord,eta;\n")
            ofile.write("\t\t unsigned int pp=0*"+VAR_TILE_SZ+"[0]*"+VAR_TILE_SZ+"[1]+"+VAR_TRD_ID_Y+"*"+VAR_TILE_SZ+"[1]+"+VAR_TRD_ID_X+";\n")
            ofile.write("\t\t  for(unsigned int k=0;k<=(ijk_lm[5]-ijk_lm[4]);++k,pp+="+VAR_TILE_SZ+"[0]*"+VAR_TILE_SZ+"[1]){\n")

            ofile.write("\t\t\t  z = ptmin[2] + (k+"+VAR_TILE_LIMITS+"[4])*dz;\n")
            ofile.write("\t\t\t  y = ptmin[1] + ("+VAR_TRD_ID_Y+"+"+VAR_TILE_LIMITS+"[2])*dy;\n")
            ofile.write("\t\t\t  x = ptmin[0] + ("+VAR_TRD_ID_X+"+"+VAR_TILE_LIMITS+"[0])*dx;\n")
            ofile.write("\t\t\t  r_coord = sqrt(x*x + y*y + z*z);\n")
            ofile.write("\t\t\t  eta=ETA_CONST;\n")
            ofile.write("\t\t\t  if (r_coord >= ETA_R0) {\n")
            ofile.write("\t\t\t     eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);\n")
            ofile.write("\t\t\t  }\n\n")



            ofile.write("\t\t\t  // Dendro: {{{ \n")
            ofile.write("\t\t\t  // Dendro: original ops: "+str(sympy.count_ops(lexp))+"\n")

            rops=0
            ofile.write("\t\t\t     // Dendro: printing temp variables\n")
            for (v1, v2) in _v[0]:
                ofile.write("\t\t const double ")
                ofile.write(dendro.change_deriv_names(sympy.ccode(v2, assign_to=v1, user_functions=custom_functions))+"\n")
                rops = rops + sympy.count_ops(v2)

            ofile.write("\t\t\t      // Dendro: printing variables\n\n")
            for i, e in enumerate(_v[1]):
                ofile.write("\t\t   "+dendro.change_deriv_names(sympy.ccode(e, assign_to=lname[i], user_functions=custom_functions))+"\n")
                rops = rops + sympy.count_ops(e)


            ofile.write("\t\t\t      // Dendro: reduced ops: "+str(rops)+"\n")
            ofile.write("\t\t\t      // Dendro: }}} \n")
            ofile.write("\t\t\t     } //loop z end \n")

            ofile.write("\t}// end of the if for the thread idx \n")

            ofile.write("\t\t\t__syncthreads();\n\n")

            ofile.write("\t\t\t// sotre computed variables\n\n")
            for var in bssnOutputVars:
                ofile.write("\t\t"+FUNC_STORE_VAR+"("+var+", &"+VAR_UNZIP_OUT+"["+VAR_ENUM_TO_OUTPUT_SYM[var]+"][offset],(const unsigned int *) "+VAR_TILE_LIMITS+",(const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+", (const unsigned int *) "+VAR_TILE_LIMITS_STORE+",(const unsigned int *) "+VAR_TILE_SZ+");\n")

            ofile.write("\t\t   __syncthreads();\n")

            ofile.write("\t  } // end of block assigned to gpu block loop x \n\n")
            ofile.write("\t } // end of block assigned to gpu block loop y \n\n")
            ofile.write("\t} // end of block assigned to gpu block loop z \n\n")

            ofile.write("} // end of function" +FUNC_COMP_RHS_PRE+"_"+varnames[var_id]+" \n\n")


        memManager.deallocAll()
        memManager.clearScopeVariables()

        ofile.write("/**@brief apply KO dissipation \n")
        ofile.write(" @param[in] "+VAR_UNZIP_IN+": unzipped input array (global memory) \n")
        ofile.write(" @param[in] "+TYPE_DERIV_STRUCT+": allocated workspace for derivative computations \n")
        ofile.write(" @param[in] "+VAR_DENDRO_BLK_LIST+": dendro block list \n")
        ofile.write(" @param[in] "+VAR_GPU_BLOCK_MAP+": gpu block map  \n")
        ofile.write(" @param[in] "+VAR_CUDA_DEVICE+": cuda device properties  \n")
        ofile.write(" @param[out] "+VAR_UNZIP_OUT+": unzipped output computed rhs  \n")
        ofile.write("*/ \n")
        ofile.write("__device__ void "+FUNC_KO_DISS+"(double **"+VAR_UNZIP_OUT+", const double**"+VAR_UNZIP_IN+","+TYPE_DERIV_STRUCT+"* __derivWorkspace, const "+ TYPE_BLK_CU+ "* dblock, const unsigned int * "+VAR_GPU_BLOCK_MAP+",const "+TYPE_BSSN_PARAMS+" * "+VAR_BSSN_PARAMS+",const cudaDeviceProp* "+VAR_CUDA_DEVICE+", "+TYPE_SHARED_MEM+"* "+VAR_SHARED_MEM+", unsigned int "+VAR_STREAM_ID+"){\n")
        ofile.write("\n")

        ofile.write("\t// bssn compute parameters \n")

        ofile.write("\tconst double lambda[4]={"+VAR_BSSN_PARAMS+"->BSSN_LAMBDA[0],"+VAR_BSSN_PARAMS+"->BSSN_LAMBDA[1],"+VAR_BSSN_PARAMS+"->BSSN_LAMBDA[2],"+VAR_BSSN_PARAMS+"->BSSN_LAMBDA[3]};\n")
        ofile.write("\tconst double lambda_f[2]={"+VAR_BSSN_PARAMS+"->BSSN_LAMBDA_F[0],"+VAR_BSSN_PARAMS+"->BSSN_LAMBDA_F[1]};\n")
        ofile.write("\tconst double kosigma="+VAR_BSSN_PARAMS+"->KO_DISS_SIGMA;\n")
        ofile.write("\tconst double ETA_R0="+VAR_BSSN_PARAMS+"->ETA_R0;\n")
        ofile.write("\tconst double R0="+VAR_BSSN_PARAMS+"->ETA_R0;\n")
        ofile.write("\tconst double ETA_DAMPING="+VAR_BSSN_PARAMS+"->ETA_DAMPING;\n")
        ofile.write("\tconst double ETA_DAMPING_EXP="+VAR_BSSN_PARAMS+"->ETA_DAMPING_EXP;\n")
        ofile.write("\tconst double ETA_CONST="+VAR_BSSN_PARAMS+"->ETA_CONST;\n")
        ofile.write("\tconst double eta_power[2]={"+VAR_BSSN_PARAMS+"->BSSN_ETA_POWER[0],"+VAR_BSSN_PARAMS+"->BSSN_ETA_POWER[1]};\n")


        ofile.write("\tconst unsigned int NUM_SM_UNITS="+VAR_CUDA_DEVICE+"->multiProcessorCount;\n")
        ofile.write("\tconst unsigned int SM_ID=get_smid();//"+VAR_BLK_ID_X+"%NUM_SM_UNITS;\n")
        ofile.write("\tconst unsigned int offset=dblock->getOffset();\n")
        ofile.write("\tconst unsigned int *sz=dblock->getSz();\n")
        ofile.write("\tconst unsigned int *"+VAR_DENDRO_BLK_ALIGNED_SZ+"=dblock->getAlignedSz();\n")
        ofile.write("\tconst double* hx=dblock->getDx();\n")

        ofile.write("\tconst double dx=hx[0];\n")
        ofile.write("\tconst double dy=hx[1];\n")
        ofile.write("\tconst double dz=hx[2];\n")

        ofile.write("\tconst double* ptmin=dblock->getPtMin();\n")
        ofile.write("\tconst double* ptmax=dblock->getPtMax();\n")
        ofile.write("\tconst unsigned int bflag=dblock->getBFlag();\n")

        total_dep=4
        rhs_req_pad=0
        rhs_tile_size_1d=math.floor(((memManager.getMemUsable())/(total_dep*8))**(1.0/3.0))



        ofile.write("\tconst unsigned int "+VAR_TILE_SZ+"[3]={"+str(rhs_tile_size_1d)+","+str(rhs_tile_size_1d)+","+str(rhs_tile_size_1d)+"};\n")
        rhs_tile_size=rhs_tile_size_1d**3

        VAR_KO_TEMP=["kograd_0","kograd_1","kograd_2"]
        VAR_KO_TEMP_RHS=["unZipSharedOut"]

        for var in VAR_KO_TEMP:
            memManager.malloc(var,rhs_tile_size,ofile,prefix="\t")

        for var in VAR_KO_TEMP_RHS:
            memManager.malloc(var,rhs_tile_size,ofile,prefix="\t")


        ofile.write("\tconst unsigned int Lb = "+str(deriv_max_pad-rhs_req_pad)+";// load begin bound\n")
        ofile.write("\tconst unsigned int Le = sz[0]-"+str(deriv_max_pad-rhs_req_pad)+";// load end bound\n")

        # !! Note that we assume tile size are cubic.
        ofile.write("//!! Note that we assume tile size are cubic.\n")
        ofile.write("\tconst unsigned int BLK_ITERATIONS_X = ((Le-Lb)<"+VAR_TILE_SZ+"[0])? 1: ((int)ceil((double)(Le-Lb-"+VAR_TILE_SZ+"[0])/("+VAR_TILE_SZ+"[0]-2*" +str(rhs_req_pad)+")))+1;\n")
        ofile.write("\tconst unsigned int BLK_ITERATIONS_Y = BLK_ITERATIONS_X;\n")
        ofile.write("\tconst unsigned int BLK_ITERATIONS_Z = BLK_ITERATIONS_X;\n")
        ofile.write("\n")

        ofile.write("\tunsigned int "+VAR_TILE_LIMITS+"[3*2];\n")
        ofile.write("\tunsigned int "+VAR_TILE_LIMITS_STORE+"[3*2];\n")

        ofile.write("\tfor(unsigned int iter_z=0;iter_z<BLK_ITERATIONS_Z;iter_z++){\n\n")
        ofile.write("\t\t "+VAR_TILE_LIMITS+"[2*2+0]=max("+str(deriv_max_pad-rhs_req_pad)+",(int)("+str(deriv_max_pad-rhs_req_pad)+" + "+VAR_TILE_SZ+"[2]*iter_z -2*iter_z*"+str(rhs_req_pad)+"));\n")
        ofile.write("\t\t "+VAR_TILE_LIMITS+"[2*2+1]=min("+VAR_TILE_LIMITS+"[2*2+0]+"+VAR_TILE_SZ+"[2]-1,sz[2]-"+str(deriv_max_pad-rhs_req_pad)+"-1);\n")
        ofile.write("\n")

        ofile.write("\t  for(unsigned int iter_y=0;iter_y<BLK_ITERATIONS_Y;iter_y++){\n\n")
        ofile.write("\t\t "+VAR_TILE_LIMITS+"[2*1+0]=max("+str(deriv_max_pad-rhs_req_pad)+",(int)("+str(deriv_max_pad-rhs_req_pad)+" + "+VAR_TILE_SZ+"[1]*iter_y -2*iter_y*"+str(rhs_req_pad)+"));\n")
        ofile.write("\t\t "+VAR_TILE_LIMITS+"[2*1+1]=min("+VAR_TILE_LIMITS+"[2*1+0]+"+VAR_TILE_SZ+"[1]-1,sz[1]-"+str(deriv_max_pad-rhs_req_pad)+"-1);\n")
        ofile.write("\n")

        ofile.write("\t    for(unsigned int iter_x=0;iter_x<BLK_ITERATIONS_X;iter_x++){\n")

        ofile.write("\t\t "+VAR_TILE_LIMITS+"[2*0+0]=max("+str(deriv_max_pad-rhs_req_pad)+",(int)("+str(deriv_max_pad-rhs_req_pad)+" + "+VAR_TILE_SZ+"[0]*iter_x -2*iter_x*"+str(rhs_req_pad)+"));\n")
        ofile.write("\t\t "+VAR_TILE_LIMITS+"[2*0+1]=min("+VAR_TILE_LIMITS+"[2*0+0]+"+VAR_TILE_SZ+"[0]-1,sz[0]-"+str(deriv_max_pad-rhs_req_pad)+"-1);\n")
        ofile.write("\n")
        ofile.write("\n")

        ofile.write("\t\t"+VAR_TILE_LIMITS_STORE+"[0]=0;\n")
        ofile.write("\t\t"+VAR_TILE_LIMITS_STORE+"[1]="+VAR_TILE_LIMITS+"[1] - "+VAR_TILE_LIMITS+"[0];\n")
        ofile.write("\n")
        ofile.write("\t\t"+VAR_TILE_LIMITS_STORE+"[2]=0;\n")
        ofile.write("\t\t"+VAR_TILE_LIMITS_STORE+"[3]="+VAR_TILE_LIMITS+"[3] - "+VAR_TILE_LIMITS+"[2];\n")
        ofile.write("\n")
        ofile.write("\t\t"+VAR_TILE_LIMITS_STORE+"[4]=0;\n")
        ofile.write("\t\t"+VAR_TILE_LIMITS_STORE+"[5]="+VAR_TILE_LIMITS+"[5] - "+VAR_TILE_LIMITS+"[4];\n")

        ofile.write("\t\t //if(threadIdx.x ==0 && threadIdx.y==0 && threadIdx.z==0)\n")
        ofile.write("\t\t //printf(\" iter %d %d %d : threadid (%d,%d,%d) tile begin: (%d,%d,%d) tile end: (%d,%d,%d) \\n\",iter_x,iter_y,iter_z, threadIdx.x,threadIdx.y,threadIdx.z,ijk_lm[0],ijk_lm[2],ijk_lm[4],ijk_lm[1],ijk_lm[3],ijk_lm[5]);\n\n")

        ofile.write("\n\n")

        ofile.write("\t\t unsigned int pp;\n")
        for var_id in range(0,len(varnames)):

            mi = [0, 1, 2, 4, 5, 8]
            midx = ['00', '01', '02', '11', '12', '22']

            varOut=varnames[var_id]
            exp=outs[var_id]

            num_e = 0
            lexp = []
            lname = []
            idx=""

            if type(exp) == list:
                num_e = num_e + len(exp)
                for j, ev in enumerate(exp):
                    lexp.append(ev)
                    lname.append(varOut+repr(j)+idx)
            elif type(exp) == sympy.Matrix:
                num_e = num_e + len(exp)
                for j, k in enumerate(mi):
                    lexp.append(exp[k])
                    lname.append(varOut+midx[j]+idx)
            else:
                num_e = num_e + 1
                lexp.append(exp)
                lname.append(varOut+idx)

            for rhs in lname:
                var_d=list(VAR_ENUM_TO_INPUT_SYM.keys())[list(VAR_ENUM_TO_INPUT_SYM.values()).index(VAR_ENUM_TO_OUTPUT_SYM[rhs])]

                ofile.write("\t\t //ko dissipation for variable "+var_d+"\n\n")
                for var in VAR_KO_TEMP:
                    ofile.write("\t\t "+FUNC_LOAD_VAR+"(&("+VAR_DERIV_WORKSPACE+"->__"+var+"_"+var_d+"[("+VAR_DERIV_WORKSPACE_OFFSET+")]),(double *) "+var+",(const unsigned int *) "+VAR_TILE_LIMITS+",(const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+",(const unsigned int *) "+VAR_TILE_SZ+");\n")

                #ofile.write("\t\t "+FUNC_LOAD_VAR+"(&"+VAR_UNZIP_IN+"["+VAR_ENUM_TO_INPUT_SYM[var_d]+"][offset],(double *) "+VAR_KO_TEMP_RHS[0]+",(const unsigned int *) "+VAR_TILE_LIMITS+",(const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+",(const unsigned int *) "+VAR_TILE_SZ+");\n")
                ofile.write("\t\t "+FUNC_LOAD_VAR+"(&"+VAR_UNZIP_OUT+"["+VAR_ENUM_TO_OUTPUT_SYM[rhs]+"][offset],(double *) "+VAR_KO_TEMP_RHS[0]+",(const unsigned int *) "+VAR_TILE_LIMITS+",(const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+",(const unsigned int *) "+VAR_TILE_SZ+");\n")
                ofile.write("\t\t __syncthreads();\n\n")

                ofile.write("\t\tif(!(("+VAR_TRD_ID_X+">("+VAR_TILE_LIMITS+"[1]-"+VAR_TILE_LIMITS+"[0])) || ("+VAR_TRD_ID_Y+">("+VAR_TILE_LIMITS+"[3]-"+VAR_TILE_LIMITS+"[2]))) ){ \n\n")

                ofile.write("\t\t pp=0*"+VAR_TILE_SZ+"[0]*"+VAR_TILE_SZ+"[1]+"+VAR_TRD_ID_Y+"*"+VAR_TILE_SZ+"[1]+"+VAR_TRD_ID_X+";\n")
                ofile.write("\t\t  for(unsigned int k=0;k<=(ijk_lm[5]-ijk_lm[4]);++k,pp+="+VAR_TILE_SZ+"[0]*"+VAR_TILE_SZ+"[1]){\n")

                ofile.write("\t\t  "+VAR_KO_TEMP_RHS[0]+"[pp]  += kosigma * ("+VAR_KO_TEMP[0]+"[pp] +"+VAR_KO_TEMP[1]+"[pp] + "+VAR_KO_TEMP[2]+"[pp]);\n")

                ofile.write("\t\t  } //loop z end \n")
                ofile.write("\t\t}// end of the if for the thread idx \n")
                ofile.write("\t\t__syncthreads();\n\n")

                ofile.write("\t\t// sotre computed variables\n\n")
                ofile.write("\t\t"+FUNC_STORE_VAR+"("+VAR_KO_TEMP_RHS[0]+", &"+VAR_UNZIP_OUT+"["+VAR_ENUM_TO_OUTPUT_SYM[rhs]+"][offset],(const unsigned int *) "+VAR_TILE_LIMITS+",(const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+", (const unsigned int *) "+VAR_TILE_LIMITS_STORE+",(const unsigned int *) "+VAR_TILE_SZ+");\n")
                ofile.write("\t\t__syncthreads();\n\n")

        ofile.write("\t  } // end of block assigned to gpu block loop x \n\n")
        ofile.write("\t } // end of block assigned to gpu block loop y \n\n")
        ofile.write("\t} // end of block assigned to gpu block loop z \n\n")


        ofile.write("}// end of function "+FUNC_KO_DISS+"\n")


        ofile.write("}// end of namespace cuda\n")
        ofile.close()


        






def main():


    # shared derivs

    '''
    dxn  = "grad_0"
    dxxn = "grad2_0_0"

    dyn  = "grad_1"
    dyyn = "grad2_1_1"

    dzn  = "grad_2"
    dzzn = "grad2_2_2"

    dxyn = "grad2_0_1"
    dxzn = "grad2_0_2"
    dyzn = "grad2_1_2"

    adxn = "agrad_0"
    adyn = "agrad_1"
    adzn = "agrad_2"

    kodxn = "kograd_0"
    kodyn = "kograd_1"
    kodzn = "kograd_2"
    '''
    dxn  =  VAR_OUT_SHARED_0
    dxxn =  VAR_OUT_SHARED_0

    dyn  = VAR_OUT_SHARED_0
    dyyn = VAR_OUT_SHARED_0

    dzn  = VAR_OUT_SHARED_0
    dzzn = VAR_OUT_SHARED_0

    dxyn = VAR_OUT_SHARED_1
    dxzn = VAR_OUT_SHARED_1
    dyzn = VAR_OUT_SHARED_1

    adxn = VAR_OUT_SHARED_0
    adyn = VAR_OUT_SHARED_0
    adzn = VAR_OUT_SHARED_0

    kodxn = VAR_OUT_SHARED_0
    kodyn = VAR_OUT_SHARED_0
    kodzn = VAR_OUT_SHARED_0


    func_dx="deriv42_x((double *) "+dxn+",(const double *) "+VAR_IN_SHARED+",dx, (const unsigned int *) "+VAR_TILE_LIMITS+" , (const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+" , (const unsigned int *) "+VAR_TILE_SZ+", 3, bflag);"
    func_dy="deriv42_y((double *) "+dyn+",(const double *) "+VAR_IN_SHARED+",dy, (const unsigned int *) "+VAR_TILE_LIMITS+" , (const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+" , (const unsigned int *) "+VAR_TILE_SZ+", 3, bflag);"
    func_dz="deriv42_z((double *) "+dzn+",(const double *) "+VAR_IN_SHARED+",dz, (const unsigned int *) "+VAR_TILE_LIMITS+" , (const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+" , (const unsigned int *) "+VAR_TILE_SZ+", 3, bflag);"

    func_dxx="deriv42_xx((double *) "+dxxn+",(const double *) "+VAR_IN_SHARED+",dx, (const unsigned int *) "+VAR_TILE_LIMITS+" , (const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+" , (const unsigned int *) "+VAR_TILE_SZ+", 3, bflag);"
    func_dxy="deriv42_y((double *) "+dxyn+",(const double *) "+dxn+",dy, (const unsigned int *) "+VAR_TILE_LIMITS+" , (const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+" , (const unsigned int *) "+VAR_TILE_SZ+", 3, bflag);"
    func_dxz="deriv42_z((double *) "+dxzn+",(const double *) "+dxn+",dz, (const unsigned int *) "+VAR_TILE_LIMITS+" , (const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+" , (const unsigned int *) "+VAR_TILE_SZ+", 3, bflag);"

    func_dyy="deriv42_yy((double *) "+dyyn+",(const double *) "+VAR_IN_SHARED+",dy, (const unsigned int *) "+VAR_TILE_LIMITS+" , (const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+" , (const unsigned int *) "+VAR_TILE_SZ+", 3, bflag);"
    func_dyz="deriv42_z((double *) " +dyzn+",(const double *) "+dyn+",dz, (const unsigned int *) "+VAR_TILE_LIMITS+" , (const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+" , (const unsigned int *) "+VAR_TILE_SZ+", 3, bflag);"
    func_dzz="deriv42_zz((double *) "+dzzn+",(const double *) "+VAR_IN_SHARED+",dz, (const unsigned int *) "+VAR_TILE_LIMITS+" , (const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+" , (const unsigned int *) "+VAR_TILE_SZ+", 3, bflag);"

    func_adx="deriv42adv_x((double *) "+adxn+",(const double *) "+VAR_IN_SHARED+",dx, (const unsigned int *) "+VAR_TILE_LIMITS+" , (const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+" , (const unsigned int *) "+VAR_TILE_SZ+", (const bool*) "+ VAR_BETA0_BOOL+" , 3, bflag);"
    func_ady="deriv42adv_y((double *) "+adyn+",(const double *) "+VAR_IN_SHARED+",dy, (const unsigned int *) "+VAR_TILE_LIMITS+" , (const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+" , (const unsigned int *) "+VAR_TILE_SZ+", (const bool*) "+ VAR_BETA1_BOOL+" , 3, bflag);"
    func_adz="deriv42adv_z((double *) "+adzn+",(const double *) "+VAR_IN_SHARED+",dz, (const unsigned int *) "+VAR_TILE_LIMITS+" , (const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+" , (const unsigned int *) "+VAR_TILE_SZ+", (const bool*) "+ VAR_BETA2_BOOL+" , 3, bflag);"

    func_kodx="ko_deriv42_x((double *) "+kodxn+",(const double *) "+VAR_IN_SHARED+",dx,(const unsigned int *) "+VAR_TILE_LIMITS+" , (const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+" , (const unsigned int *) "+VAR_TILE_SZ+", 3, bflag);"
    func_kody="ko_deriv42_y((double *) "+kodyn+",(const double *) "+VAR_IN_SHARED+",dy,(const unsigned int *) "+VAR_TILE_LIMITS+" , (const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+" , (const unsigned int *) "+VAR_TILE_SZ+", 3, bflag);"
    func_kodz="ko_deriv42_z((double *) "+kodzn+",(const double *) "+VAR_IN_SHARED+",dz,(const unsigned int *) "+VAR_TILE_LIMITS+" , (const unsigned int *) "+VAR_DENDRO_BLK_ALIGNED_SZ+" , (const unsigned int *) " +VAR_TILE_SZ+", 3, bflag);"


    ## number of passes for cuda derivatives.
    Derivative = namedtuple("Derivative", "DerivType DerivDir DerivName DerivTile1D DerivInput DerivOutput IB IE JB JE KB KE padWidth DerivFuncCall")

    ####
    ## Since the block shared memory is not enough to compute the all the derivs (15) for a given variable,
    ## we use multiple passes of deriv computations.
    ##
    ## We assume that the deriv TILE is cubic, for simplicity !!!
    ##


    ### !!!!!!! NOTE: WHEN SPECIFYING THE TILE SZ MAKE SURE YOU HAVE 5 POINTS FOR ONE SIDED DERIVS, WHEN THE TILE LOAD THE BLOCK IN THE ITERATIONS
    TileSz1D=12
    bssn_derivs=[
            Derivative(DerivType="d",DerivDir="x",DerivName="deriv_x",DerivTile1D=TileSz1D,DerivInput=VAR_IN_SHARED,DerivOutput="grad_0",IB=3,IE=-3,JB=1,JE=-1,KB=1,KE=-1,padWidth=3,DerivFuncCall="_RSWS_"+func_dx),
            Derivative(DerivType="d",DerivDir="y",DerivName="deriv_y",DerivTile1D=TileSz1D,DerivInput=VAR_IN_SHARED,DerivOutput="grad_1",IB=3,IE=-3,JB=3,JE=-3,KB=1,KE=-1,padWidth=3,DerivFuncCall="_RSWS_"+func_dy),
            Derivative(DerivType="d",DerivDir="z",DerivName="deriv_z",DerivTile1D=TileSz1D,DerivInput=VAR_IN_SHARED,DerivOutput="grad_2",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall="_RSWS_"+func_dz),
            Derivative(DerivType="dd",DerivDir="xx",DerivName="deriv_xx",DerivTile1D=TileSz1D,DerivInput=VAR_IN_SHARED,DerivOutput="grad2_0_0",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall="_RSWS_"+func_dxx),
            Derivative(DerivType="dd",DerivDir="yy",DerivName="deriv_yy",DerivTile1D=TileSz1D,DerivInput=VAR_IN_SHARED,DerivOutput="grad2_1_1",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall="_RSWS_"+func_dyy),
            Derivative(DerivType="dd",DerivDir="zz",DerivName="deriv_zz",DerivTile1D=TileSz1D,DerivInput=VAR_IN_SHARED,DerivOutput="grad2_2_2",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall="_RSWS_"+func_dzz),
            Derivative(DerivType="ko",DerivDir="x",DerivName="ko_deriv_x",DerivTile1D=TileSz1D,DerivInput=VAR_IN_SHARED,DerivOutput="kograd_0",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall="_RSWS_"+func_kodx),
            Derivative(DerivType="ko",DerivDir="y",DerivName="ko_deriv_y",DerivTile1D=TileSz1D,DerivInput=VAR_IN_SHARED,DerivOutput="kograd_1",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall="_RSWS_"+func_kody),
            Derivative(DerivType="ko",DerivDir="z",DerivName="ko_deriv_z",DerivTile1D=TileSz1D,DerivInput=VAR_IN_SHARED,DerivOutput="kograd_2",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall="_RSWS_"+func_kodz),
            Derivative(DerivType="ad",DerivDir="x",DerivName="adv_deriv_x",DerivTile1D=TileSz1D,DerivInput=VAR_IN_SHARED,DerivOutput="agrad_0",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall="_RSWS_"+func_adx),
            Derivative(DerivType="ad",DerivDir="y",DerivName="adv_deriv_y",DerivTile1D=TileSz1D,DerivInput=VAR_IN_SHARED,DerivOutput="agrad_1",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall="_RSWS_"+func_ady),
            Derivative(DerivType="ad",DerivDir="z",DerivName="adv_deriv_z",DerivTile1D=TileSz1D,DerivInput=VAR_IN_SHARED,DerivOutput="agrad_2",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall="_RSWS_"+func_adz),
            Derivative(DerivType="dd",DerivDir="xy",DerivName="deriv_xy",DerivTile1D=TileSz1D,DerivInput=dxn,DerivOutput="grad2_0_1",IB=3,IE=-3,JB=3,JE=-3,KB=1,KE=-1,padWidth=3,DerivFuncCall="_RSWS_"+func_dxy),
            Derivative(DerivType="dd",DerivDir="xz",DerivName="deriv_xz",DerivTile1D=TileSz1D,DerivInput=dxn,DerivOutput="grad2_0_2",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall="_RSWS_"+func_dxz),
            Derivative(DerivType="dd",DerivDir="yz",DerivName="deriv_yz",DerivTile1D=TileSz1D,DerivInput=dyn,DerivOutput="grad2_1_2",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall="_RSWS_"+func_dyz)
        ]


    

    #cudaDerivAllocDeallocHeader("../bssn/cuda_gr/include/bssn_rhs_deriv_mem_cuda.h")
    #cudaDerivAllocDeallocSource("../bssn/cuda_gr/src/bssn_rhs_deriv_mem_cuda.cpp",["bssn_rhs_deriv_mem_cuda.h"])

    
    subset_exp=bssn.outs#[0:4]
    subset_var=bssn.vnames#[0:4]

    cudaCompute("../bssn/cuda_gr/include/rhs_bssn.cuh","../bssn/cuda_gr/src/rhs_bssn.cu",bssn_derivs,subset_exp,subset_var,"__computeBSSNRHS",["block_cu.h","params_cu.h","bssn_rhs_deriv_mem_cuda.h","cudaUtils.cuh","derivs.cuh","cudaUtils.h"])


if __name__ == "__main__":
    main()

