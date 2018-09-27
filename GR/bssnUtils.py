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


## ==== BSSN GPU code generation paramerters

d = ["alpha", "beta0", "beta1", "beta2",
      "B0", "B1", "B2",
      "chi", "Gt0", "Gt1", "Gt2", "K",
      "gt0", "gt1", "gt2", "gt3", "gt4", "gt5",
      "At0", "At1", "At2", "At3", "At4", "At5" ]

# variable names, to access the 2D array.
varEnum=["cuda::VAR::U_ALPHA","cuda::VAR::U_BETA0","cuda::VAR::U_BETA1","cuda::VAR::U_BETA2","cuda::VAR::U_B0","cuda::VAR::U_B1","cuda::VAR::U_B2","cuda::VAR::U_CHI","cuda::VAR::U_GT0","cuda::VAR::U_GT1","cuda::VAR::U_GT2","cuda::VAR::U_K","cuda::VAR::U_SYMGT0","cuda::VAR::U_SYMGT1","cuda::VAR::U_SYMGT2","cuda::VAR::U_SYMGT3","cuda::VAR::U_SYMGT4","cuda::VAR::U_SYMGT5","cuda::VAR::U_SYMAT0","cuda::VAR::U_SYMAT1","cuda::VAR::U_SYMAT2","cuda::VAR::U_SYMAT3","cuda::VAR::U_SYMAT4","cuda::VAR::U_SYMAT5"]

# enum to symbolic input vars dictionary
varEnumToInputSymbol={ "alpha" : "cuda::VAR::U_ALPHA",
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
varEnumToOutputSymbol={    "a_rhs" : "cuda::VAR::U_ALPHA",
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
dd = ["gt0", "gt1", "gt2", "gt3", "gt4", "gt5", "chi",
       "alpha", "beta0", "beta1", "beta2" ]

# advective derivatives
ad = ["gt0", "gt1", "gt2", "gt3", "gt4", "gt5",
       "At0", "At1", "At2", "At3", "At4", "At5",
       "alpha", "beta0", "beta1", "beta2", "chi", "Gt0", "Gt1", "Gt2", "K",
       "B0", "B1", "B2"] 

# first derivs required for constraints--no gauge variables
con_d = [ "chi", "Gt0", "Gt1", "Gt2", "K",
           "gt0", "gt1", "gt2", "gt3", "gt4", "gt5",
           "At0", "At1", "At2", "At3", "At4", "At5" ]

# second derivs required for constraints--no gauge variables
con_dd = ["gt0", "gt1", "gt2", "gt3", "gt4", "gt5", "chi"]


pd = ["grad_0_", "grad_1_", "grad_2_"]
pad = ["agrad_0_", "agrad_1_", "agrad_2_"]
pkod = ["kograd_0_", "kograd_1_", "kograd_2_"]
pdd = ["grad2_0_0_", "grad2_0_1_", "grad2_0_2_", "grad2_1_1_", "grad2_1_2_", "grad2_2_2_"]

# first derivative in i direction
func_i=[]
for f in d:
    for p in pd:
        func_i.append(p+f)

# second derivative in ij direction
func_ij=[]
for f in dd:
    for p in pdd:
        func_ij.append(p+f)

#advective derivative in i direction
afunc_i=[]
for f in ad:
    for p in pad:
        afunc_i.append(p+f)


#Kriess-Oliger derivative in i direction
kofunc_i=[]
for f in d:
    for p in pkod:
        kofunc_i.append(p+f)


# cuda utility functions 
## Note all the device vars which is global starts with __
loadVar="cuda::__loadGlobalToShared3D<double>"
storeVar="cuda::__storeSharedToGlobal3D<double>"
signExtract="cuda::__extractSign3D<double>"

unzipIn="__unzipInVar"
unzipout="__unzipOutVar"


## shift vector block shared variables to compute advective derivs
beta0="beta0"
beta1="beta1"
beta2="beta2"

# shared input variable name for derivative kernels
varInShared="unzipVarInShared"
# shared output variable name for derivative kernels
varOutShared="unzipVarOutShared"
varOutShared1="unzipVarOutShared1"


# block ids
blockId_x="blockIdx.x"
blockId_y="blockIdx.y"
blockId_z="blockIdx.z"

# thread ids
threadId_x="threadIdx.x"
threadId_y="threadIdx.y"
threadId_z="threadIdx.x"

# block dim
blockDim_x="blockDim.x"
blockDim_y="blockDim.y"
blockDim_z="blockDim.z"


# x,y,z bounds of the time i_lm[0] is the min and i_lm[1] is the max. 
tile_sz="tile_sz"
dendro_block_sz="alignedSz"
tile_limits="ijk_lm"


#max number of pad provided by the unzip operation in dendro
dendro_blk_pad=3
# required pad to compute the deriv in the current pass. 
dendro_blk_req_pad=0

## tile size of the gpu load functions.
dendro_tile_sz=0;


DerivStruct="MemoryDerivs"
Block_CU="cuda::_Block"
BSSNComputePars="BSSNComputeParams"


##
# generate the code to allocate derivative memory variables (allocated size unzip_dof)
##
def cudaDerivAllocDeallocHeader(fname,headers=[]):
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


        ofile.write("\tstruct "+DerivStruct+"{\n\n")

        ofile.write("\t unsigned int __maxBlkSz;\n")

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
        ofile.write("\tvoid allocateDerivMemory(unsigned int maxBlkSz, unsigned int numSM); \n")

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
        ofile.write("\tvoid cuda::"+DerivStruct+"::allocateDerivMemory(unsigned int maxBlkSz, unsigned int numSM){ \n")

        ofile.write("\t\t __maxBlkSz=maxBlkSz;\n")
        ofile.write("\t\t const size_t bytes=sizeof(double)*numSM*maxBlkSz;\n")

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
        ofile.write("\tvoid cuda::"+DerivStruct+"::deallocateDerivMemory(){ \n")

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

##
# generate the code to compute derivs  (header)
# derivatives are written indirectly to the global memory.
##
def cudaComputeDerivKernelHeader(fname,cuda_deriv_kernel_names, cuda_deriv_passes,headers=[]):
    with open(fname, 'w') as ofile:
        ofile.write("// generated by Dendro-GR SymPyGR code gernation framework\n")
        ofile.write("//date: "+str(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+"\n")
        fileName,fileExt=os.path.splitext(os.path.basename(fname))
        ofile.write("#ifndef "+fileName.upper()+"_"+fileExt[1:].upper()+" \n")
        ofile.write("#define "+fileName.upper()+"_"+fileExt[1:].upper()+" \n")
        ofile.write("#include<iostream>\n")
        ofile.write("#include\"cuda_runtime.h\"\n")
        ofile.write("#include<device_launch_parameters.h>\n")
        for header in headers:
            ofile.write("#include \""+header+"\"\n")

        ofile.write("namespace cuda {\n")
        ofile.write("\n")
        for deriv_kernel in cuda_deriv_kernel_names:
            ofile.write("/**compute derivative kernel "+deriv_kernel+"*/\n")
            ofile.write("__global__ void "+deriv_kernel+"(const double**"+unzipIn+","+DerivStruct+"* __derivWorkspace, const "+ Block_CU+ "* __dendroBlkList, const cudaDeviceProp* __deviceProperties);\n")
            ofile.write("\n")
        
        ofile.write("}// end of namespace cuda\n")

        ofile.write("\n")
        ofile.write("\n")
        ofile.write("#endif\n")
    
    ofile.close()



##
# generate the code to compute derivs (source)
# # derivatives are written indirectly (use shared mem then to global memory) to the global memory.
##
def cudaComputeDerivKernelSource(fname,cuda_deriv_kernel_names,cuda_deriv_passes,headers=[]):

    # cuda device properties
    cuda_device="__deviceProperties"
    # dendro block list parameters
    dendro_blkList="__dendroBlkList"
    numBlocks="cuda::__DENDRO_NUM_BLOCKS"
    #x% of block shared memory utilised for the bssn computations
    sharedMemUtil="cuda::__GPU_BLOCK_SHARED_MEM_UTIL"
    derivWorkSpace="__derivWorkspace"
    max_dendro_blk_sz=derivWorkSpace+"->__maxBlkSz"


    with open(fname, 'w') as ofile:
        ofile.write("// generated by Dendro-GR SymPyGR code gernation framework\n")
        ofile.write("//date: "+str(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+"\n")
        for header in headers:
            ofile.write("#include \""+header+"\"\n")
        ofile.write("\n")
        ofile.write("namespace cuda {\n")

        ofile.write("\n")

        for deriv_pass in range(0,len(cuda_deriv_kernel_names)):
            ofile.write("/**compute derivative kernel "+cuda_deriv_kernel_names[deriv_pass]+" */\n")
            ofile.write("__global__ void "+cuda_deriv_kernel_names[deriv_pass]+"(const double**"+unzipIn+","+DerivStruct+"* __derivWorkspace, const "+ Block_CU+ "* __dendroBlkList, const cudaDeviceProp*__deviceProperties){\n")
            ofile.write("\n")

            ofile.write("const _Block dblock="+dendro_blkList+"["+blockId_x+"];\n")
            ofile.write("const unsigned int NUM_SM_UNITS="+cuda_device+"->multiProcessorCount;\n")
            ofile.write("const unsigned int SM_ID="+blockId_x+"%NUM_SM_UNITS;\n")
            ofile.write("const unsigned int offset=dblock.getOffset();\n")
            ofile.write("const unsigned int *sz=dblock.getSz();\n")
            ofile.write("const unsigned int *alignedSz=dblock.getAlignedSz();\n")
            ofile.write("const double* hx=dblock.getDx();\n")

            ofile.write("const double dx=hx[0];\n")
            ofile.write("const double dy=hx[1];\n")
            ofile.write("const double dz=hx[2];\n")

            ofile.write("const double* ptmin=dblock.getPtMin();\n")
            ofile.write("const double* ptmax=dblock.getPtMax();\n")
            ofile.write("const unsigned int bflag=dblock.getBFlag();\n")
            ofile.write("\n")
            #ofile.write("const unsigned int blkSz=sz[0]*sz[1]*sz[2];\n")

            ofile.write("\n")

            ofile.write("const unsigned int "+tile_sz+"[3]={"+str(cuda_deriv_passes[deriv_pass][0].DerivTile1D)+","+str(cuda_deriv_passes[deriv_pass][0].DerivTile1D)+","+str(cuda_deriv_passes[deriv_pass][0].DerivTile1D)+"};\n")
            #ofile.write("const unsigned int TILE_SZ="+tile_sz+"[0]*"+tile_sz+"[1]*"+tile_sz+"[2]"+";\n")

            dendro_tile_sz=(cuda_deriv_passes[deriv_pass][0].DerivTile1D)**3
            ofile.write("\n")
            ofile.write("\n")

            # allocate memory for shared deriv variables.
            ofile.write("//allocate memory for shared deriv variables. \n")

            dendro_blk_req_pad=0
            for deriv in cuda_deriv_passes[deriv_pass]:
                if(dendro_blk_req_pad<deriv.padWidth):
                    dendro_blk_req_pad=deriv.padWidth


            if(dendro_blk_req_pad>dendro_blk_pad):
                print("code generation error : maxPadwith for derivatives is larger than the dendro block pad width\n")

            ofile.write("const unsigned int Lb = "+str(dendro_blk_pad-dendro_blk_req_pad)+";// load begin bound\n")
            ofile.write("const unsigned int Le = sz[0]-"+str(dendro_blk_pad-dendro_blk_req_pad)+";// load end bound\n")

            #ofile.write("const unsigned int BLK_INTERATIONS = max(1,(int)ceil((double)(Le-Lb-"+tile_sz+"[0])/("+tile_sz+"[0]-2*" +str(dendro_blk_req_pad)+")));\n")
            ofile.write("\tconst unsigned int BLK_INTERATIONS = ((Le-Lb)<"+tile_sz+"[0])? 1: ((int)ceil((double)(Le-Lb-"+tile_sz+"[0])/("+tile_sz+"[0]-2*" +str(dendro_blk_req_pad)+")))+1;\n")
            ofile.write("\n")

            ofile.write("unsigned int ijk_lm[3*2];")
            ofile.write("\n")

            # allocate memory for shared deriv variables. 
            ofile.write("//allocate memory for shared deriv variables. \n")

            ofile.write("\n")
            ofile.write("\n")
            
            ofile.write("//allocate memory for shared unzip input. \n")
            ofile.write("__shared__ double "+varInShared+"["+str(dendro_tile_sz)+"];\n")
            ofile.write("__shared__ double "+varOutShared+"["+str(dendro_tile_sz)+"];\n")
            ofile.write("__shared__ double "+varOutShared1+"["+str(dendro_tile_sz)+"];\n")


            '''for deriv in cuda_deriv_passes[deriv_pass]:
                ofile.write("double * "+deriv.DerivOutput+"=&derivOut[0];\n")'''


            ## beta0, beta1, beta2 for advective deriv computations

            advecAlloc=False
            for deriv in cuda_deriv_passes[deriv_pass]:
                if(deriv.DerivType=="ad"):
                    advecAlloc=True

            if(advecAlloc):
                ofile.write("__shared__ bool "+beta0+ "["+str(dendro_tile_sz)+"];\n")
                ofile.write("__shared__ bool "+beta1+ "["+str(dendro_tile_sz)+"];\n")
                ofile.write("__shared__ bool "+beta2+ "["+str(dendro_tile_sz)+"];\n\n")


            ofile.write("for(unsigned int iter=0;iter<BLK_INTERATIONS;iter++){\n\n")

            ofile.write("\t\t "+tile_limits+"[2*0+0]=max("+str(dendro_blk_pad-dendro_blk_req_pad)+",(int)("+str(dendro_blk_pad-dendro_blk_req_pad)+" + "+tile_sz+"[0]*iter -2*iter*"+str(dendro_blk_req_pad)+"));\n")
            ofile.write("\t\t "+tile_limits+"[2*0+1]=min("+tile_limits+"[2*0+0]+"+tile_sz+"[0],sz[0]-"+str(dendro_blk_pad-dendro_blk_req_pad)+ ");\n")

            ofile.write("\n")

            ofile.write("\t\t "+tile_limits+"[2*1+0]=max("+str(dendro_blk_pad-dendro_blk_req_pad)+",(int)("+str(dendro_blk_pad-dendro_blk_req_pad)+" + "+tile_sz+"[1]*iter -2*iter*"+str(dendro_blk_req_pad)+"));\n")
            ofile.write("\t\t "+tile_limits+"[2*1+1]=min("+tile_limits+"[2*1+0]+"+tile_sz+"[1],sz[1]-"+str(dendro_blk_pad-dendro_blk_req_pad)+ ");\n")

            ofile.write("\n")

            ofile.write("\t\t "+tile_limits+"[2*2+0]=max("+str(dendro_blk_pad-dendro_blk_req_pad)+",(int)("+str(dendro_blk_pad-dendro_blk_req_pad)+" + "+tile_sz+"[2]*iter -2*iter*"+str(dendro_blk_req_pad)+"));\n")
            ofile.write("\t\t "+tile_limits+"[2*2+1]=min("+tile_limits+"[2*2+0]+"+tile_sz+"[2],sz[2]-"+str(dendro_blk_pad-dendro_blk_req_pad)+ ");\n")

            ofile.write("\t\t //printf(\" iter : %d threadid (%d,%d,%d) tile begin: (%d,%d,%d) tile end: (%d,%d,%d) \\n\",iter, threadIdx.x,threadIdx.y,threadIdx.z,ijk_lm[0],ijk_lm[2],ijk_lm[4],ijk_lm[1],ijk_lm[3],ijk_lm[5]);\n\n")
            ofile.write("\n")

            if(advecAlloc):
                ofile.write("\t\t"+loadVar+"(&"+unzipIn+"[cuda::VAR::U_BETA0][offset],(double *) "+varInShared+",(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                ofile.write("\t\t"+loadVar+"(&"+unzipIn+"[cuda::VAR::U_BETA1][offset],(double *) "+varOutShared+",(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                ofile.write("\t\t"+loadVar+"(&"+unzipIn+"[cuda::VAR::U_BETA2][offset],(double *) "+varOutShared1+",(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                ofile.write("\t\t__syncthreads();\n")
                ofile.write("\t\t"+signExtract+"((double *)"+varInShared+",(bool *) "+beta0+",(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                ofile.write("\t\t"+signExtract+"((double *)"+varOutShared+",(bool *) "+beta1+",(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                ofile.write("\t\t"+signExtract+"((double *)"+varOutShared1+",(bool *) "+beta2+",(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                ofile.write("\t\t__syncthreads();\n")

            for e in d:
                enumStr=varEnum[d.index(e)]
                ofile.write("\n")
                ofile.write("\t\t//load input data from global to shared memory\n")

                ofile.write("\t\t"+loadVar+"(&"+unzipIn+"["+enumStr+"][offset],(double *) "+varInShared+",(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                ofile.write("\t\t__syncthreads();\n")
                ofile.write("\t\t//sync to make sure all the data is loaded\n")

                for deriv in cuda_deriv_passes[deriv_pass]:
                    if((deriv.DerivType=="d") and (deriv.DerivDir=="x")):
                        ofile.write("\t\t// computing deriv "+deriv.DerivDir+" for variable "+e+" \n")
                        ofile.write("\t\t"+deriv.DerivFuncCall+"\n")
                        ofile.write("\t\t__syncthreads();\n") # not essential if each thread writes only the points it has computed in the block.
                        for deriv1 in cuda_deriv_passes[deriv_pass]:
                            if((e in dd) and (deriv1.DerivType=="dd") and ((deriv1.DerivDir=="xy") or (deriv1.DerivDir=="xz"))):
                                ofile.write("\t\t// computing deriv "+deriv1.DerivDir+" for variable "+e+" \n")
                                #ofile.write("\t\t"+deriv1.DerivOutput+" = &"+beta0+"[0];\n")
                                ofile.write("\t\t"+deriv1.DerivFuncCall+"\n")
                                ofile.write("\t\t__syncthreads();\n") # not essential if each thread writes only the points it has computed in the block.
                                # write the mixed derivatives
                                #ofile.write("\t\t"+storeVar+"((double *) "+deriv1.DerivOutput+",&("+derivWorkSpace+"->__"+deriv1.DerivOutput+"_"+ e +"[SM_ID*("+max_dendro_blk_sz+")]),(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                                ofile.write("\t\t"+storeVar+"((double *) "+varOutShared1+",&("+derivWorkSpace+"->__"+deriv1.DerivOutput+"_"+ e +"[SM_ID*("+max_dendro_blk_sz+")]),(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                                if((deriv1.DerivDir=="xy")):
                                    ofile.write("\t\t__syncthreads();\n")
                                ofile.write("\n")
                        #write the x, y,z derivs.
                        #ofile.write("\t\t"+storeVar+"((double *) "+deriv.DerivOutput+",&("+derivWorkSpace+"->__"+deriv.DerivOutput+"_"+ e +"[SM_ID*("+max_dendro_blk_sz+")]),(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                        ofile.write("\t\t"+storeVar+"((double *) "+varOutShared+",&("+derivWorkSpace+"->__"+deriv.DerivOutput+"_"+ e +"[SM_ID*("+max_dendro_blk_sz+")]),(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                        ofile.write("\t\t__syncthreads();\n")
                        ofile.write("\n")

                    if((deriv.DerivType=="d") and (deriv.DerivDir=="y")):
                        ofile.write("\t\t// computing deriv "+deriv.DerivDir+" for variable "+e+" \n")
                        ofile.write("\t\t"+deriv.DerivFuncCall+"\n")
                        ofile.write("\t\t__syncthreads();\n") # not essential if each thread writes only the points it has computed in the block.
                        for deriv1 in cuda_deriv_passes[deriv_pass]:
                            if((e in dd) and (deriv1.DerivType=="dd") and (deriv1.DerivDir=="yz")):
                                ofile.write("\t\t// computing deriv "+deriv1.DerivDir+" for variable "+e+" \n")
                                #ofile.write("\t\t"+deriv1.DerivOutput+" = &"+beta0+"[0];\n")
                                ofile.write("\t\t"+deriv1.DerivFuncCall+"\n")
                                ofile.write("\t\t__syncthreads();\n") # not essential if each thread writes only the points it has computed in the block.
                                # write the mixed derivatives
                                #ofile.write("\t\t"+storeVar+"((double *) "+deriv1.DerivOutput+",&("+derivWorkSpace+"->__"+deriv1.DerivOutput+"_"+ e +"[SM_ID*("+max_dendro_blk_sz+")]),(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                                ofile.write("\t\t"+storeVar+"((double *) "+varOutShared1+",&("+derivWorkSpace+"->__"+deriv1.DerivOutput+"_"+ e +"[SM_ID*("+max_dendro_blk_sz+")]),(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                                #ofile.write("\t\t__syncthreads();\n")
                                ofile.write("\n")
                        #write the x, y,z derivs.
                        #ofile.write("\t\t"+storeVar+"((double *) "+deriv.DerivOutput+",&("+derivWorkSpace+"->__"+deriv.DerivOutput+"_"+ e +"[SM_ID*("+max_dendro_blk_sz+")]),(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                        ofile.write("\t\t"+storeVar+"((double *) "+varOutShared+",&("+derivWorkSpace+"->__"+deriv.DerivOutput+"_"+ e +"[SM_ID*("+max_dendro_blk_sz+")]),(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                        ofile.write("\t\t__syncthreads();\n")
                        ofile.write("\n")

                    if((deriv.DerivType=="d") and (deriv.DerivDir=="z")):
                        ofile.write("\t\t// computing deriv "+deriv.DerivDir+" for variable "+e+" \n")
                        ofile.write("\t\t"+deriv.DerivFuncCall+"\n")
                        ofile.write("\t\t__syncthreads();\n") # not essential if each thread writes only the points it has computed in the block.
                        ofile.write("\t\t"+storeVar+"((double *) "+varOutShared+",&("+derivWorkSpace+"->__"+deriv.DerivOutput+"_"+ e +"[SM_ID*("+max_dendro_blk_sz+")]),(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                        ofile.write("\t\t__syncthreads();\n")

                ofile.write("\n")

                derivCount=0

                for deriv in cuda_deriv_passes[deriv_pass]:
                    if((e in dd) and (deriv.DerivType=="dd")  and ((deriv.DerivDir=="xx") or (deriv.DerivDir=="yy") or (deriv.DerivDir=="zz"))):
                        ofile.write("\t\t"+deriv.DerivFuncCall+"\n")
                        derivCount=derivCount+1
                        if(derivCount==2):
                            derivCount=0
                            ofile.write("\t\t__syncthreads();\n") # not essential if each thread writes only the points it has computed in the block.
                            #ofile.write("\t\t"+storeVar+"((double *) "+deriv.DerivOutput+",&("+derivWorkSpace+"->__"+deriv.DerivOutput+"_"+ e +"[SM_ID*("+max_dendro_blk_sz+")]),(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                            ofile.write("\t\t"+storeVar+"((double *) "+varOutShared+",&("+derivWorkSpace+"->__"+deriv.DerivOutput+"_"+ e +"[SM_ID*("+max_dendro_blk_sz+")]),(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                            ofile.write("\t\t"+storeVar+"((double *) "+varOutShared1+",&("+derivWorkSpace+"->__"+deriv.DerivOutput+"_"+ e +"[SM_ID*("+max_dendro_blk_sz+")]),(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                            ofile.write("\t\t__syncthreads();\n")
                        if(deriv.DerivDir=="zz"):
                            ofile.write("\t\t__syncthreads();\n")
                            ofile.write("\t\t"+storeVar+"((double *) "+varOutShared+",&("+derivWorkSpace+"->__"+deriv.DerivOutput+"_"+ e +"[SM_ID*("+max_dendro_blk_sz+")]),(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                            ofile.write("\t\t__syncthreads();\n")

                ofile.write("\n")

                for deriv in cuda_deriv_passes[deriv_pass]:
                    if( (e in ad) and (deriv.DerivType=="ko")):
                        ofile.write("\t\t"+deriv.DerivFuncCall+"\n")
                        derivCount=derivCount+1
                        if(derivCount==2):
                            derivCount=0
                            ofile.write("\t\t__syncthreads();\n") # not essential if each thread writes only the points it has computed in the block.
                            #ofile.write("\t\t"+storeVar+"((double *) "+deriv.DerivOutput+",&("+derivWorkSpace+"->__"+deriv.DerivOutput+"_"+ e +"[SM_ID*("+max_dendro_blk_sz+")]),(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                            ofile.write("\t\t"+storeVar+"((double *) "+varOutShared+",&("+derivWorkSpace+"->__"+deriv.DerivOutput+"_"+ e +"[SM_ID*("+max_dendro_blk_sz+")]),(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                            ofile.write("\t\t"+storeVar+"((double *) "+varOutShared1+",&("+derivWorkSpace+"->__"+deriv.DerivOutput+"_"+ e +"[SM_ID*("+max_dendro_blk_sz+")]),(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                            ofile.write("\t\t__syncthreads();\n")
                        if(deriv.DerivDir=="z"):
                            ofile.write("\t\t__syncthreads();\n")
                            ofile.write("\t\t"+storeVar+"((double *) "+varOutShared+",&("+derivWorkSpace+"->__"+deriv.DerivOutput+"_"+ e +"[SM_ID*("+max_dendro_blk_sz+")]),(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                            ofile.write("\t\t__syncthreads();\n")

                for deriv in cuda_deriv_passes[deriv_pass]:
                    if( (e in ad) and (deriv.DerivType=="ad")):
                        ofile.write("\t\t"+deriv.DerivFuncCall+"\n")
                        derivCount=derivCount+1
                        if(derivCount==2):
                            derivCount=0
                            ofile.write("\t\t__syncthreads();\n") # not essential if each thread writes only the points it has computed in the block.
                            #ofile.write("\t\t"+storeVar+"((double *) "+deriv.DerivOutput+",&("+derivWorkSpace+"->__"+deriv.DerivOutput+"_"+ e +"[SM_ID*("+max_dendro_blk_sz+")]),(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                            ofile.write("\t\t"+storeVar+"((double *) "+varOutShared+",&("+derivWorkSpace+"->__"+deriv.DerivOutput+"_"+ e +"[SM_ID*("+max_dendro_blk_sz+")]),(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                            ofile.write("\t\t"+storeVar+"((double *) "+varOutShared1+",&("+derivWorkSpace+"->__"+deriv.DerivOutput+"_"+ e +"[SM_ID*("+max_dendro_blk_sz+")]),(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                            ofile.write("\t\t__syncthreads();\n")
                        if(deriv.DerivDir=="z"):
                            ofile.write("\t\t__syncthreads();\n")
                            ofile.write("\t\t"+storeVar+"((double *) "+varOutShared+",&("+derivWorkSpace+"->__"+deriv.DerivOutput+"_"+ e +"[SM_ID*("+max_dendro_blk_sz+")]),(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                            ofile.write("\t\t__syncthreads();\n")


            ofile.write("\t\t} // end of block tile loop\n\n")
            #ofile.write("\t\t__syncthreads();\n")

            ofile.write("} // end of kernel "+cuda_deriv_kernel_names[deriv_pass]+"\n\n")
            ofile.write("\n")

        ofile.write("}// end of namespace cuda\n")
        ofile.close()

# generates the header file for bssn equation computations.
# fname: file name to be written to
# outs: symbolic expressions for bssn equations in sympy dendro frameworl
# varnames: variable names.


##
# RHS code generation unstaged
##
def cudaComputeRHSHeaderUnStaged(fname,outs,varnames,headers=[]):
    with open(fname, 'w') as ofile:
        ofile.write("// generated by Dendro-GR SymPyGR code gernation framework\n")
        ofile.write("//date: "+str(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+"\n")
        fileName,fileExt=os.path.splitext(os.path.basename(fname))
        ofile.write("#ifndef "+fileName.upper()+"_"+fileExt[1:].upper()+" \n")
        ofile.write("#define "+fileName.upper()+"_"+fileExt[1:].upper()+" \n")
        ofile.write("#include<iostream>\n")
        ofile.write("#include\"cuda_runtime.h\"\n")
        ofile.write("#include<device_launch_parameters.h>\n")
        for header in headers:
            ofile.write("#include \""+header+"\"\n")

        ofile.write("namespace cuda {\n")
        ofile.write("\n")

        for var in varnames:
            ofile.write("/** computes rhs "+var+"*/\n")
            ofile.write("__global__ void __compute_"+var+"(double** "+unzipout+", const double**"+unzipIn+", "+DerivStruct+"* __derivWorkspace, const "+ Block_CU+ "* __dendroBlkList, const "+BSSNComputePars+"* __bssnPar, const cudaDeviceProp*__deviceProperties);\n")

        #ofile.write("/**computes bssn equations using gpu. */\n")
        #ofile.write("void rhs_bssn(double** "+unzipout+", const double**"+unzipIn+", "+DerivStruct+"* __derivWorkspace, const "+ Block_CU+ "* __dendroBlkList ,unsigned int numBlocks, const "+BSSNComputePars+"* __bssnPar, const cudaDeviceProp*__deviceProperties);\n")
        #ofile.write("\n")

        ofile.write("}// end of namespace cuda\n")
        ofile.write("\n")
        ofile.write("\n")
        ofile.write("#endif\n")

##
# # RHS code generation unstaged source.
##
def cudaComputeRHSSourceUnStaged(fname,outs,varnames,headers=[],sharedMemSz=48*1024):

    # cuda device properties
    cuda_device="__deviceProperties"
    # dendro block list parameters
    dendro_blkList="__dendroBlkList"
    bssnParams="__bssnPar"
    numBlocks="cuda::__DENDRO_NUM_BLOCKS"
    #x% of block shared memory utilised for the bssn computations
    sharedMemUtil="cuda::__GPU_BLOCK_SHARED_MEM_UTIL"
    derivWorkSpace="__derivWorkspace"
    max_dendro_blk_sz=derivWorkSpace+"->__maxBlkSz"

    sharedMemUtilFactor=0.8

    with open(fname, 'w') as ofile:
        ofile.write("// generated by Dendro-GR SymPyGR code gernation framework\n")
        ofile.write("//date: "+str(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+"\n")
        ofile.write("\n")
        ofile.write("\n")
        for header in headers:
            ofile.write("#include \""+header+"\"\n")

        ofile.write("namespace cuda { \n\n\n\n")

        mi = [0, 1, 2, 4, 5, 8]
        midx = ['00', '01', '02', '11', '12', '22']

        ofile.write("\n")
        idx="[pp]"

        for var_id in range(0,len(varnames)):
            varOut=varnames[var_id]
            exp=outs[var_id]

            print("code generation for : "+varOut)

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
                        if varDep[0:-4] in varEnumToInputSymbol.keys():
                            bssnInputVars.append(varDep[0:-4])
                        elif varDep[0:-4] in varEnumToOutputSymbol.keys():
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
                    if varDep[0:-4] in varEnumToInputSymbol.keys():
                        bssnInputVars.append(varDep[0:-4])
                    elif varDep[0:-4] in varEnumToOutputSymbol.keys():
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
                    if varDep[0:-4] in varEnumToInputSymbol.keys():
                        bssnInputVars.append(varDep[0:-4])
                    elif varDep[0:-4] in varEnumToOutputSymbol.keys():
                        bssnOutputVars.append(varDep[0:-4])
                    else:
                        for key,value in custom_functions.items():
                            if value in varDep[0:-4]:
                                derivVars.append(varDep[0:-4])
                                break


            for lvar in lname:
                if lvar[0:-4] in varEnumToOutputSymbol.keys():
                    bssnOutputVars.append(lvar[0:-4])
                else:
                    bssnStagedVars.append(lvar[0:-4])

            bssnInputVars=list(set(bssnInputVars))
            bssnOutputVars=list(set(bssnOutputVars))
            bssnStagedVars=list(set(bssnStagedVars))
            derivVars=list(set(derivVars))

            total_dep=len(bssnInputVars)+len(bssnStagedVars)+len(derivVars)+len(bssnOutputVars);
            rhs_tile_size=math.floor(((sharedMemUtilFactor*sharedMemSz)/(total_dep*8))**(1.0/3.0))

            print("dependenacy computation completed\n")


            ofile.write("/** computes rhs "+varOut+"*/\n")
            ofile.write("__global__ void __compute_"+varOut+"(double** "+unzipout+", const double**"+unzipIn+", "+DerivStruct+"* __derivWorkspace, const "+ Block_CU+ "* __dendroBlkList, const "+BSSNComputePars+"* __bssnPar, const cudaDeviceProp*__deviceProperties){\n")

            ofile.write("\tconst _Block dblock="+dendro_blkList+"["+blockId_x+"];\n")
            ofile.write("\tconst unsigned int NUM_SM_UNITS="+cuda_device+"->multiProcessorCount;\n")
            ofile.write("\tconst unsigned int SM_ID="+blockId_x+"%NUM_SM_UNITS;\n")
            ofile.write("\tconst unsigned int offset=dblock.getOffset();\n")
            ofile.write("\tconst unsigned int *sz=dblock.getSz();\n")
            ofile.write("\tconst double* hx=dblock.getDx();\n")
            ofile.write("\tconst double* ptmin=dblock.getPtMin();\n")
            ofile.write("\tconst double* ptmax=dblock.getPtMax();\n")


            ofile.write("\t// bssn compute parameters \n")

            ofile.write("\tconst double lambda[4]={__bssnPar->BSSN_LAMBDA[0],__bssnPar->BSSN_LAMBDA[1],__bssnPar->BSSN_LAMBDA[2],__bssnPar->BSSN_LAMBDA[3]};\n")
            ofile.write("\tconst double lambda_f[2]={__bssnPar->BSSN_LAMBDA_F[0],__bssnPar->BSSN_LAMBDA_F[1]};\n")
            ofile.write("\tconst double kosigma=__bssnPar->KO_DISS_SIGMA;\n")
            ofile.write("\tconst double ETA_R0=__bssnPar->ETA_R0;\n")
            ofile.write("\tconst double R0=__bssnPar->ETA_R0;\n")
            ofile.write("\tconst double ETA_DAMPING=__bssnPar->ETA_DAMPING;\n")
            ofile.write("\tconst double ETA_DAMPING_EXP=__bssnPar->ETA_DAMPING_EXP;\n")
            ofile.write("\tconst double ETA_CONST=__bssnPar->ETA_CONST;\n")
            ofile.write("\tconst double eta_power[2]={__bssnPar->BSSN_ETA_POWER[0],__bssnPar->BSSN_ETA_POWER[1]};\n")

            ofile.write("\tconst double dx=hx[0];\n")
            ofile.write("\tconst double dy=hx[1];\n")
            ofile.write("\tconst double dz=hx[2];\n")

            ofile.write("\tconst unsigned int bflag=dblock.getBFlag();\n")
            ofile.write("\n")
            ofile.write("\tconst unsigned int blkSz=sz[0]*sz[1]*sz[2];\n")
            ofile.write("\n")

            ofile.write("\tconst unsigned int "+tile_sz+"[3]={"+str(rhs_tile_size)+","+str(rhs_tile_size)+","+str(rhs_tile_size)+"};\n")
            ofile.write("\tconst unsigned int TILE_SZ="+tile_sz+"[0]*"+tile_sz+"[1]*"+tile_sz+"[2]"+";\n")

            ofile.write("\t\n")
            dendro_blk_req_pad=0;

            ofile.write("\tconst unsigned int Lb = "+str(dendro_blk_pad-dendro_blk_req_pad)+";// load begin bound\n")
            ofile.write("\tconst unsigned int Le = sz[0]-"+str(dendro_blk_pad-dendro_blk_req_pad)+";// load end bound\n")

            #ofile.write("const unsigned int BLK_INTERATIONS = max(1,(int)ceil((double)(Le-Lb-"+tile_sz+"[0])/("+tile_sz+"[0]-2*" +str(dendro_blk_req_pad)+")));\n")
            ofile.write("\tconst unsigned int BLK_INTERATIONS = ((Le-Lb)<"+tile_sz+"[0])? 1: ((int)ceil((double)(Le-Lb-"+tile_sz+"[0])/("+tile_sz+"[0]-2*" +str(dendro_blk_req_pad)+")))+1;\n")
            ofile.write("\n")

            ofile.write("\tunsigned int ijk_lm[3*2];")
            ofile.write("\n")

            # allocate memory for shared deriv variables.
            ofile.write("\t//allocate memory for shared deriv variables. \n")

            ofile.write("\n")
            ofile.write("\n")
            ofile.write("\t //input vars begin\n")
            for var in bssnInputVars:
                ofile.write("\t__shared__ double "+var+"[" + str(rhs_tile_size**3) + "];\n")
            ofile.write("\t //input vars end\n")

            ofile.write("\t // staged vars begin\n")
            for var in bssnStagedVars:
                ofile.write("\t__shared__ double "+var+"[" + str(rhs_tile_size**3) + "];\n")

            ofile.write("\t // staged vars end\n")

            ofile.write("\t // deriv vars begin\n")
            for var in derivVars:
                ofile.write("\t__shared__ double "+var+"[" + str(rhs_tile_size**3) + "];\n")
            ofile.write("\t // deriv vars end\n")

            ofile.write("\t // output vars begin\n")
            for var in bssnOutputVars:
                ofile.write("\t__shared__ double "+var+"[" + str(rhs_tile_size**3) + "];\n")
            ofile.write("\t // output vars end\n")


            ofile.write("\tfor(unsigned int iter=0;iter<BLK_INTERATIONS;iter++){\n\n")
            ofile.write("\t\t "+tile_limits+"[2*0+0]=max("+str(dendro_blk_pad-dendro_blk_req_pad)+",(int)("+str(dendro_blk_pad-dendro_blk_req_pad)+" + "+tile_sz+"[0]*iter -2*iter*"+str(dendro_blk_req_pad)+"));\n")
            ofile.write("\t\t "+tile_limits+"[2*0+1]=min("+tile_limits+"[2*0+0]+"+tile_sz+"[0],sz[0]-"+ str(dendro_blk_pad) +");\n")
            ofile.write("\t\t "+tile_limits+"[2*1+0]=max("+str(dendro_blk_pad-dendro_blk_req_pad)+",(int)("+str(dendro_blk_pad-dendro_blk_req_pad)+" + "+tile_sz+"[1]*iter -2*iter*"+str(dendro_blk_req_pad)+"));\n")
            ofile.write("\t\t "+tile_limits+"[2*1+1]=min("+tile_limits+"[2*1+0]+"+tile_sz+"[1],sz[1]-"+ str(dendro_blk_pad) +");\n")
            ofile.write("\t\t "+tile_limits+"[2*2+0]=max("+str(dendro_blk_pad-dendro_blk_req_pad)+",(int)("+str(dendro_blk_pad-dendro_blk_req_pad)+" + "+tile_sz+"[2]*iter -2*iter*"+str(dendro_blk_req_pad)+"));\n")
            ofile.write("\t\t "+tile_limits+"[2*2+1]=min("+tile_limits+"[2*2+0]+"+tile_sz+"[2],sz[2]-"+ str(dendro_blk_pad) +");\n")

            ofile.write("\t\t //printf(\" iter : %d threadid (%d,%d,%d) tile begin: (%d,%d,%d) tile end: (%d,%d,%d) \\n\",iter, threadIdx.x,threadIdx.y,threadIdx.z,ijk_lm[0],ijk_lm[2],ijk_lm[4],ijk_lm[1],ijk_lm[3],ijk_lm[5]);\n\n")
            ofile.write("\n")

            ofile.write("\t\t//load data from global to shared memory\n")
            for var in bssnInputVars:
                ofile.write("\t\t"+loadVar+"(&"+unzipIn+"["+varEnumToInputSymbol[var]+"][offset],(double *) "+var+",(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")

            for var in derivVars:
                ofile.write("\t\t"+loadVar+"(&("+derivWorkSpace+"->__"+var+"[offset]),(double *) "+var+",(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")

            ofile.write("\t__syncthreads();\n\n")

            ofile.write("\t\tconst unsigned int i_b=ijk_lm[2*0+0];\n")
            ofile.write("\t\tconst unsigned int i_e=ijk_lm[2*0+1];\n")

            ofile.write("\t\tconst unsigned int j_b=ijk_lm[2*1+0];\n")
            ofile.write("\t\tconst unsigned int j_e=ijk_lm[2*1+1];\n")

            ofile.write("\t\tconst unsigned int k_b=ijk_lm[2*2+0];\n")
            ofile.write("\t\tconst unsigned int k_e=ijk_lm[2*2+1];\n")

            ofile.write("\t\tunsigned int l_x=i_e-i_b;\n")
            ofile.write("\t\tunsigned int l_y=j_e-j_b;\n")
            ofile.write("\t\tunsigned int l_z=k_e-k_b;\n")

            ofile.write("\t\tif(threadIdx.x>=l_x || threadIdx.y >= l_y || threadIdx.z>=l_z) return;\n")

            ofile.write("\t\tif(l_x<blockDim.x) l_x=blockDim.x;\n")
            ofile.write("\t\tif(l_y<blockDim.y) l_y=blockDim.y;\n")
            ofile.write("\t\tif(l_z<blockDim.z) l_z=blockDim.z;\n")

            ofile.write("\t\tconst unsigned int ix_b= (i_b + ("+threadId_x+" * l_x)/"+blockDim_x+")-ijk_lm[0];\n")
            ofile.write("\t\tconst unsigned int ix_e= (i_b + (("+threadId_x+" +1)* l_x)/"+blockDim_x+")-ijk_lm[0];\n\n")

            ofile.write("\t\tconst unsigned int jy_b= (j_b + ("+threadId_y+" * l_y)/"+blockDim_y+")-ijk_lm[2];\n")
            ofile.write("\t\tconst unsigned int jy_e= (j_b + (("+threadId_y+"+1)* l_y)/"+blockDim_y+")-ijk_lm[2];\n")

            ofile.write("\t\tconst unsigned int kz_b= (k_b + ("+threadId_z+" * l_z)/"+blockDim_z+")-ijk_lm[4];\n")
            ofile.write("\t\tconst unsigned int kz_e= (k_b + (("+threadId_z+"+1)* l_z)/"+blockDim_z+")-ijk_lm[4];\n")

            ofile.write("\n\n")

            ofile.write("\t\t double x,y,z,r_coord,eta;\n")

            ofile.write("\t\t for (unsigned int k=kz_b;k<kz_e;k++){\n")
            ofile.write("\t\t  z = ptmin[2] + (k+ijk_lm[4])*dz;\n")
            ofile.write("\t\t    for (unsigned int j=jy_b;j<jy_e;j++){\n")
            ofile.write("\t\t     y = ptmin[1] + (j+ijk_lm[2])*dy;\n")
            ofile.write("\t\t      for (unsigned int i=ix_b;i<ix_e;i++){\n")
            ofile.write("\t\t       x = ptmin[0] + (i+ijk_lm[0])*dx;\n")

            ofile.write("\t\t\n\n")
            ofile.write("\t\t       r_coord = sqrt(x*x + y*y + z*z);\n")
            ofile.write("\t\t       eta=ETA_CONST;\n")
            ofile.write("\t\t       if (r_coord >= ETA_R0) {\n")
            ofile.write("\t\t          eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);\n")
            ofile.write("\t\t       }\n\n")


            ofile.write("\t\t      const unsigned int pp=k*"+tile_sz+"[0]*"+tile_sz+"[1]+j*"+tile_sz+"[1]+i;\n")

            ofile.write("\t\t      // Dendro: {{{ \n")
            ofile.write("\t\t      // Dendro: original ops: "+str(sympy.count_ops(lexp))+"\n")

            rops=0
            ofile.write("\t\t      // Dendro: printing temp variables\n")
            for (v1, v2) in _v[0]:
                ofile.write('\t\t   double ')
                ofile.write("\t\t"+dendro.change_deriv_names(sympy.ccode(v2, assign_to=v1, user_functions=custom_functions))+"\n")
                rops = rops + sympy.count_ops(v2)

            ofile.write("\t\t      // Dendro: printing variables\n\n")
            for i, e in enumerate(_v[1]):
                ofile.write("\t\t      "+dendro.change_deriv_names(sympy.ccode(e, assign_to=lname[i], user_functions=custom_functions))+"\n")
                rops = rops + sympy.count_ops(e)


            ofile.write("\t\t      // Dendro: reduced ops: "+str(rops)+"\n")
            ofile.write("\t\t      // Dendro: }}} \n")


            ofile.write("\t\t    }\n")
            ofile.write("\t\t  }\n")
            ofile.write("\t\t}\n")

            ofile.write("\t__syncthreads();\n\n")

            ofile.write("\t// sotre computed variables\n\n")
            for var in bssnOutputVars:
                ofile.write("\t\t"+storeVar+"("+var+", &"+unzipout+"["+varEnumToOutputSymbol[var]+"][offset],(const unsigned int *) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")

            ofile.write("\t__syncthreads();\n")

            ofile.write("\t}\n\n")
            ofile.write("}\n")

        '''ofile.write("/**computes bssn equations using gpu. This is a generated code by Dendro-GR symbolic framework */\n")
        ofile.write("void rhs_bssn(double** "+unzipout+", const double**"+unzipIn+", "+DerivStruct+"* __derivWorkspace, const "+ Block_CU+ "* __dendroBlkList, unsigned int numBlocks, const "+BSSNComputePars+"* __bssnPar, const cudaDeviceProp*__deviceProperties){\n")
        ofile.write("\n")

        ofile.write("\tconst _Block dblock="+dendro_blkList+"["+blockId_x+"];\n")
        ofile.write("\tconst unsigned int NUM_SM_UNITS="+cuda_device+"->multiProcessorCount;\n")
        ofile.write("\tconst unsigned int SM_ID="+blockId_x+"%NUM_SM_UNITS;\n")
        ofile.write("\tconst unsigned int offset=dblock.getOffset();\n")
        ofile.write("\tconst unsigned int *sz=dblock.getSz();\n")
        ofile.write("\tconst double* hx=dblock.getDx();\n")

        ofile.write("\tconst double dx=hx[0];\n")
        ofile.write("\tconst double dy=hx[1];\n")
        ofile.write("\tconst double dz=hx[2];\n")

        ofile.write("\tconst double* ptmin=dblock.getPtMin();\n")
        ofile.write("\tconst double* ptmax=dblock.getPtMax();\n")
        ofile.write("\tconst unsigned int bflag=dblock.getBFlag();\n")
        ofile.write("\n")
        ofile.write("\tconst unsigned int blkSz=sz[0]*sz[1]*sz[2];\n")
        ofile.write("\n")

        ofile.write("//compute all the derivs for the all bssn variables begin============================================================================================================\n")
        for deriv_pass in range(0,len(cuda_deriv_kernel_names)):
            ofile.write("\t// performing deriv pass "+str(deriv_pass)+"\n")
            ofile.write("\tcuda::"+cuda_deriv_kernel_names[deriv_pass]+"<<<dim3(numBlocks,1),dim3(10,10,10)>>>("+unzipIn+", "+derivWorkSpace+", "+dendro_blkList+", "+cuda_device+");\n")
            ofile.write("\tCUDA_CHECK_ERROR();\n")
            ofile.write("\tcudaDeviceSynchronize();\n")
            ofile.write("\tCUDA_CHECK_ERROR();\n")


        ofile.write("\t//sync threads after deriv computation \n")
        #ofile.write("\t__syncthreads();\n")

        ofile.write("\t//compute all the derivs for the all bssn variables end ============================================================================================================\n")

        ofile.write("\n\n\n")
        ofile.write("\t//call for rhs compute begin============================================================================================================\n")

        for var_id in range(0,len(varnames)):
            ofile.write("\tcuda::__compute_"+varnames[var_id]+"<<<dim3(numBlocks,1),dim3(10,10,10)>>>("+unzipout+","+unzipIn+", "+derivWorkSpace+", "+dendro_blkList+", "+bssnParams+", "+cuda_device+");\n")
            ofile.write("\tCUDA_CHECK_ERROR();\n")
            ofile.write("\tcudaDeviceSynchronize();\n")
            ofile.write("\tCUDA_CHECK_ERROR();\n")


        #ofile.write("cudaDeviceSynchronize();\n")
        #ofile.write("CUDA_CHECK_ERROR();\n")

        ofile.write("//call for rhs compute end==============================================================================================================\n")
        ofile.write("\n\n")
        ofile.write("}// end of main kernel\n\n")'''
        ofile.write("} // end of namespace cuda \n\n\n\n")



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
    dxn  =  varOutShared
    dxxn =  varOutShared

    dyn  = varOutShared
    dyyn = varOutShared1

    dzn  = varOutShared
    dzzn = varOutShared

    dxyn = varOutShared1
    dxzn = varOutShared1
    dyzn = varOutShared1

    adxn = varOutShared
    adyn = varOutShared1
    adzn = varOutShared

    kodxn = varOutShared
    kodyn = varOutShared1
    kodzn = varOutShared


    func_dx="deriv42_x((double *) "+dxn+",(const double *) "+varInShared+",dx, (const unsigned int *) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", 2, bflag);"
    func_dy="deriv42_y((double *) "+dyn+",(const double *) "+varInShared+",dy, (const unsigned int *) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", 2, bflag);"
    func_dz="deriv42_z((double *) "+dzn+",(const double *) "+varInShared+",dz, (const unsigned int *) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", 2, bflag);"

    func_dxx="deriv42_xx((double *) "+dxxn+",(const double *) "+varInShared+",dx, (const unsigned int *) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", 2, bflag);"
    func_dxy="deriv42_y((double *) "+dxyn+",(const double *) "+dxn+",dy, (const unsigned int *) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", 2, bflag);"
    func_dxz="deriv42_z((double *) "+dxzn+",(const double *) "+dxn+",dz, (const unsigned int *) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", 2, bflag);"

    func_dyy="deriv42_yy((double *) "+dyyn+",(const double *) "+varInShared+",dy, (const unsigned int *) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", 2, bflag);"
    func_dyz="deriv42_z((double *) " +dyzn+",(const double *) "+dyn+",dz, (const unsigned int *) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", 2, bflag);"
    func_dzz="deriv42_zz((double *) "+dzzn+",(const double *) "+varInShared+",dz, (const unsigned int *) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", 2, bflag);"

    func_adx="deriv42adv_x((double *) "+adxn+",(const double *) "+varInShared+",dx, (const unsigned int *) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", (const bool*) "+ beta0+" , 3, bflag);"
    func_ady="deriv42adv_y((double *) "+adyn+",(const double *) "+varInShared+",dy, (const unsigned int *) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", (const bool*) "+ beta1+" , 3, bflag);"
    func_adz="deriv42adv_z((double *) "+adzn+",(const double *) "+varInShared+",dz, (const unsigned int *) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", (const bool*) "+ beta2+" , 3, bflag);"

    func_kodx="ko_deriv42_x((double *) "+kodxn+",(const double *) "+varInShared+",dx,(const unsigned int *) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", 3, bflag);"
    func_kody="ko_deriv42_y((double *) "+kodyn+",(const double *) "+varInShared+",dy,(const unsigned int *) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", 3, bflag);"
    func_kodz="ko_deriv42_z((double *) "+kodzn+",(const double *) "+varInShared+",dz,(const unsigned int *) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) " +tile_sz+", 3, bflag);"


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
    cuda_deriv_passes=[

        # deriv pass 1
        [
            Derivative(DerivType="d",DerivDir="x",DerivName="deriv_x",DerivTile1D=TileSz1D,DerivInput=varInShared,DerivOutput="grad_0",IB=3,IE=-3,JB=1,JE=-1,KB=1,KE=-1,padWidth=3,DerivFuncCall="_RSWS_"+func_dx),
            Derivative(DerivType="d",DerivDir="y",DerivName="deriv_y",DerivTile1D=TileSz1D,DerivInput=varInShared,DerivOutput="grad_1",IB=3,IE=-3,JB=3,JE=-3,KB=1,KE=-1,padWidth=3,DerivFuncCall="_RSWS_"+func_dy),
            Derivative(DerivType="d",DerivDir="z",DerivName="deriv_z",DerivTile1D=TileSz1D,DerivInput=varInShared,DerivOutput="grad_2",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall="_RSWS_"+func_dz),
            Derivative(DerivType="dd",DerivDir="xx",DerivName="deriv_xx",DerivTile1D=TileSz1D,DerivInput=varInShared,DerivOutput="grad2_0_0",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall="_RSWS_"+func_dxx),
            Derivative(DerivType="dd",DerivDir="yy",DerivName="deriv_yy",DerivTile1D=TileSz1D,DerivInput=varInShared,DerivOutput="grad2_1_1",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall="_RSWS_"+func_dyy),
            Derivative(DerivType="dd",DerivDir="zz",DerivName="deriv_zz",DerivTile1D=TileSz1D,DerivInput=varInShared,DerivOutput="grad2_2_2",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall="_RSWS_"+func_dzz),
            Derivative(DerivType="ko",DerivDir="x",DerivName="ko_deriv_x",DerivTile1D=TileSz1D,DerivInput=varInShared,DerivOutput="kograd_0",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall="_RSWS_"+func_kodx),
            Derivative(DerivType="ko",DerivDir="y",DerivName="ko_deriv_y",DerivTile1D=TileSz1D,DerivInput=varInShared,DerivOutput="kograd_1",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall="_RSWS_"+func_kody),
            Derivative(DerivType="ko",DerivDir="z",DerivName="ko_deriv_z",DerivTile1D=TileSz1D,DerivInput=varInShared,DerivOutput="kograd_2",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall="_RSWS_"+func_kodz),
            Derivative(DerivType="ad",DerivDir="x",DerivName="adv_deriv_x",DerivTile1D=TileSz1D,DerivInput=varInShared,DerivOutput="agrad_0",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall="_RSWS_"+func_adx),
            Derivative(DerivType="ad",DerivDir="y",DerivName="adv_deriv_y",DerivTile1D=TileSz1D,DerivInput=varInShared,DerivOutput="agrad_1",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall="_RSWS_"+func_ady),
            Derivative(DerivType="ad",DerivDir="z",DerivName="adv_deriv_z",DerivTile1D=TileSz1D,DerivInput=varInShared,DerivOutput="agrad_2",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall="_RSWS_"+func_adz),
            Derivative(DerivType="dd",DerivDir="xy",DerivName="deriv_xy",DerivTile1D=TileSz1D,DerivInput=dxn,DerivOutput="grad2_0_1",IB=3,IE=-3,JB=3,JE=-3,KB=1,KE=-1,padWidth=3,DerivFuncCall="_RSWS_"+func_dxy),
            Derivative(DerivType="dd",DerivDir="xz",DerivName="deriv_xz",DerivTile1D=TileSz1D,DerivInput=dxn,DerivOutput="grad2_0_2",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall="_RSWS_"+func_dxz),
            Derivative(DerivType="dd",DerivDir="yz",DerivName="deriv_yz",DerivTile1D=TileSz1D,DerivInput=dyn,DerivOutput="grad2_1_2",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall="_RSWS_"+func_dyz)


        ]


    ]


    cuda_deriv_kernel_names=["__RSWS_computeDerivs"]


    cudaDerivAllocDeallocHeader("../bssn/cuda_gr/include/bssn_rhs_deriv_mem_cuda.h")
    cudaDerivAllocDeallocSource("../bssn/cuda_gr/src/bssn_rhs_deriv_mem_cuda.cpp",["bssn_rhs_deriv_mem_cuda.h"])

    cudaComputeDerivKernelHeader("../bssn/cuda_gr/include/derivs_bssn.cuh",cuda_deriv_kernel_names,cuda_deriv_passes,["block_cu.h","params_cu.h","bssn_rhs_deriv_mem_cuda.h","cudaUtils.cuh","derivs.cuh"])
    cudaComputeDerivKernelSource("../bssn/cuda_gr/src/derivs_bssn.cu",cuda_deriv_kernel_names,cuda_deriv_passes,["derivs_bssn.cuh"])


    #subset_exp=bssn_stages.outs[0:1]
    #subset_var=bssn_stages.vnames[0:1]
    #subset_exp=bssn_stages.outs[5:len(bssn_stages.outs)]
    #subset_var=bssn_stages.vnames[5:len(bssn_stages.outs)]
    #subset_exp=bssn.outs[0:len(bssn.outs)]
    #subset_var=bssn.vnames[0:len(bssn.vnames)]
    #cudaComputeRHSHeaderUnStaged("../../../BSSN_GR/cuda_gr/include/rhs_bssn.cuh",subset_exp,subset_var,["block_cu.h","params_cu.h","bssn_rhs_deriv_mem_cuda.h","cudaUtils.cuh","derivs.cuh","derivs_bssn.cuh","cudaUtils.h"])
    #cudaComputeRHSSourceUnStaged("../../../BSSN_GR/cuda_gr/src/rhs_bssn.cu",subset_exp,subset_var,["rhs_bssn.cuh"])


if __name__ == "__main__":
    main()

