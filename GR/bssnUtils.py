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

d = ["alpha", "beta0", "beta1", "beta2",
      "B0", "B1", "B2",
      "chi", "Gt0", "Gt1", "Gt2", "K",
      "gt0", "gt1", "gt2", "gt3", "gt4", "gt5",
      "At0", "At1", "At2", "At3", "At4", "At5" ]

# variable names, to access the 2D array. 
varEnum=["cuda::VAR::U_ALPHA","cuda::VAR::U_BETA0","cuda::VAR::U_BETA1","cuda::VAR::U_BETA2","cuda::VAR::U_B0","cuda::VAR::U_B1","cuda::VAR::U_B2","cuda::VAR::U_CHI","cuda::VAR::U_GT0","cuda::VAR::U_GT1","cuda::VAR::U_GT2","cuda::VAR::U_K","cuda::VAR::U_SYMGT0","cuda::VAR::U_SYMGT1","cuda::VAR::U_SYMGT2","cuda::VAR::U_SYMGT3","cuda::VAR::U_SYMGT4","cuda::VAR::U_SYMGT5","cuda::VAR::U_SYMAT0","cuda::VAR::U_SYMAT1","cuda::VAR::U_SYMAT2","cuda::VAR::U_SYMAT3","cuda::VAR::U_SYMAT4","cuda::VAR::U_SYMAT5"]

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
loadVar="cuda::__loadGlobalToShared<double>"
storeVar="cuda::__storeSharedToGlobal<double>"

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

# cuda device properties
cuda_device="__CUDA_DEVICE_PROPERTIES"
# dendro block list parameters 
dendro_blkList="__DENDRO_BLOCK_LIST"
numBlocks="__DENDRO_NUM_BLOCKS"
#x% of block shared memory utilised for the bssn computations
sharedMemUtil="__GPU_BLOCK_SHARED_MEM_UTIL"
derivWorkSpace="__BSSN_DERIV_WORKSPACE"
max_dendro_blk_sz="__BSSN_DERIV_WORKSPACE.__maxBlkSz"

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
dendro_block_sz="sz"
tile_limits="ijk_lm"


#max number of pad provided by the unzip operation in dendro
dendro_blk_pad=3
# required pad to compute the deriv in the current pass. 
dendro_blk_req_pad=0

## tile size of the gpu load functions.
dendro_tile_sz=0;


DerivStruct="MemoryDerivs"
Block_CU="cuda::_Block"


# shared derivs
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

func_dx="deriv42_x((double *) "+dxn+",(const double *) "+varInShared+",dx, (const unsigned int **) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", 2, bflag);"
func_dy="deriv42_y((double *) "+dyn+",(const double *) "+varInShared+",dy, (const unsigned int **) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", 2, bflag);"
func_dz="deriv42_z((double *) "+dzn+",(const double *) "+varInShared+",dz, (const unsigned int **) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", 2, bflag);"

func_dxx="deriv42_xx((double *) "+dxxn+",(const double *) "+varInShared+",dx, (const unsigned int **) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", 2, bflag);"
func_dxy="deriv42_y((double *) "+dxyn+",(const double *) "+dxn+",dy, (const unsigned int **) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", 2, bflag);"
func_dxz="deriv42_z((double *) "+dxzn+",(const double *) "+dxn+",dz, (const unsigned int **) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", 2, bflag);"

func_dyy="deriv42_yy((double *) "+dyyn+",(const double *) "+varInShared+",dy, (const unsigned int **) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", 2, bflag);"
func_dyz="deriv42_z((double *) " +dyzn+",(const double *) "+dyn+",dz, (const unsigned int **) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", 2, bflag);"
func_dzz="deriv42_zz((double *) "+dzzn+",(const double *) "+varInShared+",dz, (const unsigned int **) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", 2, bflag);"

func_adx="deriv42adv_x((double *) "+adxn+",(const double *) "+varInShared+",dx, (const unsigned int **) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", (const double*) "+ beta0+" , 3, bflag);"
func_ady="deriv42adv_y((double *) "+adyn+",(const double *) "+varInShared+",dy, (const unsigned int **) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", (const double*) "+ beta0+" , 3, bflag);"
func_adz="deriv42adv_z((double *) "+adzn+",(const double *) "+varInShared+",dz, (const unsigned int **) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", (const double*) "+ beta0+" , 3, bflag);"

func_kodx="ko_deriv42_x((double *) "+kodxn+",(const double *) "+varInShared+",dx,(const unsigned int **) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", 2, bflag);"
func_kody="ko_deriv42_y((double *) "+kodyn+",(const double *) "+varInShared+",dy,(const unsigned int **) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) "+tile_sz+", 2, bflag);"
func_kodz="ko_deriv42_z((double *) "+kodzn+",(const double *) "+varInShared+",dz,(const unsigned int **) "+tile_limits+" , (const unsigned int *) "+dendro_block_sz+" , (const unsigned int *) " +tile_sz+", 2, bflag);"

## number of passes for cuda derivatives. 
Derivative = namedtuple("Derivative", "DerivType DerivName DerivTile1D DerivInput DerivOutput IB IE JB JE KB KE padWidth DerivFuncCall")

####
## Since the block shared memory is not enough to compute the all the derivs (15) for a given variable, 
## we use multiple passes of deriv computations. 
## 
## We assume that the deriv TILE is cubic, for simplicity !!!
##

cuda_deriv_passes=[
    # deriv pass 1    
    [
     Derivative(DerivType="d",DerivName="deriv_x",DerivTile1D=9,DerivInput=varInShared,DerivOutput="grad_0",IB=3,IE=-3,JB=1,JE=-1,KB=1,KE=-1,padWidth=2,DerivFuncCall=func_dx),
     Derivative(DerivType="d",DerivName="deriv_y",DerivTile1D=9,DerivInput=varInShared,DerivOutput="grad_1",IB=3,IE=-3,JB=3,JE=-3,KB=1,KE=-1,padWidth=2,DerivFuncCall=func_dy),
     Derivative(DerivType="d",DerivName="deriv_z",DerivTile1D=9,DerivInput=varInShared,DerivOutput="grad_2",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=2,DerivFuncCall=func_dz),
     Derivative(DerivType="dd",DerivName="deriv_xy",DerivTile1D=9,DerivInput=dxn,DerivOutput="grad2_0_1",IB=3,IE=-3,JB=3,JE=-3,KB=1,KE=-1,padWidth=2,DerivFuncCall=func_dxy),
     Derivative(DerivType="dd",DerivName="deriv_xz",DerivTile1D=9,DerivInput=dxn,DerivOutput="grad2_0_2",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=2,DerivFuncCall=func_dxz),
     Derivative(DerivType="dd",DerivName="deriv_yz",DerivTile1D=9,DerivInput=dyn,DerivOutput="grad2_1_2",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=2,DerivFuncCall=func_dyz)],
     
    #deriv pass 2  
    [
    Derivative(DerivType="dd",DerivName="deriv_xx",DerivTile1D=9,DerivInput=varInShared,DerivOutput="grad2_0_0",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=2,DerivFuncCall=func_dxx),
    Derivative(DerivType="dd",DerivName="deriv_yy",DerivTile1D=9,DerivInput=varInShared,DerivOutput="grad2_1_1",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=2,DerivFuncCall=func_dyy),
    Derivative(DerivType="dd",DerivName="deriv_zz",DerivTile1D=9,DerivInput=varInShared,DerivOutput="grad2_2_2",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=2,DerivFuncCall=func_dzz),
    Derivative(DerivType="ko",DerivName="ko_deriv_x",DerivTile1D=9,DerivInput=varInShared,DerivOutput="kograd_0",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=2,DerivFuncCall=func_kodx),
    Derivative(DerivType="ko",DerivName="ko_deriv_y",DerivTile1D=9,DerivInput=varInShared,DerivOutput="kograd_1",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=2,DerivFuncCall=func_kody),
    Derivative(DerivType="ko",DerivName="ko_deriv_z",DerivTile1D=9,DerivInput=varInShared,DerivOutput="kograd_2",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=2,DerivFuncCall=func_kodz)],

    #deriv pass 3

    [Derivative(DerivType="ad",DerivName="adv_deriv_x",DerivTile1D=11,DerivInput=varInShared,DerivOutput="agrad_0",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall=func_adx)],

    #deriv pass 4
    [Derivative(DerivType="ad",DerivName="adv_deriv_y",DerivTile1D=11,DerivInput=varInShared,DerivOutput="agrad_1",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall=func_ady)],

    #deriv pass 5
    [Derivative(DerivType="ad",DerivName="adv_deriv_z",DerivTile1D=11,DerivInput=varInShared,DerivOutput="agrad_2",IB=3,IE=-3,JB=3,JE=-3,KB=3,KE=-3,padWidth=3,DerivFuncCall=func_adz)]

    
    ]

cuda_deriv_kernel_names=["__computeDerivPass1","__computeDerivPass2","__computeDerivPass3","__computeDerivPass4","__computeDerivPass5"]





            



def cudaDerivAllocDeallocHeader(fname,headers=[]):
    with open(fname, 'w') as ofile:
        
        ofile.write("// generated by bssnUtils.py code for derivative allocations and deallocations\n")
        ofile.write("//date: "+str(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+"\n")

        ofile.write("#ifndef BSSN_RHS_DERIV_ALLOC \n")
        ofile.write("#define BSSN_RHS_DERIV_ALLOC\n")
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
        ofile.write("\tvoid allocateDerivMemory(unsigned int maxBlkSz); \n")

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


def cudaDerivAllocDeallocSource(fname,headers=[]):
    with open(fname, 'w') as ofile:
        ofile.write("// generated code for derivative allocations and deallocations \n")
        ofile.write("//date: "+str(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+"\n")
        for header in headers:
            ofile.write("#include \""+header+"\"\n")
        ofile.write("\n")
        ofile.write("namespace cuda {\n")

        ofile.write("\n")
        ofile.write("\n")

        ofile.write("\t/**@brief memory allocation for deriv variables */\n")
        ofile.write("\tvoid cuda::"+DerivStruct+"::allocateDerivMemory(unsigned int maxBlkSz){ \n")

        ofile.write("\t\t __maxBlkSz=maxBlkSz;\n")
        ofile.write("\t\t const size_t bytes=sizeof(double)*__maxBlkSz;\n")

        for deriv in func_i:
            ofile.write("\t\t cudaMalloc((void **)&__"+deriv+",bytes);\n")

        for deriv in func_ij:
            ofile.write("\t\t cudaMalloc((void **)&__"+deriv+",bytes);\n")

        for deriv in afunc_i:
            ofile.write("\t\t cudaMalloc((void **)&__"+deriv+",bytes);\n")

        for deriv in kofunc_i:
            ofile.write("\t\t cudaMalloc((void **)&__"+deriv+",bytes);\n")

        ofile.write("\n")
        ofile.write("\n")

        ofile.write("} \n")

        ofile.write("\n")
        ofile.write("\n")

        ofile.write("\t/**@brief memory deallocation for deriv variables */\n")
        ofile.write("\tvoid cuda::"+DerivStruct+"::deallocateDerivMemory(){ \n")

        for deriv in func_i:
            ofile.write("\t\t cudaFree((void **)&__"+deriv+");\n")

        for deriv in func_ij:
            ofile.write("\t\t cudaFree((void **)&__"+deriv+");\n")

        for deriv in afunc_i:
            ofile.write("\t\t cudaFree((void **)&__"+deriv+");\n")

        for deriv in kofunc_i:
            ofile.write("\t\t cudaFree((void **)&__"+deriv+");\n")
        
        ofile.write("\n")
        ofile.write("\n")

        ofile.write("} \n")

        ofile.write("\n")
        ofile.write("\n")

        ofile.write("}// end of namespace cuda\n")

    ofile.close()

def cudaComputeDerivKernelHeader(fname,headers=[]):
    with open(fname, 'w') as ofile:
        ofile.write("// generated by bssnUtils.py code for derivative computation\n")
        ofile.write("//date: "+str(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+"\n")
        ofile.write("#ifndef BSSN_RHS_DERIV_COMP \n")
        ofile.write("#define BSSN_RHS_DERIV_COMP\n")
        ofile.write("#include<iostream>\n")
        ofile.write("#include\"cuda_runtime.h\"\n")
        ofile.write("#include<device_launch_parameters.h>\n")
        for header in headers:
            ofile.write("#include \""+header+"\"\n")

        ofile.write("namespace cuda {\n")

        ofile.write("\n")
        for deriv_kernel in cuda_deriv_kernel_names:
            ofile.write("/**compute derivative kernel "+deriv_kernel+"*/\n")
            ofile.write("__global__ void "+deriv_kernel+"(const double**"+unzipIn+","+DerivStruct+"& __BSSN_DERIV_WORKSPACE, const "+ Block_CU+ "* __DENDRO_BLOCK_LIST, const cudaDeviceProp*__CUDA_DEVICE_PROPERTIES);\n")
            ofile.write("\n")
        
        ofile.write("}// end of namespace cuda\n")

        ofile.write("\n")
        ofile.write("\n")
        ofile.write("#endif\n")
    
    ofile.close()


def cudaComputeDerivKernelSource(fname,headers=[]):
    
    with open(fname, 'w') as ofile:
        ofile.write("// generated by bssnUtils.py code for derivative computation\n")
        ofile.write("//date: "+str(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+"\n")
        for header in headers:
            ofile.write("#include \""+header+"\"\n")
        ofile.write("\n")
        ofile.write("namespace cuda {\n")

        ofile.write("\n")

        for deriv_pass in range(0,len(cuda_deriv_kernel_names)):
            ofile.write("/**compute derivative kernel "+cuda_deriv_kernel_names[deriv_pass]+" */\n")
            ofile.write("__global__ void "+cuda_deriv_kernel_names[deriv_pass]+"(const double**"+unzipIn+","+DerivStruct+"& __BSSN_DERIV_WORKSPACE, const "+ Block_CU+ "* __DENDRO_BLOCK_LIST, const cudaDeviceProp*__CUDA_DEVICE_PROPERTIES){\n")
            ofile.write("\n")

            ofile.write("const _Block dblock="+dendro_blkList+"["+blockId_x+"];\n")
            ofile.write("const unsigned int NUM_SM_UNITS="+cuda_device+"->multiProcessorCount;\n")
            ofile.write("const unsigned int SM_ID="+blockId_x+"%NUM_SM_UNITS;\n")
            ofile.write("const unsigned int offset=dblock.getOffset();\n")
            ofile.write("const unsigned int *sz=dblock.getSz();\n")
            ofile.write("const double* hx=dblock.getDx();\n")

            ofile.write("const double dx=hx[0];\n")
            ofile.write("const double dy=hx[1];\n")
            ofile.write("const double dz=hx[2];\n")

            ofile.write("const double* ptmin=dblock.getPtMin();\n")
            ofile.write("const double* ptmax=dblock.getPtMax();\n")
            ofile.write("const unsigned int bflag=dblock.getBFlag();\n")
            ofile.write("\n")
            ofile.write("const unsigned int blkSz=sz[0]*sz[1]*sz[2];\n")


            ofile.write("\n")

            ofile.write("const unsigned int "+tile_sz+"[3]={"+str(cuda_deriv_passes[deriv_pass][0].DerivTile1D)+","+str(cuda_deriv_passes[deriv_pass][0].DerivTile1D)+","+str(cuda_deriv_passes[deriv_pass][0].DerivTile1D)+"};\n")
            ofile.write("const unsigned int TILE_SZ="+tile_sz+"[0]*"+tile_sz+"[1]*"+tile_sz+"[2]"+";\n")

            dendro_tile_sz=(cuda_deriv_passes[deriv_pass][0].DerivTile1D)**3
            ofile.write("\n")
            ofile.write("\n")

            ofile.write("const unsigned int BLK_INTERATIONS = std::ceil((double)sz[0]/TILE_SZ);")
            ofile.write("\n")

            ofile.write("unsigned int ijk_lm[3][2];")
            ofile.write("\n")

            # allocate memory for shared deriv variables. 
            ofile.write("//allocate memory for shared deriv variables. \n")

            dendro_blk_req_pad=0
            for deriv in cuda_deriv_passes[deriv_pass]:
                ofile.write("__shared__ double "+deriv.DerivOutput+"["+str(dendro_tile_sz)+"];\n")
                if(dendro_blk_req_pad<deriv.padWidth):
                    dendro_blk_req_pad=deriv.padWidth
            
            if(dendro_blk_req_pad>dendro_blk_pad):
                print("code generation error : maxPadwith for derivatives is larger than the dendro block pad width\n")
            
            ofile.write("\n")
            ofile.write("\n")
            
            ofile.write("//allocate memory for shared unzip input. \n")
            ofile.write("__shared__ double "+varInShared+"["+str(dendro_tile_sz)+"];\n")

            ## beta0, beta1, beta2 for advective deriv computations

            advecAlloc=False
            for deriv in cuda_deriv_passes[deriv_pass]:
                if(deriv.DerivType=="ad"):
                    advecAlloc=True

            if(advecAlloc):
                ofile.write("__shared__ double "+beta0+ "["+str(dendro_tile_sz)+"];\n")
                #ofile.write("__shared__ double "+beta1+ "["+str(dendro_tile_sz)+"];\n")
                #ofile.write("__shared__ double "+beta2+ "["+str(dendro_tile_sz)+"];\n\n")


            ofile.write("for(unsigned int iter=0;iter<BLK_INTERATIONS;iter++){\n\n")
            
            ofile.write("\t\t "+tile_limits+"[0][0]=max("+str(dendro_blk_pad-dendro_blk_req_pad)+","+str(dendro_blk_pad-dendro_blk_req_pad)+" + "+tile_sz+"[0]*iter -"+str(2*dendro_blk_req_pad)+");\n")
            if(dendro_blk_req_pad==2):
                ofile.write("\t\t "+tile_limits+"[0][1]=min("+tile_limits+"[0][0]+"+tile_sz+"[0],sz[0]-1);\n")
            else:
                ofile.write("\t\t "+tile_limits+"[0][1]=min("+tile_limits+"[0][0]+"+tile_sz+"[0],sz[0]);\n")
            
            ofile.write("\n")

            ofile.write("\t\t "+tile_limits+"[1][0]=max("+str(dendro_blk_pad-dendro_blk_req_pad)+","+str(dendro_blk_pad-dendro_blk_req_pad)+" + "+tile_sz+"[1]*iter -"+str(2*dendro_blk_req_pad)+");\n")

            if(dendro_blk_req_pad==2):
                ofile.write("\t\t "+tile_limits+"[1][1]=min("+tile_limits+"[1][0]+"+tile_sz+"[1],sz[1]-1);\n")
            else:
                ofile.write("\t\t "+tile_limits+"[1][1]=min("+tile_limits+"[1][0]+"+tile_sz+"[1],sz[1]);\n")
            
            ofile.write("\n")

            ofile.write("\t\t "+tile_limits+"[2][0]=max("+str(dendro_blk_pad-dendro_blk_req_pad)+","+str(dendro_blk_pad-dendro_blk_req_pad)+" + "+tile_sz+"[2]*iter -"+str(2*dendro_blk_req_pad)+");\n")

            if(dendro_blk_req_pad==2):
                ofile.write("\t\t "+tile_limits+"[2][1]=min("+tile_limits+"[2][0]+"+tile_sz+"[2],sz[2]-1);\n")
            else:
                ofile.write("\t\t "+tile_limits+"[2][1]=min("+tile_limits+"[2][0]+"+tile_sz+"[2],sz[2]);\n")
            
            ofile.write("\n")

            for e in d:
                enumStr=varEnum[d.index(e)]
                ofile.write("\n")
                ofile.write("\t\t//load input data from global to shared memory\n")
                ofile.write("\t\t"+loadVar+"(&"+unzipIn+"["+enumStr+"][offset],(double *) "+varInShared+",(const unsigned int **) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")

                #if we are computing advective derivative we need to load the shift vectors. 
                if(e in ad and deriv.DerivType=="ad"):
                    ofile.write("\t\t"+loadVar+"(&"+unzipIn+"["+varEnum[1+(deriv_pass-2)]+"][offset], (double *) beta0, (const unsigned int **) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                    #ofile.write("\t\t"+loadVar+"(&"+unzipIn+"["+varEnum[2]+"][offset], (double *) beta1, (const unsigned int **) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
                    #ofile.write("\t\t"+loadVar+"(&"+unzipIn+"["+varEnum[3]+"][offset], (double *) beta2, (const unsigned int **) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")


                ofile.write("\t\t__syncthreads();\n")
                
                for deriv in cuda_deriv_passes[deriv_pass]:
                    if(deriv.DerivType=="d"):
                        ofile.write("\t\t"+deriv.DerivFuncCall+"\n")

                if(deriv_pass==0): # we need sync only when we compute the mixed 2nd order derivatives
                    ofile.write("\t\t__syncthreads();\n\n")

                for deriv in cuda_deriv_passes[deriv_pass]:
                    
                    if(e in dd and deriv.DerivType=="dd"):
                        ofile.write("\t\t"+deriv.DerivFuncCall+"\n")

                    if(e in ad and deriv.DerivType=="ad"):
                        ofile.write("\t\t"+deriv.DerivFuncCall+"\n")

                    if(e in ad and deriv.DerivType=="ko"):
                        ofile.write("\t\t"+deriv.DerivFuncCall+"\n")
                
                ofile.write("\t\t__syncthreads();\n")

                ofile.write("\n")
                ofile.write("\t\t//store derivs from shared to global memory\n")
                for deriv in cuda_deriv_passes[deriv_pass]:
                    if(deriv.DerivType=="d"):
                        ofile.write("\t\t"+storeVar+"((double *) "+deriv.DerivOutput+",&("+derivWorkSpace+".__"+deriv.DerivOutput+"_"+ e +"[SM_ID * ("+max_dendro_blk_sz+")]),(const unsigned int **) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")

                    if(e in dd and deriv.DerivType=="dd"):
                        ofile.write("\t\t"+storeVar+"((double *) "+deriv.DerivOutput+",&("+derivWorkSpace+".__"+deriv.DerivOutput+"_"+ e +"[SM_ID * ("+max_dendro_blk_sz+")]),(const unsigned int **) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")

                    if(e in ad and deriv.DerivType=="ad"):
                        ofile.write("\t\t"+storeVar+"((double *) "+deriv.DerivOutput+",&("+derivWorkSpace+".__"+deriv.DerivOutput+"_"+ e +"[SM_ID * ("+max_dendro_blk_sz+")]),(const unsigned int **) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")

                    if(e in ad and deriv.DerivType=="ko"):
                        ofile.write("\t\t"+storeVar+"((double *) "+deriv.DerivOutput+",&("+derivWorkSpace+".__"+deriv.DerivOutput+"_"+ e +"[SM_ID * ("+max_dendro_blk_sz+")]),(const unsigned int **) "+tile_limits+",(const unsigned int *) "+dendro_block_sz+",(const unsigned int *) "+tile_sz+");\n")
            
            ofile.write("\t\t} // end of block tile loop\n\n")


            ofile.write("} // end of kernel "+cuda_deriv_kernel_names[deriv_pass]+"\n\n")
            ofile.write("\n")

        ofile.write("}// end of namespace cuda\n")
        ofile.close()



def main():
    cudaDerivAllocDeallocHeader("../bssn/cuda_gr/include/bssn_rhs_deriv_mem_cuda.h")
    cudaDerivAllocDeallocSource("../bssn/cuda_gr/src/bssn_rhs_deriv_mem_cuda.cpp",["bssn_rhs_deriv_mem_cuda.h"])

    cudaComputeDerivKernelHeader("../bssn/cuda_gr/include/rhs.cuh",["block_cu.h","params_cu.h","bssn_rhs_deriv_mem_cuda.h","cudaUtils.cuh","derivs.cuh"])
    cudaComputeDerivKernelSource("../bssn/cuda_gr/src/rhs.cu",["rhs.cuh"])

if __name__ == "__main__":
    main()

