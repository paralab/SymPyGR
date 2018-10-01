//
// Created by milinda on 8/20/18.
//


#include "gpuBSSNExample.h"


int main (int argc, char** argv)
{
    if(argc<2)
        std::cout<<"Usage: "<<argv[0]<<" low_level high_level numBlocks threadX threadY threadZ numStreams"<<std::endl;


    const unsigned int low_level=atoi(argv[1]);
    const unsigned int high_level=atoi(argv[2]);
    const unsigned int numBlocks=atoi(argv[3]);

    unsigned int threadX=2;
    unsigned int threadY=2;
    unsigned int threadZ=2;
    unsigned int streamCount=2;

    if(argc>4)
    {
        threadX=atoi(argv[4]);
        threadY=atoi(argv[5]);
        threadZ=atoi(argv[6]);
    }

    bool useAsync=false;

    if(argc>7)
    {
        useAsync=true;
        streamCount=atoi(argv[7]);
    }



    const double mean = 0.5*(high_level-low_level);
    const double sd=(high_level-mean);

    const unsigned int MAX_STARS=50;
    const unsigned int eleOrder=4;

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(mean,sd);



    unsigned int* blkLevs=new unsigned int[numBlocks];
    unsigned int count=0;
    while(count<numBlocks){
        double number = distribution(generator);
        if ((number>=low_level)&&(number<high_level))
        {
            blkLevs[count]=(int)number;
            count++;
        }


    }

    std::vector<ot::Block> blkList;
    blkList.resize(numBlocks);

    unsigned int unzipSz=0;
    const unsigned int BLOCK_ALIGNMENT_FACTOR=32;

    for(unsigned int i=0;i<numBlocks;i++)
    {
        blkList[i]=ot::Block(2,2,2,blkLevs[i],eleOrder);
        blkList[i].setBlkNodeFlag(0);
        blkList[i].setOffset(unzipSz);
        unzipSz+=(std::ceil((blkList[i].getAllocationSzX()*blkList[i].getAllocationSzY()*blkList[i].getAllocationSzZ()/(double)BLOCK_ALIGNMENT_FACTOR))*BLOCK_ALIGNMENT_FACTOR);
    }

    delete [] blkLevs;

    const unsigned int UNZIP_DOF=unzipSz;


    // variable input
    double ** varUnzipIn = new double*[bssn::BSSN_NUM_VARS];
    double ** varUnzipOutGPU = new double*[bssn::BSSN_NUM_VARS];
    double ** varUnzipOutCPU = new double*[bssn::BSSN_NUM_VARS];


    for(unsigned int var=0;var<bssn::BSSN_NUM_VARS;var++)
    {
        varUnzipIn[var]=new double[UNZIP_DOF];
        varUnzipOutGPU[var]=new double[UNZIP_DOF];
        varUnzipOutCPU[var]=new double[UNZIP_DOF];
    }


    double u_var[bssn::BSSN_NUM_VARS];

    std::cout<<YLW<<" ================================"<<NRM<<std::endl;
    std::cout<<YLW<<"     data init begin             "<<NRM<<std::endl;
    std::cout<<YLW<<" ================================"<<NRM<<std::endl;

    // initialize the input data
    unsigned int offset,bflag;
    unsigned int sz[3];
    double dx,dy,dz;
    double x,y,z;
    const double ptmin[3]={1.0,1.0,1.0};
    double ptmax[3];
    for(unsigned int blk=0;blk<blkList.size();blk++)
    {
        offset=blkList[blk].getOffset();
        sz[0]=blkList[blk].getAllocationSzX();
        sz[1]=blkList[blk].getAllocationSzY();
        sz[2]=blkList[blk].getAllocationSzZ();

        dx=blkList[blk].computeGridDx();
        dy=blkList[blk].computeGridDy();
        dz=blkList[blk].computeGridDz();

        ptmax[0]=ptmin[0]+dx*(sz[0]);
        ptmax[1]=ptmin[1]+dx*(sz[1]);
        ptmax[2]=ptmin[2]+dx*(sz[2]);

        bflag=blkList[blk].getBlkNodeFlag();

        for (unsigned int k = 0; k < sz[2]; k++) {
           z = ptmin[2] + k*dz;

            for (unsigned int j = 0; j < sz[1]; j++) {
                y = ptmin[1] + j*dy;

                for (unsigned int i = 0; i < sz[0]; i++) {
                    x = ptmin[0] + i*dx;
                    bssn::fake_initial_data(x,y,z,u_var);
                    for(unsigned int var=0;var<bssn::BSSN_NUM_VARS;var++)
                        varUnzipIn[var][offset+k*(sz[1]*sz[0])+j*(sz[0])+i]=u_var[var];
                }
            }
        }

    }

    std::cout<<YLW<<" ================================"<<NRM<<std::endl;
    std::cout<<YLW<<"     data init end             "<<NRM<<std::endl;
    std::cout<<YLW<<" ================================"<<NRM<<std::endl;

    std::cout<<"\n\n"<<std::endl;


    std::cout<<YLW<<" ================================"<<NRM<<std::endl;
    std::cout<<YLW<<"     CPU begin             "<<NRM<<std::endl;
    std::cout<<YLW<<" ================================"<<NRM<<std::endl;

    bssn::timer::initialize();

    bssn::timer::total_runtime.start();


    /*for(unsigned int blk=0;blk<blkList.size();blk++) {
        offset = blkList[blk].getOffset();
        sz[0] = blkList[blk].getAllocationSzX();
        sz[1] = blkList[blk].getAllocationSzY();
        sz[2] = blkList[blk].getAllocationSzZ();

        dx = blkList[blk].computeGridDx();
        dy = blkList[blk].computeGridDy();
        dz = blkList[blk].computeGridDz();

        bflag=blkList[blk].getBlkNodeFlag();

        ptmax[0] = ptmin[0] + dx * (sz[0]);
        ptmax[1] = ptmin[1] + dx * (sz[1]);
        ptmax[2] = ptmin[2] + dx * (sz[2]);

        bssnrhs_sep(varUnzipOutCPU, (const double **)varUnzipIn, offset, ptmin, ptmax, sz, bflag);

    }*/

    bssn::timer::total_runtime.stop();

    std::cout<<YLW<<" ================================"<<NRM<<std::endl;
    std::cout<<YLW<<"     CPU end             "<<NRM<<std::endl;
    std::cout<<YLW<<" ================================"<<NRM<<std::endl;

    bssn::timer::profileInfo();


#ifdef BSSN_ENABLE_CUDA


    std::cout<<YLW<<" ================================"<<NRM<<std::endl;
    std::cout<<YLW<<"     GPU begin             "<<NRM<<std::endl;
    std::cout<<YLW<<" ================================"<<NRM<<std::endl;

    cuda::BSSNComputeParams bssnParams;
    bssnParams.BSSN_LAMBDA[0]=bssn::BSSN_LAMBDA[0];
    bssnParams.BSSN_LAMBDA[1]=bssn::BSSN_LAMBDA[1];
    bssnParams.BSSN_LAMBDA[2]=bssn::BSSN_LAMBDA[2];
    bssnParams.BSSN_LAMBDA[3]=bssn::BSSN_LAMBDA[3];

    bssnParams.BSSN_LAMBDA_F[0]=bssn::BSSN_LAMBDA_F[0];
    bssnParams.BSSN_LAMBDA_F[1]=bssn::BSSN_LAMBDA_F[1];

    bssnParams.BSSN_ETA_POWER[0]=bssn::BSSN_ETA_POWER[0];
    bssnParams.BSSN_ETA_POWER[1]=bssn::BSSN_ETA_POWER[1];

    bssnParams.BSSN_ETA_R0=bssn::ETA_R0;
    bssnParams.ETA_CONST=bssn::ETA_CONST;
    bssnParams.ETA_DAMPING=bssn::ETA_DAMPING;
    bssnParams.ETA_DAMPING_EXP=bssn::ETA_DAMPING_EXP;
    bssnParams.KO_DISS_SIGMA=bssn::KO_DISS_SIGMA;
    

    dim3 gpuGrid;
    std::vector<unsigned int > gpuBlockMap;
    dim3 threadBlock(threadX,threadY,threadZ);


    cuda::profile::initialize();
    double hx[3];
    if(useAsync){
       

        unsigned int numberOfTotalBlks = blkList.size();
        cuda::computeDendroBlockToGPUMap(&(*(blkList.begin())),blkList.size(),gpuBlockMap,gpuGrid);
        unsigned int maxBlkSz=0;
        for(unsigned int i=0;i<numberOfTotalBlks;i++)
        {
            int blkSize = (std::ceil((blkList[i].getAllocationSzX()*blkList[i].getAllocationSzY()
                        *blkList[i].getAllocationSzZ()/(double)BLOCK_ALIGNMENT_FACTOR))*BLOCK_ALIGNMENT_FACTOR);
            if(maxBlkSz<blkSize)
                maxBlkSz=blkSize;
        }
         
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp,0);
        const unsigned int numSM=deviceProp.multiProcessorCount;
        streamCount = 2;
        cudaStream_t stream;
        cudaStream_t streams[streamCount];

        // data copy to the GPU    
        cuda::__CUDA_DEVICE_PROPERTIES = cuda::getGPUDeviceInfo(0);

         
        
        double*** referencesTo2DInputArray = new double**[streamCount];
        double*** blkOutput = new double**[streamCount];
        double*** cpuAllocationsOn2DInputArray = new double**[streamCount];
        cuda::_Block** blkLists = new cuda::_Block* [streamCount];
        int copySizes[streamCount];
        int copyOffsets[streamCount];

        for (int index=0; index<streamCount; index++) {
            cudaStreamCreate(&streams[index]);
            referencesTo2DInputArray[index] = cuda::getReferenceTo2DArray<double>(bssn::BSSN_NUM_VARS);
            cpuAllocationsOn2DInputArray[index] = cuda::getReferencesToArrays<double>(bssn::BSSN_NUM_VARS, 
                                                maxBlkSz * numSM);
            blkLists[index] = cuda::allocateMemoryForArray<cuda::_Block>(numberOfTotalBlks);
            blkOutput[index] = cuda::alloc2DGPUArray<double>(bssn::BSSN_NUM_VARS, maxBlkSz * numSM);
        }

        cuda::__BSSN_COMPUTE_PARMS = cuda::copyValueToDevice(&bssnParams);
        cuda::profile::t_cudaMalloc_derivs.start();

                cuda::MemoryDerivs derivWorkSpace;
                derivWorkSpace.allocateDerivMemory(maxBlkSz,numSM);
                CUDA_CHECK_ERROR();

                cuda::__BSSN_DERIV_WORKSPACE=cuda::copyValueToDevice(&derivWorkSpace);
                CUDA_CHECK_ERROR();

        cuda::profile::t_cudaMalloc_derivs.stop();


        unsigned int sendCount[streamCount];
        unsigned int sendOffset[streamCount];
        const unsigned int NUM_GPU_BLOCKS=(gpuBlockMap.size()>>1u);
        cuda::__GPU_BLOCK_MAP=cuda::copyArrayToDeviceAsync(&(*(gpuBlockMap.begin())),NUM_GPU_BLOCKS,(const cudaStream_t*)streams,(unsigned int *)sendCount,(unsigned int *)sendOffset,streamCount);
        
        int counter = 0;
        int globalUnzipDOF = 0;
        bool isStarted = false;
        for(unsigned int blk=0; blk<numberOfTotalBlks; blk+=numSM) {

            int numberOfBlksForNextIteration = (blk+numSM > numberOfTotalBlks?(numberOfTotalBlks - blk) : numSM);
            int temp_unzip_dof = 0;
            int streamSelected = counter%streamCount;
            std::vector<cuda::_Block> blkCudaList;
            blkCudaList.resize(numberOfBlksForNextIteration);
            
            for(unsigned int m= 0; m < numberOfBlksForNextIteration; m++)
            {
                int blkOffset=temp_unzip_dof;
                sz[0]=blkList[blk + m].getAllocationSzX();
                sz[1]=blkList[blk + m].getAllocationSzY();
                sz[2]=blkList[blk + m].getAllocationSzZ();

                bflag=blkList[blk + m].getBlkNodeFlag();

                dx = blkList[blk + m].computeGridDx();
                dy = blkList[blk + m].computeGridDy();
                dz = blkList[blk + m].computeGridDz();

                hx[0]=dx;
                hx[1]=dy;
                hx[2]=dz;

                ptmax[0] = ptmin[0] + dx * (sz[0]);
                ptmax[1] = ptmin[1] + dx * (sz[1]);
                ptmax[2] = ptmin[2] + dx * (sz[2]);

                blkCudaList[m]=cuda::_Block((const double *)ptmin, (const double *)ptmax, blkOffset, bflag,
                                (const unsigned int*)sz, (const double *)hx);
                temp_unzip_dof += (std::ceil((sz[0]*sz[1]*sz[2]/(double)BLOCK_ALIGNMENT_FACTOR))*BLOCK_ALIGNMENT_FACTOR);
            }

            copySizes[streamSelected] = temp_unzip_dof;
            copyOffsets[streamSelected] = globalUnzipDOF;

            cuda::copy2DCudaArray<double>((const double **)varUnzipIn, 
                    cpuAllocationsOn2DInputArray[streamSelected], bssn::BSSN_NUM_VARS,
                    copySizes[streamSelected], referencesTo2DInputArray[streamSelected], globalUnzipDOF,
                    streams[streamSelected]);
            cuda::__DENDRO_BLOCK_LIST = cuda::copyArrayToDevice<cuda::_Block>(blkLists[streamSelected], 
                                    &(*(blkCudaList.begin())), blkCudaList.size(), streams[streamSelected]);
            if(isStarted) {
                //wait on other one to finish - first one will not weight on this
                cudaStreamSynchronize(streams[(streamSelected+1)%streamCount]);
            }
            cuda::computeRHSAsync(blkOutput[streamSelected], referencesTo2DInputArray[streamSelected], 
                cuda::__DENDRO_BLOCK_LIST, blkCudaList.size(), cuda::__BSSN_COMPUTE_PARMS, 
                gpuBlockMap,gpuGrid,threadBlock, streamCount, streams[streamSelected]);

            //data copy back 
            if(isStarted) {
                //wait on other one to finish
                cuda::copy2DArrayToHost<double>(blkOutput[(streamSelected+1)%streamCount], varUnzipOutGPU, bssn::BSSN_NUM_VARS, 
                    copySizes[(streamSelected+1)%streamCount], copyOffsets[(streamSelected+1)%streamCount], streams[(streamSelected+1)%streamCount]);
            }

            if(!isStarted)isStarted=true;
            
            globalUnzipDOF+=temp_unzip_dof;
            counter++;
        }

        cuda::profile::t_cudaMalloc_derivs.start();
            derivWorkSpace.deallocateDerivMemory();
            CUDA_CHECK_ERROR();
        cuda::profile::t_cudaMalloc_derivs.stop();
        for (int index=0; index<streamCount; index++) {
            cudaStreamDestroy(streams[index]);
            cudaFree(blkLists[index]);
            cuda::dealloc2DCudaArray( referencesTo2DInputArray[index],bssn::BSSN_NUM_VARS);
            cuda::dealloc2DCudaArray(blkOutput[index], bssn::BSSN_NUM_VARS);
        }
        cudaFree(cuda::__CUDA_DEVICE_PROPERTIES);
        
        cudaDeviceSynchronize();
        CUDA_CHECK_ERROR();

    }
    else{

        std::vector<cuda::_Block> blkCudaList;
        blkCudaList.resize(blkList.size());

        for(unsigned int blk=0;blk<blkList.size();blk++)
        {
            offset=blkList[blk].getOffset();
            sz[0]=blkList[blk].getAllocationSzX();
            sz[1]=blkList[blk].getAllocationSzY();
            sz[2]=blkList[blk].getAllocationSzZ();

            bflag=blkList[blk].getBlkNodeFlag();

            dx = blkList[blk].computeGridDx();
            dy = blkList[blk].computeGridDy();
            dz = blkList[blk].computeGridDz();

            hx[0]=dx;
            hx[1]=dy;
            hx[2]=dz;

            ptmax[0] = ptmin[0] + dx * (sz[0]);
            ptmax[1] = ptmin[1] + dx * (sz[1]);
            ptmax[2] = ptmin[2] + dx * (sz[2]);

            blkCudaList[blk]=cuda::_Block((const double *)ptmin,(const double *)ptmax,offset,bflag,(const unsigned int*)sz, (const double *)hx);

        }
        cuda::computeDendroBlockToGPUMap(&(*(blkList.begin())),blkList.size(),gpuBlockMap,gpuGrid);
        cuda::computeRHS(varUnzipOutGPU,(const double **)varUnzipIn,&(*(blkCudaList.begin())),blkCudaList.size(),(const cuda::BSSNComputeParams*) &bssnParams,gpuBlockMap,gpuGrid,threadBlock);
    }
        
    cuda::profile::printOutput(blkList);

    std::cout<<YLW<<" ================================"<<NRM<<std::endl;
    std::cout<<YLW<<"     GPU end             "<<NRM<<std::endl;
    std::cout<<YLW<<" ================================"<<NRM<<std::endl;

#endif






    for(unsigned int var=0;var<bssn::BSSN_NUM_VARS;var++)
    {
        delete [] varUnzipIn[var];
        delete [] varUnzipOutCPU[var];
        delete [] varUnzipOutGPU[var];
    }


    delete [] varUnzipIn;
    delete [] varUnzipOutCPU;
    delete [] varUnzipOutGPU;



    return 0;


}



