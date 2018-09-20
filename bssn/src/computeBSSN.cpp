//
// Created by milinda on 1/15/18.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief
*/
//
#include "computeBSSN.h"
#include "rhs.h"

int main (int argc, char** argv){
    /**
     *
     * parameters:
     * lowerLevel: block element 1d lower bound (int)
     * upperLevel: block element 1d upper bound (int) (blk_up>=blk_lb)
     * numberOfBlocks : total blocks going to process 
     * mean(double optional) 
     * std(double optional) 
     * isRandom(1/0; optional): default=0
     *
     * */

    int lower_bound;
    int upper_bound;
    int numberOfBlocks;
    double mean;
    double std = 1;
    bool isRandom = 0;

    if(argc>=4){
        lower_bound=atoi(argv[1]);
        upper_bound=atoi(argv[2]);
        numberOfBlocks=atoi(argv[3]);

        mean = lower_bound + (upper_bound-lower_bound)/2;

        if (argc>=5) mean=atoi(argv[4]);
        if (argc>=6) std=atoi(argv[5]);
        if (argc>=7) isRandom=atoi(argv[6]);
    }else{
        std::cerr << "Correct usage: " << argv[0] << " lowerLevel(int) upperLevel(int) numberOfBlocks(int) mean(double optional) std(double optional) isRandom(1/0; optional)" << std::endl;
        exit(0);
    }

    bssn::timer::total_runtime.start();

    Block * blkList = new Block[numberOfBlocks];
    double ** var_in = new double*[BSSN_NUM_VARS];;
    double ** var_out = new double*[BSSN_NUM_VARS];;

    const unsigned int maxDepth=12;

    // Create distribution
    int seed = 100;
    if (isRandom) seed=time(0);
    std::default_random_engine generator(seed);
    std::normal_distribution<double> distributionN(mean, std);

    int p[upper_bound-lower_bound+1]; // distribution repr purpose
    for (int i=0; i<upper_bound-lower_bound+1; i++) p[i] = 0;

    int level;
    int block_no = 0;
    unsigned long unzipSz=0;
    while (block_no<numberOfBlocks){
        level = int(distributionN(generator)); // generate a number

        if ((level>=lower_bound)&&(level<=upper_bound)) {
            Block & blk=blkList[block_no];
            blk=Block(0, 0, 0, 2*level, level, maxDepth);
            blk.block_no = block_no;
            blk.offset=unzipSz;
            unzipSz+=blk.blkSize;
            block_no++;
            p[level-lower_bound]++; // distribution repr purpose
        }
    }

    // distribution representation requirement
    std::cout << "---Level distribution---" << std::endl;

    for (int i=0; i<=upper_bound-lower_bound; i++) {
        std::cout << i+lower_bound << ": ";
        std::cout << std::setw(3) << p[i]*100/numberOfBlocks << "% " << std::string(p[i]*100/numberOfBlocks, '*') << std::endl;
    }

    const unsigned long unzip_dof=unzipSz;

    // Data generation for CPU version
    for(int i=0;i<BSSN_NUM_VARS;i++){
        var_in[i] = new double[unzip_dof];
        var_out[i] = new double[unzip_dof];
    }
        
    #pragma omp parallel for num_threads(20)
    for(unsigned int blk=0; blk<numberOfBlocks; blk++){
        
        double coord[3];
        double u[BSSN_NUM_VARS];

        Block tmpBlock=blkList[blk];
        double x=(double)tmpBlock.x;
        double y=(double)tmpBlock.y;
        double z=(double)tmpBlock.z;

        double hx=0.001;
        double hy=0.001;
        double hz=0.001;

        unsigned int offset=tmpBlock.offset;
        unsigned int size_x=tmpBlock.node1D_x;
        unsigned int size_y=tmpBlock.node1D_y;
        unsigned int size_z=tmpBlock.node1D_z;

        for(unsigned int k=0;k<tmpBlock.node1D_z;k++){
            for(unsigned int j=0;j<tmpBlock.node1D_y;j++){
                for(unsigned int i=0;i<tmpBlock.node1D_x;i++){
                    coord[0]=x+i*hx;
                    coord[1]=y+j*hy;
                    coord[2]=z+k*hz;

                    initial_data(u,coord);

                    for(unsigned int var=0;var<BSSN_NUM_VARS;var++){
                        var_in[var][offset+k*size_y*size_x+j*size_y+i]=u[var];
                        var_out[var][offset+k*size_y*size_x+j*size_y+i]=0;
                    }
                }
            }
        }
    }

    double ptmin[3], ptmax[3];
    unsigned int sz[3];
    unsigned int bflag;
    unsigned int offset;
    double dx, dy, dz;

    for(unsigned int blk=0; blk<numberOfBlocks; blk++){

        offset=blkList[blk].offset;
        sz[0]=blkList[blk].node1D_x; 
        sz[1]=blkList[blk].node1D_y;
        sz[2]=blkList[blk].node1D_z;

        bflag=0; // indicates if the block is bdy block.

        dx=0.1;
        dy=0.1;
        dz=0.1;

        ptmin[0]=0.0;
        ptmin[1]=0.0;
        ptmin[2]=0.0;

        ptmax[0]=1.0;
        ptmax[1]=1.0;
        ptmax[2]=1.0;

        bssnrhs(var_out, (const double **)var_in, offset, ptmin, ptmax, sz, bflag);
    }

    for (int var=0; var<BSSN_NUM_VARS; var++){
        delete [] var_in[var];
        delete [] var_out[var];
    }

    bssn::timer::total_runtime.stop();

    bssn::timer::profileInfo();
    return 0;
}
