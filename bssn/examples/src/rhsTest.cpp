
//
// Created by milinda on 8/20/18.
//
#include "rhsTest.h"

int main(int argc, char **argv)
{

    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " low_level high_level numBlocks mode(unstaged=0 staged=1)" << std::endl;
        return 0;
    }

    int rank = 0, npes = 1;

    const unsigned int lLev = atoi(argv[1]);
    const unsigned int hLev = atoi(argv[2]);

    const unsigned int numBlocks = atoi(argv[3]);

    const unsigned int dim = 3;
    const unsigned int mode = atoi(argv[4]);

    if (!rank)
    {
        printf("===================================================================================================================\n");
        printf("low : %d high %d numBlocks %d \n", lLev, hLev, numBlocks);
        printf("===================================================================================================================\n");
    }

    const double mean = 0.5 * (hLev - lLev);
    const double sd = (hLev - mean);

    const unsigned int MAX_STARS = 50;
    const unsigned int eleOrder = 4;

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(mean, sd);

    unsigned int *blkLevs = new unsigned int[numBlocks];
    unsigned int count = 0;
    while (count < numBlocks)
    {
        double number = distribution(generator);
        if ((number >= lLev) && (number < hLev))
        {
            blkLevs[count] = (int)number;
            count++;
        }
    }

    std::vector<ot::Block> blkList;
    blkList.resize(numBlocks);

    unsigned int unzipSz = 0;
    const unsigned int BLOCK_ALIGNMENT_FACTOR = 32;

    for (unsigned int i = 0; i < numBlocks; i++)
    {
        blkList[i] = ot::Block(0, 0, 0, blkLevs[i], eleOrder);
        blkList[i].setBlkNodeFlag(0);
        blkList[i].setOffset(unzipSz);
        unzipSz += blkList[i].getAlignedBlockSz();
    }

    delete[] blkLevs;

    const unsigned int UNZIP_DOF = unzipSz;

    // variable input
    double **varUnzipIn = new double *[bssn::BSSN_NUM_VARS];
    double **varUnzipOutCPU0 = new double *[bssn::BSSN_NUM_VARS]; // staged
    double **varUnzipOutCPU1 = new double *[bssn::BSSN_NUM_VARS]; // unstaged

    for (unsigned int var = 0; var < bssn::BSSN_NUM_VARS; var++)
    {
        varUnzipIn[var] = new double[UNZIP_DOF];
        varUnzipOutCPU0[var] = new double[UNZIP_DOF];
        varUnzipOutCPU1[var] = new double[UNZIP_DOF];
    }

    double u_var[bssn::BSSN_NUM_VARS];
    m_uiMaxDepth = 8;

    std::cout << YLW << " ================================" << NRM << std::endl;
    std::cout << YLW << "     data init begin             " << NRM << std::endl;
    std::cout << YLW << " ================================" << NRM << std::endl;

    bssn::BSSN_COMPD_MIN[0] = -1e-3;
    bssn::BSSN_COMPD_MIN[1] = -1e-3;
    bssn::BSSN_COMPD_MIN[2] = -1e-3;

    bssn::BSSN_COMPD_MAX[0] = 1e-3;
    bssn::BSSN_COMPD_MAX[1] = 1e-3;
    bssn::BSSN_COMPD_MAX[2] = 1e-3;

    bssn::BSSN_OCTREE_MIN[0] = 0;
    bssn::BSSN_OCTREE_MIN[1] = 0;
    bssn::BSSN_OCTREE_MIN[2] = 0;

    bssn::BSSN_OCTREE_MAX[0] = 1u << m_uiMaxDepth;
    bssn::BSSN_OCTREE_MAX[1] = 1u << m_uiMaxDepth;
    bssn::BSSN_OCTREE_MAX[2] = 1u << m_uiMaxDepth;

    bssn::KO_DISS_SIGMA = 0;
    bssn::BH1 = bssn::BH(0.48, -0.01, 0, 1e-6, 0.1, 0, 0, 0, 0, 0);
    bssn::BH2 = bssn::BH(0.52, 0.01, 0, 1e-6, -0.1, 0, 0, 0, 0, 0);

    Point pt_min(bssn::BSSN_COMPD_MIN[0], bssn::BSSN_COMPD_MIN[1], bssn::BSSN_COMPD_MIN[2]);
    Point pt_max(bssn::BSSN_COMPD_MAX[0], bssn::BSSN_COMPD_MAX[1], bssn::BSSN_COMPD_MAX[2]);

    // initialize the input data
    unsigned int offset, bflag;
    unsigned int sz[3];
    double dx, dy, dz;
    double x, y, z;
    double ptmin[3];
    double ptmax[3];
    for (unsigned int blk = 0; blk < blkList.size(); blk++)
    {
        offset = blkList[blk].getOffset();
        sz[0] = blkList[blk].getAllocationSzX();
        sz[1] = blkList[blk].getAllocationSzY();
        sz[2] = blkList[blk].getAllocationSzZ();

        dx = blkList[blk].computeDx(pt_min, pt_max);
        dy = blkList[blk].computeDy(pt_min, pt_max);
        dz = blkList[blk].computeDz(pt_min, pt_max);

        ptmin[0] = GRIDX_TO_X(0) - 3 * dx;
        ptmin[1] = GRIDY_TO_Y(0) - 3 * dy;
        ptmin[2] = GRIDZ_TO_Z(0) - 3 * dz;

        ptmax[0] = GRIDX_TO_X((1u << m_uiMaxDepth)) + 3 * dx;
        ptmax[1] = GRIDY_TO_Y((1u << m_uiMaxDepth)) + 3 * dy;
        ptmax[2] = GRIDZ_TO_Z((1u << m_uiMaxDepth)) + 3 * dz;

        bflag = blkList[blk].getBlkNodeFlag();

        for (unsigned int k = 0; k < sz[2]; k++)
        {
            z = ptmin[2] + k * dz;

            for (unsigned int j = 0; j < sz[1]; j++)
            {
                y = ptmin[1] + j * dy;

                for (unsigned int i = 0; i < sz[0]; i++)
                {
                    x = ptmin[0] + i * dx;
                    bssn::fake_initial_data(x, y, z, u_var);
                    for (unsigned int var = 0; var < bssn::BSSN_NUM_VARS; var++)
                        varUnzipIn[var][offset + k * (sz[1] * sz[0]) + j * (sz[0]) + i] = u_var[var];
                }
            }
        }
    }

    std::cout << YLW << " ================================" << NRM << std::endl;
    std::cout << YLW << "     data init end             " << NRM << std::endl;
    std::cout << YLW << " ================================" << NRM << std::endl;

    std::cout << "\n\n"<< std::endl;

    std::cout << YLW << " ================================" << NRM << std::endl;
    std::cout << YLW << "     CPU begin             " << NRM << std::endl;
    std::cout << YLW << " ================================" << NRM << std::endl;

    if (mode == 0)
    {


        double totalTime = 0.0;
        double derivTime = 0.0;
        double rhsTime = 0.0;
        int numIter = 10; 

        for(int iter = 0; iter<numIter; iter++)
        {

            bssn::timer::initialize();

            auto t1 = Time::now();

            for (unsigned int blk = 0; blk < blkList.size(); blk++)
            {

                const unsigned int offset = blkList[blk].getOffset();
                unsigned int sz[3];
                double hx[3];
                double ptmin[3];
                double ptmax[3];

                sz[0] = blkList[blk].getAllocationSzX();
                sz[1] = blkList[blk].getAllocationSzY();
                sz[2] = blkList[blk].getAllocationSzZ();

                const unsigned int bflag = blkList[blk].getBlkNodeFlag();

                hx[0] = blkList[blk].computeDx(pt_min, pt_max);
                hx[1] = blkList[blk].computeDy(pt_min, pt_max);
                hx[2] = blkList[blk].computeDz(pt_min, pt_max);

                ptmin[0] = GRIDX_TO_X(0) - 3 * hx[0];
                ptmin[1] = GRIDY_TO_Y(0) - 3 * hx[1];
                ptmin[2] = GRIDZ_TO_Z(0) - 3 * hx[2];

                ptmax[0] = GRIDX_TO_X((1u << m_uiMaxDepth)) + 3 * hx[0];
                ptmax[1] = GRIDY_TO_Y((1u << m_uiMaxDepth)) + 3 * hx[1];
                ptmax[2] = GRIDZ_TO_Z((1u << m_uiMaxDepth)) + 3 * hx[2];

                bssnrhs(varUnzipOutCPU1, (const double **)varUnzipIn, offset, ptmin, ptmax, sz, bflag);
            }
            auto t2 = Time::now();
            fsec fs = t2 - t1;
            totalTime += fs.count();
            derivTime += bssn::timer::t_deriv.seconds;
            bssn::timer::t_deriv.seconds = 0.0;
            rhsTime += bssn::timer::t_rhs.seconds;
            bssn::timer::t_rhs.seconds = 0.0;
            //std::cout<<"iteration "<<iter<<" " << derivTime<<" "<<rhsTime<<std::endl;
        }

        totalTime = totalTime/numIter;
        derivTime = derivTime/numIter;
        rhsTime = rhsTime/numIter;
     
        std::cout << "CPU compute unstaged time : " << totalTime << std::endl;
        std::cout << "Derivative  time : " << derivTime << std::endl;
        std::cout << "RHS         time : " << rhsTime << std::endl;
    }
    else if (mode == 1)
    {
        double totalTime = 0.0;
        double derivTime = 0.0;
        double rhsTime = 0.0;
        int numIter = 10; 

        
        for(int iter = 0; iter<numIter; iter++)
        {
            bssn::timer::initialize();

            auto t1 = Time::now();
            for (unsigned int blk = 0; blk < blkList.size(); blk++)
            {

                const unsigned int offset = blkList[blk].getOffset();
                unsigned int sz[3];
                double hx[3];
                double ptmin[3];
                double ptmax[3];

                sz[0] = blkList[blk].getAllocationSzX();
                sz[1] = blkList[blk].getAllocationSzY();
                sz[2] = blkList[blk].getAllocationSzZ();

                const unsigned int bflag = blkList[blk].getBlkNodeFlag();

                hx[0] = blkList[blk].computeDx(pt_min, pt_max);
                hx[1] = blkList[blk].computeDy(pt_min, pt_max);
                hx[2] = blkList[blk].computeDz(pt_min, pt_max);

                ptmin[0] = GRIDX_TO_X(0) - 3 * hx[0];
                ptmin[1] = GRIDY_TO_Y(0) - 3 * hx[1];
                ptmin[2] = GRIDZ_TO_Z(0) - 3 * hx[2];

                ptmax[0] = GRIDX_TO_X((1u << m_uiMaxDepth)) + 3 * hx[0];
                ptmax[1] = GRIDY_TO_Y((1u << m_uiMaxDepth)) + 3 * hx[1];
                ptmax[2] = GRIDZ_TO_Z((1u << m_uiMaxDepth)) + 3 * hx[2];

                bssnrhs_auto_sep(varUnzipOutCPU1, (const double **)varUnzipIn, offset, ptmin, ptmax, sz, bflag);
            }
            auto t2 = Time::now();
            fsec fs = t2 - t1;
            totalTime += fs.count();
            derivTime += bssn::timer::t_deriv.seconds;
            bssn::timer::t_deriv.seconds = 0.0;
            rhsTime += bssn::timer::t_rhs.seconds;
            bssn::timer::t_rhs.seconds = 0.0;
            //std::cout<<"iteration "<<iter<<" " << derivTime<<" "<<rhsTime<<std::endl;
        }

        totalTime = totalTime/numIter;
        derivTime = derivTime/numIter;
        rhsTime = rhsTime/numIter;
        
        std::cout << "CPU compute auto staged time : " << totalTime << std::endl;
        std::cout << "Derivative  time : " << derivTime << std::endl;
        std::cout << "RHS         time : " << rhsTime << std::endl;
    }
    if(mode == 3)
    {   
        //control
        double totalTime = 0.0;
        double derivTime = 0.0;
        double rhsTime = 0.0;
        int numIter = 10; 

        for(int iter = 0; iter<numIter; iter++)
        {

            bssn::timer::initialize();

            auto t1 = Time::now();

            for (unsigned int blk = 0; blk < blkList.size(); blk++)
            {

                const unsigned int offset = blkList[blk].getOffset();
                unsigned int sz[3];
                double hx[3];
                double ptmin[3];
                double ptmax[3];

                sz[0] = blkList[blk].getAllocationSzX();
                sz[1] = blkList[blk].getAllocationSzY();
                sz[2] = blkList[blk].getAllocationSzZ();

                const unsigned int bflag = blkList[blk].getBlkNodeFlag();

                hx[0] = blkList[blk].computeDx(pt_min, pt_max);
                hx[1] = blkList[blk].computeDy(pt_min, pt_max);
                hx[2] = blkList[blk].computeDz(pt_min, pt_max);

                ptmin[0] = GRIDX_TO_X(0) - 3 * hx[0];
                ptmin[1] = GRIDY_TO_Y(0) - 3 * hx[1];
                ptmin[2] = GRIDZ_TO_Z(0) - 3 * hx[2];

                ptmax[0] = GRIDX_TO_X((1u << m_uiMaxDepth)) + 3 * hx[0];
                ptmax[1] = GRIDY_TO_Y((1u << m_uiMaxDepth)) + 3 * hx[1];
                ptmax[2] = GRIDZ_TO_Z((1u << m_uiMaxDepth)) + 3 * hx[2];

                bssnrhs(varUnzipOutCPU1, (const double **)varUnzipIn, offset, ptmin, ptmax, sz, bflag);
            }
            auto t2 = Time::now();
            fsec fs = t2 - t1;
            totalTime += fs.count();
            derivTime += bssn::timer::t_deriv.seconds;
            bssn::timer::t_deriv.seconds = 0.0;
            rhsTime += bssn::timer::t_rhs.seconds;
            bssn::timer::t_rhs.seconds = 0.0;
            //std::cout<<"iteration "<<iter<<" " << derivTime<<" "<<rhsTime<<std::endl;
        }

        totalTime = totalTime/numIter;
        derivTime = derivTime/numIter;
        rhsTime = rhsTime/numIter;
     
        std::cout << "CPU compute unstaged time : " << totalTime << std::endl;
        std::cout << "Derivative  time : " << derivTime << std::endl;
        std::cout << "RHS         time : " << rhsTime << std::endl;

        //staged version
        double totalTimeStaged = 0.0;
        double derivTimeStaged = 0.0;
        double rhsTimeStaged = 0.0;
        int numIterStaged = 10; 

        
        for(int iter = 0; iter<numIterStaged; iter++)
        {
            bssn::timer::initialize();

            auto t1 = Time::now();
            for (unsigned int blk = 0; blk < blkList.size(); blk++)
            {

                const unsigned int offset = blkList[blk].getOffset();
                unsigned int sz[3];
                double hx[3];
                double ptmin[3];
                double ptmax[3];

                sz[0] = blkList[blk].getAllocationSzX();
                sz[1] = blkList[blk].getAllocationSzY();
                sz[2] = blkList[blk].getAllocationSzZ();

                const unsigned int bflag = blkList[blk].getBlkNodeFlag();

                hx[0] = blkList[blk].computeDx(pt_min, pt_max);
                hx[1] = blkList[blk].computeDy(pt_min, pt_max);
                hx[2] = blkList[blk].computeDz(pt_min, pt_max);

                ptmin[0] = GRIDX_TO_X(0) - 3 * hx[0];
                ptmin[1] = GRIDY_TO_Y(0) - 3 * hx[1];
                ptmin[2] = GRIDZ_TO_Z(0) - 3 * hx[2];

                ptmax[0] = GRIDX_TO_X((1u << m_uiMaxDepth)) + 3 * hx[0];
                ptmax[1] = GRIDY_TO_Y((1u << m_uiMaxDepth)) + 3 * hx[1];
                ptmax[2] = GRIDZ_TO_Z((1u << m_uiMaxDepth)) + 3 * hx[2];

                bssnrhs_auto_sep(varUnzipOutCPU1, (const double **)varUnzipIn, offset, ptmin, ptmax, sz, bflag);
            }
            auto t2 = Time::now();
            fsec fs = t2 - t1;
            totalTimeStaged += fs.count();
            derivTimeStaged += bssn::timer::t_deriv.seconds;
            bssn::timer::t_deriv.seconds = 0.0;
            rhsTimeStaged += bssn::timer::t_rhs.seconds;
            bssn::timer::t_rhs.seconds = 0.0;
            //std::cout<<"iteration "<<iter<<" " << derivTime<<" "<<rhsTime<<std::endl;
        }

        totalTimeStaged = totalTimeStaged/numIter;
        derivTimeStaged = derivTimeStaged/numIter;
        rhsTimeStaged = rhsTimeStaged/numIter;
        
        std::cout << "CPU compute auto staged time : " << totalTimeStaged << std::endl;
        std::cout << "Derivative  time : " << derivTimeStaged << std::endl;
        std::cout << "RHS         time : " << rhsTimeStaged << std::endl;

        //check for equality
        double l_inf;
        for (unsigned int var = 0; var < bssn::BSSN_NUM_VARS; var++)
        {
            l_inf = 0;
            for (unsigned int blk = 0; blk < blkList.size(); blk++)
            {

                const unsigned int offset = blkList[blk].getOffset();
                unsigned int sz[3];

                sz[0] = blkList[blk].getAllocationSzX();
                sz[1] = blkList[blk].getAllocationSzY();
                sz[2] = blkList[blk].getAllocationSzZ();

                for (unsigned int k = 3; k < sz[2] - 3; k++)
                    for (unsigned int j = 3; j < sz[1] - 3; j++)
                        for (unsigned int i = 3; i < sz[0] - 3; i++)
                            if (l_inf < fabs(varUnzipOutCPU0[var][offset + k * sz[0] * sz[1] + j * sz[0] + i] - varUnzipOutCPU1[var][offset + k * sz[0] * sz[1] + j * sz[0] + i])){
                                printf("\n");
                                printf("\n");
                                printf(bssn::BSSN_VAR_NAMES[var]);
                                printf("\n");
                                printf("%f",varUnzipOutCPU0[var][offset + k * sz[0] * sz[1] + j * sz[0] + i]);
                                printf("\n");
                                printf("%f", varUnzipOutCPU1[var][offset + k * sz[0] * sz[1] + j * sz[0] + i]);
                                printf("\n");
                                printf("\n");
                                l_inf = fabs(varUnzipOutCPU0[var][offset + k * sz[0] * sz[1] + j * sz[0] + i] - varUnzipOutCPU1[var][offset + k * sz[0] * sz[1] + j * sz[0] + i]);
                            }
            }

            std::cout << "comparison for var: " << var << bssn::BSSN_VAR_NAMES[var] << " l_inf : " << l_inf << std::endl;
        }    
    }

    for (unsigned int var = 0; var < bssn::BSSN_NUM_VARS; var++)
    {
        delete[] varUnzipIn[var];
        delete[] varUnzipOutCPU0[var];
        delete[] varUnzipOutCPU1[var];
    }

    delete[] varUnzipIn;
    delete[] varUnzipOutCPU0;
    delete[] varUnzipOutCPU1;

    return 0;
}
