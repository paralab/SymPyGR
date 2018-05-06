//
// Created by milinda on 1/15/18.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Computes BSSN (Generated bssn code) for a given block list.
*/
//

#ifndef BSSN_COMPUTEBSSN_H
#define BSSN_COMPUTEBSSN_H

#define PI 3.14159265

#include <cmath>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "utils.h"

#include "profile_param.h"

using namespace std;

const unsigned int ELE_ORDER=4;
const unsigned int BSSN_NUM_VARS=24;
const unsigned int PAD_WIDTH=3;


struct Block 
{

    unsigned int x;  // x coord
    unsigned int y;  // y coord
    unsigned int z;  // z coord

    unsigned int regLevel; // regular grid level
    unsigned int blkLevel; // block level
    unsigned int maxDepth; // max depth . note that blkLevel < regLevel < maxDepth

    unsigned int node1D_x; // unzip size in x direction
    unsigned int node1D_y; // unzip size in y direction
    unsigned int node1D_z; // unzip size in z direction
    unsigned int offset;  // offset


    Block()
    {
        x=0;y=0;z=0;

        regLevel=0;
        blkLevel=0;
        maxDepth=8;

        regLevel=0;
        node1D_x=0;
        node1D_y=0;
        node1D_z=0;
        offset=0;

    }

    Block(unsigned int px,unsigned int py,unsigned int pz,unsigned int pregLevel,unsigned int pblkLevel,unsigned int pmaxDepth)
    {

        x=px;
        y=py;
        z=pz;

        regLevel=pregLevel;
        blkLevel=pblkLevel;
        maxDepth=pmaxDepth;

        const unsigned int ele1D=(1u<<(regLevel-blkLevel));

        // 3 sizes kept independent if we need ot experiment with mem. alignment (especially for vectorized code)
        node1D_x=ELE_ORDER*ele1D+1+2*PAD_WIDTH;
        node1D_y=ELE_ORDER*ele1D+1+2*PAD_WIDTH;
        node1D_z=ELE_ORDER*ele1D+1+2*PAD_WIDTH;
    }


};




#endif //BSSN_COMPUTEBSSN_H
