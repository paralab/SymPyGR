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

#include "rhs.h"
#include "profile_param.h"
#include "kernal.h"
const unsigned int ELE_ORDER=4;
const unsigned int BSSN_NUM_VARS=24;
const unsigned int PAD_WIDTH=3;




struct Block
{

    unsigned int ele1D;
    unsigned int node1D_x;
    unsigned int node1D_y;
    unsigned int node1D_z;
    unsigned int offset;

    Block()
    {
        ele1D=0;
        node1D_x=0;
        node1D_y=0;
        node1D_z=0;
        offset=0;

    }

    Block(unsigned int pEle1D)
    {
        ele1D=pEle1D; // should be a power of 2
        // 3 sizes kept independent if we need ot experiment with mem. alignment (especially for vectorized code)
        node1D_x=ELE_ORDER*ele1D+1+2*PAD_WIDTH;
        node1D_y=ELE_ORDER*ele1D+1+2*PAD_WIDTH;
        node1D_z=ELE_ORDER*ele1D+1+2*PAD_WIDTH;
    }


};


#endif //BSSN_COMPUTEBSSN_H
