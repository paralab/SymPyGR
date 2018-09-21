/**
 * Create on: Sep 21, 2018
 *      Adapted by: Akila, Eranga, Eminda, Ruwan
 **/
 
#ifndef BLOCK_H
#define BLOCK_H

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
    unsigned int blkSize;  // size
    unsigned int offset;  // offset
    unsigned int block_no;  // offset


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
        block_no=0;

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
        blkSize=node1D_x*node1D_y*node1D_z;    
    }
};

#endif