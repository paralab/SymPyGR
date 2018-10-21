//
// Created by milinda on 4/21/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief contains the class for the block. In order to perform a stencil on the adaptive grid, we treat the adaptive grid as a collection of regular blocks.
 *  Let \f$ \tau \f$ be a 2:1 balance complete sorted (according SFC ordering) octree then we can decompose \f$\tau =\{b_i^l\}\f$   sequence of finite number of regular
 *  blocks.
 *
 *  This class contains the block coordinates level of the regular grid embedded, stencil ghost width and number of local points available.
*/

#ifndef SFCSORTBENCH_BLOCK_H
#define SFCSORTBENCH_BLOCK_H

#include "dendro.h"
#include <assert.h>
#include "point.h"

extern unsigned int m_uiMaxDepth;

#define GHOST_WIDTH 3
#define DX 1e-3
#define DY 1e-3
#define DZ 1e-3
namespace ot
{

   class Block
   {



   private:
     /**coords of the block*/
     unsigned int m_uiX;

     unsigned int m_uiY;

     unsigned int m_uiZ;

     /** size of the regular grid inside the block. */
     unsigned int m_uiRegGridLev;

     /**num elements per one side*/
     unsigned int m_uiBlkElem_1D;

     /** offset used for local memory allocation.*/
     DendroIntL m_uiOffset;

     /** array size (1D for the current block. */
     unsigned int m_uiSize1D;

     /**padding width (1D) used for pad the block for neighbour blocks. */
     unsigned int m_uiPaddingWidth;

     /**element order */
     unsigned int m_uiEleOrder;

     /** allocation length on X direction*/
     unsigned int m_uiSzX;

     /** allocation length on Y direction*/
     unsigned int m_uiSzY;

     /** allocation length on Z direction*/
     unsigned int m_uiSzZ;

    /**boundary */
    unsigned int m_uiBflag;




   public:
    /**@brief Default constructor*/
    Block();

    /**
     * @brief constructor to initialize and create a block.
     * @param [in] pNode ot::TreeNode for the block.
     * @param [in] rotID rotation ID for the block.
     * @param [in] regLev level of the regular grid embedded by the block.
     * @param [in] regEleBegin Local element begin location for the for the octree embedded by the block.
     * @param [in] regEleEnd Local element end location for the octree embedded by the block .
     * @param [in] eleorder: element order of the mesh.
     * */
     Block(unsigned int px,unsigned int py,unsigned int pz,unsigned int regLev,unsigned int  eleOrder);

     ~Block();


    /**
     * @brief returns the regular grid lev (m_uiRegGridLev) value.
     * note: In octree2BlockDecomposition m_uiRegGridLev is used to store the rotation id of the block.
     *  */
     inline unsigned int getRegularGridLev()const {return m_uiRegGridLev;}

     /** @brief returns 1D padding width */
     inline unsigned int get1DPadWidth()const {return m_uiPaddingWidth;}

     /**@brief returns the element order*/
     inline unsigned int getElementOrder() const {return m_uiEleOrder;}

     /**@brief set the block offset*/
     void setOffset(DendroIntL offset);

     /**@brief set the blkFlag with the correct bdy*/
     inline void setBlkNodeFlag(unsigned int flag){ m_uiBflag=flag;};

     /**@brief set the blkFlag with the correct bdy*/
     inline unsigned int getBlkNodeFlag() const { return m_uiBflag;};

     /** @brief get offset*/
     inline DendroIntL getOffset() const {return m_uiOffset;}

     /** @brief returns the 1D array size*/
     inline unsigned int get1DArraySize() const {return m_uiSize1D;}

     /**@brief allocation length on X direction*/
     inline unsigned int getAllocationSzX() const {return m_uiSzX;}

     /**@brief allocation length on Y direction*/
     inline unsigned int getAllocationSzY() const {return m_uiSzY;}

     /**@brief allocation length on Z direction*/
     inline unsigned int getAllocationSzZ() const {return m_uiSzZ;}



     /**@brief computes and returns the space discretization (grid domain) */
     double computeGridDx () const;

     /**@brief computes and returns the space discretization (grid domain) */
     double computeGridDy () const;

     /**@brief computes and returns the space discretization (grid domain) */
     double computeGridDz () const;

       /**@brief computes and returns the space discretization on x direction (problem domain)*/
     double computeDx(const Point & d_min,const Point & d_max) const ;
     /**@brief computes and returns the space discretization on x direction (problem domain)*/
     double computeDy(const Point & d_min,const Point & d_max) const ;
     /**@brief computes and returns the space discretization on x direction (problem domain)*/
     double computeDz(const Point & d_min,const Point & d_max) const ;

    inline unsigned int getAlignedBlockSz() const
     {
       unsigned int tmp;
      ((m_uiSzX & ((1u<<DENDRO_BLOCK_ALIGN_FACTOR_LOG)-1))==0)? tmp=m_uiSzX : tmp=((m_uiSzX/(1u<<DENDRO_BLOCK_ALIGN_FACTOR_LOG))+1)*(1u<<DENDRO_BLOCK_ALIGN_FACTOR_LOG);
      return tmp*m_uiSzY*m_uiSzZ;
     }


   };

} // end of namespace ot





#endif //SFCSORTBENCH_BLOCK_H
