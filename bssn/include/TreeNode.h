/**
 * @file TreeNode.h
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief Basic treenode form Dendro-5
 * @version 0.1
 * @date 2020-11-23
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#pragma once 
#include <ostream>
#include "dendro.h"
#include "point.h"
extern unsigned int m_uiMaxDepth;

namespace ot
{

    class TreeNode
    {
        protected:
            unsigned int m_uiX;
            unsigned int m_uiY;
            unsigned int m_uiZ;
            unsigned int m_uiLevel;

            /**
            @brief The type of boundary
            */
            enum BoundaryType1 { NEGATIVE= 2, POSITIVE= 4};

            /**
             @brief The type of boundary
            */
            enum BoundaryType2 {
                X_NEG_BDY=1, Y_NEG_BDY=2, Z_NEG_BDY=4, NEG_POS_DEMARCATION=8,
                EXTERNAL_BDY=16, X_POS_BDY=32, Y_POS_BDY=64, Z_POS_BDY=128
            };

            /**
             @brief  The type of boundary
            */
            enum BoundaryType3 {
                FACE_BDY=1, EDGE_BDY=2, CORNER_BDY=3
            };

            enum OctantFlagType {
                MAX_LEVEL=31, BOUNDARY=64, NODE=128
            };

        public: 

            TreeNode();

            TreeNode(const unsigned int dim,const unsigned int maxDepth);

            TreeNode(const unsigned int x, const unsigned int y,const unsigned int z, const unsigned int lev, const unsigned int dim, const unsigned int maxDepth);

            inline unsigned int getX() const {
                return m_uiX;
            }

            inline unsigned int getY() const {
                return m_uiY;
            }

            inline unsigned int getZ() const {
                return m_uiZ;
            }


            inline unsigned int getDim() const {
                return m_uiDim;
            }

            inline unsigned int getMaxDepth() const {
                return m_uiMaxDepth;
            }

            inline unsigned int getLevel() const {
                return (m_uiLevel & MAX_LEVEL);
            }

            inline bool isRoot() const {
                return ((this->getLevel()==0) & (m_uiX==0) & (m_uiY==0) & (m_uiZ==0));
            }

            inline unsigned int minX() const {
                return getX();
            } //end fn.

            inline unsigned int minY() const {
                if (m_uiDim < 2) { return 0; }
                return getY();
            } //end fn.

            inline unsigned int minZ() const {
                if (m_uiDim < 3) { return 0; }
                return getZ();
            } //end fn.

            inline unsigned int maxX() const {
                unsigned int len = (1u << (m_uiMaxDepth - getLevel()));
                return (minX() + len);
            } //end fn.

            inline unsigned int maxY() const {
                if (m_uiDim < 2) { return 1; }
                unsigned int len = (1u << (m_uiMaxDepth - getLevel()));
                return (minY() + len);
            } //end fn.

            inline unsigned int maxZ() const {
                if (m_uiDim < 3) { return 1; }
                unsigned int len = (1u << (m_uiMaxDepth - getLevel()));
                return (minZ() + len);
            } //end fn.

            inline void  setFlag(unsigned int flag) {  m_uiLevel=flag;  }
            inline unsigned int getFlag() const { return m_uiLevel; }

           
    };

     std::ostream & operator << (std::ostream & os,TreeNode const & node);

}// end of namespace ot