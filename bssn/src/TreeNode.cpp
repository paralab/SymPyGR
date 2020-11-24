#include "TreeNode.h"
unsigned int m_uiMaxDepth=30;

namespace ot
{

    TreeNode::TreeNode(const unsigned int x, const unsigned int y,
                          const unsigned int z, const unsigned int lev, const unsigned int dim, const unsigned int maxDepth) {

        //m_uiMaxDepth = maxDepth;
        m_uiX = x;
        if (dim > 1) {m_uiY = y; } else {m_uiY = 0; }
        if (dim > 2) {m_uiZ = z; } else {m_uiZ = 0; }

        m_uiLevel = lev;
        //m_uiMaxDepth=maxDepth;

    } //end function

    TreeNode::TreeNode()
    {

      //m_uiMaxDepth=31;
      m_uiX=0; m_uiY=0; m_uiZ=0;
      m_uiLevel=0;

    }

    TreeNode:: TreeNode(const unsigned int dim,const unsigned int maxDepth)
    {
        m_uiX=0;
        m_uiY=0;
        m_uiZ=0;
        m_uiLevel=0;
        //m_uiMaxDepth=maxDepth;
    }

    std::ostream& operator<<(std::ostream& os, TreeNode const& other) {
            return (os << other.getX() << " " << other.getY() << " " << other.getZ() << " " << other.getLevel());
    } //end fn.

}
    