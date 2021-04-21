/*
 * cl_MTK_Block_Sets_Info.hpp
 *
 *  Created on: Sep 17, 2018
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_BLOCK_SETS_INFO_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_BLOCK_SETS_INFO_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_Mesh_Enums.hpp"
namespace moris
{
namespace mtk
{
//////////////////////////
// STRUC FOR BLOCK SET  //
//////////////////////////
/*
 * A block set contains the following information:
 *
 * Cell indices in a block - vector mCellIndicesInSet
 * Block set name          - name of block set
 *
 * A block set contains cells of the same topology, this
 * is mostly a restriction for visualization.
 *
 */

struct MtkBlockSetInfo
{
    Matrix< IdMat >*  mCellIdsInSet;
    std::string       mBlockSetName;
    enum CellTopology mBlockSetTopo;
    bool              mParallelConsistencyReq = true;

    MtkBlockSetInfo():
        mCellIdsInSet(),
        mBlockSetName(){}

    bool
    blockset_has_name()
    {
        return !mBlockSetName.empty();
    }


};

inline
std::ostream &
operator<<(std::ostream & os, mtk::MtkBlockSetInfo const * const & dt)
{
    os<<"Block Set Name: "<< dt->mBlockSetName << " | Block Cell Topology: "<<get_enum_str(dt->mBlockSetTopo)<<" | Number of Cells: "<<dt->mCellIdsInSet->numel();
    return os;
}

}
}



#endif /* PROJECTS_MTK_SRC_CL_MTK_BLOCK_SETS_INFO_HPP_ */
