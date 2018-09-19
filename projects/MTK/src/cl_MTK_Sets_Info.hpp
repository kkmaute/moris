/*
 * cl_MTK_Sets_Info.hpp
 *
 *  Created on: Sep 17, 2018
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_SETS_INFO_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_SETS_INFO_HPP_

#include "cl_MTK_Node_Sets_Info.hpp"
#include "cl_MTK_Side_Sets_Info.hpp"
#include "cl_MTK_Block_Sets_Info.hpp"

namespace moris
{

namespace mtk
{
///////////////////////////////
// STRUC FOR SETS CONTAINER  //
///////////////////////////////

struct MtkSetsInfo
{
    MtkNodeSetsInfo* NodeSetsInfo;
    MtkSideSetsInfo* SideSetsInfo;
    MtkBlockSetsInfo* BlockSetsInfo;

    MtkSetsInfo():
        NodeSetsInfo(),
        SideSetsInfo(),
        BlockSetsInfo(){}
};
}
}


#endif /* PROJECTS_MTK_SRC_CL_MTK_SETS_INFO_HPP_ */
