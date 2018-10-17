/*
 * cl_MTK_Node_Sets_Info.hpp
 *
 *  Created on: Sep 17, 2018
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_NODE_SETS_INFO_HPP_
#define PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_NODE_SETS_INFO_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Cell.hpp"
namespace moris
{
namespace mtk
{
//////////////////////////
// STRUC FOR NODE SET  //
//////////////////////////
struct MtkNodeSetsInfo
{
    moris::Cell< Matrix< IdMat > >*    EntIds;
    moris::Cell< std::string >     NSetNames;

    MtkNodeSetsInfo():
        EntIds(),
        NSetNames(){}
};
}
}




#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_NODE_SETS_INFO_HPP_ */
