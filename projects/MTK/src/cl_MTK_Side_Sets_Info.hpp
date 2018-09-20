/*
 * cl_MTK_Side_Sets_Info.hpp
 *
 *  Created on: Sep 17, 2018
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_SIDE_SETS_INFO_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_SIDE_SETS_INFO_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Cell.hpp"

namespace moris
{
namespace mtk
{
/////////////////////////
// STRUC FOR SIDE SET  //
/////////////////////////
struct MtkSideSetsInfo
{
    moris::Cell< Matrix< IdMat > >*    ElemIdsAndSideOrds;
    moris::Cell< std::string >     SSetNames;

    MtkSideSetsInfo():
        ElemIdsAndSideOrds(),
        SSetNames(){}
};
}
}



#endif /* PROJECTS_MTK_SRC_CL_MTK_SIDE_SETS_INFO_HPP_ */
