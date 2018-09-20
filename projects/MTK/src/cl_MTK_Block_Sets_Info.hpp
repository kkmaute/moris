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

namespace moris
{
namespace mtk
{
//////////////////////////
// STRUC FOR BLOCK SET  //
//////////////////////////
struct MtkBlockSetsInfo
{
    Matrix< DDUMat >*          BSetInds;
    moris::Cell< std::string >   BSetNames;

    MtkBlockSetsInfo():
        BSetInds(),
        BSetNames(){}
};
}
}



#endif /* PROJECTS_MTK_SRC_CL_MTK_BLOCK_SETS_INFO_HPP_ */
