/*
 * cl_MTK_Fields_Info.hpp
 *
 *  Created on: Sep 17, 2018
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_FIELDS_INFO_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_FIELDS_INFO_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Mesh_Enums.hpp"

namespace moris
{
namespace mtk
{
///////////////////////
// STRUC FOR FIELDS  //
///////////////////////
//template <typename T>
//template< typename Variant = boost::variant< bool, sint, real, const char* > >
struct MtkFieldsInfo
{
    moris::Cell< Matrix< DDRMat > >*    FieldsData;
    moris::Cell< std::string >     FieldsName;
    moris::Cell< enum EntityRank > FieldsRank;
    moris::Cell< std::string >*    SetsOwner;

    MtkFieldsInfo():
        FieldsData(),
        FieldsName(),
        FieldsRank(),
        SetsOwner() {}
};
}
}



#endif /* PROJECTS_MTK_SRC_CL_MTK_FIELDS_INFO_HPP_ */
