/*
 * cl_MTK_Integration_Mesh.hpp
 *
 *  Created on: Apr 15, 2019
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_INTEGRATION_MESH_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_INTEGRATION_MESH_HPP_

#include "cl_MTK_Mesh_Core.hpp"

#include "assert.hpp"
#include "cl_Matrix.hpp"

namespace moris
{
namespace mtk
{
class Integration_Mesh : public virtual Mesh
{
    // Functions only valid for integration meshes

    //##############################################
    // Mesh Sets Access
    //##############################################


//    virtual
//    moris::Cell<std::string>
//    get_set_names(enum EntityRank aSetEntityRank) const
//    {
//        MORIS_ERROR(0," get_set_names has no base implementation");
//        return moris::Cell<std::string>(0);
//    }
//
//    virtual
//    Matrix< IndexMat >
//    get_set_entity_loc_inds( enum EntityRank aSetEntityRank,
//                             std::string     aSetName) const
//                             {
//        MORIS_ERROR(0," get_set_entity_ids has no base implementation");
//        return Matrix< IndexMat >(0,0);
//                             }
//
//    virtual
//    void
//    get_sideset_elems_loc_inds_and_ords(
//            const  std::string     & aSetName,
//            Matrix< IndexMat >     & aElemIndices,
//            Matrix< IndexMat >     & aSidesetOrdinals ) const
//    {
//        MORIS_ERROR(0," get_sideset_elems_loc_inds_and_ords has no base implementation");
//    }
//
//
//    virtual
//    void
//    get_sideset_cells_and_ords(
//            const  std::string & aSetName,
//            moris::Cell< mtk::Cell const * > & aCells,
//            Matrix< IndexMat > &       aSidesetOrdinals ) const
//    {
//        MORIS_ERROR(0,"get_sideset_cells_and_ords not implemented");
//    }

};
}
}


#endif /* PROJECTS_MTK_SRC_CL_MTK_INTEGRATION_MESH_HPP_ */
