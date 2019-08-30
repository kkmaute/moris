/*
 * cl_MTK_Interpolation_Mesh_STK.hpp
 *
 *  Created on: Apr 15, 2019
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_INTERPOLATION_MESH_STK_HPP_
#define PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_INTERPOLATION_MESH_STK_HPP_

#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Mesh_Data_STK.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"
namespace moris
{
namespace mtk
{
class Interpolation_Mesh_STK : public Mesh_Core_STK, public Interpolation_Mesh
{
    // Functions only valid for interpolation mIntegrationMeshes

public:
    Interpolation_Mesh_STK(std::shared_ptr<Mesh_Data_STK> aSTKMeshData):
        Mesh_Core_STK(aSTKMeshData)
    {

    }

    Interpolation_Mesh_STK(
            std::string    aFileName,
            MtkMeshData*   aSuppMeshData,
            const bool     aCreateFacesAndEdges = true ):
                Mesh_Core_STK(aFileName,aSuppMeshData,aCreateFacesAndEdges)

    {

    }

    Interpolation_Mesh_STK(MtkMeshData & aMeshData ):
                Mesh_Core_STK(aMeshData)

    {

    }


    std::shared_ptr<Mesh_Data_STK>
    get_stk_data_shared_pointer()
    {
        return mSTKMeshData;
    }

    Matrix< IndexMat >
    get_elements_in_support_of_basis(moris_index aBasisIndex,
                                     moris_index aInterpIndex = 0)
    {
        return get_entity_connected_to_entity_loc_inds(aBasisIndex, EntityRank::NODE, EntityRank::ELEMENT);
    }

};
}
}



#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_INTERPOLATION_MESH_STK_HPP_ */
