/* cl_MTK_Integration_Mesh_STK.hpp
 *
 *  Created on: Apr 15, 2019
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_INTEGRATION_MESH_STK_HPP_
#define PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_INTEGRATION_MESH_STK_HPP_

#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Mesh_Data_STK.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"
namespace moris
{
namespace mtk
{
class Integration_Mesh_STK : public Mesh_Core_STK, public Integration_Mesh
{
    // Functions only valid for interpolation mIntegrationMeshes

public:
    Integration_Mesh_STK(std::shared_ptr<Mesh_Data_STK> aSTKMeshData):
        Mesh_Core_STK(aSTKMeshData)
    {

    }

    Integration_Mesh_STK(
            std::string    aFileName,
            MtkMeshData*   aSuppMeshData,
            const bool     aCreateFacesAndEdges = true ):
                Mesh_Core_STK(aFileName,aSuppMeshData,aCreateFacesAndEdges)

    {

    }

    Integration_Mesh_STK(MtkMeshData & aMeshData ):
                Mesh_Core_STK(aMeshData)

    {

    }

    explicit
    Integration_Mesh_STK(Interpolation_Mesh & aInterpMesh)
    {
        MORIS_ERROR(aInterpMesh.get_mesh_type() == MeshType::STK,"operator= between an interpolation and integration mesh only valid between stk meshes");

        Interpolation_Mesh_STK* tInterpolationSTK = dynamic_cast<Interpolation_Mesh_STK*>(&aInterpMesh);

        // get the shared
        mSTKMeshData = tInterpolationSTK->get_stk_data_shared_pointer();
    }
};
}
}



#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_INTEGRATION_MESH_STK_HPP_ */
