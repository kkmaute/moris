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
        class Interpolation_Mesh_STK : public Mesh_Core_STK
                , public Interpolation_Mesh
        {
            // Functions only valid for interpolation mIntegrationMeshes

          public:
            Interpolation_Mesh_STK( std::shared_ptr< Mesh_Data_STK > aSTKMeshData );

            Interpolation_Mesh_STK(
                    std::string  aFileName,
                    MtkMeshData* aSuppMeshData,
                    const bool   aCreateFacesAndEdges = true );

            Interpolation_Mesh_STK( MtkMeshData& aMeshData );

            ~Interpolation_Mesh_STK();

            std::shared_ptr< Mesh_Data_STK >
            get_stk_data_shared_pointer();
        };
    }    // namespace mtk
}    // namespace moris
#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_INTERPOLATION_MESH_STK_HPP_ */
