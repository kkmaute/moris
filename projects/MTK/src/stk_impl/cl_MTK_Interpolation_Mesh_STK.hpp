/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Interpolation_Mesh_STK.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_INTERPOLATION_MESH_STK_HPP_
#define PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_INTERPOLATION_MESH_STK_HPP_

#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Mesh_Data_STK.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"

namespace moris::mtk
{
    class Interpolation_Mesh_STK : public Mesh_Core_STK
            , public Interpolation_Mesh
    {
        /// dummy variable need to avoid segfault when compiling gcc-8 and higher
        uint tDummy = 0;

      public:
        Interpolation_Mesh_STK( std::shared_ptr< Mesh_Data_STK > aSTKMeshData );

        Interpolation_Mesh_STK(
                std::string  aFileName,
                MtkMeshData* aSuppMeshData,
                const bool   aCreateFacesAndEdges = true );

        Interpolation_Mesh_STK( MtkMeshData& aMeshData );

        ~Interpolation_Mesh_STK() override;

        std::shared_ptr< Mesh_Data_STK >
        get_stk_data_shared_pointer();

        int get_dummy() { return tDummy; }
    };
}    // namespace moris::mtk

#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_INTERPOLATION_MESH_STK_HPP_ */
