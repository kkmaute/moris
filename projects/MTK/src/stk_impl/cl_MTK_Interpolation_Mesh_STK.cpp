/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Interpolation_Mesh_STK.cpp
 *
 */

#include "cl_MTK_Interpolation_Mesh_STK.hpp"

#include <utility>
#include "cl_MTK_Mesh_Data_STK.hpp"

namespace moris::mtk
{
    Interpolation_Mesh_STK::Interpolation_Mesh_STK( std::shared_ptr< Mesh_Data_STK > aSTKMeshData )
            : Mesh_Core_STK( std::move( aSTKMeshData ) )
    {
    }

    //-----------------------------------------------------------------------------------------------
    Interpolation_Mesh_STK::Interpolation_Mesh_STK(
            std::string  aFileName,
            MtkMeshData* aSuppMeshData,
            const bool   aCreateFacesAndEdge )
            : Mesh_Core_STK( std::move( aFileName ), aSuppMeshData, aCreateFacesAndEdge )
    {
    }

    //-----------------------------------------------------------------------------------------------
    Interpolation_Mesh_STK::Interpolation_Mesh_STK( MtkMeshData& aMeshData )
            : Mesh_Core_STK( aMeshData )
    {
    }

    //-----------------------------------------------------------------------------------------------
    Interpolation_Mesh_STK::~Interpolation_Mesh_STK()
    {
    }

    //-----------------------------------------------------------------------------------------------
    std::shared_ptr< Mesh_Data_STK >
    Interpolation_Mesh_STK::get_stk_data_shared_pointer()
    {
        return mSTKMeshData;
    }

}    // namespace moris::mtk
