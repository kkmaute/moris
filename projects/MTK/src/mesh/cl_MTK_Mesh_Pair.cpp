/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Mesh_Pair.cpp
 *
 */

#include "cl_MTK_Mesh_Pair.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"

namespace moris::mtk
{

    //--------------------------------------------------------------------------------------------------------------

    Mesh_Pair::Mesh_Pair(
            Interpolation_Mesh* aInterpolationMesh,
            Integration_Mesh*   aIntegrationMesh,
            bool                aIsOwned )
            : mInterpolationMesh( aInterpolationMesh )
            , mIntegrationMesh( aIntegrationMesh )
            , mIsOwned( aIsOwned )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Mesh_Pair::Mesh_Pair( const Mesh_Pair& aMeshPair )
            : mInterpolationMesh( aMeshPair.mInterpolationMesh )
            , mIntegrationMesh( aMeshPair.mIntegrationMesh )
            , mIsOwned( false )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Mesh_Pair::~Mesh_Pair()
    {
        if ( mIsOwned )
        {
            delete mInterpolationMesh;
            delete mIntegrationMesh;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Interpolation_Mesh* Mesh_Pair::get_interpolation_mesh() const
    {
        MORIS_ASSERT( mInterpolationMesh, "Interpolation mesh does not exist." );
        return mInterpolationMesh;
    }

    //--------------------------------------------------------------------------------------------------------------

    Integration_Mesh* Mesh_Pair::get_integration_mesh() const
    {
        MORIS_ASSERT( mIntegrationMesh, "Integration mesh does not exist." );
        return mIntegrationMesh;
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::mtk
