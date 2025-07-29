/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Mesh_Manager.cpp
 *
 */

#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Mesh_Pair.hpp"
#include "cl_MTK_Field.hpp"

#include "cl_HMR.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"

namespace moris::mtk
{
    //-------------------------------------------------------------------------

    Mesh_Manager::Mesh_Manager()
    {
    }

    //-------------------------------------------------------------------------

    Mesh_Manager::~Mesh_Manager()
    {
    }

    //-------------------------------------------------------------------------

    uint
    Mesh_Manager::register_mesh_pair(
            Interpolation_Mesh* aInterpolationMesh,
            Integration_Mesh*   aIntegrationMesh,
            bool                aIsOwned,
            std::string const & aMeshPairName )
    {
        // Create new mesh pair
        Mesh_Pair tMeshPair( aInterpolationMesh, aIntegrationMesh, false );

        // Push back new pair
        mMeshPairs.push_back( tMeshPair );

        moris_index tMeshPairIndex = mMeshPairs.size() - 1;

        // Ownership begins after adding
        mMeshPairs( tMeshPairIndex ).mIsOwned = aIsOwned;

        if ( not aMeshPairName.empty() )
        {
            mMeshPairNameToIndexMap[ aMeshPairName ] = tMeshPairIndex;
        }

        // report on this operation
        MORIS_LOG_INFO( "Registered mesh pair #%i", tMeshPairIndex );

        mMeshPairToFieldIndexMap.resize( tMeshPairIndex );

        return tMeshPairIndex;
    }

    //-------------------------------------------------------------------------

    uint
    Mesh_Manager::register_mesh_pair( Mesh_Pair& aMeshPair )
    {
        // Push back new pair
        mMeshPairs.push_back( aMeshPair );

        // Ensure correct ownership
        mMeshPairs( mMeshPairs.size() - 1 ).mIsOwned = aMeshPair.mIsOwned;
        aMeshPair.mIsOwned                           = false;

        return mMeshPairs.size() - 1;
    }

    //-------------------------------------------------------------------------

    const Mesh_Pair&
    Mesh_Manager::get_mesh_pair( moris_index aPairIndex )
    {
        MORIS_ASSERT( aPairIndex < (moris_index)mMeshPairs.size(),
                "Mesh_Manager::get_mesh_pair: requested mesh pair does not exist." );
        return mMeshPairs( aPairIndex );
    }

    //-------------------------------------------------------------------------

    const Mesh_Pair&
    Mesh_Manager::get_mesh_pair( const std::string& aMeshPairName )
    {
        MORIS_ASSERT( mMeshPairNameToIndexMap.key_exists( aMeshPairName ), "Mesh pair with name %s does not exist", aMeshPairName.c_str() );

        return mMeshPairs( mMeshPairNameToIndexMap.find( aMeshPairName ) );
    }

    //-------------------------------------------------------------------------

    void
    Mesh_Manager::remove_mesh_pair( const std::string& aMeshPairName )
    {
        if ( mMeshPairNameToIndexMap.key_exists( aMeshPairName ) )
        {
            moris_index tMeshPairIndex = mMeshPairNameToIndexMap.find( aMeshPairName );

            mMeshPairs.erase( tMeshPairIndex );

            mMeshPairNameToIndexMap.erase( aMeshPairName );
        }
    }

    //-------------------------------------------------------------------------

    void
    Mesh_Manager::get_mesh_pair(
            moris_index          aPairIndex,
            Interpolation_Mesh*& aInterpolationMesh,
            Integration_Mesh*&   aIntegrationMesh )
    {
        MORIS_ASSERT( aPairIndex < (moris_index)mMeshPairs.size(),
                "Mesh_Manager::get_mesh_pair: requested mesh pair does not exist." );

        aInterpolationMesh = mMeshPairs( aPairIndex ).get_interpolation_mesh();
        aIntegrationMesh   = mMeshPairs( aPairIndex ).get_integration_mesh();
    }

    //-------------------------------------------------------------------------

    Interpolation_Mesh*
    Mesh_Manager::get_interpolation_mesh( moris_index aMeshIndex )
    {
        MORIS_ASSERT( aMeshIndex < (moris_index)mMeshPairs.size(),
                "Mesh_Manager::get_interpolation_mesh: requested mesh does not exist." );

        return mMeshPairs( aMeshIndex ).get_interpolation_mesh();
    }

    //-------------------------------------------------------------------------

    Integration_Mesh*
    Mesh_Manager::get_integration_mesh( moris_index aMeshIndex )
    {
        MORIS_ASSERT( aMeshIndex < (moris_index)mMeshPairs.size(),
                "Mesh_Manager::get_integration_mesh: requested mesh does not exist." );

        return mMeshPairs( aMeshIndex ).get_integration_mesh();
    }

    //--------------------------------------------------------------------

    void
    Mesh_Manager::register_field( const std::shared_ptr< mtk::Field >& aField )
    {
        const std::string& tFieldLabel = aField->get_label();

        if ( not mFieldLabelToIndexMap.key_exists( tFieldLabel ) )
        {
            mFields.push_back( aField );

            // get field index in cell
            moris_index tFieldIndex = mFields.size() - 1;

            // register field index in label to index map
            mFieldLabelToIndexMap[ tFieldLabel ] = tFieldIndex;

            // get mesh pair label of field
            const std::string& tMeshPairLabel = aField->get_mesh_pair_label();

            MORIS_ERROR( mMeshPairNameToIndexMap.key_exists( tMeshPairLabel ),
                    "Mesh_Manager::register_field(), Mesh Pair does not exist" );

            // get mesh pair index of mesh pair label
            moris_index tMeshPairIndex = mMeshPairNameToIndexMap.find( tMeshPairLabel );

            // put field index in mesh pair to field index map
            mMeshPairToFieldIndexMap( tMeshPairIndex ).push_back( tFieldIndex );
        }
        else
        {
            MORIS_ASSERT( false,
                    "Mesh_Manager::register_field(): Field name was registered earlier" );
        }
    }

    //-------------------------------------------------------------------------
}    // namespace moris::mtk
