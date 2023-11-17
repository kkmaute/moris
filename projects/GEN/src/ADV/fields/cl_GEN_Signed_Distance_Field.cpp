/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Geometry_Signed_Distance_Field.cpp
 *
 */

#include "cl_GEN_Signed_Distance_Field.hpp"
#include "cl_SDF_Generator.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"

namespace moris::ge
{
    //--------------------------------------------------------------------------------------------------------------

    Signed_Distance_Field::Signed_Distance_Field(
            std::string               aObjectPath,
            Matrix< DDRMat >          aObjectOffset,
            real                      aSDFShift )
            : Field_Discrete_Integration( Matrix< DDRMat >( 1, 1, 0.0 ), 0 )
    {
        mObjectPath = aObjectPath;

        if ( aObjectOffset.numel() == 3 )
        {
            mObjectOffset = aObjectOffset;
        }

        mShift = aSDFShift;
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Signed_Distance_Field::get_field_value( uint aNodeIndex )
    {
        MORIS_ASSERT( mValues.numel() > aNodeIndex,
                "Geometry_SDF::get_field_value() - "
                "Node index out of bounds, mValues numel: %zu , aNodeIndex %i",
                mValues.numel(),
                aNodeIndex );
        return mValues( aNodeIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    //        const Matrix<DDRMat>& Geometry_SDF::get_dfield_dadvs(uint aNodeIndex)
    //        {
    //            MORIS_ERROR(false, "get_dfield_dadvs function is not implemented for a mesh field geometry.");
    //            return mSensitivities;
    //        }

    //--------------------------------------------------------------------------------------------------------------

    void
    Signed_Distance_Field::evaluate_nodal_values()
    {
        mtk::Mesh* tMesh = mMeshPair.get_interpolation_mesh();

        sdf::SDF_Generator tSDFGenerator( mObjectPath, mObjectOffset, true );

        tSDFGenerator.calculate_sdf( tMesh, mValues );

        for ( uint Ik = 0; Ik < mValues.numel(); Ik++ )
        {
            mValues( Ik ) = mValues( Ik ) + mShift;
        }

        //------------------------------------------------------
        // communicate owned to shared values
        uint tNumberOfNodes = mValues.n_rows();

        uint tCount = 0;
        // loop over all nodes and build owned map
        for ( uint tNodeIndex = 0; tNodeIndex < mValues.n_rows(); ++tNodeIndex )
        {
            // get node owner
            sint tNodeOwner = tMesh->get_entity_owner( tNodeIndex, mtk::EntityRank::NODE );

            // process only owned nodes
            if ( par_rank() == tNodeOwner )
            {
                tCount++;
            }
        }

        // set size nodal vectors
        Matrix< DDSMat > tOwnedNodeIDs( tCount, 1 );
        Matrix< DDSMat > tSharedNodeIDs( tNumberOfNodes, 1 );

        tCount = 0;
        // loop over all nodes and build owned map
        for ( uint tNodeIndex = 0; tNodeIndex < mValues.n_rows(); ++tNodeIndex )
        {
            // get node owner
            sint tNodeOwner = tMesh->get_entity_owner( tNodeIndex, mtk::EntityRank::NODE );

            // get node ID
            sint tNodeID = tMesh->get_glb_entity_id_from_entity_loc_index(
                    tNodeIndex,
                    mtk::EntityRank::NODE );

            tSharedNodeIDs( tNodeIndex ) = tNodeID;

            // process only owned nodes
            if ( par_rank() == tNodeOwner )
            {
                tOwnedNodeIDs( tCount++ ) = tNodeID;
            }
        }

        // create distributed vectors for node value computation
        sol::Matrix_Vector_Factory tDistributedFactory;

        sol::Dist_Map* tOwnedNodeMap  = tDistributedFactory.create_map( tOwnedNodeIDs );
        sol::Dist_Map* tSharedNodeMap = tDistributedFactory.create_map( tSharedNodeIDs );

        sol::Dist_Vector* tOwnedNodalValues  = tDistributedFactory.create_vector( tOwnedNodeMap, 1, false, true );
        sol::Dist_Vector* tSharedNodalValues = tDistributedFactory.create_vector( tSharedNodeMap, 1, false, true );

        tOwnedNodalValues->vec_put_scalar( 0.0 );
        tSharedNodalValues->vec_put_scalar( 0.0 );

        // loop over all nodes
        for ( uint tNodeIndex = 0; tNodeIndex < mValues.n_rows(); ++tNodeIndex )
        {
            // get node owner
            sint tNodeOwner = tMesh->get_entity_owner( tNodeIndex, mtk::EntityRank::NODE );

            // process only owned nodes
            if ( par_rank() == tNodeOwner )
            {
                // get node ID
                sint tNodeID = tMesh->get_glb_entity_id_from_entity_loc_index(
                        tNodeIndex,
                        mtk::EntityRank::NODE );

                // copy nodal value on distributed vector
                ( *tOwnedNodalValues )( tNodeID ) = mValues( tNodeIndex );
            }
        }

        // get values for shared and owned nodes
        tSharedNodalValues->import_local_to_global( *tOwnedNodalValues );

        // copy shared and own nodal values onto local nodal field
        for ( uint tNodeIndex = 0; tNodeIndex < mValues.n_rows(); ++tNodeIndex )
        {
            // get node ID
            sint tNodeID = tMesh->get_glb_entity_id_from_entity_loc_index(
                    tNodeIndex,
                    mtk::EntityRank::NODE );

            // extract nodal value
            real tValue = ( *tSharedNodalValues )( tNodeID );

            // apply nodal value to field
            mValues( tNodeIndex ) = tValue;
        }

        delete tOwnedNodalValues;
        delete tSharedNodalValues;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Signed_Distance_Field::reset_nodal_data()
    {
        // Reset child nodes
        Field_Discrete_Integration::reset_nodal_data();

        // Re-evaluate field values
        this->evaluate_nodal_values();
    }

    //--------------------------------------------------------------------------------------------------------------

}
