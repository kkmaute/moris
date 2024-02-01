/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_BSpline_Field.cpp
 *
 */

#include "cl_GEN_BSpline_Field.hpp"
#include "cl_MTK_Mesh_Pair.hpp"
#include "cl_MTK_Field.hpp"
#include "cl_MTK_Field_Discrete.hpp"
#include "cl_MTK_Mapper.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Dist_Map.hpp"
#include "fn_trans.hpp"

// Logging package
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

namespace moris::gen
{
    //--------------------------------------------------------------------------------------------------------------

    BSpline_Field::BSpline_Field(
            mtk::Mesh_Pair           aMeshPair,
            sol::Dist_Vector*        aOwnedADVs,
            const Matrix< DDSMat >&  aSharedADVIds,
            uint                     aADVOffsetID,
            uint                     aDiscretizationIndex,
            std::shared_ptr< Field > aField )
            : Field_Discrete_Integration( aSharedADVIds, aMeshPair, aField->get_name() )
            , mADVOffsetID( aADVOffsetID )
            , mMeshPair( aMeshPair )
            , mDiscretizationIndex( aDiscretizationIndex )
    {
        // Map to B-splines
        Matrix< DDRMat > tTargetField = this->map_to_bsplines( aField );

        // Distribute coefficients
        this->distribute_coeffs(
                tTargetField,
                aOwnedADVs,
                aSharedADVIds.length() );
    }

    //--------------------------------------------------------------------------------------------------------------

    BSpline_Field::BSpline_Field(
            sol::Dist_Vector*             aOwnedADVs,
            const Matrix< DDSMat >&       aSharedADVIds,
            uint                          aADVOffsetID,
            uint                          aDiscretizationIndex,
            std::shared_ptr< mtk::Field > aMTKField,
            mtk::Mesh_Pair                aMeshPair )
            : Field_Discrete_Integration( aSharedADVIds, aMeshPair, aMTKField->get_label() )
            , mADVOffsetID( aADVOffsetID )
            , mMeshPair( aMeshPair )
            , mDiscretizationIndex( aDiscretizationIndex )
    {
        // Map to B-splines
        const Matrix< DDRMat >& tTargetField = aMTKField->get_coefficients();

        // Distribute coefficients
        this->distribute_coeffs(
                tTargetField,
                aOwnedADVs,
                aSharedADVIds.length() );
    }

    //--------------------------------------------------------------------------------------------------------------

    BSpline_Field::~BSpline_Field()
    {
        delete mOwnedNodalValues;
        delete mSharedNodalValues;
    }
    //--------------------------------------------------------------------------------------------------------------

    void
    BSpline_Field::distribute_coeffs(
            const Matrix< DDRMat >& aTargetField,
            sol::Dist_Vector*       aOwnedADVs,
            uint                    aNumberOfCoefficients )
    {
        // Get B-spline mesh index
        mtk::Mesh* tMesh                    = mMeshPair.get_interpolation_mesh();
        uint       tDiscretizationMeshIndex = mDiscretizationIndex;

        // Assign ADVs
        for ( uint iCoefficientIndex = 0; iCoefficientIndex < aNumberOfCoefficients; iCoefficientIndex++ )
        {
            if ( (uint)par_rank() == tMesh->get_entity_owner( iCoefficientIndex, mtk::EntityRank::BSPLINE, tDiscretizationMeshIndex ) )
            {
                // Calculate ADV ID using offset
                sint tADVId = mADVOffsetID +                                     //
                              tMesh->get_glb_entity_id_from_entity_loc_index(    //
                                      iCoefficientIndex,
                                      mtk::EntityRank::BSPLINE,
                                      tDiscretizationMeshIndex );

                // Assign distributed vector element based on ID
                ( *aOwnedADVs )( tADVId ) = aTargetField( iCoefficientIndex );
            }
        }

        // Determine number of owned and shared node IDs
        uint tOwnedNodeCount = 0;
        for ( uint tNodeIndex = 0; tNodeIndex < tMesh->get_num_nodes(); tNodeIndex++ )
        {
            if ( (uint)par_rank() == tMesh->get_entity_owner( tNodeIndex, mtk::EntityRank::NODE, tDiscretizationMeshIndex ) )
            {
                tOwnedNodeCount++;
            }
        }
        Matrix< DDSMat > tOwnedNodeIDs( tOwnedNodeCount, 1 );
        Matrix< DDSMat > tSharedNodeIDs( tMesh->get_num_nodes(), 1 );

        // Assign owned and shared node IDs
        tOwnedNodeCount = 0;
        for ( uint tNodeIndex = 0; tNodeIndex < tMesh->get_num_nodes(); tNodeIndex++ )
        {
            sint tNodeID = tMesh->get_glb_entity_id_from_entity_loc_index(
                    tNodeIndex,
                    mtk::EntityRank::NODE,
                    tDiscretizationMeshIndex );
            tSharedNodeIDs( tNodeIndex ) = tNodeID;
            if ( (uint)par_rank() == tMesh->get_entity_owner( tNodeIndex, mtk::EntityRank::NODE, tDiscretizationMeshIndex ) )
            {
                tOwnedNodeIDs( tOwnedNodeCount++ ) = tNodeID;
            }
        }

        // Create owned and shared distributed vectors
        sol::Matrix_Vector_Factory tDistributedFactory;
        sol::Dist_Map*             tOwnedNodeMap  = tDistributedFactory.create_map( tOwnedNodeIDs );
        sol::Dist_Map*             tSharedNodeMap = tDistributedFactory.create_map( tSharedNodeIDs );
        mOwnedNodalValues                         = tDistributedFactory.create_vector( tOwnedNodeMap, 1, false, true );
        mSharedNodalValues                        = tDistributedFactory.create_vector( tSharedNodeMap, 1, false, true );

        // Import ADVs and assign nodal values
        this->import_advs( aOwnedADVs );
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    BSpline_Field::get_field_value( uint aNodeIndex )
    {
        mtk::Interpolation_Mesh* tMesh = mMeshPair.get_interpolation_mesh();
        sint tNodeID = tMesh->get_glb_entity_id_from_entity_loc_index(
                (moris_index) aNodeIndex,
                mtk::EntityRank::NODE,
                (moris_index) mDiscretizationIndex );

        return ( *mSharedNodalValues )( tNodeID );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    BSpline_Field::get_dfield_dadvs( uint aNodeIndex )
    {
        mtk::Mesh* tMesh = mMeshPair.get_interpolation_mesh();

        moris::mtk::Vertex& tVertex = tMesh->get_mtk_vertex( aNodeIndex );

        if ( tVertex.has_interpolation( mDiscretizationIndex ) )
        {
            mSensitivities = trans( tMesh->get_t_matrix_of_node_loc_ind( aNodeIndex, mDiscretizationIndex ) );
        }
        else
        {
            // return empty matrix when no T-matrix has been found
            // FIXME better solution: ensure that all T-matrices are available
            mSensitivities = Matrix< DDRMat >( 0, 0 );
        }
        return mSensitivities;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat >
    BSpline_Field::get_determining_adv_ids( uint aNodeIndex )
    {
        mtk::Mesh* tMesh = mMeshPair.get_interpolation_mesh();

        moris::mtk::Vertex& tVertex = tMesh->get_mtk_vertex( aNodeIndex );

        if ( tVertex.has_interpolation( mDiscretizationIndex ) )
        {
            return trans( tMesh->get_coefficient_IDs_of_node( aNodeIndex, mDiscretizationIndex ) ) + mADVOffsetID;
        }
        else
        {
            // return empty matrix when no T-matrix has been found
            // FIXME better solution: ensure that all T-matrices are available
            return Matrix< DDSMat >( 0, 0 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    BSpline_Field::import_advs( sol::Dist_Vector* aOwnedADVs )
    {
        // Import ADVs as usual
        mADVManager.import_advs( aOwnedADVs );

        // Reset evaluated field
        mOwnedNodalValues->vec_put_scalar( 0 );

        // Evaluate into matrix TODO just make the field accept a cell to begin with
        uint tNumberOfCoefficients = mADVManager.get_determining_adv_ids().length();
        moris::Matrix< DDRMat > tCoeff( tNumberOfCoefficients, 1 );
        for ( uint Ik = 0; Ik < tNumberOfCoefficients; Ik++ )
        {
            tCoeff( Ik ) = mADVManager.get_variable( Ik );
        }

        // Create field
        auto tField = new mtk::Field_Discrete( mMeshPair, mDiscretizationIndex );

        // Set coefficients
        tField->unlock_field();
        tField->set_coefficients( tCoeff );

        // Use mapper
        mtk::Mapper tMapper;
        tMapper.perform_mapping( tField, mtk::EntityRank::BSPLINE, mtk::EntityRank::NODE );

        // Get coefficients
        Matrix< DDRMat > tNodalValues = tField->get_values();

        // Clean up
        delete tField;

        // Evaluate field at owned nodes
        mtk::Mesh* tMesh = mMeshPair.get_interpolation_mesh();
        for ( uint tNodeIndex = 0; tNodeIndex < tMesh->get_num_nodes(); tNodeIndex++ )
        {
            if ( (uint)par_rank() == tMesh->get_entity_owner( tNodeIndex, mtk::EntityRank::NODE ) )
            {
                sint tNodeID = tMesh->get_glb_entity_id_from_entity_loc_index(
                        tNodeIndex,
                        mtk::EntityRank::NODE );
                MORIS_ASSERT( tNodalValues( tNodeIndex ) != MORIS_REAL_MAX, "value is MORIS_REAL_MAX, check mapper" );
                ( *mOwnedNodalValues )( tNodeID ) = tNodalValues( tNodeIndex );
            }
        }

        // Global assembly
        mOwnedNodalValues->vector_global_assembly();

        // Import nodal values
        mSharedNodalValues->import_local_to_global( *mOwnedNodalValues );
    }

    //--------------------------------------------------------------------------------------------------------------

    std::shared_ptr< mtk::Field > BSpline_Field::get_mtk_field()
    {
        // Output field
        auto tMTKField = std::make_shared< mtk::Field_Discrete >( mMeshPair, mDiscretizationIndex );

        // Set coefficient vector TODO vector instead of matrix, like above
        uint tNumberOfCoefficients = mADVManager.get_determining_adv_ids().length();
        moris::Matrix< DDRMat > tCoefficients( tNumberOfCoefficients, 1 );
        for ( uint iCoefficientIndex = 0; iCoefficientIndex < tNumberOfCoefficients; iCoefficientIndex++ )
        {
            tCoefficients( iCoefficientIndex ) = mADVManager.get_variable( iCoefficientIndex );
        }

        // Set coefficients
        tMTKField->unlock_field();
        tMTKField->set_coefficients( tCoefficients );

        // Set name
        tMTKField->set_label( this->get_name() );

        return tMTKField;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    BSpline_Field::map_to_bsplines( std::shared_ptr< Field > aField )
    {
        // Mapper
        mtk::Mapper tMapper;

        // New mesh
        mtk::Interpolation_Mesh* tMesh = mMeshPair.get_interpolation_mesh();

        // Output field
        auto tOutputField = new mtk::Field_Discrete( mMeshPair, mDiscretizationIndex );

        // Nodal values
        Matrix< DDRMat > tNodalValues( tMesh->get_num_nodes(), 1 );
        for ( uint tNodeIndex = 0; tNodeIndex < tMesh->get_num_nodes(); tNodeIndex++ )
        {
            tNodalValues( tNodeIndex ) =
                    aField->get_field_value( tNodeIndex, tMesh->get_node_coordinate( tNodeIndex ) );
        }
        tOutputField->unlock_field();
        tOutputField->set_values( tNodalValues );

        // Perform mapping
        tMapper.perform_mapping( tOutputField, mtk::EntityRank::NODE, mtk::EntityRank::BSPLINE );

        // Get coefficients
        Matrix< DDRMat > tCoefficients = tOutputField->get_coefficients();

        // Clean up
        delete tOutputField;

        // Return mapped field
        return tCoefficients;
    }

    //--------------------------------------------------------------------------------------------------------------

}
