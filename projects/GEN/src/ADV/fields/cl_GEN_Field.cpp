/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Field.cpp
 *
 */

#include "cl_GEN_Field.hpp"
#include "cl_GEN_Node_Manager.hpp"
#include "cl_GEN_Derived_Node.hpp"
#include "cl_GEN_Basis_Node.hpp"
#include "cl_MTK_Field_Discrete.hpp"
#include <utility>

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Field::Field(
            Matrix< DDRMat >& aADVs,
            Matrix< DDUMat >  aFieldVariableIndices,
            Matrix< DDUMat >  aADVIndices,
            Matrix< DDRMat >  aConstants,
            std::string       aName )
            : mADVManager( aADVs, aFieldVariableIndices, aADVIndices, aConstants )
            , mSensitivities( 1, aFieldVariableIndices.length() + aConstants.length() )
            , mName( std::move( aName ) )
    {
        this->verify_name();
    }

    //--------------------------------------------------------------------------------------------------------------

    Field::Field(
            Matrix< DDRMat > aConstants,
            std::string      aName )
            : mADVManager( aConstants )
            , mSensitivities( 1, aConstants.length(), 0.0 )
            , mName( aName )
    {
        this->verify_name();
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
    Field::Field(
            const Matrix< DDSMat >& aSharedADVIds,
            std::string             aName)
            : mADVManager( aSharedADVIds )
            , mSensitivities( 1, aSharedADVIds.length() )
            , mName( std::move( aName ) )
    {
        this->verify_name();
    }

    //--------------------------------------------------------------------------------------------------------------

    Field::Field( const Field& aCopy,
            const Cell< uint >& aReplaceVariables,
            const Cell< real >& aNewConstants )
            : mADVManager( aCopy.mADVManager, aReplaceVariables, aNewConstants )
            , mSensitivities( aCopy.mSensitivities )
            , mName( aCopy.mName )
    {
        // Do not verify name, as we are making a copy
    }

    //--------------------------------------------------------------------------------------------------------------

    std::shared_ptr< Field > Field::copy(
            const Cell< uint >& aReplaceVariables,
            const Cell< real >& aNewConstants )
    {
        return nullptr;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Field::import_advs( sol::Dist_Vector* aOwnedADVs )
    {
        mADVManager.import_advs( aOwnedADVs );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Field::set_node_manager( Node_Manager& aNodeManager )
    {
        mNodeManager = aNodeManager;
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Field::get_number_of_reference_coordinates()
    {
        return 0;
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
    real Field::get_interpolated_field_value(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        if ( aNodeIndex < mNodeManager.get_number_of_base_nodes() )
        {
            // Get base field value
            return this->get_field_value( aNodeIndex, aCoordinates );
        }
        else
        {
            // Initialize field value
            real tFieldValue = 0.0;
            
            // Add contributions from each locator
            for ( auto iBasisNode : mNodeManager.get_derived_node( aNodeIndex )->get_basis_nodes() )
            {
                tFieldValue += this->get_field_value( iBasisNode.get_index(), iBasisNode.get_global_coordinates() ) * iBasisNode.get_basis();
            }
            
            // Return result
            return tFieldValue;
        }
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
    const Matrix< DDRMat >& Field::get_interpolated_dfield_dadvs(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        if ( aNodeIndex < mNodeManager.get_number_of_base_nodes() )
        {
            // Get base sensitivities
            return this->get_dfield_dadvs( aNodeIndex, aCoordinates );
        }
        else
        {
            // Add contributions from each locator
            for ( auto iBasisNode : mNodeManager.get_derived_node( aNodeIndex )->get_basis_nodes() )
            {
                // Get locator sensitivities
                Matrix< DDRMat > tLocatorSensitivities = this->get_dfield_dadvs( iBasisNode.get_index(), iBasisNode.get_global_coordinates() ) * iBasisNode.get_basis();

                // Have to do a resize, since each locator can depend on different number of ADVs
                mInterpolatedSensitivities.resize( 1, mInterpolatedSensitivities.length() + tLocatorSensitivities.length() );

                // Get current joined sensitivity length
                uint tJoinedSensitivityLength = mInterpolatedSensitivities.length();

                // Append to current list
                for ( uint iBasisNodeSensitivity = 0; iBasisNodeSensitivity < tLocatorSensitivities.length(); iBasisNodeSensitivity++ )
                {
                    mInterpolatedSensitivities( tJoinedSensitivityLength + iBasisNodeSensitivity ) = tLocatorSensitivities( iBasisNodeSensitivity );
                }
            }

            // Return result
            return mInterpolatedSensitivities;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat > Field::get_determining_adv_ids(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        return mADVManager.get_determining_adv_ids();
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat > Field::get_interpolated_determining_adv_ids(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        if ( aNodeIndex < mNodeManager.get_number_of_base_nodes() )
        {
            // Get base sensitivities
            return this->get_determining_adv_ids( aNodeIndex, aCoordinates );
        }
        else
        {
            // Add contributions from each locator
            for ( auto iBasisNode : mNodeManager.get_derived_node( aNodeIndex )->get_basis_nodes() )
            {
                // Get locator sensitivities
                Matrix< DDSMat > tLocatorADVIDs = this->get_determining_adv_ids( iBasisNode.get_index(), iBasisNode.get_global_coordinates() ) * iBasisNode.get_basis();

                // Have to do a resize, since each locator can depend on different number of ADVs
                mInterpolatedADVIDs.resize( 1, mInterpolatedADVIDs.length() + tLocatorADVIDs.length() );

                // Get current joined ADV ID length
                uint tJoinedADVIDLength = mInterpolatedADVIDs.length();

                // Append to current list
                for ( uint iBasisNodeADV = 0; iBasisNodeADV < tLocatorADVIDs.length(); iBasisNodeADV++ )
                {
                    mInterpolatedADVIDs( tJoinedADVIDLength + iBasisNodeADV ) = tLocatorADVIDs( iBasisNodeADV );
                }
            }

            // Return result
            return mInterpolatedADVIDs;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Field::has_advs()
    {
        return mADVManager.has_advs();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Field::set_dependencies( Cell< std::shared_ptr< Field > > aDependencyFields )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Field::add_child_node( uint aNodeIndex, std::shared_ptr< Child_Node > aChildNode )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Field::add_nodal_data( mtk::Interpolation_Mesh* aMesh )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Field::reset_nodal_data()
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Field::set_num_original_nodes( uint aNumOriginalNodes )
    {
        mNumOriginalNodes = aNumOriginalNodes;
    }

    //--------------------------------------------------------------------------------------------------------------

    std::string Field::get_name()
    {
        return mName;
    }

    //--------------------------------------------------------------------------------------------------------------

    std::shared_ptr< mtk::Field > Field::create_mtk_field( const mtk::Mesh_Pair& aMeshPair )
    {
        // Output field
        auto tMTKField = std::make_shared< mtk::Field_Discrete >( aMeshPair );

        // Get interpolation mesh
        mtk::Interpolation_Mesh* tInterpolationMesh = aMeshPair.get_interpolation_mesh();

        // Get nodal values
        uint tNumberOfNodes = tInterpolationMesh->get_num_nodes();
        Matrix< DDRMat > tNodalValues( tNumberOfNodes, 1 );
        for ( uint tNodeIndex = 0; tNodeIndex < tNumberOfNodes; tNodeIndex++ )
        {
            tNodalValues( tNodeIndex ) =
                    this->get_field_value( tNodeIndex, tInterpolationMesh->get_node_coordinate( tNodeIndex ) );
        }

        // Set nodal values
        tMTKField->unlock_field();
        tMTKField->set_values( tNodalValues );

        // Set name
        tMTKField->set_label( mName );

        return tMTKField;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Field::verify_name()
    {
        // Assign default name if needed
        if ( mName.empty() )
        {
            mName = "Field " + std::to_string( Field::mCount );
        }

        // Increment count
        Field::mCount++;
    }

    //--------------------------------------------------------------------------------------------------------------
}
