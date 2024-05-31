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
#include "cl_GEN_Derived_Node.hpp"
#include "cl_GEN_Basis_Node.hpp"
#include "cl_MTK_Field_Discrete.hpp"

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    Field::Field( Vector< real >& aADVs,
            const Vector< uint >& aFieldVariableIndices,
            const Vector< uint >& aADVIndices,
            const Vector< real >& aConstants,
            std::string           aName )
            : mADVHandler( aADVs, aFieldVariableIndices, aADVIndices, aConstants )
            , mSensitivities( 1, aFieldVariableIndices.size() + aConstants.size() )
            , mName( std::move( aName ) )
    {
        this->verify_name();
    }

    //--------------------------------------------------------------------------------------------------------------

    Field::Field(
            const Vector< ADV >& aADVs,
            std::string          aName )
            : mADVHandler( aADVs )
            , mSensitivities( 1, aADVs.size(), 0.0 )
            , mName( std::move( aName ) )
    {
        this->verify_name();
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
    Field::Field(
            const Vector< sint >& aSharedADVIds,
            std::string           aName)
            : mADVHandler( aSharedADVIds )
            , mSensitivities( 1, aSharedADVIds.size() )
            , mName( std::move( aName ) )
    {
        this->verify_name();
    }

    //--------------------------------------------------------------------------------------------------------------

    Field::Field( const Field& aCopy,
            const Vector< uint >& aReplaceVariables,
            const Vector< real >& aNewConstants )
            : mADVHandler( aCopy.mADVHandler, aReplaceVariables, aNewConstants )
            , mSensitivities( aCopy.mSensitivities )
            , mName( aCopy.mName )
    {
        // Do not verify name, as we are making a copy
    }

    //--------------------------------------------------------------------------------------------------------------

    std::shared_ptr< Field > Field::copy(
            const Vector< uint >& aReplaceVariables,
            const Vector< real >& aNewConstants )
    {
        return nullptr;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Field::import_advs( sol::Dist_Vector* aOwnedADVs )
    {
        mADVHandler.import_advs( aOwnedADVs );
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Field::get_number_of_reference_coordinates()
    {
        return 0;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Field::get_interpolated_field_value(
            const Vector< Basis_Node >& aBasisNodes,
            const Node_Manager&       aNodeManager )
    {
        // Initialize field value
        real tFieldValue = 0.0;

        // Add contributions from each basis node
        for ( auto iBasisNode : aBasisNodes )
        {
            // Check if we have a background node
            if ( aNodeManager.is_background_node( iBasisNode.get_index() ) )
            {
                // Background node, we can directly add its contribution
                tFieldValue += this->get_field_value( iBasisNode.get_index(), iBasisNode.get_global_coordinates() ) * iBasisNode.get_basis();
            }
            else
            {
                // Not a background node, add contribution from its locators
                tFieldValue += this->get_interpolated_field_value( iBasisNode.get_locator_nodes(), aNodeManager ) * iBasisNode.get_basis();
            }
        }

        // Return result
        return tFieldValue;
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
    void Field::append_interpolated_dfield_dadvs(
            Matrix< DDRMat >&         aInterpolatedSensitivities,
            const Vector< Basis_Node >& aBasisNodes,
            const Node_Manager&       aNodeManager,
            real                      aBasisFactor )
    {
        // Add contributions from each basis node
        for ( auto iBasisNode : aBasisNodes )
        {
            // Check if we have a background node
            if ( aNodeManager.is_background_node( iBasisNode.get_index() ) )
            {
                // Get locator sensitivities
                Matrix< DDRMat > tBasisNodeSensitivities = this->get_dfield_dadvs( iBasisNode.get_index(), iBasisNode.get_global_coordinates() ) * iBasisNode.get_basis() * aBasisFactor;

                // Get current joined sensitivity length
                uint tJoinedSensitivityLength = aInterpolatedSensitivities.length();

                // Have to do a resize, since each basis node can depend on different number of ADVs
                aInterpolatedSensitivities.resize( 1, aInterpolatedSensitivities.length() + tBasisNodeSensitivities.length() );

                // Append to current list
                for ( uint iBasisNodeSensitivity = 0; iBasisNodeSensitivity < tBasisNodeSensitivities.length(); iBasisNodeSensitivity++ )
                {
                    aInterpolatedSensitivities( tJoinedSensitivityLength + iBasisNodeSensitivity ) = tBasisNodeSensitivities( iBasisNodeSensitivity );
                }
            }
            else
            {
                // Not a background node, so we add contribution from its locators
                this->append_interpolated_dfield_dadvs(
                        aInterpolatedSensitivities,
                        iBasisNode.get_locator_nodes(),
                        aNodeManager,
                        iBasisNode.get_basis() * aBasisFactor );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< sint > Field::get_determining_adv_ids(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        return mADVHandler.get_determining_adv_ids();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Field::get_determining_adv_ids(
            Vector< sint >&     aDeterminingADVIDs,
            const Derived_Node& aDerivedNode,
            const Node_Manager& aNodeManager )
    {
        aDeterminingADVIDs = this->get_determining_adv_ids( aDerivedNode.get_index(), aDerivedNode.get_global_coordinates() );
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
    void Field::append_interpolated_determining_adv_ids(
            Vector< sint >&             aInterpolatedADVIDs,
            const Vector< Basis_Node >& aBasisNodes,
            const Node_Manager&         aNodeManager )
    {
        // Add contributions from each basis node
        for ( auto iBasisNode : aBasisNodes )
        {
            // Check if we have a background node
            if ( aNodeManager.is_background_node( iBasisNode.get_index() ) )
            {
                // Get locator sensitivities
                Vector< sint > tBasisNodeADVIDs = this->get_determining_adv_ids( iBasisNode.get_index(), iBasisNode.get_global_coordinates() );

                // Get current joined ADV ID length
                uint tJoinedADVIDLength = aInterpolatedADVIDs.size();

                // Have to do a resize, since each basis node can depend on different number of ADVs
                aInterpolatedADVIDs.resize( tJoinedADVIDLength + tBasisNodeADVIDs.size() );

                // Append to current list
                for ( uint iBasisNodeADV = 0; iBasisNodeADV < tBasisNodeADVIDs.size(); iBasisNodeADV++ )
                {
                    aInterpolatedADVIDs( tJoinedADVIDLength + iBasisNodeADV ) = tBasisNodeADVIDs( iBasisNodeADV );
                }
            }
            else
            {
                // Not a background node, so we add contribution from its locators
                this->append_interpolated_determining_adv_ids(
                        aInterpolatedADVIDs,
                        iBasisNode.get_locator_nodes(),
                        aNodeManager );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Field::has_advs()
    {
        return mADVHandler.has_advs();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Field::update_dependencies( Vector< std::shared_ptr< Field > > aUpdatedFields )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Field::reset_nodal_data( mtk::Interpolation_Mesh* aMesh )
    {
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
        auto tMTKField = std::make_shared< mtk::Field_Discrete >( aMeshPair, gDiscretizationIndex );

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
