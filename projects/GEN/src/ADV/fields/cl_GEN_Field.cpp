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
#include <utility>

namespace moris::gen
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

    uint Field::get_number_of_reference_coordinates()
    {
        return 0;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Field::get_interpolated_field_value(
            const Cell< Basis_Node >& aBasisNodes,
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
            const Cell< Basis_Node >& aBasisNodes,
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

    Matrix< DDSMat > Field::get_determining_adv_ids(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        return mADVManager.get_determining_adv_ids();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Field::get_determining_adv_ids(
            Matrix< DDSMat >& aDeterminingADVIDs,
            const Derived_Node& aDerivedNode,
            const Node_Manager& aNodeManager )
    {
        aDeterminingADVIDs = this->get_determining_adv_ids( aDerivedNode.get_index(), aDerivedNode.get_global_coordinates() );
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
    void Field::append_interpolated_determining_adv_ids(
            Matrix< DDSMat >&         aInterpolatedADVIDs,
            const Cell< Basis_Node >& aBasisNodes,
            const Node_Manager&       aNodeManager )
    {
        // Add contributions from each basis node
        for ( auto iBasisNode : aBasisNodes )
        {
            // Check if we have a background node
            if ( aNodeManager.is_background_node( iBasisNode.get_index() ) )
            {
                // Get locator sensitivities
                Matrix< DDSMat > tBasisNodeADVIDs = this->get_determining_adv_ids( iBasisNode.get_index(), iBasisNode.get_global_coordinates() );

                // Get current joined ADV ID length
                uint tJoinedADVIDLength = aInterpolatedADVIDs.length();

                // Have to do a resize, since each basis node can depend on different number of ADVs
                aInterpolatedADVIDs.resize( 1, tJoinedADVIDLength + tBasisNodeADVIDs.length() );

                // Append to current list
                for ( uint iBasisNodeADV = 0; iBasisNodeADV < tBasisNodeADVIDs.length(); iBasisNodeADV++ )
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
        return mADVManager.has_advs();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Field::set_dependencies( Cell< std::shared_ptr< Field > > aDependencyFields )
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
