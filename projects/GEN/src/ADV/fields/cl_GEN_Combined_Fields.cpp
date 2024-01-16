/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Combined_Fields.cpp
 *
 */

#include "cl_GEN_Combined_Fields.hpp"
#include "cl_MTK_Field_Discrete.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Combined_Fields::Combined_Fields(
            Cell< std::shared_ptr< Field > > aFields,
            bool                             aUseMinimum,
            std::string                      aName )
            : Field( Matrix< DDRMat >{{}}, aName )
            , mFields( aFields )
            , mScale( 2 * aUseMinimum - 1 )
    {
        MORIS_ERROR( mFields.size() > 0, "A GEN Combined_Fields must be created with at least one field." );
    }

    //--------------------------------------------------------------------------------------------------------------

    real Combined_Fields::get_field_value(
            uint                  aNodeIndex,
            const Matrix<DDRMat>& aCoordinates )
    {
        real tResult = mScale * mFields( 0 )->get_field_value( aNodeIndex, aCoordinates );
        for ( const auto& iField : mFields )
        {
            tResult = std::min( tResult, mScale * iField->get_field_value( aNodeIndex, aCoordinates ) );
        }
        return mScale * tResult;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Combined_Fields::get_field_value( Derived_Node* aDerivedNode, const Node_Manager& aNodeManager )
    {
        real tResult = mScale * mFields( 0 )->get_field_value( aDerivedNode, aNodeManager );
        for ( const auto& iField : mFields )
        {
            tResult = std::min( tResult, mScale * iField->get_field_value( aDerivedNode, aNodeManager ) );
        }
        return mScale * tResult;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >& Combined_Fields::get_dfield_dadvs(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates)
    {
        // Find which field is the minimum
        real tMin = mScale * mFields( 0 )->get_field_value( aNodeIndex, aCoordinates );
        uint tMinFieldIndex = 0;
        for ( uint iFieldIndex = 1; iFieldIndex < mFields.size(); iFieldIndex++ )
        {
            real tResult = mScale * mFields( iFieldIndex )->get_field_value( aNodeIndex, aCoordinates );
            if ( tResult < tMin )
            {
                tMin = tResult;
                tMinFieldIndex = iFieldIndex;
            }
        }

        // Return relevant sensitivity
        return mFields( tMinFieldIndex )->get_dfield_dadvs( aNodeIndex, aCoordinates );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Combined_Fields::get_dfield_dadvs( Matrix< DDRMat >& aSensitivities, Derived_Node* aDerivedNode, const Node_Manager& aNodeManager )
    {
        // Find which field is the minimum
        real tMin = mScale * mFields( 0 )->get_field_value( aDerivedNode, aNodeManager );
        uint tMinFieldIndex = 0;
        for ( uint iFieldIndex = 1; iFieldIndex < mFields.size(); iFieldIndex++ )
        {
            real tResult = mScale * mFields( iFieldIndex )->get_field_value( aDerivedNode, aNodeManager );
            if ( tResult < tMin )
            {
                tMin = tResult;
                tMinFieldIndex = iFieldIndex;
            }
        }

        // Return relevant sensitivity
        mFields( tMinFieldIndex )->get_dfield_dadvs( aSensitivities, aDerivedNode, aNodeManager );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat > Combined_Fields::get_determining_adv_ids(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        // Find which field is the minimum
        real tMin = mScale * mFields( 0 )->get_field_value( aNodeIndex, aCoordinates );
        uint tMinFieldIndex = 0;
        for ( uint iFieldIndex = 1; iFieldIndex < mFields.size(); iFieldIndex++ )
        {
            real tResult = mScale * mFields( iFieldIndex )->get_field_value( aNodeIndex, aCoordinates );
            if ( tResult < tMin )
            {
                tMin = tResult;
                tMinFieldIndex = iFieldIndex;
            }
        }

        // Return relevant sensitivity
        return mFields( tMinFieldIndex )->get_determining_adv_ids( aNodeIndex, aCoordinates );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Combined_Fields::get_determining_adv_ids( Matrix< DDSMat >& aDeterminingADVIDs, Derived_Node* aDerivedNode, const Node_Manager& aNodeManager )
    {
        // Find which field is the minimum
        real tMin = mScale * mFields( 0 )->get_field_value( aDerivedNode, aNodeManager );
        uint tMinFieldIndex = 0;
        for ( uint iFieldIndex = 1; iFieldIndex < mFields.size(); iFieldIndex++ )
        {
            real tResult = mScale * mFields( iFieldIndex )->get_field_value( aDerivedNode, aNodeManager );
            if ( tResult < tMin )
            {
                tMin = tResult;
                tMinFieldIndex = iFieldIndex;
            }
        }

        // Return relevant sensitivity
        mFields( tMinFieldIndex )->get_determining_adv_ids( aDeterminingADVIDs, aDerivedNode, aNodeManager );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Combined_Fields::get_dfield_dcoordinates(
            uint                  aNodeIndex,
            const Matrix<DDRMat>& aCoordinates,
            Matrix<DDRMat>&       aSensitivities)
    {
        // Find which field is the minimum
        real tMin = mScale * mFields( 0 )->get_field_value( aNodeIndex, aCoordinates );
        uint tMinFieldIndex = 0;
        for ( uint iFieldIndex = 1; iFieldIndex < mFields.size(); iFieldIndex++ )
        {
            real tResult = mScale * mFields( iFieldIndex )->get_field_value( aNodeIndex, aCoordinates );
            if ( tResult < tMin )
            {
                tMin = tResult;
                tMinFieldIndex = iFieldIndex;
            }
        }

        // Get relevant sensitivity
        mFields( tMinFieldIndex )->get_dfield_dcoordinates( aNodeIndex, aCoordinates, aSensitivities );
    }

    //--------------------------------------------------------------------------------------------------------------
    
    std::shared_ptr< mtk::Field > Combined_Fields::get_mtk_field()
    {
        // Trivial mesh pair, since it may or may not be needed
        mtk::Mesh_Pair tMeshPair( nullptr, nullptr );

        // Test if MTK field is needed
        bool tNeedMTKField = false;
        for ( auto iField : mFields )
        {
            // Get MTK field
            auto tMTKField = iField->get_mtk_field();

            // If field exists, we have at least one discrete field and need a full MTK field
            if ( tMTKField )
            {
                // Set that we need an MTK field
                tNeedMTKField = true;

                // Grab the relevant mesh pair
                tMeshPair = tMTKField->get_mesh_pair();

                // No use in testing the rest of the fields
                break;
            }
        }

        // Create MTK field if needed
        if ( tNeedMTKField )
        {
            return this->create_mtk_field( tMeshPair );
        }
        else
        {
            return nullptr;
        }
    }
    
    //--------------------------------------------------------------------------------------------------------------

}
