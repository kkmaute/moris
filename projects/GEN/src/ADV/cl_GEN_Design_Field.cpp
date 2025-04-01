/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Field.cpp
 *
 */

#include "cl_GEN_Design_Field.hpp"

#include <utility>
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_GEN_BSpline_Field.hpp"
#include "cl_GEN_Stored_Field.hpp"

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    Field_Parameters::Field_Parameters( const Parameter_List& aParameterList )
            : mDiscretizationIndex( aParameterList.get< sint >( "discretization_mesh_index" ) )
            , mDiscretizationLowerBound( aParameterList.get< real >( "discretization_lower_bound" ) )
            , mDiscretizationUpperBound( aParameterList.get< real >( "discretization_upper_bound" ) )
            , mUseMultilinearInterpolation( aParameterList.get< bool >( "use_multilinear_interpolation" ) )
            , mDelaunay( aParameterList.get< bool >( "delaunay" ) )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Design_Field::Design_Field(
            std::shared_ptr< Field > aField,
            const Field_Parameters&  aParameters,
            Node_Manager&            aNodeManager )
            : mField( std::move( aField ) )
            , mNodeManager( &aNodeManager )
            , mParameters( aParameters )
    {
        // Check for a provided field
        MORIS_ERROR( mField, "A design must be provided a field for computing values." );

        // Override linear interpolation if multilinear intersections will be used
        if ( mParameters.mUseMultilinearInterpolation )
        {
            Derived_Node::set_override_linear_interpolation();
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Design_Field::discretize(
            const mtk::Mesh_Pair& aMeshPair,
            sol::Dist_Vector*     aOwnedADVs,
            const Vector< sint >& aSharedADVIds,
            uint                  aADVOffsetID )
    {
        if ( mParameters.mDiscretizationIndex >= 0 )
        {
            // Create a B-spline field
            mField = std::make_shared< BSpline_Field >(
                    aMeshPair,
                    aOwnedADVs,
                    aSharedADVIds,
                    aADVOffsetID,
                    mParameters.mDiscretizationIndex,
                    mField );

            // Set analytic field index, for now
            Field::gDiscretizationIndex = mParameters.mDiscretizationIndex;
        }
        else if ( mParameters.mDiscretizationIndex == -1 )
        {
            // Just store nodal values
            mField = std::make_shared< Stored_Field >(
                    aMeshPair.get_interpolation_mesh(),
                    mField );
        }
        mField->mMeshPairForAnalytic = aMeshPair;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Design_Field::discretize(
            const std::shared_ptr< mtk::Field >& aMTKField,
            const mtk::Mesh_Pair&                aMeshPair,
            sol::Dist_Vector*                    aOwnedADVs,
            const Vector< sint >&                aSharedADVIds,
            uint                                 aADVOffsetID )
    {
        if ( mParameters.mDiscretizationIndex >= 0 )
        {
            // Create a B-spline field
            mField = std::make_shared< BSpline_Field >(
                    aOwnedADVs,
                    aSharedADVIds,
                    aADVOffsetID,
                    mParameters.mDiscretizationIndex,
                    aMTKField,
                    aMeshPair );
        }
        else if ( mParameters.mDiscretizationIndex == -1 )
        {
            // TODO
            MORIS_ERROR( false, "Stored field cannot be remeshed for now" );
        }
        mField->mMeshPairForAnalytic = aMeshPair;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Design_Field::get_field_value(
            const uint              aNodeIndex,
            const Matrix< DDRMat >& aCoordinates ) const
    {
        if ( mNodeManager->is_background_node( aNodeIndex ) )
        {
            // Get background node field value
            return mField->get_field_value( aNodeIndex, aCoordinates );
        }
        else
        {
            // Get derived node
            const Node_Manager& tNodeManager = *mNodeManager;
            const Derived_Node& tDerivedNode = tNodeManager.get_derived_node( aNodeIndex );

            // Determine how to perform interpolation
            if ( mParameters.mUseMultilinearInterpolation )
            {
                // If we use multilinear interpolation, it is needed for all derived nodes
                return mField->get_interpolated_field_value(
                        tDerivedNode.get_background_nodes(),
                        tNodeManager );
            }
            else
            {
                // Let field decide the value
                return mField->get_field_value(
                        tDerivedNode,
                        tNodeManager );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >& Design_Field::get_dfield_dadvs(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        if ( mNodeManager->is_background_node( aNodeIndex ) )
        {
            // Get background node sensitivities
            return mField->get_dfield_dadvs( aNodeIndex, aCoordinates );
        }
        else
        {
            // Derived node, so we use interpolated vector
            mInterpolatedSensitivities.set_size( 0, 0 );

            // Get derived node
            const Node_Manager& tNodeManager = *mNodeManager;
            const Derived_Node& tDerivedNode = tNodeManager.get_derived_node( aNodeIndex );

            // Determine how to perform interpolation
            if ( mParameters.mUseMultilinearInterpolation )
            {
                // If we use multilinear interpolation, it is needed for all derived nodes
                mField->append_interpolated_dfield_dadvs(
                        mInterpolatedSensitivities,
                        tDerivedNode.get_background_nodes(),
                        tNodeManager );
            }
            else
            {
                // Let field decide the sensitivities
                mField->get_dfield_dadvs(
                        mInterpolatedSensitivities,
                        tDerivedNode,
                        tNodeManager );
            }

            // Return result
            return mInterpolatedSensitivities;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< sint > Design_Field::get_determining_adv_ids(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        if ( mNodeManager->is_background_node( aNodeIndex ) )
        {
            // Get background node determining ADV IDs
            return mField->get_determining_adv_ids( aNodeIndex, aCoordinates );
        }
        else
        {
            // Derived node, so we might need interpolation
            mInterpolatedADVIDs.clear();

            // Get derived node
            const Node_Manager& tNodeManager = *mNodeManager;
            const Derived_Node& tDerivedNode = tNodeManager.get_derived_node( aNodeIndex );

            // Determine how to perform interpolation
            if ( mParameters.mUseMultilinearInterpolation )
            {
                // If we use multilinear interpolation, it is needed for all derived nodes
                mField->append_interpolated_determining_adv_ids(
                        mInterpolatedADVIDs,
                        tDerivedNode.get_background_nodes(),
                        tNodeManager );
            }
            else
            {
                // Let field provide the IDs
                mField->get_determining_adv_ids(
                        mInterpolatedADVIDs,
                        tDerivedNode,
                        tNodeManager );
            }

            // Return result
            return mInterpolatedADVIDs;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Design_Field::import_advs( sol::Dist_Vector* aOwnedADVs )
    {
        mField->import_advs( aOwnedADVs );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Design_Field::reset_nodal_data( mtk::Interpolation_Mesh* aMesh )
    {
        mField->reset_nodal_data( aMesh );
    }

    //--------------------------------------------------------------------------------------------------------------

    std::string
    Design_Field::get_name()
    {
        return mField->get_name();
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Design_Field::use_multilinear_interpolation() const
    {
        return mParameters.mUseMultilinearInterpolation;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Design_Field::update_dependencies( const Vector< std::shared_ptr< Field > >& aUpdatedFields )
    {
        mField->update_dependencies( aUpdatedFields );
    }

    //--------------------------------------------------------------------------------------------------------------

    std::shared_ptr< mtk::Field > Design_Field::get_mtk_field()
    {
        return mField->get_mtk_field();
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::gen
