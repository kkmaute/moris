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

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Field_Parameters::Field_Parameters( const ParameterList& aParameterList )
            : mDiscretizationIndex( aParameterList.get< sint >( "discretization_mesh_index" ) )
            , mDiscretizationLowerBound( aParameterList.get< real >( "discretization_lower_bound" ) )
            , mDiscretizationUpperBound( aParameterList.get< real >( "discretization_upper_bound" ) )
            , mUseMultilinearInterpolation( aParameterList.get< bool >( "use_multilinear_interpolation" ) )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Design_Field::Design_Field(
            std::shared_ptr< Field > aField,
            Field_Parameters         aParameters,
            Node_Manager&            aNodeManager )
            : mField( std::move( aField ) )
            , mNodeManager( &aNodeManager )
            , mParameters( std::move( aParameters ) )
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
            mtk::Mesh_Pair          aMeshPair,
            sol::Dist_Vector*       aOwnedADVs,
            const Matrix< DDSMat >& aSharedADVIds,
            uint                    aADVOffsetID )
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
            std::shared_ptr< mtk::Field > aMTKField,
            mtk::Mesh_Pair                aMeshPair,
            sol::Dist_Vector*             aOwnedADVs,
            const Matrix< DDSMat >&       aSharedADVIds,
            uint                          aADVOffsetID )
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
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        if ( mNodeManager->is_background_node( aNodeIndex ) )
        {
            // Get background node field value
            return mField->get_field_value( aNodeIndex, aCoordinates );
        }
        else
        {
            if ( mParameters.mUseMultilinearInterpolation )
            {
                // If we use multilinear interpolation, it is needed for all derived nodes
                return mField->get_interpolated_field_value(
                        mNodeManager->get_derived_node( aNodeIndex )->get_background_nodes(),
                        *mNodeManager );
            }
            else
            {
                // Let field decide the value
                return mField->get_field_value(
                        mNodeManager->get_derived_node( aNodeIndex ),
                        *mNodeManager );
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

            // Check for multilinear
            if ( mParameters.mUseMultilinearInterpolation )
            {
                // If we use multilinear interpolation, it is needed for all derived nodes
                mField->append_interpolated_dfield_dadvs(
                        mInterpolatedSensitivities,
                        mNodeManager->get_derived_node( aNodeIndex )->get_background_nodes(),
                        *mNodeManager );
            }
            else
            {
                // Let field decide the sensitivities
                mField->get_dfield_dadvs(
                        mInterpolatedSensitivities,
                        mNodeManager->get_derived_node( aNodeIndex ),
                        *mNodeManager );
            }

            // Return result
            return mInterpolatedSensitivities;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat > Design_Field::get_determining_adv_ids(
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
            mInterpolatedADVIDs.set_size( 0, 0 );

            // Check for multilinear
            if ( mParameters.mUseMultilinearInterpolation )
            {
                // If we use multilinear interpolation, it is needed for all derived nodes
                mField->append_interpolated_determining_adv_ids(
                        mInterpolatedADVIDs,
                        mNodeManager->get_derived_node( aNodeIndex )->get_background_nodes(),
                        *mNodeManager );
            }
            else
            {
                // Let field provide the IDs
                mField->get_determining_adv_ids(
                        mInterpolatedADVIDs,
                        mNodeManager->get_derived_node( aNodeIndex ),
                        *mNodeManager );
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

    bool Design_Field::use_multilinear_interpolation()
    {
        return mParameters.mUseMultilinearInterpolation;
    }

    //--------------------------------------------------------------------------------------------------------------

    std::shared_ptr< mtk::Field > Design_Field::get_mtk_field()
    {
        return mField->get_mtk_field();
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::ge
