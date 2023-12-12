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
            : mNumberOfRefinements( aParameterList.get_cell< uint >( "number_of_refinements" ) )
            , mRefinementMeshIndices( aParameterList.get_cell< uint >( "refinement_mesh_index" ) )
            , mRefinementFunctionIndex( aParameterList.get< sint >( "refinement_function_index" ) )
            , mDiscretizationIndex( aParameterList.get< sint >( "discretization_mesh_index" ) )
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
        MORIS_ERROR( mField, "A design must be provided a field for computing values." );
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Design_Field::intended_discretization()
    {
        return ( mParameters.mDiscretizationIndex >= 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Design_Field::discretize(
            mtk::Mesh_Pair        aMeshPair,
            sol::Dist_Vector*     aOwnedADVs,
            const Matrix<DDSMat>& aSharedADVIds,
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
            const Matrix<DDSMat>&         aSharedADVIds,
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
        if ( mNodeManager->is_base_node( aNodeIndex ) )
        {
            // Get base node field value
            return mField->get_field_value( aNodeIndex, aCoordinates );
        }
        else
        {
            if ( mParameters.mUseMultilinearInterpolation )
            {
                // If we use multilinear interpolation, it is needed for all derived nodes
                return mField->get_interpolated_field_value( mNodeManager->get_derived_node( aNodeIndex )->get_background_nodes() );
            }
            else
            {
                // Let field decide the value
                return mField->get_field_value( mNodeManager->get_derived_node( aNodeIndex ) );
            }
        }
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
    const Matrix< DDRMat >& Design_Field::get_dfield_dadvs(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        if ( mNodeManager->is_base_node( aNodeIndex ) )
        {
            // Get base node sensitivities
            return mField->get_dfield_dadvs( aNodeIndex, aCoordinates );
        }
        else
        {
            if ( mParameters.mUseMultilinearInterpolation )
            {
                // If we use multilinear interpolation, it is needed for all derived nodes
                return mField->get_interpolated_dfield_dadvs( mNodeManager->get_derived_node( aNodeIndex )->get_background_nodes() );
            }
            else
            {
                // Let field decide the sensitivities
                return mField->get_dfield_dadvs( mNodeManager->get_derived_node( aNodeIndex ) );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat > Design_Field::get_determining_adv_ids(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        if ( mNodeManager->is_base_node( aNodeIndex ) )
        {
            // Get base node determining ADV IDs
            return mField->get_determining_adv_ids( aNodeIndex, aCoordinates );
        }
        else
        {
            if ( mParameters.mUseMultilinearInterpolation )
            {
                // If we use multilinear interpolation, it is needed for all derived nodes
                return mField->get_interpolated_determining_adv_ids( mNodeManager->get_derived_node( aNodeIndex )->get_background_nodes() );
            }
            else
            {
                // Let field provide the IDs
                return mField->get_determining_adv_ids( mNodeManager->get_derived_node( aNodeIndex ) );
            }
        }
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
    void Design_Field::get_dfield_dcoordinates(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates,
            Matrix< DDRMat >&       aSensitivities )
    {
        mField->get_dfield_dcoordinates( aNodeIndex, aCoordinates, aSensitivities );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Design_Field::import_advs( sol::Dist_Vector* aOwnedADVs )
    {
        mField->import_advs( aOwnedADVs );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Design_Field::add_child_node( uint aNodeIndex, std::shared_ptr< Child_Node > aChildNode )
    {
        mField->add_child_node( aNodeIndex, aChildNode );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Design_Field::add_nodal_data( mtk::Interpolation_Mesh* aMesh )
    {
        mField->add_nodal_data( aMesh );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Design_Field::reset_nodal_data()
    {
        mField->reset_nodal_data();
    }

    //--------------------------------------------------------------------------------------------------------------

    std::string Design_Field::get_name()
    {
        return mField->get_name();
    }

    //--------------------------------------------------------------------------------------------------------------

    const Cell< uint >&
    Design_Field::get_num_refinements()
    {
        return mParameters.mNumberOfRefinements;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Cell< uint >&
    Design_Field::get_refinement_mesh_indices()
    {
        return mParameters.mRefinementMeshIndices;
    }

    //--------------------------------------------------------------------------------------------------------------

    sint
    Design_Field::get_refinement_function_index()
    {
        return mParameters.mRefinementFunctionIndex;
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_index
    Design_Field::get_discretization_mesh_index() const
    {
        MORIS_ASSERT( mParameters.mDiscretizationIndex >= 0,
                "A discretization is not intended for this field. Check this with intended_discretization() first." );

        return mParameters.mDiscretizationIndex;
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Design_Field::get_discretization_lower_bound()
    {
        return mParameters.mDiscretizationLowerBound;
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Design_Field::get_discretization_upper_bound()
    {
        return mParameters.mDiscretizationUpperBound;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Design_Field::use_multilinear_interpolation()
    {
        return mParameters.mUseMultilinearInterpolation;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Design_Field::set_num_original_nodes( uint aNumOriginalNodes )
    {
        mField->set_num_original_nodes( aNumOriginalNodes );
    }

    //--------------------------------------------------------------------------------------------------------------

    std::shared_ptr< mtk::Field > Design_Field::get_mtk_field()
    {
        return mField->get_mtk_field();
    }

    //--------------------------------------------------------------------------------------------------------------

}
