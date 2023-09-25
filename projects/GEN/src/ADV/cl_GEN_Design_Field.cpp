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

    Design_Field::Design_Field(
            std::shared_ptr< Field > aField,
            Field_Parameters aParameters )
            : mField( std::move( aField ) )
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
            const Matrix<DDUMat>& aCoefficientIndices,
            const Matrix<DDSMat>& aSharedADVIds,
            uint                  aADVOffsetID )
    {
        if ( mParameters.mDiscretizationIndex >= 0 )
        {
            mField = std::make_shared< BSpline_Field >(
                    aMeshPair,
                    aOwnedADVs,
                    aCoefficientIndices,
                    aSharedADVIds,
                    aADVOffsetID,
                    mParameters.mDiscretizationIndex,
                    mField );
        }
        else if ( mParameters.mDiscretizationIndex == -1 )
        {
            mField = std::make_shared< Stored_Field >(
                    aMeshPair.get_interpolation_mesh(),
                    mField );
        }
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
    real Design_Field::get_field_value(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        return mField->get_field_value( aNodeIndex, aCoordinates );
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
    const Matrix< DDRMat >& Design_Field::get_dfield_dadvs(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        return mField->get_dfield_dadvs( aNodeIndex, aCoordinates );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat > Design_Field::get_determining_adv_ids(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        return mField->get_determining_adv_ids( aNodeIndex, aCoordinates );
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Design_Field::depends_on_advs()
    {
        return mField->has_advs();
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
        return mParameters.mName;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDSMat >&
    Design_Field::get_num_refinements()
    {
        return mParameters.mNumRefinements;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDSMat >&
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

    void
    Design_Field::set_num_original_nodes( uint aNumOriginalNodes )
    {
        mField->set_num_original_nodes( aNumOriginalNodes );
    }

    //--------------------------------------------------------------------------------------------------------------

}
