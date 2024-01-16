/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Design.cpp
 *
 */

#include "cl_GEN_Design.hpp"

namespace moris::ge
{
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

    Design::Design( Field_Parameters aParameters )
            : mParameters( std::move( aParameters ) )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    const Cell< uint >& Design::get_num_refinements()
    {
        return mParameters.mNumberOfRefinements;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Cell< uint >&
    Design::get_refinement_mesh_indices(){
        return mParameters.mRefinementMeshIndices;
    }

    //--------------------------------------------------------------------------------------------------------------

    sint
    Design::get_refinement_function_index()
    {
        return mParameters.mRefinementFunctionIndex;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Design::intended_discretization()
    {
        return ( mParameters.mDiscretizationIndex >= 0 );
    }

    //--------------------------------------------------------------------------------------------------------------
    
    moris_index
    Design::get_discretization_mesh_index() const
    {
        MORIS_ASSERT( mParameters.mDiscretizationIndex >= 0,
                "A discretization is not intended for this field. Check this with intended_discretization() first." );

        return mParameters.mDiscretizationIndex;
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Design::get_discretization_lower_bound()
    {
        return mParameters.mDiscretizationLowerBound;
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Design::get_discretization_upper_bound()
    {
        return mParameters.mDiscretizationUpperBound;
    }
}    // namespace moris::ge
