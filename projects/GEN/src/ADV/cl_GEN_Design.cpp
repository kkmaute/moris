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

    //--------------------------------------------------------------------------------------------------------------

    Design_Parameters::Design_Parameters( const ParameterList& aParameterList )
            : mNumberOfRefinements( aParameterList.get_cell< uint >( "number_of_refinements" ) )
            , mRefinementMeshIndices( aParameterList.get_cell< uint >( "refinement_mesh_index" ) )
            , mRefinementFunctionIndex( aParameterList.get< sint >( "refinement_function_index" ) )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Design_Parameters::Design_Parameters()
            : mNumberOfRefinements( {} )
            , mRefinementMeshIndices( {} )
            , mRefinementFunctionIndex( -1 )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Design::Design( Design_Parameters aParameters )
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
    Design::get_refinement_mesh_indices()
    {
        return mParameters.mRefinementMeshIndices;
    }

    //--------------------------------------------------------------------------------------------------------------

    sint
    Design::get_refinement_function_index()
    {
        return mParameters.mRefinementFunctionIndex;
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::ge
