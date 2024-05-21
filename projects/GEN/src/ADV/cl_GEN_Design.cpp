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

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    Design_Parameters::Design_Parameters( const Parameter_List& aParameterList )
            : mNumberOfRefinements( aParameterList.get< Vector< uint > >( "number_of_refinements" ) )
            , mRefinementMeshIndices( aParameterList.get< Vector< uint > >( "refinement_mesh_index" ) )
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

    const Vector< uint >& Design::get_num_refinements()
    {
        return mParameters.mNumberOfRefinements;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Vector< uint >&
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

}    // namespace moris::gen
