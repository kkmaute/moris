/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_GEN_Performer.cpp
 *
 */

#include "cl_WRK_GEN_Performer.hpp"

#include <utility>
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_Matrix.hpp"

namespace moris::wrk
{

    // ----------------------------------------------------------------------------

    Gen_Performer::Gen_Performer( std::shared_ptr< moris::gen::Geometry_Engine > aGeometryEngine )
            : mGeometryEngine( std::move( aGeometryEngine ) )
    {
    }

    // ----------------------------------------------------------------------------

    uint
    Gen_Performer::get_num_refinement_fields()
    {
        return mGeometryEngine->get_num_refinement_fields();
    }

    // ----------------------------------------------------------------------------

    gen::Geometric_Region Gen_Performer::get_geometric_region(
            uint                    aFieldIndex,
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        return mGeometryEngine->get_geometric_region( aFieldIndex, aNodeIndex, aCoordinates );
    }

    // ----------------------------------------------------------------------------

    const Vector< uint >&
    Gen_Performer::get_num_refinements( uint aFieldIndex )
    {
        return mGeometryEngine->get_num_refinements( aFieldIndex );
    }

    // ----------------------------------------------------------------------------

    const Vector< uint >&
    Gen_Performer::get_refinement_mesh_indices( uint aFieldIndex )
    {
        return mGeometryEngine->get_refinement_mesh_indices( aFieldIndex );
    }

    // ----------------------------------------------------------------------------

    sint
    Gen_Performer::get_refinement_function_index(
            uint aFieldIndex,
            uint aRefinementIndex )
    {
        return mGeometryEngine->get_refinement_function_index( aFieldIndex, aRefinementIndex );
    }

    // ----------------------------------------------------------------------------

}
