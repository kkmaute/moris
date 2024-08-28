/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_SDF_Object.cpp
 *
 */

#include <catch.hpp>
#include <algorithm>
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "paths.hpp"

#include "linalg_typedefs.hpp"
#include "op_minus.hpp"
#include "fn_norm.hpp"
#include "fn_dot.hpp"
#include "fn_all_true.hpp"

#include "cl_SDF_Facet_Vertex.hpp"
#include "cl_SDF_Line.hpp"
#include "cl_SDF_Object.hpp"

namespace moris::sdf
{
    TEST_CASE( "SDF::Object()", "[gen], [sdf], [object transformation test]" )
    {
        real tEpsilon = 1e-9;

        // create a 2D object
        std::string tObjectPath = get_base_moris_dir() + "projects/GEN/SDF/test/data/rhombus.obj";
        Object      tObject( tObjectPath );

        // expected results
        Matrix< DDRMat > tRotationCoordsExpected = { { 0.3535533906, 0.3535533906 }, { -0.3535533906, 0.3535533906 } };
        Matrix< DDRMat > tResetCoordsExpected    = { { 0.5, 0.0 }, { 0.0, 0.5 } };
        Matrix< DDRMat > tShiftCoordsExpected    = { { -0.25, 0.75 }, { -0.75, 0.25 } };
        Matrix< DDRMat > tScaleCoordsExpected    = { { -1.0, 0.0 }, { 0.0, -0.25 } };

        // rotate the object by 45 degrees and get the first facet
        Matrix< DDRMat > tRotationMatrix = { { 0.7071067812, -0.7071067812 }, { 0.7071067812, 0.7071067812 } };
        tObject.rotate( tRotationMatrix );
        Matrix< DDRMat > tFacetCoords = tObject.get_facet( 0 ).get_vertex_coords();

        // check rotation
        CHECK( all_true( abs( tFacetCoords - tRotationCoordsExpected ) < tEpsilon ) );

        // reset and get first facet coordinates
        tObject.reset_coordinates();
        tFacetCoords = tObject.get_facet( 0 ).get_vertex_coords();

        // check reset
        CHECK( all_true( abs( tFacetCoords - tResetCoordsExpected ) < tEpsilon ) );

        // shift the object and get the second facet
        Vector< real > tShift = { -0.25, 0.25 };
        tObject.shift( tShift );
        tFacetCoords = tObject.get_facet( 1 ).get_vertex_coords();

        // check shift
        CHECK( all_true( abs( tFacetCoords - tShiftCoordsExpected ) < tEpsilon ) );

        // reset
        tObject.reset_coordinates();

        // scale the object and get the third facet
        Vector< real > tScale = { 2.0, 0.5 };
        tObject.scale( tScale );
        tFacetCoords = tObject.get_facet( 2 ).get_vertex_coords();

        // check scale
        CHECK( all_true( abs( tFacetCoords - tScaleCoordsExpected ) < tEpsilon ) );

        // reset
        tObject.reset_coordinates();
    }
}    // namespace moris::sdf
