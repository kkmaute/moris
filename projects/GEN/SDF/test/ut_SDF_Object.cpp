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
        // Matrix< DDRMat > tRotationCoordsExpected = { { 0.3535533906, 0.3535533906 }, { -0.3535533906, 0.3535533906 } };
        Matrix< DDRMat > tResetCoordsExpected = { { 0.5, 0.0 }, { 0.0, 0.5 } };
        Matrix< DDRMat > tShiftCoordsExpected = { { -0.25, -0.5 }, { 0.75, 0.0 } };
        // Matrix< DDRMat > tScaleCoordsExpected    = { { -1.0, 0.0 }, { 0.0, -0.25 } };

        // rotate the object by 45 degrees and get the first facet
        // Matrix< DDRMat > tRotationMatrix = { { 0.7071067812, -0.7071067812 }, { 0.7071067812, 0.7071067812 } };
        // tObject.set_rotation( tRotationMatrix );
        // Matrix< DDRMat > tFacetCoords = tObject.get_all_vertex_coordinates_of_facet( 0 );

        // check rotation
        // CHECK( all_true( abs( tFacetCoords - tRotationCoordsExpected ) < tEpsilon ) );

        // reset and get first facet coordinates
        tObject.reset_coordinates();
        Matrix< DDRMat > tFacetCoords = tObject.get_all_vertex_coordinates_of_facet( 0 );

        // check reset
        CHECK( all_true( abs( tFacetCoords - tResetCoordsExpected ) < tEpsilon ) );

        // shift the object and get the second facet
        Matrix< DDRMat > tShift = { { -0.25 }, { 0.25 } };
        tObject.set_vertex_displacement( 1, tShift );
        tFacetCoords = tObject.get_all_vertex_coordinates_of_facet( 1 );

        // check shift
        CHECK( all_true( abs( tFacetCoords - tShiftCoordsExpected ) < tEpsilon ) );

        // reset and get first facet coordinates
        tObject.reset_coordinates();
        tFacetCoords = tObject.get_all_vertex_coordinates_of_facet( 0 );

        // check reset
        CHECK( all_true( abs( tFacetCoords - tResetCoordsExpected ) < tEpsilon ) );

    }
}    // namespace moris::sdf
