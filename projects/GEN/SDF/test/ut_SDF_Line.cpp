/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_SDF_Line.cpp
 *
 */

#include <catch.hpp>
#include <algorithm>
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"

#include "linalg_typedefs.hpp"
#include "op_minus.hpp"
#include "fn_norm.hpp"
#include "fn_dot.hpp"

#include "cl_Vector.hpp"
#include "cl_SDF_Facet_Vertex.hpp"
#include "cl_SDF_Line.hpp"

namespace moris::sdf
{
    TEST_CASE( "SDF::Line()", "[gen], [sdf], [line]" )
    {
        real tEpsilon = 1e-9;

        // example coordinartes for the line
        // create list of vertices
        Vector< std::shared_ptr< Facet_Vertex > > tVertices;
        tVertices.resize( 2, nullptr );
        tVertices( 0 ) = std::make_shared< Facet_Vertex >( 0, Matrix< DDRMat >( { { 2.0 }, { 2.0 } } ) );
        tVertices( 1 ) = std::make_shared< Facet_Vertex >( 1, Matrix< DDRMat >( { { -2.0 }, { -1.0 } } ) );

        // Create the line with the given vertices
        Line tLine( 0, tVertices );

        // Check that the coordinates were applied correctly
        Matrix< DDRMat > tCoordinates = tLine.get_vertex_coords();
        REQUIRE( abs( tCoordinates( 0, 0 ) - 2.0 ) < tEpsilon );
        REQUIRE( abs( tCoordinates( 0, 1 ) - 2.0 ) < tEpsilon );
        REQUIRE( abs( tCoordinates( 1, 0 ) + 2.0 ) < tEpsilon );
        REQUIRE( abs( tCoordinates( 1, 1 ) + 1.0 ) < tEpsilon );

        // Check the center vector is correct
        Matrix< DDRMat > tCenter = { { 0.0 }, { 0.5 } };
        CHECK( abs( tCenter( 0 ) - tLine.get_center()( 0 ) ) < tEpsilon );
        CHECK( abs( tCenter( 1 ) - tLine.get_center()( 1 ) ) < tEpsilon );

        Matrix< DDRMat > tNormal = { { 0.6 }, { -0.8 } };
        CHECK( abs( tNormal( 0 ) - tLine.get_normal()( 0 ) ) < tEpsilon );
        CHECK( abs( tNormal( 1 ) - tLine.get_normal()( 1 ) ) < tEpsilon );

        real tMinX = -2.0;
        real tMinY = -1.0;

        real tMaxX = 2.0;
        real tMaxY = 2.0;

        CHECK( abs( tLine.get_min_coord( 0 ) - tMinX ) < tEpsilon );
        CHECK( abs( tLine.get_min_coord( 1 ) - tMinY ) < tEpsilon );

        CHECK( abs( tLine.get_max_coord( 0 ) - tMaxX ) < tEpsilon );
        CHECK( abs( tLine.get_max_coord( 1 ) - tMaxY ) < tEpsilon );

        real tHesseActual = -0.4;
        CHECK( abs( tHesseActual - tLine.get_hesse() ) < tEpsilon );

        SECTION( "SDF Line - intersect with X axis test" )
        {
            // Check intersection for point left of the line
            Matrix< DDRMat > tPoint        = { { -1.0 }, { 1.0 } };
            real             tIntersection = 2.0 / 3.0;
            real             tComputedIntersection;
            bool             tError;
            tLine.intersect_with_coordinate_axis( tPoint, 0, tComputedIntersection, tError );
            CHECK( !tError );
            CHECK( abs( tComputedIntersection - tIntersection ) < tEpsilon );

            // reset error flag
            tError = true;

            // Check intersection for point right of the line
            tPoint        = { { 4.0 }, { 1.5 } };
            tIntersection = 4.0 / 3.0;
            tLine.intersect_with_coordinate_axis( tPoint, 0, tComputedIntersection, tError );
            CHECK( !tError );
            CHECK( abs( tComputedIntersection - tIntersection ) < tEpsilon );

            // reset error flag
            tError = true;

            // Check that a proper intersection is computed if the vertices are flipped
            std::reverse( tVertices.begin(), tVertices.end() );
            Line tLineReversed( 1, tVertices );
            tLineReversed.intersect_with_coordinate_axis( tPoint, 0, tComputedIntersection, tError );
            CHECK( !tError );
            CHECK( abs( tComputedIntersection - tIntersection ) < tEpsilon );

            // reset error flag
            tError = true;

            // Check intersection for point that directly hits a vertex
            tPoint        = { { 0.0 }, { 2.0 } };
            tIntersection = 2.0;
            tLine.intersect_with_coordinate_axis( tPoint, 0, tComputedIntersection, tError );
            CHECK( !tError );
            CHECK( abs( tComputedIntersection - tIntersection ) < tEpsilon );

            // reset error flag
            tError = true;

            // Check intersection for a point that is above the line
            // NOTE: the computed intersection lies on the infinite line created by the two points, ie. not within the Line object
            // There is no check to guarantee the intersection lies in between the two points
            tPoint        = { { 6.0 }, { 3.5 } };
            tIntersection = 4.0;
            tLine.intersect_with_coordinate_axis( tPoint, 0, tComputedIntersection, tError );
            CHECK( !tError );
            CHECK( abs( tComputedIntersection - tIntersection ) < tEpsilon );

            // Check to ensure proper error for a horizontal line
            Vector< std::shared_ptr< Facet_Vertex > > tHorizontalVertices;
            tHorizontalVertices.resize( 2, nullptr );
            tHorizontalVertices( 0 ) = std::make_shared< Facet_Vertex >( 2, Matrix< DDRMat >( { { -4.0 }, { 2.0 } } ) );
            tHorizontalVertices( 1 ) = std::make_shared< Facet_Vertex >( 3, Matrix< DDRMat >( { { 4.0 }, { 2.0 } } ) );
            Line tHorizontalLine( 2, tHorizontalVertices );

            tHorizontalLine.intersect_with_coordinate_axis( tPoint, 0, tComputedIntersection, tError );
            CHECK( tError );
            CHECK( std::isnan( tComputedIntersection ) );

            // extra check of the min/max coords for a vertical line
            real tMinX = -4.0;
            real tMinY = 2.0;

            real tMaxX = 4.0;
            real tMaxY = 2.0;

            CHECK( abs( tHorizontalLine.get_min_coord( 0 ) - tMinX ) < tEpsilon );
            CHECK( abs( tHorizontalLine.get_min_coord( 1 ) - tMinY ) < tEpsilon );

            CHECK( abs( tHorizontalLine.get_max_coord( 0 ) - tMaxX ) < tEpsilon );
            CHECK( abs( tHorizontalLine.get_max_coord( 1 ) - tMaxY ) < tEpsilon );
        }
        SECTION( "SDF Line - intersect with Y axis test" )
        {
            // Check intersection for point below the line
            Matrix< DDRMat > tPoint        = { { -1.0 }, { 1.0 } };
            real             tIntersection = -1.0 / 4.0;
            real             tComputedIntersection;
            bool             tError;
            tLine.intersect_with_coordinate_axis( tPoint, 1, tComputedIntersection, tError );
            CHECK( !tError );
            CHECK( abs( tComputedIntersection - tIntersection ) < tEpsilon );

            // reset error flag
            tError = true;

            // Check intersection for point above the line
            tPoint        = { { 4.0 }, { 1.5 } };
            tIntersection = 7.0 / 2.0;
            tLine.intersect_with_coordinate_axis( tPoint, 1, tComputedIntersection, tError );
            CHECK( !tError );
            CHECK( abs( tComputedIntersection - tIntersection ) < tEpsilon );

            // reset error flag
            tError = true;

            // Check intersection for point that directly hits a vertex
            tPoint        = { { 2.0 }, { 0.0 } };
            tIntersection = 2.0;
            tLine.intersect_with_coordinate_axis( tPoint, 1, tComputedIntersection, tError );
            CHECK( !tError );
            CHECK( abs( tComputedIntersection - tIntersection ) < tEpsilon );

            // reset error flag
            tError = false;

            // Check to ensure proper error for a vertical line
            Vector< std::shared_ptr< Facet_Vertex > > tVerticalVertices;
            tVerticalVertices.resize( 2, nullptr );
            tVerticalVertices( 0 ) = std::make_shared< Facet_Vertex >( 2, Matrix< DDRMat >( { { 2.0 }, { 2.0 } } ) );
            tVerticalVertices( 1 ) = std::make_shared< Facet_Vertex >( 3, Matrix< DDRMat >( { { 2.0 }, { 6.0 } } ) );
            Line tVerticalLine( 3, tVerticalVertices );

            tVerticalLine.intersect_with_coordinate_axis( tPoint, 1, tComputedIntersection, tError );
            CHECK( tError );
            CHECK( std::isnan( tComputedIntersection ) );

            // extra check of the min/max coords for a vertical line
            real tMinX = 2.0;
            real tMinY = 2.0;

            real tMaxX = 2.0;
            real tMaxY = 6.0;

            CHECK( abs( tVerticalLine.get_min_coord( 0 ) - tMinX ) < tEpsilon );
            CHECK( abs( tVerticalLine.get_min_coord( 1 ) - tMinY ) < tEpsilon );

            CHECK( abs( tVerticalLine.get_max_coord( 0 ) - tMaxX ) < tEpsilon );
            CHECK( abs( tVerticalLine.get_max_coord( 1 ) - tMaxY ) < tEpsilon );
        }
    }
}    // namespace moris::sdf
