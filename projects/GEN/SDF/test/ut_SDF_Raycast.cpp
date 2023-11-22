/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_SDF_Core.cpp
 *
 */

#include <string>
#include <catch.hpp>

// core
#include "typedefs.hpp"
#include "paths.hpp"

// comm
#include "cl_Communication_Manager.hpp"    // COM/src
#include "cl_Communication_Tools.hpp"      // COM/src

// linalg
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_minus.hpp"
#include "op_equal_equal.hpp"
#include "fn_norm.hpp"
#include "fn_all_true.hpp"

#include "cl_MTK_Mesh_Factory.hpp"

// SDF
#include "cl_SDF_Object.hpp"
#include "cl_SDF_Raycast.hpp"
#include "fn_SDF_Raycast.hpp"
#include "SDF_Tools.hpp"

namespace moris::sdf
{

    TEST_CASE(
            "ge::sdf::Raycast",
            "[geomeng],[sdf],[Raycaster]" )
    {
        // get root from environment
        std::string tMorisRoot = moris::get_base_moris_dir();

        if ( par_size() == 1 )
        {
            SECTION( "SDF: Raycast Free Function 3D Test" )
            {
                // create triangle object from object file
                std::string    tObjectPath = tMorisRoot + "projects/GEN/SDF/test/data/tetrahedron.obj";
                Object         tObject( tObjectPath );
                Cell< Facet* > tFacets = tObject.get_facets();

                // define test point
                Matrix< DDRMat > tTestPoint = {
                    { 0.9 },
                    { 0.6 },
                    { 0.7 }
                };

                // preselect in x direction and ensure they are correct
                Cell< uint > tCandidatesExpected = { 0, 1, 2 };
                Cell< uint > tCandidateTriangles = preselect_triangles( tObject, tTestPoint, 0 );

                REQUIRE( tCandidateTriangles.size() == 3 );
                CHECK( tCandidatesExpected( 0 ) == tCandidateTriangles( 0 ) );
                CHECK( tCandidatesExpected( 1 ) == tCandidateTriangles( 1 ) );
                CHECK( tCandidatesExpected( 2 ) == tCandidateTriangles( 2 ) );

                // repeat for y direction
                tCandidatesExpected = { 0, 1 };
                tCandidateTriangles = preselect_triangles( tObject, tTestPoint, 1 );

                REQUIRE( tCandidateTriangles.size() == 2 );
                CHECK( tCandidatesExpected( 0 ) == tCandidateTriangles( 0 ) );
                CHECK( tCandidatesExpected( 1 ) == tCandidateTriangles( 1 ) );

                // repeat for z direction
                tCandidatesExpected = { 0, 1, 3 };
                tCandidateTriangles = preselect_triangles( tObject, tTestPoint, 2 );

                REQUIRE( tCandidateTriangles.size() == 3 );
                CHECK( tCandidatesExpected( 0 ) == tCandidateTriangles( 0 ) );
                CHECK( tCandidatesExpected( 1 ) == tCandidateTriangles( 1 ) );
                CHECK( tCandidatesExpected( 2 ) == tCandidateTriangles( 2 ) );

                // check for intersection with facets and ensure they are correct
                Cell< Facet* > tIntersectedTriangles = intersect_triangles( tCandidateTriangles, tObject, tTestPoint, 2 );
                Cell< Facet* > tIntersectedFacetsExpected( 2 );
                tIntersectedFacetsExpected( 0 ) = tFacets( 0 );
                tIntersectedFacetsExpected( 1 ) = tFacets( 3 );

                REQUIRE( tIntersectedTriangles.size() == 2 );
                CHECK( tIntersectedFacetsExpected( 0 ) == tIntersectedTriangles( 0 ) );
                CHECK( tIntersectedFacetsExpected( 1 ) == tIntersectedTriangles( 1 ) );

                // compute the intersection locations and ensure they are correct
                Cell< real > tIntersectionCoordinatesExpected = { 0.4718, 0.9024 };
                Cell< real > tIntersectionCoordinates         = intersect_ray_with_facets( tIntersectedTriangles, tTestPoint, 2 );
                REQUIRE( tIntersectionCoordinates.size() == 2 );
                CHECK( tIntersectionCoordinates( 0 ) - tIntersectionCoordinatesExpected( 0 ) < gSDFepsilon );
                CHECK( tIntersectionCoordinates( 1 ) - tIntersectionCoordinatesExpected( 1 ) < gSDFepsilon );

                // check if the point is inside and compare it to expectations
                Object_Region tPointIsInside         = check_if_node_is_inside_triangles( tIntersectionCoordinates, tTestPoint, 2 );
                Object_Region tPointIsInsideExpected = INSIDE;

                CHECK( tPointIsInside == tPointIsInsideExpected );

                // check with the full raycast algorithm and see if they match
                tPointIsInside = raycast_point( tObject, tTestPoint );

                CHECK( tPointIsInside == tPointIsInsideExpected );

                // repeat test for point that is outside
                tTestPoint = {
                    { 0.2 },
                    { 0.6 },
                    { 0.7 }
                };

                tCandidateTriangles = preselect_triangles( tObject, tTestPoint, 2 );

                REQUIRE( tCandidateTriangles.size() == 0 );

                tCandidatesExpected = { 0, 1, 2 };
                tCandidateTriangles = preselect_triangles( tObject, tTestPoint, 0 );


                REQUIRE( tCandidateTriangles.size() == 3 );
                CHECK( tCandidateTriangles( 0 ) == tCandidatesExpected( 0 ) );
                CHECK( tCandidateTriangles( 1 ) == tCandidatesExpected( 1 ) );
                CHECK( tCandidateTriangles( 2 ) == tCandidatesExpected( 2 ) );

                tIntersectedTriangles           = intersect_triangles( tCandidateTriangles, tObject, tTestPoint, 0 );
                tIntersectedFacetsExpected( 0 ) = tFacets( 1 );
                tIntersectedFacetsExpected( 1 ) = tFacets( 2 );

                REQUIRE( tIntersectedTriangles.size() == 2 );
                CHECK( tIntersectedTriangles( 0 ) == tIntersectedFacetsExpected( 0 ) );
                CHECK( tIntersectedTriangles( 1 ) == tIntersectedFacetsExpected( 1 ) );

                // although two facets are intersected, one of them should produce an error and be removed
                tIntersectionCoordinatesExpected = { 0.715517872727636, 1.284482127272365 };
                tIntersectionCoordinates         = intersect_ray_with_facets( tIntersectedTriangles, tTestPoint, 0 );


                REQUIRE( tIntersectionCoordinates.size() == 2 );
                CHECK( tIntersectionCoordinates( 0 ) - tIntersectionCoordinatesExpected( 0 ) < gSDFepsilon );
                CHECK( tIntersectionCoordinates( 1 ) - tIntersectionCoordinatesExpected( 1 ) < gSDFepsilon );


                tPointIsInside = check_if_node_is_inside_triangles( tIntersectionCoordinates, tTestPoint, 0 );
                tPointIsInsideExpected = OUTSIDE;

                CHECK( tPointIsInside == tPointIsInsideExpected );

                tPointIsInside = raycast_point( tObject, tTestPoint );

                CHECK( tPointIsInside == tPointIsInsideExpected );
            }
            SECTION( "SDF: Raycast Free Function 2D Test" )
            {
                // create triangle object from object file
                std::string    tObjectPath = tMorisRoot + "projects/GEN/SDF/test/data/rhombus.obj";
                Object         tObject( tObjectPath );
                Cell< Facet* > tFacets = tObject.get_facets();

                // define test point
                Matrix< DDRMat > tTestPoint = { { -.25 }, { -0.3 } };

                // preselect in y direction and ensure the candidates and intersected facets are marked
                Cell< uint >   tIntersectedLines;
                Cell< Facet* > tCandidateLines;
                preselect_lines( tObject, tTestPoint, 1, tIntersectedLines, tCandidateLines );
                uint   tIntersectedLinesExpected = 1;
                Facet* tCandidateLinesExpected   = tFacets( 2 );

                REQUIRE( tIntersectedLines.size() == 1 );
                REQUIRE( tCandidateLines.size() == 1 );
                CHECK( tIntersectedLines( 0 ) == tIntersectedLinesExpected );
                CHECK( tCandidateLines( 0 ) == tCandidateLinesExpected );

                // preselect in x direction and ensure the candidates and intersected facets are marked
                preselect_lines( tObject, tTestPoint, 0, tIntersectedLines, tCandidateLines );
                tIntersectedLinesExpected = 3;

                REQUIRE( tIntersectedLines.size() == 1 );
                REQUIRE( tCandidateLines.size() == 1 );
                CHECK( tIntersectedLines( 0 ) == tIntersectedLinesExpected );
                CHECK( tCandidateLines( 0 ) == tCandidateLinesExpected );

                // intersect the candidate facets and determine the intersection location
                Cell< real > tIntersectionCoordinates        = intersect_ray_with_facets( tCandidateLines, tTestPoint, 0 );
                real         tIntersectionCoordinateExpected = -0.2;

                REQUIRE( tIntersectionCoordinates.size() == 1 );
                CHECK( tIntersectionCoordinates( 0 ) - tIntersectionCoordinateExpected < gSDFepsilon );

                // determine if the point is inside/outside and check expectation
                Object_Region tRegion              = check_if_node_is_inside_lines( tIntersectionCoordinates, tIntersectedLines, tTestPoint, 0 );
                Object_Region tPointInsideExpected = OUTSIDE;

                CHECK( tRegion == tPointInsideExpected );

                // ensure the test gives the same result when the entire algorithm is called at once
                tRegion = raycast_point( tObject, tTestPoint );

                CHECK( tRegion == tPointInsideExpected );

                // repeat for a point inside the surface
                tTestPoint = { { -.25 }, { 0.2 } };

                preselect_lines( tObject, tTestPoint, 1, tIntersectedLines, tCandidateLines );
                tCandidateLinesExpected = tFacets( 1 );

                REQUIRE( tCandidateLines.size() == 1 );
                REQUIRE( tIntersectedLines.size() == 0 );
                CHECK( tCandidateLines( 0 ) == tCandidateLinesExpected );

                tIntersectionCoordinates        = intersect_ray_with_facets( tCandidateLines, tTestPoint, 1 );
                tIntersectionCoordinateExpected = 0.25;

                REQUIRE( tIntersectionCoordinates.size() == 1 );
                CHECK( tIntersectionCoordinates( 0 ) - tIntersectionCoordinateExpected < gSDFepsilon );

                tRegion              = check_if_node_is_inside_lines( tIntersectionCoordinates, tIntersectedLines, tTestPoint, 1 );
                tPointInsideExpected = INSIDE;

                CHECK( tRegion == tPointInsideExpected );

                tRegion = raycast_point( tObject, tTestPoint );

                CHECK( tRegion == tPointInsideExpected );
            }
        }
    }
};    // namespace moris::sdf