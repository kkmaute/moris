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
#include "cl_SDF_Data.hpp"
#include "cl_SDF_Object.hpp"
#include "cl_SDF_Raycast.hpp"
#include "SDF_Tools.hpp"

using namespace moris;

TEST_CASE(
        "ge::sdf::Raycast",
        "[geomeng],[sdf],[Raycaster]" )
{
    // get root from environment
    std::string tMorisRoot = moris::get_base_moris_dir();

    if ( par_size() == 1 )
    {
        SECTION( "SDF: Raycast 3D Test" )
        {
            // create triangle object from object file
            std::string         tObjectPath = tMorisRoot + "/projects/GEN/GEN_MAIN/SDF/test/data/tetrahedron.obj";
            sdf::Object         tObject( tObjectPath );
            Cell< sdf::Facet* > tFacets = tObject.get_facets();

            // create raycaster
            sdf::Raycast tRaycaster( tObject );

            // define test point
            Matrix< DDRMat > tTestPoint = {
                { 0.9 },
                { 0.6 },
                { 0.7 }
            };

            tRaycaster.set_point( tTestPoint );

            // preselect in x direction and ensure they are correct
            tRaycaster.preselect_triangles_x();
            Cell< uint > tCandidatesExpected = { 0, 1, 2 };
            Cell< uint > tCandidateTriangles = tRaycaster.get_candidate_facets();

            REQUIRE( tCandidateTriangles.size() == 3 );
            CHECK( tCandidatesExpected( 0 ) == tCandidateTriangles( 0 ) );
            CHECK( tCandidatesExpected( 1 ) == tCandidateTriangles( 1 ) );
            CHECK( tCandidatesExpected( 2 ) == tCandidateTriangles( 2 ) );

            // repeat for y direction
            tRaycaster.preselect_triangles_y();
            tCandidatesExpected = { 0, 1 };
            tCandidateTriangles = tRaycaster.get_candidate_facets();

            REQUIRE( tCandidateTriangles.size() == 2 );
            CHECK( tCandidatesExpected( 0 ) == tCandidateTriangles( 0 ) );
            CHECK( tCandidatesExpected( 1 ) == tCandidateTriangles( 1 ) );

            // repeat for z direction
            tRaycaster.preselect_triangles_z();
            tCandidatesExpected = { 0, 1, 3 };
            tCandidateTriangles = tRaycaster.get_candidate_facets();

            REQUIRE( tCandidateTriangles.size() == 3 );
            CHECK( tCandidatesExpected( 0 ) == tRaycaster.get_candidate_facets()( 0 ) );
            CHECK( tCandidatesExpected( 1 ) == tRaycaster.get_candidate_facets()( 1 ) );
            CHECK( tCandidatesExpected( 2 ) == tRaycaster.get_candidate_facets()( 2 ) );

            // check for intersection with facets and ensure they are correct
            tRaycaster.intersect_triangles( 2 );
            Cell< sdf::Facet* > tIntersectedTriangles = tRaycaster.get_intersected_facets();
            Cell< sdf::Facet* > tIntersectedFacetsExpected( 2 );
            tIntersectedFacetsExpected( 0 ) = tFacets( 0 );
            tIntersectedFacetsExpected( 1 ) = tFacets( 3 );

            REQUIRE( tIntersectedTriangles.size() == 2 );
            CHECK( tIntersectedFacetsExpected( 0 ) == tIntersectedTriangles( 0 ) );
            CHECK( tIntersectedFacetsExpected( 1 ) == tIntersectedTriangles( 1 ) );

            // compute the intersection locations and ensure they are correct BRENDAN FINISH UNIT TEST
            // Cell< real > tIntersectionCoordinates = 
        }
        SECTION( "2D SDF Raycast Test" )
        {
            // create triangle object from object file
            std::string         tObjectPath = tMorisRoot + "/projects/GEN/GEN_MAIN/SDF/test/data/rhombus.obj";
            sdf::Object         tObject( tObjectPath );
            Cell< sdf::Facet* > tFacets = tObject.get_facets();

            // create raycaster
            sdf::Raycast tRaycaster( tObject );

            // define test point
            Matrix< DDRMat > tTestPoint = { { -.25 }, { -0.3 } };

            tRaycaster.set_point( tTestPoint );

            // preselect in y direction and ensure the candidates and intersected facets are marked
            tRaycaster.preselect_lines( 1 );
            uint                tIntersectedLinesExpected = 1;
            sdf::Facet*         tCandidateLinesExpected   = tFacets( 2 );
            Cell< uint >        tIntersectedLines         = tRaycaster.get_candidate_facets();
            Cell< sdf::Facet* > tCandidateLines           = tRaycaster.get_intersected_facets();

            REQUIRE( tIntersectedLines.size() == 1 );
            REQUIRE( tCandidateLines.size() == 1 );
            CHECK( tIntersectedLines( 0 ) == tIntersectedLinesExpected );
            CHECK( tCandidateLines( 0 ) == tCandidateLinesExpected );

            // preselect in x direction and ensure the candidates and intersected facets are marked
            tRaycaster.preselect_lines( 0 );
            tIntersectedLinesExpected = 3;
            tIntersectedLines         = tRaycaster.get_candidate_facets();
            tCandidateLines           = tRaycaster.get_intersected_facets();

            REQUIRE( tIntersectedLines.size() == 1 );
            REQUIRE( tCandidateLines.size() == 1 );
            CHECK( tIntersectedLines( 0 ) == tIntersectedLinesExpected );
            CHECK( tCandidateLines( 0 ) == tCandidateLinesExpected );

            // intersect the candidate facets and determine the intersection location
            tRaycaster.intersect_ray_with_facets( 0 );
            real         tIntersectionCoordinateExpected = -0.2;
            Cell< real > tIntersectionCoordinates        = tRaycaster.get_intersection_coordinates();

            REQUIRE( tIntersectionCoordinates.size() == 1 );
            CHECK( tIntersectionCoordinates( 0 ) - tIntersectionCoordinateExpected < sdf::gSDFepsilon);

            // determine if the point is inside/outside and check expectation
            tRaycaster.check_if_node_is_inside_lines( 0 );
            uint tPointInsideExpected = 0;
            uint tPointInside         = tRaycaster.is_point_inside();

            CHECK( tPointInside == tPointInsideExpected );

            // ensure the test gives the same result when the entire algorithm is called at once
            tRaycaster.calculate_raycast( tTestPoint );
            tPointInside = tRaycaster.is_point_inside();

            CHECK( tPointInside == tPointInsideExpected );

            // repeat for a point inside the surface
            tTestPoint = { { -.25 }, { 0.2 } };
            tRaycaster.set_point( tTestPoint );

            tRaycaster.preselect_lines( 1 );
            tCandidateLinesExpected = tFacets( 1 );
            tCandidateLines         = tRaycaster.get_intersected_facets();
            tIntersectedLines       = tRaycaster.get_candidate_facets();

            REQUIRE( tCandidateLines.size() == 1 );
            REQUIRE( tIntersectedLines.size() == 0 );
            CHECK( tCandidateLines( 0 ) == tCandidateLinesExpected );

            tRaycaster.intersect_ray_with_facets( 1 );
            tIntersectionCoordinateExpected = 0.25;
            tIntersectionCoordinates        = tRaycaster.get_intersection_coordinates();

            REQUIRE( tIntersectionCoordinates.size() == 1 );
            CHECK( tIntersectionCoordinates( 0 ) - tIntersectionCoordinateExpected  < sdf::gSDFepsilon );

            tRaycaster.check_if_node_is_inside_lines( 1 );
            tPointInsideExpected = 1;
            tPointInside = tRaycaster.is_point_inside();

            CHECK( tPointInside == tPointInsideExpected );

            tRaycaster.calculate_raycast( tTestPoint );

            CHECK( tPointInside == tPointInsideExpected );
        }
    }
}