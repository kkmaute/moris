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

            // compute the intersection locations and ensure they are correct
            tRaycaster.intersect_ray_with_facets( 2 );
            Cell< real > tIntersectionCoordinatesExpected = { 0.4718, 0.9024 };
            Cell< real > tIntersectionCoordinates         = tRaycaster.get_intersection_coordinates();
            REQUIRE( tIntersectionCoordinates.size() == 2 );
            CHECK( tIntersectionCoordinates( 0 ) - tIntersectionCoordinatesExpected( 0 ) < sdf::gSDFepsilon );
            CHECK( tIntersectionCoordinates( 1 ) - tIntersectionCoordinatesExpected( 1 ) < sdf::gSDFepsilon );

            // check if the point is insde and compare it to expectations
            tRaycaster.check_if_node_is_inside_triangles( 2 );
            uint tPointIsInsideExpected = 1;

            CHECK( tRaycaster.is_point_inside() == tPointIsInsideExpected );

            // check with the full raycast algorithm and see if they match
            tRaycaster.raycast_point( tTestPoint );

            CHECK( tRaycaster.is_point_inside() == tPointIsInsideExpected );

            // repeat test for point that is outside
            tTestPoint = {
                { 0.2 },
                { 0.6 },
                { 0.7 }
            };
            tRaycaster.set_point( tTestPoint );

            tRaycaster.preselect_triangles_z();

            REQUIRE( tRaycaster.get_candidate_facets().size() == 0 );

            tRaycaster.preselect_triangles_x();
            tCandidatesExpected = { 0, 1, 2 };
            tCandidateTriangles = tRaycaster.get_candidate_facets();

            REQUIRE( tRaycaster.get_candidate_facets().size() == 3 );
            CHECK( tCandidateTriangles( 0 ) == tCandidatesExpected( 0 ) );
            CHECK( tCandidateTriangles( 1 ) == tCandidatesExpected( 1 ) );
            CHECK( tCandidateTriangles( 2 ) == tCandidatesExpected( 2 ) );

            tRaycaster.intersect_triangles( 0 );
            tIntersectedTriangles           = tRaycaster.get_intersected_facets();
            tIntersectedFacetsExpected( 0 ) = tFacets( 1 );
            tIntersectedFacetsExpected( 1 ) = tFacets( 2 );

            REQUIRE( tIntersectedTriangles.size() == 2 );
            CHECK( tIntersectedTriangles( 0 ) == tIntersectedFacetsExpected( 0 ) );
            CHECK( tIntersectedTriangles( 1 ) == tIntersectedFacetsExpected( 1 ) );

            // although two facets are intersected, one of them should produce an error and be removed
            tRaycaster.intersect_ray_with_facets( 0 );
            tIntersectionCoordinatesExpected = { 0.715517872727636, 1.284482127272365 };
            tIntersectionCoordinates = tRaycaster.get_intersection_coordinates();

            REQUIRE( tIntersectionCoordinates.size() == 2 );
            CHECK( tIntersectionCoordinates( 0 ) - tIntersectionCoordinatesExpected( 0 ) < sdf::gSDFepsilon );
            CHECK( tIntersectionCoordinates( 1 ) - tIntersectionCoordinatesExpected( 1 ) < sdf::gSDFepsilon );


            tRaycaster.check_if_node_is_inside_triangles( 0 );
            tPointIsInsideExpected = 0;

            CHECK( tRaycaster.is_point_inside() == tPointIsInsideExpected );

            tRaycaster.raycast_point( tTestPoint );
            
            CHECK( tRaycaster.is_point_inside() == tPointIsInsideExpected );
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
            CHECK( tIntersectionCoordinates( 0 ) - tIntersectionCoordinateExpected < sdf::gSDFepsilon );

            // determine if the point is inside/outside and check expectation
            tRaycaster.check_if_node_is_inside_lines( 0 );
            uint tPointInsideExpected = 0;
            uint tPointInside         = tRaycaster.is_point_inside();

            CHECK( tPointInside == tPointInsideExpected );

            // ensure the test gives the same result when the entire algorithm is called at once
            tRaycaster.raycast_point( tTestPoint );
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
            CHECK( tIntersectionCoordinates( 0 ) - tIntersectionCoordinateExpected < sdf::gSDFepsilon );

            tRaycaster.check_if_node_is_inside_lines( 1 );
            tPointInsideExpected = 1;
            tPointInside         = tRaycaster.is_point_inside();

            CHECK( tPointInside == tPointInsideExpected );

            tRaycaster.raycast_point( tTestPoint );

            CHECK( tPointInside == tPointInsideExpected );
        }
    }
}