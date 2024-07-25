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
#include "moris_typedefs.hpp"
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
#include "fn_SDF_Raycast.hpp"
#include "SDF_Tools.hpp"

namespace moris::sdf
{

    TEST_CASE(
            "gen::sdf::Raycast",
            "[geomeng],[sdf],[Raycaster]" )
    {
        // get root from environment
        std::string tMorisRoot = moris::get_base_moris_dir();

        if ( par_size() == 1 )
        {
            SECTION( "SDF: Raycast Free Function Test - 3D" )
            {
                // create triangle object from object file
                std::string tObjectPath = tMorisRoot + "projects/GEN/SDF/test/data/tetrahedron.obj";
                Object      tObject( tObjectPath );

                // define test point that is inside the object
                Matrix< DDRMat > tTestPoint = { { 0.9, 0.6, 0.7 } };

                // preselect in x direction and ensure they are correct
                Vector< uint >      tCandidatesExpected = { 0, 1, 2 };
                Vector< uint >      tCandidateTriangles;
                Preselection_Result tPreselection = preselect_triangles( tObject, tTestPoint, 0, tCandidateTriangles );

                REQUIRE( tCandidateTriangles.size() == 3 );
                REQUIRE( tPreselection == Preselection_Result::SUCCESS );
                CHECK( tCandidatesExpected( 0 ) == tCandidateTriangles( 0 ) );
                CHECK( tCandidatesExpected( 1 ) == tCandidateTriangles( 1 ) );
                CHECK( tCandidatesExpected( 2 ) == tCandidateTriangles( 2 ) );

                // repeat for y direction
                tCandidatesExpected = { 0, 1 };
                tPreselection       = preselect_triangles( tObject, tTestPoint, 1, tCandidateTriangles );

                REQUIRE( tCandidateTriangles.size() == 2 );
                REQUIRE( tPreselection == Preselection_Result::SUCCESS );
                CHECK( tCandidatesExpected( 0 ) == tCandidateTriangles( 0 ) );
                CHECK( tCandidatesExpected( 1 ) == tCandidateTriangles( 1 ) );

                // repeat for z direction
                tCandidatesExpected = { 0, 1, 3 };
                tPreselection       = preselect_triangles( tObject, tTestPoint, 2, tCandidateTriangles );

                // Check the preselection results
                REQUIRE( tCandidateTriangles.size() == 3 );
                REQUIRE( tPreselection == Preselection_Result::SUCCESS );
                CHECK( tCandidatesExpected( 0 ) == tCandidateTriangles( 0 ) );
                CHECK( tCandidatesExpected( 1 ) == tCandidateTriangles( 1 ) );
                CHECK( tCandidatesExpected( 2 ) == tCandidateTriangles( 2 ) );

                // check moller trumbore algorithm for proper intersection location and parent facet determination
                Vector< real > tIntersectionCoordinatesExpected = { 0.408248, 0.88055613061054816 };
                Facet&         tFirstIntersectedFacetExpected   = tObject.get_facet( 3 );
                Facet&         tSecondIntersectedFacetExpected  = tObject.get_facet( 0 );

                Vector< Facet* > tIntersectedFacets;
                Vector< real >   tIntersections = intersect_triangles_moller_trumbore( tCandidateTriangles, tObject, tTestPoint, 2, tIntersectedFacets );
                REQUIRE( tIntersectedFacets.size() == 2 );
                CHECK( &tFirstIntersectedFacetExpected == tIntersectedFacets( 0 ) );
                CHECK( &tSecondIntersectedFacetExpected == tIntersectedFacets( 1 ) );
                CHECK( std::abs( tIntersections( 0 ) - tIntersectionCoordinatesExpected( 0 ) ) < tObject.get_intersection_tolerance() );
                CHECK( std::abs( tIntersections( 1 ) - tIntersectionCoordinatesExpected( 1 ) ) < tObject.get_intersection_tolerance() );

                // check if the point is inside and compare it to expectations
                Object_Region tPointIsInside = check_if_node_is_inside_triangles( tIntersections, tTestPoint, 2, 1e-8 );

                CHECK( tPointIsInside == INSIDE );

                // check with the full raycast algorithm and see if they match
                tPointIsInside = raycast_point( tObject, tTestPoint );

                CHECK( tPointIsInside == INSIDE );

                // repeat test for point that is outside
                tTestPoint = { { 0.2, 0.6, 0.7 } };

                // Expected results
                tIntersectionCoordinatesExpected = { 0.715517872727636, 1.284482127272365 };

                tPreselection = preselect_triangles( tObject, tTestPoint, 2, tCandidateTriangles );

                REQUIRE( tCandidateTriangles.size() == 0 );
                REQUIRE( tPreselection == Preselection_Result::SUCCESS );

                tCandidatesExpected = { 0, 1, 2 };
                tPreselection       = preselect_triangles( tObject, tTestPoint, 0, tCandidateTriangles );
                REQUIRE( tPreselection == Preselection_Result::SUCCESS );

                REQUIRE( tCandidateTriangles.size() == 3 );
                CHECK( tCandidateTriangles( 0 ) == tCandidatesExpected( 0 ) );
                CHECK( tCandidateTriangles( 1 ) == tCandidatesExpected( 1 ) );
                CHECK( tCandidateTriangles( 2 ) == tCandidatesExpected( 2 ) );

                tIntersections = intersect_triangles_moller_trumbore( tCandidateTriangles, tObject, tTestPoint, 0, tIntersectedFacets );

                // tIntersectedTriangles           = intersect_triangles( tCandidateTriangles, tObject, tTestPoint, 0 );

                REQUIRE( tIntersectedFacets.size() == 2 );
                CHECK( *( tIntersectedFacets( 0 ) ) == tObject.get_facet( 1 ) );
                CHECK( *( tIntersectedFacets( 1 ) ) == tObject.get_facet( 2 ) );

                // although two facets are intersected, one of them should produce an error and be removed
                // tIntersectionCoordinates         = intersect_ray_with_facets( tIntersectedTriangles, tTestPoint, 0, Preselection_Result::SUCCESS );

                REQUIRE( tIntersections.size() == 2 );
                CHECK( std::abs( tIntersections( 0 ) - tIntersectionCoordinatesExpected( 0 ) ) < tObject.get_intersection_tolerance() );
                CHECK( std::abs( tIntersections( 1 ) - tIntersectionCoordinatesExpected( 1 ) ) < tObject.get_intersection_tolerance() );

                tPointIsInside = check_if_node_is_inside_triangles( tIntersections, tTestPoint, 0, 1e-8 );

                CHECK( tPointIsInside == OUTSIDE );

                tPointIsInside = raycast_point( tObject, tTestPoint );

                CHECK( tPointIsInside == OUTSIDE );
            }
            SECTION( "SDF: Raycast Free Function Test - 2D" )
            {
                // create object from object file
                std::string tObjectPath = tMorisRoot + "projects/GEN/SDF/test/data/rhombus.obj";
                Object      tObject( tObjectPath );

                // define test point
                Matrix< DDRMat > tTestPoint = { { -.25 }, { -0.3 } };

                // preselect in y direction and ensure the candidates and intersected facets are marked
                Vector< uint >   tIntersectedLines;
                Vector< Facet* > tCandidateLines;
                preselect_lines( tObject, tTestPoint, 1, tIntersectedLines, tCandidateLines );
                uint tIntersectedLinesExpected = 1;

                REQUIRE( tIntersectedLines.size() == 1 );
                REQUIRE( tCandidateLines.size() == 1 );
                CHECK( tIntersectedLines( 0 ) == tIntersectedLinesExpected );
                CHECK( tCandidateLines( 0 ) == &tObject.get_facet( 2 ) );

                // preselect in x direction and ensure the candidates and intersected facets are marked
                preselect_lines( tObject, tTestPoint, 0, tIntersectedLines, tCandidateLines );
                tIntersectedLinesExpected = 3;

                REQUIRE( tIntersectedLines.size() == 1 );
                REQUIRE( tCandidateLines.size() == 1 );
                CHECK( tIntersectedLines( 0 ) == tIntersectedLinesExpected );
                CHECK( tCandidateLines( 0 ) == &tObject.get_facet( 2 ) );

                // intersect the candidate facets and determine the intersection location
                Vector< real > tIntersectionCoordinates        = intersect_ray_with_facets( tCandidateLines, tTestPoint, 0, Preselection_Result::SUCCESS );
                real           tIntersectionCoordinateExpected = -0.2;

                REQUIRE( tIntersectionCoordinates.size() == 1 );
                CHECK( std::abs( tIntersectionCoordinates( 0 ) - tIntersectionCoordinateExpected ) < tObject.get_intersection_tolerance() );

                // determine if the point is inside/outside and check expectation
                Object_Region tRegion = check_if_node_is_inside_lines( tObject, tIntersectionCoordinates, tIntersectedLines, tTestPoint, 0 );

                CHECK( tRegion == OUTSIDE );

                // ensure the test gives the same result when the entire algorithm is called at once
                tRegion = raycast_point( tObject, tTestPoint );

                CHECK( tRegion == OUTSIDE );

                // repeat for a point inside the surface
                tTestPoint = { { -.25 }, { 0.2 } };

                preselect_lines( tObject, tTestPoint, 1, tIntersectedLines, tCandidateLines );

                REQUIRE( tCandidateLines.size() == 1 );
                REQUIRE( tIntersectedLines.size() == 0 );
                CHECK( *tCandidateLines( 0 ) == tObject.get_facet( 1 ) );

                tIntersectionCoordinates        = intersect_ray_with_facets( tCandidateLines, tTestPoint, 1, Preselection_Result::SUCCESS );
                tIntersectionCoordinateExpected = 0.25;

                REQUIRE( tIntersectionCoordinates.size() == 1 );
                CHECK( std::abs( tIntersectionCoordinates( 0 ) - tIntersectionCoordinateExpected ) < tObject.get_intersection_tolerance() );

                tRegion = check_if_node_is_inside_lines( tObject, tIntersectionCoordinates, tIntersectedLines, tTestPoint, 1 );

                CHECK( tRegion == INSIDE );

                tRegion = raycast_point( tObject, tTestPoint );

                CHECK( tRegion == INSIDE );

                // Repeat with a point that is on a facet
                tTestPoint( 0, 0 ) = 0.25;
                tTestPoint( 1, 0 ) = 0.25;

                tRegion = raycast_point( tObject, tTestPoint );

                CHECK( tRegion == INTERFACE );

                // Repeat with a point that is on a vertex
                tTestPoint( 0, 0 ) = 0.0;
                tTestPoint( 1, 0 ) = 0.5;

                tRegion = raycast_point( tObject, tTestPoint );

                CHECK( tRegion == INTERFACE );
            }
            SECTION( "SDF: Compute distance to facets test - 3D" )
            {
                // create triangle object from object file
                std::string tObjectPath = tMorisRoot + "projects/GEN/SDF/test/data/tetrahedron.obj";
                Object      tObject( tObjectPath );

                // define test point that is inside the object
                Matrix< DDRMat > tTestPoint = { { 0.9, 0.6, 0.7 } };

                Vector< Facet* > tIntersectionIndices;

                real tLineDistanceXExpected = 1.2844821272723648;     // facet index = 2
                real tLineDistanceYExpected = 0.91953300647060843;    // facet index = 1
                real tLineDistanceZExpected = 0.88055613061054816;    // facet index = 0

                // compute with raycast function
                Vector< real > tLineDistanceX = compute_distance_to_facets( tObject, tTestPoint, 0, tIntersectionIndices );
                Vector< real > tLineDistanceY = compute_distance_to_facets( tObject, tTestPoint, 1, tIntersectionIndices );
                Vector< real > tLineDistanceZ = compute_distance_to_facets( tObject, tTestPoint, 2, tIntersectionIndices );

                // compare
                REQUIRE( tLineDistanceX.size() == 1 );
                REQUIRE( tLineDistanceY.size() == 1 );
                REQUIRE( tLineDistanceZ.size() == 1 );
                CHECK( std::abs( tLineDistanceX( 0 ) - tLineDistanceXExpected ) < tObject.get_intersection_tolerance() );
                CHECK( std::abs( tLineDistanceY( 0 ) - tLineDistanceYExpected ) < tObject.get_intersection_tolerance() );
                CHECK( std::abs( tLineDistanceZ( 0 ) - tLineDistanceZExpected ) < tObject.get_intersection_tolerance() );
            }
            SECTION( "SDF: Compute distance to facets test - 2D" )
            {
                // create triangle object from object file
                std::string tObjectPath = tMorisRoot + "projects/GEN/SDF/test/data/rhombus.obj";
                Object      tObject( tObjectPath );

                // define test point
                Matrix< DDRMat > tTestPoint = { { -.25 }, { -0.3 } };

                // vector of facets that intersect the ray
                Vector< Facet* > tIntersectionIndices;

                // expected results
                Vector< real > tLineDistanceXExpected = { -0.2, 0.2 };
                Vector< real > tLineDistanceYExpected = { -0.25, 0.25 };

                // compute with raycast
                Vector< real > tLineDistanceX = compute_distance_to_facets( tObject, tTestPoint, 0, tIntersectionIndices );
                Vector< real > tLineDistanceY = compute_distance_to_facets( tObject, tTestPoint, 1, tIntersectionIndices );

                // compare
                REQUIRE( tLineDistanceX.size() == 2 );
                REQUIRE( tLineDistanceY.size() == 2 );
                CHECK( std::abs( tLineDistanceX( 0 ) - tLineDistanceXExpected( 0 ) ) < tObject.get_intersection_tolerance() );
                CHECK( std::abs( tLineDistanceX( 1 ) - tLineDistanceXExpected( 1 ) ) < tObject.get_intersection_tolerance() );
                CHECK( std::abs( tLineDistanceY( 0 ) - tLineDistanceYExpected( 0 ) ) < tObject.get_intersection_tolerance() );
                CHECK( std::abs( tLineDistanceY( 1 ) - tLineDistanceYExpected( 1 ) ) < tObject.get_intersection_tolerance() );
            }
        }
    }
}    // namespace moris::sdf
