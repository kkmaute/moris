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

// mtk
#include "cl_MTK_Enums.hpp"

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
#include "SDF_Tools.hpp"

namespace moris::sdf
{

#ifdef MORIS_HAVE_ARBORX
    // initialize Kokkos for the use in the spatial tree library ArborX
    std::unique_ptr< Kokkos::ScopeGuard > guard = !Kokkos::is_initialized() && !Kokkos::is_finalized() ? std::make_unique< Kokkos::ScopeGuard >() : nullptr;
#endif

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
                // Tolerance for results
                // real tEpsilon = 1e-8;

                // number of random rays to cast to check result
                // uint tNumRays = 100;

                // create triangle object from object file
                std::string tObjectPath = tMorisRoot + "projects/GEN/SDF/test/data/tetrahedron.obj";
                Object      tObject( tObjectPath );

                // define test point that is inside the object
                Matrix< DDRMat > tTestPoint = { { 0.9, 0.6, 0.7 } };

                // Define ray direction
                // Matrix< DDRMat > tDirection = { { 1.0 }, { 0.0 }, { 0.0 } };

                // // preselect in x direction and ensure they are correct
                // Vector< uint > tCandidatesExpected = { 1, 0, 2 };
                // Vector< uint > tCandidateTriangles = tObject.preselect_with_arborx( tTestPoint, tDirection );

                // REQUIRE( tCandidateTriangles.size() == 3 );
                // // REQUIRE( tPreselection == Preselection_Result::SUCCESS );
                // CHECK( tCandidatesExpected( 0 ) == tCandidateTriangles( 0 ) );
                // CHECK( tCandidatesExpected( 1 ) == tCandidateTriangles( 1 ) );
                // CHECK( tCandidatesExpected( 2 ) == tCandidateTriangles( 2 ) );

                // // repeat for y direction
                // tDirection          = { { 0.0 }, { 1.0 }, { 0.0 } };
                // tCandidatesExpected = { 1, 0 };
                // tCandidateTriangles = tObject.preselect_with_arborx( tTestPoint, tDirection );

                // REQUIRE( tCandidateTriangles.size() == 2 );
                // // REQUIRE( tPreselection == Preselection_Result::SUCCESS );
                // CHECK( tCandidatesExpected( 0 ) == tCandidateTriangles( 0 ) );
                // CHECK( tCandidatesExpected( 1 ) == tCandidateTriangles( 1 ) );

                // // repeat for z direction
                // tDirection          = { { 0.0 }, { 0.0 }, { 1.0 } };
                // tCandidatesExpected = { 1, 0 };
                // tCandidateTriangles = tObject.preselect_with_arborx( tTestPoint, tDirection );
                // // tPreselection       = preselect_triangles( tObject, tTestPoint, 2, tCandidateTriangles ); brendan

                // // Check the preselection results
                // REQUIRE( tCandidateTriangles.size() == 2 );
                // // REQUIRE( tPreselection == Preselection_Result::SUCCESS ); brendan
                // CHECK( tCandidatesExpected( 0 ) == tCandidateTriangles( 0 ) );
                // CHECK( tCandidatesExpected( 1 ) == tCandidateTriangles( 1 ) );

                // // check moller trumbore algorithm for proper intersection location and parent facet determination
                // Vector< real > tIntersectionCoordinatesExpected = { 0.1805561306105482 };

                // Vector< real > tIntersections( 2 );
                // for ( uint iCandidate = 0; iCandidate < 2; ++iCandidate )
                // {
                //     tIntersections( iCandidate ) = tObject.moller_trumbore( tCandidateTriangles( iCandidate ), tTestPoint, tDirection );
                // }
                // CHECK( std::isnan( tIntersections( 0 ) ) );
                // CHECK( std::abs( tIntersections( 1 ) - tIntersectionCoordinatesExpected( 0 ) ) < tEpsilon );

                // check if the point is inside and compare it to expectations
                // mtk::mtk::Mesh_Region tPointIsInside = check_if_node_is_inside_triangles( tIntersections, tTestPoint, 2, 1e-8 );

                // CHECK( tPointIsInside == INSIDE );

                // cast a bunch of random rays and ensure they all return the correct result
                mtk::Mesh_Region tPointIsInside;
                tPointIsInside = tObject.get_region_from_raycast( tTestPoint );
                REQUIRE( tPointIsInside == mtk::Mesh_Region::INSIDE );

                // repeat test for point that is outside
                tTestPoint = { { 0.2, 0.6, 0.7 } };

                // Expected results
                // tIntersectionCoordinatesExpected = { 0.715517872727636, 1.284482127272365 };

                // tCandidateTriangles = tObject.preselect_with_arborx( tTestPoint, tDirection );

                // REQUIRE( tCandidateTriangles.size() == 0 );

                // tCandidatesExpected = { 1, 0, 2 };
                // tDirection          = { { 1.0 }, { 0.0 }, { 0.0 } };
                // tCandidateTriangles = tObject.preselect_with_arborx( tTestPoint, tDirection );

                // REQUIRE( tCandidateTriangles.size() == 3 );
                // CHECK( tCandidateTriangles( 0 ) == tCandidatesExpected( 0 ) );
                // CHECK( tCandidateTriangles( 1 ) == tCandidatesExpected( 1 ) );
                // CHECK( tCandidateTriangles( 2 ) == tCandidatesExpected( 2 ) );

                // tIntersections.resize( 3 );
                // for ( uint iCandidate = 0; iCandidate < 3; ++iCandidate )
                // {
                //     tIntersections( iCandidate ) = tObject.moller_trumbore( tCandidateTriangles( iCandidate ), tTestPoint, tDirection );
                // }

                // tIntersectedTriangles           = intersect_triangles( tCandidateTriangles, tObject, tTestPoint, 0 );

                // REQUIRE( tIntersectedFacets.size() == 2 );
                // CHECK( *( tIntersectedFacets( 0 ) ) == tObject.get_facet( 1 ) );
                // CHECK( *( tIntersectedFacets( 1 ) ) == tObject.get_facet( 2 ) );

                // although two facets are intersected, one of them should produce an error and be removed
                // tIntersectionCoordinates         = intersect_ray_with_facets( tIntersectedTriangles, tTestPoint, 0, Preselection_Result::SUCCESS );

                // CHECK( std::abs( tIntersections( 0 ) - tIntersectionCoordinatesExpected( 0 ) ) < tEpsilon );
                // CHECK( std::isnan( tIntersections( 1 ) ) );
                // CHECK( std::abs( tIntersections( 2 ) - tIntersectionCoordinatesExpected( 1 ) ) < tEpsilon );

                // tPointIsInside = check_if_node_is_inside_triangles( tIntersections, tTestPoint, 0, 1e-8 );

                // CHECK( tPointIsInside == OUTSIDE );
                tPointIsInside = tObject.get_region_from_raycast( tTestPoint );
                REQUIRE( tPointIsInside == mtk::Mesh_Region::OUTSIDE );
            }
            SECTION( "SDF: Raycast Free Function Test - 2D" )
            {
                // Tolerance for results
                // real tEpsilon = 1e-8;

                // create object from object file
                std::string tObjectPath = tMorisRoot + "projects/GEN/SDF/test/data/rhombus.obj";
                Object      tObject( tObjectPath );

                // define test point
                Matrix< DDRMat > tTestPoint = { { -.25 }, { -0.3 } };

                // define test direction
                // Matrix< DDRMat > tDirection = { { 1.0, 0.0 } };

                // preselect in y direction and ensure the candidates and intersected facets are marked
                // Vector< uint > tCandidateLines =tObject.preselect_with_arborx( tObject tTestPoint, tDirection ); brendan
                // uint tIntersectedLinesExpected = 1;

                // REQUIRE( tIntersectedLines.size() == 1 );
                // REQUIRE( tCandidateLines.size() == 1 );
                // CHECK( tIntersectedLines( 0 ) == tIntersectedLinesExpected );
                // CHECK( tCandidateLines( 0 ) == &tObject.get_facet( 2 ) );

                // // preselect in x direction and ensure the candidates and intersected facets are marked
                // Vector< uint > tCandidateLines           = tObject.preselect_with_arborx( tTestPoint, tDirection );
                // Vector< uint > tIntersectedLinesExpected = { 2, 3 };

                // REQUIRE( tCandidateLines.size() == 2 );
                // CHECK( tCandidateLines( 0 ) == tIntersectedLinesExpected( 0 ) );
                // CHECK( tCandidateLines( 0 ) == tIntersectedLinesExpected( 1 ) );

                // // intersect the candidate facets and determine the intersection location
                // Vector< real > tIntersections( tCandidateLines.size() );
                // for ( uint iCandidate : tCandidateLines )
                // {
                //     tIntersections( iCandidate ) = tObject.moller_trumbore( iCandidate, tTestPoint, tDirection );
                // }
                // real tIntersectionCoordinateExpected = -0.2;

                // REQUIRE( tIntersections.size() == 2 );
                // CHECK( std::abs( tIntersections( 1 ) - tIntersectionCoordinateExpected ) < tEpsilon );

                // // determine if the point is inside/outside and check expectation
                // mtk::mtk::Mesh_Region tRegion = check_if_node_is_inside_lines( tObject, tIntersectionCoordinates, tIntersectedLines, tTestPoint, 0 ); brendan

                // CHECK( tRegion == OUTSIDE );

                // ensure the test gives the same result when the entire algorithm is called at once
                mtk::Mesh_Region tRegion = tObject.get_region_from_raycast( tTestPoint );

                CHECK( tRegion == mtk::Mesh_Region::OUTSIDE );

                // repeat for a point inside the surface
                tTestPoint = { { -.25 }, { 0.2 } };

                // tDirection      = { { 0.0, 1.0 } };
                // tCandidateLines = tObject.preselect_with_arborx( tTestPoint, tDirection );

                // REQUIRE( tCandidateLines.size() == 1 );
                // CHECK( tCandidateLines( 0 ) == 1 );

                // tIntersections.resize( tCandidateLines.size() );
                // for ( uint iCandidate : tCandidateLines )
                // {
                //     tIntersections( iCandidate ) = tObject.moller_trumbore( iCandidate, tTestPoint, tDirection );
                // }

                // tIntersectionCoordinateExpected = 0.25;

                // REQUIRE( tIntersections.size() == 1 );
                // CHECK( std::abs( tIntersections( 0 ) - tIntersectionCoordinateExpected ) < tEpsilon );

                // tRegion = check_if_node_is_inside_lines( tObject, tIntersectionCoordinates, tIntersectedLines, tTestPoint, 1 );

                // CHECK( tRegion == INSIDE );

                tRegion = tObject.get_region_from_raycast( tTestPoint );

                CHECK( tRegion == mtk::Mesh_Region::INSIDE );

                // Repeat with a point that is on a facet
                tTestPoint( 0, 0 ) = 0.25;
                tTestPoint( 1, 0 ) = 0.25;

                tRegion = tObject.get_region_from_raycast( tTestPoint );

                CHECK( tRegion == mtk::Mesh_Region::INTERFACE );

                // Repeat with a point that is on a vertex
                tTestPoint( 0, 0 ) = 0.0;
                tTestPoint( 1, 0 ) = 0.5;

                tRegion = tObject.get_region_from_raycast( tTestPoint );

                CHECK( tRegion == mtk::Mesh_Region::INTERFACE );
            }
            SECTION( "SDF: Compute distance to facets test - 3D" )
            {
                // Tolerance for results
                real tEpsilon = 1e-8;

                // create triangle object from object file
                std::string tObjectPath = tMorisRoot + "projects/GEN/SDF/test/data/tetrahedron.obj";
                Object      tObject( tObjectPath );

                // define test point that is inside the object
                Matrix< DDRMat > tTestPoint = { { 0.9, 0.6, 0.7 } };

                // Define the direction of the ray
                Matrix< DDRMat > tDirection = { { 1.0 }, { 0.0 }, { 0.0 } };

                real tLineDistanceXExpected = 1.2844821272723648;     // facet index = 2
                real tLineDistanceYExpected = 0.91953300647060843;    // facet index = 1
                real tLineDistanceZExpected = 0.88055613061054816;    // facet index = 0

                // compute with raycast function
                mtk::Intersection_Vector tLineDistanceX = tObject.cast_single_ray( tTestPoint, tDirection );
                tDirection                              = { { 0.0 }, { 1.0 }, { 0.0 } };
                mtk::Intersection_Vector tLineDistanceY = tObject.cast_single_ray( tTestPoint, tDirection );
                tDirection                              = { { 0.0 }, { 0.0 }, { 1.0 } };
                mtk::Intersection_Vector tLineDistanceZ = tObject.cast_single_ray( tTestPoint, tDirection );

                // compare
                REQUIRE( tLineDistanceX.size() == 1 );
                REQUIRE( tLineDistanceY.size() == 1 );
                REQUIRE( tLineDistanceZ.size() == 1 );
                CHECK( std::abs( tLineDistanceX( 0 ).second - tLineDistanceXExpected ) < tEpsilon );
                CHECK( std::abs( tLineDistanceY( 0 ).second - tLineDistanceYExpected ) < tEpsilon );
                CHECK( std::abs( tLineDistanceZ( 0 ).second - tLineDistanceZExpected ) < tEpsilon );
            }
            SECTION( "SDF: Compute distance to facets test - 2D" )
            {
                // Tolerance for results
                real tEpsilon = 1e-8;

                // create triangle object from object file
                std::string tObjectPath = tMorisRoot + "projects/GEN/SDF/test/data/rhombus.obj";
                Object      tObject( tObjectPath );

                // define test point
                Matrix< DDRMat > tTestPoint = { { -.25 }, { -0.3 } };

                // Define the direction of the ray
                Matrix< DDRMat > tDirection = { { 1.0 }, { 0.0 } };

                // expected results
                Vector< real > tLineDistanceXExpected = { -0.2, 0.2 };
                Vector< real > tLineDistanceYExpected = { -0.25, 0.25 };

                // compute with raycast
                mtk::Intersection_Vector tLineDistanceX = tObject.cast_single_ray( tTestPoint, tDirection );
                tDirection                              = { { 0.0 }, { 1.0 } };
                mtk::Intersection_Vector tLineDistanceY = tObject.cast_single_ray( tTestPoint, tDirection );

                // compare
                REQUIRE( tLineDistanceX.size() == 2 );
                REQUIRE( tLineDistanceY.size() == 2 );
                CHECK( std::abs( tLineDistanceX( 0 ).second - tLineDistanceXExpected( 0 ) ) < tEpsilon );
                CHECK( std::abs( tLineDistanceX( 1 ).second - tLineDistanceXExpected( 1 ) ) < tEpsilon );
                CHECK( std::abs( tLineDistanceY( 0 ).second - tLineDistanceYExpected( 0 ) ) < tEpsilon );
                CHECK( std::abs( tLineDistanceY( 1 ).second - tLineDistanceYExpected( 1 ) ) < tEpsilon );
            }
        }
    }
}    // namespace moris::sdf
