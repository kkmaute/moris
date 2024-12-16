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

#if MORIS_HAVE_ARBORX
    // initialize Kokkos for the use in the spatial tree library ArborX
    std::unique_ptr< Kokkos::ScopeGuard > guard = !Kokkos::is_initialized() && !Kokkos::is_finalized() ? std::make_unique< Kokkos::ScopeGuard >() : nullptr;


    // number of random rays to cast to check result
    uint tNumRays = 1;

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

                // // Check the preselection results
                // REQUIRE( tCandidateTriangles.size() == 2 );
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
                mtk::Mesh_Region tPointIsInside = tObject.get_region_from_raycast( tTestPoint );
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

                // CHECK( std::abs( tIntersections( 0 ) - tIntersectionCoordinatesExpected( 0 ) ) < tEpsilon );
                // CHECK( std::isnan( tIntersections( 1 ) ) );
                // CHECK( std::abs( tIntersections( 2 ) - tIntersectionCoordinatesExpected( 1 ) ) < tEpsilon );

                // Repeat for all of them at the same time using batching
                tTestPoint = { { 0.9, 0.2 }, { 0.6, 0.6 }, { 0.7, 0.7 } };

                tPointIsInside = tObject.get_region_from_raycast( tTestPoint );
                REQUIRE( tPointIsInside == mtk::Mesh_Region::OUTSIDE );

                // Get the regions
                Vector< mtk::Mesh_Region > tRegions  = tObject.batch_get_region_from_raycast( tTestPoint );
                Vector< mtk::Mesh_Region > tExpected = { mtk::Mesh_Region::INSIDE, mtk::Mesh_Region::OUTSIDE };

                // Check each match
                REQUIRE( tRegions.size() == 2 );
                for ( uint iRegion = 0; iRegion < 2; ++iRegion )
                {
                    CHECK( tRegions( iRegion ) == tExpected( iRegion ) );
                }
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

                // // preselect in x direction and ensure the candidates and intersected facets are marked
                // Vector< uint > tCandidateLines           = tObject.preselect_with_arborx( tTestPoint, tDirection );
                // Vector< uint > tIntersectedLinesExpected = { 2, 3 };

                // REQUIRE( tCandidateLines.size() == 2 );
                // CHECK( tCandidateLines( 0 ) == tIntersectedLinesExpected( 0 ) );
                // CHECK( tCandidateLines( 1 ) == tIntersectedLinesExpected( 1 ) );

                // // intersect the candidate facets and determine the intersection location
                // Vector< real > tIntersections( tCandidateLines.size() );
                // for ( uint iCandidate = 0; iCandidate < tCandidateLines.size(); ++iCandidate )
                // {
                //     tIntersections( iCandidate ) = tObject.moller_trumbore( tCandidateLines( iCandidate ), tTestPoint, tDirection );
                // }
                // Vector< real > tIntersectionCoordinatesExpected = { 0.05, 0.45 };

                // REQUIRE( tIntersections.size() == 2 );
                // CHECK( std::abs( tIntersections( 0 ) - tIntersectionCoordinatesExpected( 0 ) ) < tEpsilon );
                // CHECK( std::abs( tIntersections( 1 ) - tIntersectionCoordinatesExpected( 1 ) ) < tEpsilon );

                // ensure the test gives the same result when the entire algorithm is called, repeat many times for many different random vectors

                mtk::Mesh_Region tRegion = tObject.get_region_from_raycast( tTestPoint );
                REQUIRE( tRegion == mtk::Mesh_Region::OUTSIDE );


                // repeat for a point inside the surface
                tTestPoint = { { -.25 }, { 0.2 } };

                // tDirection      = { { 0.0, 1.0 } };
                // tCandidateLines = tObject.preselect_with_arborx( tTestPoint, tDirection );

                // REQUIRE( tCandidateLines.size() == 1 );
                // CHECK( tCandidateLines( 0 ) == 1 );

                // tIntersections.resize( tCandidateLines.size() );
                // for ( uint iCandidate = 0; iCandidate < tCandidateLines.size(); ++iCandidate )
                // {
                //     tIntersections( iCandidate ) = tObject.moller_trumbore( tCandidateLines( iCandidate ), tTestPoint, tDirection );
                // }

                // tIntersectionCoordinatesExpected = { 0.05 };

                // REQUIRE( tIntersections.size() == 1 );
                // CHECK( std::abs( tIntersections( 0 ) - tIntersectionCoordinatesExpected( 0 ) ) < tEpsilon );

                tRegion = tObject.get_region_from_raycast( tTestPoint );
                REQUIRE( tRegion == mtk::Mesh_Region::INSIDE );

                // Repeat with a point that is on a facet
                tTestPoint( 0, 0 ) = 0.25;
                tTestPoint( 1, 0 ) = 0.25;

                tRegion = tObject.get_region_from_raycast( tTestPoint );
                REQUIRE( tRegion == mtk::Mesh_Region::INTERFACE );


                // Repeat with a point that is on a vertex
                tTestPoint( 0, 0 ) = 0.0;
                tTestPoint( 1, 0 ) = 0.5;

                tRegion = tObject.get_region_from_raycast( tTestPoint );
                REQUIRE( tRegion == mtk::Mesh_Region::INTERFACE );

                // Repeat for all of them at the same time using batching
                tTestPoint = { { -0.25, -0.25, 0.25, 0.0 }, { -0.3, 0.2, 0.25, 0.5 } };

                // Get the regions
                Vector< mtk::Mesh_Region > tRegions  = tObject.batch_get_region_from_raycast( tTestPoint );
                Vector< mtk::Mesh_Region > tExpected = { mtk::Mesh_Region::OUTSIDE, mtk::Mesh_Region::INSIDE, mtk::Mesh_Region::INTERFACE, mtk::Mesh_Region::INTERFACE };

                // Check each match
                REQUIRE( tRegions.size() == 4 );
                for ( uint iRegion = 0; iRegion < 4; ++iRegion )
                {
                    CHECK( tRegions( iRegion ) == tExpected( iRegion ) );
                }
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

                real tLineDistanceXExpected = 0.384482127272365;    // facet index = 2
                real tLineDistanceYExpected = 0.319533006470609;    // facet index = 1
                real tLineDistanceZExpected = 0.180556130610548;    // facet index = 0

                // compute with raycast function
                bool                              tWarning;
                Vector< std::pair< uint, real > > tLineDistanceX = tObject.cast_single_ray( tTestPoint, tDirection, tWarning );
                tDirection                                       = { { 0.0 }, { 1.0 }, { 0.0 } };
                Vector< std::pair< uint, real > > tLineDistanceY = tObject.cast_single_ray( tTestPoint, tDirection, tWarning );
                tDirection                                       = { { 0.0 }, { 0.0 }, { 1.0 } };
                Vector< std::pair< uint, real > > tLineDistanceZ = tObject.cast_single_ray( tTestPoint, tDirection, tWarning );

                // compare
                REQUIRE( tLineDistanceX.size() == 1 );
                REQUIRE( tLineDistanceY.size() == 1 );
                REQUIRE( tLineDistanceZ.size() == 1 );
                CHECK( std::abs( tLineDistanceX( 0 ).second - tLineDistanceXExpected ) < tEpsilon );
                CHECK( std::abs( tLineDistanceY( 0 ).second - tLineDistanceYExpected ) < tEpsilon );
                CHECK( std::abs( tLineDistanceZ( 0 ).second - tLineDistanceZExpected ) < tEpsilon );

                // batch all 3 rays and check that the result is correct
                Matrix< DDRMat >           tOrigins    = { { 0.9, 0.9, 0.9 }, { 0.6, 0.6, 0.6 }, { 0.7, 0.7, 0.7 } };
                Vector< Matrix< DDRMat > > tDirections = { { { 1.0, 1.0 }, { 0.0, 0.0 }, { 0.0, 0.0 } }, { { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 }, { 0.0, 0.0, 0.0 } }, { { 0.0 }, { 0.0 }, { 1.0 } } };

                Vector< Vector< bool > >                              tWarnings;
                Vector< Vector< Vector< std::pair< uint, real > > > > tLineDistances = tObject.cast_batch_of_rays( tOrigins, tDirections, tWarnings );

                REQUIRE( tLineDistances.size() == 3 );
                REQUIRE( tLineDistances( 0 ).size() == 2 );
                REQUIRE( tLineDistances( 1 ).size() == 3 );
                REQUIRE( tLineDistances( 2 ).size() == 1 );
                CHECK( std::abs( tLineDistances( 0 )( 0 )( 0 ).second - tLineDistanceXExpected ) < tEpsilon );
                CHECK( tLineDistances( 0 )( 0 )( 0 ).first == 2 );
                CHECK( std::abs( tLineDistances( 0 )( 1 )( 0 ).second - tLineDistanceXExpected ) < tEpsilon );
                CHECK( tLineDistances( 0 )( 1 )( 0 ).first == 2 );
                CHECK( std::abs( tLineDistances( 1 )( 0 )( 0 ).second - tLineDistanceYExpected ) < tEpsilon );
                CHECK( tLineDistances( 1 )( 0 )( 0 ).first == 1 );
                CHECK( std::abs( tLineDistances( 1 )( 1 )( 0 ).second - tLineDistanceYExpected ) < tEpsilon );
                CHECK( tLineDistances( 1 )( 1 )( 0 ).first == 1 );
                CHECK( std::abs( tLineDistances( 1 )( 2 )( 0 ).second - tLineDistanceYExpected ) < tEpsilon );
                CHECK( tLineDistances( 1 )( 2 )( 0 ).first == 1 );
                CHECK( std::abs( tLineDistances( 2 )( 0 )( 0 ).second - tLineDistanceZExpected ) < tEpsilon );
                CHECK( tLineDistances( 2 )( 0 )( 0 ).first == 0 );

                // batch again using the other functionality to cast the same direction on every origin
                Matrix< DDRMat > tSameDirections = { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };

                Vector< Vector< bool > > tWarningSame;
                tLineDistances = tObject.cast_batch_of_rays( tOrigins, tSameDirections, tWarningSame );

                REQUIRE( tLineDistances.size() == 3 );
                REQUIRE( tLineDistances( 0 ).size() == 3 );
                REQUIRE( tLineDistances( 1 ).size() == 3 );
                REQUIRE( tLineDistances( 2 ).size() == 3 );
                CHECK( std::abs( tLineDistances( 0 )( 0 )( 0 ).second - tLineDistanceXExpected ) < tEpsilon );
                CHECK( tLineDistances( 0 )( 0 )( 0 ).first == 2 );
                CHECK( std::abs( tLineDistances( 1 )( 1 )( 0 ).second - tLineDistanceYExpected ) < tEpsilon );
                CHECK( tLineDistances( 1 )( 1 )( 0 ).first == 1 );
                CHECK( std::abs( tLineDistances( 2 )( 2 )( 0 ).second - tLineDistanceZExpected ) < tEpsilon );
                CHECK( tLineDistances( 2 )( 2 )( 0 ).first == 0 );
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
                Vector< real > tLineDistanceXExpected = { 0.05, 0.45 };
                Vector< real > tLineDistanceYExpected = { 0.05, 0.55 };

                // compute with raycast
                bool                              tWarning;
                Vector< std::pair< uint, real > > tLineDistanceX = tObject.cast_single_ray( tTestPoint, tDirection, tWarning );
                tDirection                                       = { { 0.0 }, { 1.0 } };
                Vector< std::pair< uint, real > > tLineDistanceY = tObject.cast_single_ray( tTestPoint, tDirection, tWarning );

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

#endif    // MORIS_HAVE_ARBORX
