/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_MTK_Surface_Mesh.cpp
 *
 */

#include "paths.hpp"    // for moris root
#include "catch.hpp"
#include "fn_check_equal.hpp"
#include "cl_MTK_Surface_Mesh.hpp"
#include "fn_MTK_Load_External_Surface_Mesh.hpp"
#include "cl_Communication_Tools.hpp"

namespace moris::mtk
{
#if MORIS_HAVE_ARBORX
    // initialize Kokkos for the use in the spatial tree library ArborX
    std::unique_ptr< Kokkos::ScopeGuard > guard = !Kokkos::is_initialized() && !Kokkos::is_finalized() ? std::make_unique< Kokkos::ScopeGuard >() : nullptr;
#endif

    // get root from environment
    std::string tMorisRoot = moris::get_base_moris_dir();

    // Tests loading in the surface mesh from an obj file, and computing various quanitities of interest and their sensitivities
    TEST_CASE( "MTK Surface Mesh", "[MTK],[MTK_Surface_Mesh]" )
    {
        Matrix< DDRMat >                     tCoordsExpected     = { { 2.25, 0.25, 1.0, 1.25 }, { 1.25, 0.5, -0.25, 0.5 } };
        Vector< Vector< moris_index > >      tConnExpected       = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 } };
        Matrix< DDRMat >                     tNormalsExpected    = { { -0.35112344, -0.70710678, 0.94868330, 0.60000000 }, { 0.93632918, -0.70710678, -0.31622777, -0.80000000 } };
        Vector< real >                       tMeasureExpected    = { 2.13600094, 1.06066017, 0.79056942, 1.25000000 };
        Vector< Vector< Matrix< DDRMat > > > tNormalSensExpected = {
            { { { 0.15391713, -0.41044567 }, { 0.05771892, -0.15391713 } }, { { -0.15391713, 0.41044567 }, { -0.05771892, 0.15391713 } }, { { 0.00000000, 0.00000000 }, { 0.00000000, 0.00000000 } }, { { 0.00000000, 0.00000000 }, { 0.00000000, 0.00000000 } } },    // dN1/dv4
            { { { 0.00000000, 0.00000000 }, { 0.00000000, 0.00000000 } }, { { -0.47140452, -0.47140452 }, { 0.47140452, 0.47140452 } }, { { 0.47140452, 0.47140452 }, { -0.47140452, -0.47140452 } }, { { 0.00000000, 0.00000000 }, { 0.00000000, 0.00000000 } } },
            { { { 0.00000000, 0.00000000 }, { 0.00000000, 0.00000000 } }, { { 0.00000000, 0.00000000 }, { 0.00000000, 0.00000000 } }, { { 0.37947332, -0.12649111 }, { 1.13841996, -0.37947332 } }, { { -0.37947332, 0.12649111 }, { -1.13841996, 0.37947332 } } },    // dN3/dv4
            { { { -0.38400000, 0.51200000 }, { -0.28800000, 0.38400000 } }, { { 0.00000000, 0.00000000 }, { 0.00000000, 0.00000000 } }, { { 0.00000000, 0.00000000 }, { 0.00000000, 0.00000000 } }, { { 0.38400000, -0.51200000 }, { 0.28800000, -0.38400000 } } }
        };
        Vector< Vector< Matrix< DDRMat > > > tCenterSensExpected = {
            { { { 0.50000000, 0.00000000 }, { 0.00000000, 0.50000000 } }, { { 0.50000000, 0.00000000 }, { 0.00000000, 0.50000000 } }, { { 0.00000000, 0.00000000 }, { 0.00000000, 0.00000000 } }, { { 0.00000000, 0.00000000 }, { 0.00000000, 0.00000000 } } },
            { { { 0.00000000, 0.00000000 }, { 0.00000000, 0.00000000 } }, { { 0.50000000, 0.00000000 }, { 0.00000000, 0.50000000 } }, { { 0.50000000, 0.00000000 }, { 0.00000000, 0.50000000 } }, { { 0.00000000, 0.00000000 }, { 0.00000000, 0.00000000 } } },
            { { { 0.00000000, 0.00000000 }, { 0.00000000, 0.00000000 } }, { { 0.00000000, 0.00000000 }, { 0.00000000, 0.00000000 } }, { { 0.50000000, 0.00000000 }, { 0.00000000, 0.50000000 } }, { { 0.50000000, 0.00000000 }, { 0.00000000, 0.50000000 } } },
            { { { 0.50000000, 0.00000000 }, { 0.00000000, 0.50000000 } }, { { 0.00000000, 0.00000000 }, { 0.00000000, 0.00000000 } }, { { 0.00000000, 0.00000000 }, { 0.00000000, 0.00000000 } }, { { 0.50000000, 0.00000000 }, { 0.00000000, 0.50000000 } } }
        };
        // load a surface mesh from file
        std::string    tFilePath = tMorisRoot + "/projects/GEN/test/data/triangle_sensitivity_oblique.obj";
        Vector< real > tOffsets  = { 0.0, 0.0 };
        Vector< real > tScales   = { 1.0, 1.0 };
        Surface_Mesh   tSurfaceMesh( load_vertices_from_object_file( tFilePath, tOffsets, tScales ), load_facets_from_object_file( tFilePath ) );

        // Check number of vertices and facets
        REQUIRE( tSurfaceMesh.get_number_of_vertices() == tCoordsExpected.n_cols() );
        REQUIRE( tSurfaceMesh.get_number_of_facets() == tConnExpected.size() );

        // Check vertex coordinates
        check_equal( tSurfaceMesh.get_all_vertex_coordinates(), tCoordsExpected );

        // Compute the facet measure
        Vector< real > tFacetMeasures = tSurfaceMesh.compute_facet_measure();

        // Loop over the surface mesh and check for correct facet normals and connectivity
        for ( uint iF = 0; iF < tConnExpected.size(); ++iF )
        {
            Vector< moris_index > tFacetVertices = tSurfaceMesh.get_facets_vertex_indices( iF );
            Matrix< DDRMat >      tFacetNormal   = tSurfaceMesh.get_facet_normal( iF );
            for ( uint iV = 0; iV < tConnExpected( iF ).size(); ++iV )
            {
                // Check facet connectivity
                REQUIRE( tFacetVertices( iV ) == tConnExpected( iF )( iV ) );

                // Check facet normal component
                CHECK( tFacetNormal( iV ) == Approx( tNormalsExpected( iV, iF ) ) );
            }

            // Check facet measure
            CHECK( tFacetMeasures( iF ) == Approx( tMeasureExpected( iF ) ) );

            // Check normal vector sensitivity wrt to every vertex in the mesh
            for ( uint iV = 0; iV < tCoordsExpected.n_cols(); iV++ )
            {
                Matrix< DDRMat > tNormalSens = tSurfaceMesh.compute_dfacet_normal_dvertex( iF, iV );
                Matrix< DDRMat > tCenterSens = tSurfaceMesh.compute_dfacet_centroid_dvertex( iF, iV );
                check_equal( tNormalSens, tNormalSensExpected( iF )( iV ), 1e8 );
                check_equal( tCenterSens, tCenterSensExpected( iF )( iV ), 1e8 );
            }
        }

        SECTION( "Raycast Region, Distance, and Sensitivities - 2D" )
        {
            // FIXME: It would be better to FD the sensitivities rather than hardcoding the expected values
            Matrix< DDRMat > tdRdOExpected = { { 0.5, 0.5 } };
            Matrix< DDRMat > tdRdDExpected = { { 0.1125, 0.1125 } };
            Matrix< DDRMat > tdRdVExpected = { { -0.28333333333, -0.28333333333 }, { -0.2166666666667, -0.2166666666667 } };

            // define test point and direction
            Matrix< DDRMat > tTestPoint = { { 0.8 }, { 0.4 } };
            Matrix< DDRMat > tDirection = { { -1.0, -1.0 } };

            // Check the region of this point
            mtk::Mesh_Region tRegion = tSurfaceMesh.get_region_from_raycast( tTestPoint );
            REQUIRE( tRegion == mtk::Mesh_Region::INSIDE );

            // Get the distance and sensitivities
            bool                tWarning;
            Intersection_Vector tDistance = tSurfaceMesh.cast_single_ray( tTestPoint, tDirection, tWarning );
            REQUIRE( tDistance.size() == 1 );
            CHECK( tDistance( 0 ).second == Approx( 0.2250 ) );    // distance
            CHECK( tDistance( 0 ).first == 1 );                    // intersected facet
            CHECK( not tWarning );                                 // no warning

            // Get the sensitivity wrt to the ray origin
            Matrix< DDRMat > tdRdO = tSurfaceMesh.compute_draycast_dorigin( tTestPoint, tDirection, tDistance( 0 ).first );

            // Get the sensitivity wrt to the ray direction
            Matrix< DDRMat > tdRdD = tSurfaceMesh.compute_draycast_ddirection( tTestPoint, tDirection, tDistance( 0 ).first );

            // Get the sensitivity wrt to the vertices of the intersected facet
            Matrix< DDRMat > tdRdV = tSurfaceMesh.compute_draycast_dvertices( tTestPoint, tDirection, tDistance( 0 ).first );

            // Loop over sensitivities and check they are correct
            for ( uint iS = 0; iS < 2; iS++ )    // Loop over spatial dimensions
            {
                CHECK( tdRdO( iS ) == Approx( tdRdOExpected( iS ) ) );
                CHECK( tdRdD( iS ) == Approx( tdRdDExpected( iS ) ) );

                for ( uint iV = 0; iV < 2; iV++ )    // Loop over vertices
                {
                    CHECK( tdRdV( iV, iS ) == Approx( tdRdVExpected( iV, iS ) ) );
                }
            }
        }
    }
    // Test for raycasting
    TEST_CASE( "MTK Surface Mesh Raycast", "[MTK], [MTK_Surface_Mesh], [Raycast]" )
    {
        if ( par_size() == 1 )
        {
            SECTION( "SDF: Raycast Free Function Test - 3D" )
            {
                // create surface mesh from object file
                std::string    tFilePath = tMorisRoot + "projects/GEN/SDF/test/data/tetrahedron.obj";
                Vector< real > tOffsets  = { 0.0, 0.0 };
                Vector< real > tScales   = { 1.0, 1.0 };
                Surface_Mesh   tSurfaceMesh( load_vertices_from_object_file( tFilePath, tOffsets, tScales ), load_facets_from_object_file( tFilePath ) );

                // define test point that is inside the object
                Matrix< DDRMat > tTestPoint = { { 0.9, 0.6, 0.7 } };

                // Define ray direction
                // Matrix< DDRMat > tDirection = { { 1.0 }, { 0.0 }, { 0.0 } };

                // // preselect in x direction and ensure they are correct
                // Vector< uint > tCandidatesExpected = { 1, 0, 2 };
                // Vector< uint > tCandidateTriangles = tSurfaceMesh.preselect_with_arborx( tTestPoint, tDirection );

                // REQUIRE( tCandidateTriangles.size() == 3 );
                // CHECK( tCandidatesExpected( 0 ) == tCandidateTriangles( 0 ) );
                // CHECK( tCandidatesExpected( 1 ) == tCandidateTriangles( 1 ) );
                // CHECK( tCandidatesExpected( 2 ) == tCandidateTriangles( 2 ) );

                // // repeat for y direction
                // tDirection          = { { 0.0 }, { 1.0 }, { 0.0 } };
                // tCandidatesExpected = { 1, 0 };
                // tCandidateTriangles = tSurfaceMesh.preselect_with_arborx( tTestPoint, tDirection );

                // REQUIRE( tCandidateTriangles.size() == 2 );
                // CHECK( tCandidatesExpected( 0 ) == tCandidateTriangles( 0 ) );
                // CHECK( tCandidatesExpected( 1 ) == tCandidateTriangles( 1 ) );

                // // repeat for z direction
                // tDirection          = { { 0.0 }, { 0.0 }, { 1.0 } };
                // tCandidatesExpected = { 1, 0 };
                // tCandidateTriangles = tSurfaceMesh.preselect_with_arborx( tTestPoint, tDirection );
                // // tPreselection       = preselect_triangles( tSurfaceMesh, tTestPoint, 2, tCandidateTriangles );

                // // Check the preselection results
                // REQUIRE( tCandidateTriangles.size() == 2 );
                // CHECK( tCandidatesExpected( 0 ) == tCandidateTriangles( 0 ) );
                // CHECK( tCandidatesExpected( 1 ) == tCandidateTriangles( 1 ) );

                // // check moller trumbore algorithm for proper intersection location and parent facet determination
                // Vector< real > tIntersectionCoordinatesExpected = { 0.1805561306105482 };

                // Vector< real > tIntersections( 2 );
                // for ( uint iCandidate = 0; iCandidate < 2; ++iCandidate )
                // {
                //     tIntersections( iCandidate ) = tSurfaceMesh.moller_trumbore( tCandidateTriangles( iCandidate ), tTestPoint, tDirection );
                // }
                // CHECK( std::isnan( tIntersections( 0 ) ) );
                // CHECK( std::abs( tIntersections( 1 ) - tIntersectionCoordinatesExpected( 0 ) ) < tEpsilon );

                // check if the point is inside and compare it to expectations
                // mtk::mtk::Mesh_Region tPointIsInside = check_if_node_is_inside_triangles( tIntersections, tTestPoint, 2, 1e-8 );

                // CHECK( tPointIsInside == INSIDE );

                // cast a bunch of random rays and ensure they all return the correct result
                mtk::Mesh_Region tPointIsInside = tSurfaceMesh.get_region_from_raycast( tTestPoint );
                REQUIRE( tPointIsInside == mtk::Mesh_Region::INSIDE );

                // repeat test for point that is outside
                tTestPoint = { { 0.2, 0.6, 0.7 } };

                // Expected results
                // tIntersectionCoordinatesExpected = { 0.715517872727636, 1.284482127272365 };

                // tCandidateTriangles = tSurfaceMesh.preselect_with_arborx( tTestPoint, tDirection );

                // REQUIRE( tCandidateTriangles.size() == 0 );

                // tCandidatesExpected = { 1, 0, 2 };
                // tDirection          = { { 1.0 }, { 0.0 }, { 0.0 } };
                // tCandidateTriangles = tSurfaceMesh.preselect_with_arborx( tTestPoint, tDirection );

                // REQUIRE( tCandidateTriangles.size() == 3 );
                // CHECK( tCandidateTriangles( 0 ) == tCandidatesExpected( 0 ) );
                // CHECK( tCandidateTriangles( 1 ) == tCandidatesExpected( 1 ) );
                // CHECK( tCandidateTriangles( 2 ) == tCandidatesExpected( 2 ) );

                // tIntersections.resize( 3 );
                // for ( uint iCandidate = 0; iCandidate < 3; ++iCandidate )
                // {
                //     tIntersections( iCandidate ) = tSurfaceMesh.moller_trumbore( tCandidateTriangles( iCandidate ), tTestPoint, tDirection );
                // }

                // CHECK( std::abs( tIntersections( 0 ) - tIntersectionCoordinatesExpected( 0 ) ) < tEpsilon );
                // CHECK( std::isnan( tIntersections( 1 ) ) );
                // CHECK( std::abs( tIntersections( 2 ) - tIntersectionCoordinatesExpected( 1 ) ) < tEpsilon );

                // Repeat for all of them at the same time using batching
                tTestPoint = { { 0.9, 0.2 }, { 0.6, 0.6 }, { 0.7, 0.7 } };

                tPointIsInside = tSurfaceMesh.get_region_from_raycast( tTestPoint );
                REQUIRE( tPointIsInside == mtk::Mesh_Region::OUTSIDE );

                // Get the regions
                Vector< mtk::Mesh_Region > tRegions  = tSurfaceMesh.batch_get_region_from_raycast( tTestPoint );
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
                std::string    tSurfaceMeshPath = tMorisRoot + "projects/GEN/SDF/test/data/rhombus.obj";
                Vector< real > tOffsets         = { 0.0, 0.0, 0.0 };
                Vector< real > tScales          = { 1.0, 1.0, 1.0 };
                Surface_Mesh   tSurfaceMesh( load_vertices_from_object_file( tSurfaceMeshPath, tOffsets, tScales ), load_facets_from_object_file( tSurfaceMeshPath ) );

                // define test point
                Matrix< DDRMat > tTestPoint = { { -.25 }, { -0.3 } };

                // define test direction
                // Matrix< DDRMat > tDirection = { { 1.0, 0.0 } };

                // // preselect in x direction and ensure the candidates and intersected facets are marked
                // Vector< uint > tCandidateLines           = tSurfaceMesh.preselect_with_arborx( tTestPoint, tDirection );
                // Vector< uint > tIntersectedLinesExpected = { 2, 3 };

                // REQUIRE( tCandidateLines.size() == 2 );
                // CHECK( tCandidateLines( 0 ) == tIntersectedLinesExpected( 0 ) );
                // CHECK( tCandidateLines( 1 ) == tIntersectedLinesExpected( 1 ) );

                // // intersect the candidate facets and determine the intersection location
                // Vector< real > tIntersections( tCandidateLines.size() );
                // for ( uint iCandidate = 0; iCandidate < tCandidateLines.size(); ++iCandidate )
                // {
                //     tIntersections( iCandidate ) = tSurfaceMesh.moller_trumbore( tCandidateLines( iCandidate ), tTestPoint, tDirection );
                // }
                // Vector< real > tIntersectionCoordinatesExpected = { 0.05, 0.45 };

                // REQUIRE( tIntersections.size() == 2 );
                // CHECK( std::abs( tIntersections( 0 ) - tIntersectionCoordinatesExpected( 0 ) ) < tEpsilon );
                // CHECK( std::abs( tIntersections( 1 ) - tIntersectionCoordinatesExpected( 1 ) ) < tEpsilon );

                mtk::Mesh_Region tRegion = tSurfaceMesh.get_region_from_raycast( tTestPoint );
                REQUIRE( tRegion == mtk::Mesh_Region::OUTSIDE );


                // repeat for a point inside the surface
                tTestPoint = { { -.25 }, { 0.2 } };

                // tDirection      = { { 0.0, 1.0 } };
                // tCandidateLines = tSurfaceMesh.preselect_with_arborx( tTestPoint, tDirection );

                // REQUIRE( tCandidateLines.size() == 1 );
                // CHECK( tCandidateLines( 0 ) == 1 );

                // tIntersections.resize( tCandidateLines.size() );
                // for ( uint iCandidate = 0; iCandidate < tCandidateLines.size(); ++iCandidate )
                // {
                //     tIntersections( iCandidate ) = tSurfaceMesh.moller_trumbore( tCandidateLines( iCandidate ), tTestPoint, tDirection );
                // }

                // tIntersectionCoordinatesExpected = { 0.05 };

                // REQUIRE( tIntersections.size() == 1 );
                // CHECK( std::abs( tIntersections( 0 ) - tIntersectionCoordinatesExpected( 0 ) ) < tEpsilon );

                tRegion = tSurfaceMesh.get_region_from_raycast( tTestPoint );
                REQUIRE( tRegion == mtk::Mesh_Region::INSIDE );

                // Repeat with a point that is on a facet
                tTestPoint( 0, 0 ) = 0.25;
                tTestPoint( 1, 0 ) = 0.25;

                tRegion = tSurfaceMesh.get_region_from_raycast( tTestPoint );
                REQUIRE( tRegion == mtk::Mesh_Region::INTERFACE );


                // Repeat with a point that is on a vertex
                tTestPoint( 0, 0 ) = 0.0;
                tTestPoint( 1, 0 ) = 0.5;

                tRegion = tSurfaceMesh.get_region_from_raycast( tTestPoint );
                REQUIRE( tRegion == mtk::Mesh_Region::INTERFACE );

                // Repeat for all of them at the same time using batching
                tTestPoint = { { -0.25, -0.25, 0.25, 0.0 }, { -0.3, 0.2, 0.25, 0.5 } };

                // Get the regions
                Vector< mtk::Mesh_Region > tRegions  = tSurfaceMesh.batch_get_region_from_raycast( tTestPoint );
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
                std::string    tSurfaceMeshPath = tMorisRoot + "projects/GEN/SDF/test/data/tetrahedron.obj";
                Vector< real > tOffsets         = { 0.0, 0.0, 0.0 };
                Vector< real > tScales          = { 1.0, 1.0, 1.0 };
                Surface_Mesh   tSurfaceMesh( load_vertices_from_object_file( tSurfaceMeshPath, tOffsets, tScales ), load_facets_from_object_file( tSurfaceMeshPath ) );

                // define test point that is inside the object
                Matrix< DDRMat > tTestPoint = { { 0.9, 0.6, 0.7 } };

                // Define the direction of the ray
                Matrix< DDRMat > tDirection = { { 1.0 }, { 0.0 }, { 0.0 } };

                real tLineDistanceXExpected = 0.384482127272365;    // facet index = 2
                real tLineDistanceYExpected = 0.319533006470609;    // facet index = 1
                real tLineDistanceZExpected = 0.180556130610548;    // facet index = 0

                // compute with raycast function
                bool                              tWarning;
                Vector< std::pair< uint, real > > tLineDistanceX = tSurfaceMesh.cast_single_ray( tTestPoint, tDirection, tWarning );
                tDirection                                       = { { 0.0 }, { 1.0 }, { 0.0 } };
                Vector< std::pair< uint, real > > tLineDistanceY = tSurfaceMesh.cast_single_ray( tTestPoint, tDirection, tWarning );
                tDirection                                       = { { 0.0 }, { 0.0 }, { 1.0 } };
                Vector< std::pair< uint, real > > tLineDistanceZ = tSurfaceMesh.cast_single_ray( tTestPoint, tDirection, tWarning );

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
                Vector< Vector< Vector< std::pair< uint, real > > > > tLineDistances = tSurfaceMesh.cast_batch_of_rays( tOrigins, tDirections, tWarnings );

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
                tLineDistances = tSurfaceMesh.cast_batch_of_rays( tOrigins, tSameDirections, tWarningSame );

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
                std::string    tSurfaceMeshPath = tMorisRoot + "projects/GEN/SDF/test/data/rhombus.obj";
                Vector< real > tOffsets         = { 0.0, 0.0 };
                Vector< real > tScales          = { 1.0, 1.0 };
                Surface_Mesh   tSurfaceMesh( load_vertices_from_object_file( tSurfaceMeshPath, tOffsets, tScales ), load_facets_from_object_file( tSurfaceMeshPath ) );

                // define test point
                Matrix< DDRMat > tTestPoint = { { -.25 }, { -0.3 } };

                // Define the direction of the ray
                Matrix< DDRMat > tDirection = { { 1.0 }, { 0.0 } };

                // expected results
                Vector< real > tLineDistanceXExpected = { 0.05, 0.45 };
                Vector< real > tLineDistanceYExpected = { 0.05, 0.55 };

                // compute with raycast
                bool                              tWarning;
                Vector< std::pair< uint, real > > tLineDistanceX = tSurfaceMesh.cast_single_ray( tTestPoint, tDirection, tWarning );
                tDirection                                       = { { 0.0 }, { 1.0 } };
                Vector< std::pair< uint, real > > tLineDistanceY = tSurfaceMesh.cast_single_ray( tTestPoint, tDirection, tWarning );

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
}    // namespace moris::mtk
