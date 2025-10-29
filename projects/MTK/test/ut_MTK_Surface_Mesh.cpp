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
    static std::unique_ptr< Kokkos::ScopeGuard > guard = !Kokkos::is_initialized() && !Kokkos::is_finalized() ? std::make_unique< Kokkos::ScopeGuard >() : nullptr;
#endif

    // get root from environment
    static std::string tMorisRoot = moris::get_base_moris_dir();

    // Finite difference epsilon
    static real tEps = 1e-7;

    // Tests loading in the surface mesh from an obj file, and computing various quanitities of interest and their sensitivities
    TEST_CASE( "MTK Surface Mesh", "[MTK],[MTK_Surface_Mesh]" )
    {
        // Expected values
        Matrix< DDRMat >                tCoordsExpected        = { { 2.25, 0.25, 1.0, 1.25 }, { 1.25, 0.5, -0.25, 0.5 } };
        Vector< Vector< moris_index > > tConnExpected          = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 } };
        Matrix< DDRMat >                tNormalsExpected       = { { -0.35112344, -0.70710678, 0.94868330, 0.60000000 }, { 0.93632918, -0.70710678, -0.31622777, -0.80000000 } };
        Matrix< DDRMat >                tCenterExpected        = { { 1.25000000, 0.62500000, 1.12500000, 1.75000000 }, { 0.87500000, 0.12500000, 0.12500000, 0.87500000 } };
        Vector< real >                  tMeasureExpected       = { 2.13600094, 1.06066017, 0.79056942, 1.25000000 };
        real                            tShapeDiameterExpected = 0.525359746;

        // load a surface mesh from file
        std::string tFilePath = tMorisRoot + "/projects/GEN/test/data/triangle_sensitivity_oblique.obj";
        // std::string    tFilePath = "/home/chong/work/AU25/Input_Files/QI_Test/square_corners_only.obj";
        Vector< real > tOffsets = { 0.0, 0.0 };
        Vector< real > tScales  = { 1.0, 1.0 };
        Surface_Mesh   tSurfaceMesh( load_vertices_from_object_file( tFilePath, tOffsets, tScales ), load_facets_from_object_file( tFilePath ) );

        // Config for shape diameter computation
        uint tNumRays   = 120;
        real tConeAngle = 15.0;    // degrees

        // Check number of vertices and facets
        REQUIRE( tSurfaceMesh.get_number_of_vertices() == tCoordsExpected.n_cols() );
        REQUIRE( tSurfaceMesh.get_number_of_facets() == tConnExpected.size() );

        // Check vertex coordinates
        check_equal( tSurfaceMesh.get_all_vertex_coordinates(), tCoordsExpected );

        // Compute the actual facet measure and centers
        Vector< real >   tFacetMeasures = tSurfaceMesh.compute_facet_measure();
        Matrix< DDRMat > tFacetCenters  = tSurfaceMesh.compute_facet_centroids();

        // Compute the nodal and global shape diameter
        Vector< real > tNodalShapeDiameter  = tSurfaceMesh.compute_nodal_shape_diameter( tConeAngle, tNumRays );
        real           tGlobalShapeDiameter = tSurfaceMesh.compute_global_shape_diameter( tConeAngle, tNumRays );

        // Nodal shape diameter sensitivities
        const Matrix< DDRMat >& tNodalShapeDiameterSensitivities = tSurfaceMesh.get_nodal_shape_diameter_sensitivities();

        // Check global shape diameter
        CHECK( tGlobalShapeDiameter == Approx( tShapeDiameterExpected ) );

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

                // Check facet center component
                CHECK( tFacetCenters( iV, iF ) == Approx( tCenterExpected( iV, iF ) ) );
            }

            // Check facet measure
            CHECK( tFacetMeasures( iF ) == Approx( tMeasureExpected( iF ) ) );

            // Setup perturbation matrix for FD
            Matrix< DDRMat > tPerturbation( 2, 1, 0.0 );

            // Loop over vertices and finite difference the quanitities
            for ( uint iV = 0; iV < tSurfaceMesh.get_number_of_vertices(); iV++ )
            {
                // Compute analytic sensitivities
                Matrix< DDRMat > tNormalSens        = tSurfaceMesh.compute_dfacet_normal_dvertex( iF, iV );
                Matrix< DDRMat > tCenterSens        = tSurfaceMesh.compute_dfacet_centroid_dvertex( iF, iV );
                Matrix< DDRMat > tMeasureSens       = tSurfaceMesh.compute_dfacet_measure_dvertex( iF, iV );
                Matrix< DDRMat > tShapeDiameterSens = tSurfaceMesh.compute_ddiameter_dvertex( iV );

                // Loop over dimensions
                for ( uint iDim = 0; iDim < 2; iDim++ )
                {
                    // Perturb positively
                    tPerturbation( iDim ) = tEps;
                    tSurfaceMesh.set_vertex_displacement( iV, tPerturbation );
                    Matrix< DDRMat > tNormalPlus             = tSurfaceMesh.get_facet_normal( iF );
                    Matrix< DDRMat > tCenterPlus             = tSurfaceMesh.compute_facet_centroid( iF );
                    real             tMeasurePlus            = tSurfaceMesh.compute_facet_measure( iF );
                    real             tShapeDiameterPlus      = tSurfaceMesh.compute_global_shape_diameter( tConeAngle, tNumRays );
                    Vector< real >   tNodalShapeDiameterPlus = tSurfaceMesh.compute_nodal_shape_diameter( tConeAngle, tNumRays );

                    // Perturb negatively
                    tPerturbation( iDim ) = -tEps;
                    tSurfaceMesh.set_vertex_displacement( iV, tPerturbation );
                    Matrix< DDRMat > tNormalMinus             = tSurfaceMesh.get_facet_normal( iF );
                    Matrix< DDRMat > tCenterMinus             = tSurfaceMesh.compute_facet_centroid( iF );
                    real             tMeasureMinus            = tSurfaceMesh.compute_facet_measure( iF );
                    real             tShapeDiameterMinus      = tSurfaceMesh.compute_global_shape_diameter( tConeAngle, tNumRays );
                    Vector< real >   tNodalShapeDiameterMinus = tSurfaceMesh.compute_nodal_shape_diameter( tConeAngle, tNumRays );

                    // Reset perturbation
                    tPerturbation( iDim ) = 0.0;
                    tSurfaceMesh.set_vertex_displacement( iV, tPerturbation );

                    // Compute finite difference results - nodal shape diameter
                    real tNSDForward  = 0.0;
                    real tNSDBackward = 0.0;
                    real tNSDCentral  = 0.0;
                    for ( uint iF = 0; iF < tSurfaceMesh.get_number_of_facets(); iF++ )
                    {
                        tNSDForward += tFacetMeasures( iF ) * ( tNodalShapeDiameterPlus( iF ) - tNodalShapeDiameter( iF ) ) / tEps;
                        tNSDBackward += tFacetMeasures( iF ) * ( tNodalShapeDiameter( iF ) - tNodalShapeDiameterMinus( iF ) ) / tEps;
                        tNSDCentral += tFacetMeasures( iF ) * ( tNodalShapeDiameterPlus( iF ) - tNodalShapeDiameterMinus( iF ) ) / ( 2.0 * tEps );
                    }


                    // Compute finite difference results - normal
                    Matrix< DDRMat > tNForward  = ( tNormalPlus - tFacetNormal ) / tEps;
                    Matrix< DDRMat > tNBackward = ( tFacetNormal - tNormalMinus ) / tEps;
                    Matrix< DDRMat > tNCentral  = ( tNormalPlus - tNormalMinus ) / ( 2.0 * tEps );

                    // Compute finite difference results - center
                    Matrix< DDRMat > tCForward  = ( tCenterPlus - tFacetCenters.get_column( iF ) ) / tEps;
                    Matrix< DDRMat > tCBackward = ( tFacetCenters.get_column( iF ) - tCenterMinus ) / tEps;
                    Matrix< DDRMat > tCCentral  = ( tCenterPlus - tCenterMinus ) / ( 2.0 * tEps );

                    // Compute finite difference results - measure
                    real tMForward  = ( tMeasurePlus - tFacetMeasures( iF ) ) / tEps;
                    real tMBackward = ( tFacetMeasures( iF ) - tMeasureMinus ) / tEps;
                    real tMCentral  = ( tMeasurePlus - tMeasureMinus ) / ( 2.0 * tEps );

                    // Check sensitivities for measure
                    CHECK( tMeasureSens( iDim ) == Approx( tMForward ) );
                    CHECK( tMeasureSens( iDim ) == Approx( tMBackward ) );
                    CHECK( tMeasureSens( iDim ) == Approx( tMCentral ) );


                    // Check sensitivities for normal and center
                    for ( uint iComp = 0; iComp < 2; iComp++ )
                    {
                        // Check sensitivities for normal
                        CHECK( tNormalSens( iComp, iDim ) == Approx( tNForward( iComp ) ) );
                        CHECK( tNormalSens( iComp, iDim ) == Approx( tNBackward( iComp ) ) );
                        CHECK( tNormalSens( iComp, iDim ) == Approx( tNCentral( iComp ) ) );

                        // Check sensitivities for center
                        CHECK( tCenterSens( iComp, iDim ) == Approx( tCForward( iComp ) ) );
                        CHECK( tCenterSens( iComp, iDim ) == Approx( tCBackward( iComp ) ) );
                        CHECK( tCenterSens( iComp, iDim ) == Approx( tCCentral( iComp ) ) );
                    }

                    // Check nodal shape diameter sensitivities
                    CHECK( tNodalShapeDiameterSensitivities( iV, iDim ) == Approx( tNSDForward ).epsilon( 1e-2 ) );
                    CHECK( tNodalShapeDiameterSensitivities( iV, iDim ) == Approx( tNSDBackward ).epsilon( 1e-2 ) );
                    CHECK( tNodalShapeDiameterSensitivities( iV, iDim ) == Approx( tNSDCentral ).epsilon( 1e-2 ) );

                    // Compute finite difference results - global shape diameter
                    real tSDForward  = ( tShapeDiameterPlus - tGlobalShapeDiameter ) / tEps;
                    real tSDBackward = ( tGlobalShapeDiameter - tShapeDiameterMinus ) / tEps;
                    real tSDCentral  = ( tShapeDiameterPlus - tShapeDiameterMinus ) / ( 2.0 * tEps );

                    // Check sensitivities shape diameter
                    CHECK( tShapeDiameterSens( iDim ) == Approx( tSDForward ).epsilon( 1e-2 ) );
                    CHECK( tShapeDiameterSens( iDim ) == Approx( tSDBackward ).epsilon( 1e-2 ) );
                    CHECK( tShapeDiameterSens( iDim ) == Approx( tSDCentral ).epsilon( 1e-2 ) );
                }
            }
        }
        SECTION( "Raycast Region, Distance, and Sensitivities - 2D" )
        {
            // define test point and direction
            Matrix< DDRMat > tTestPoint = { { 0.8 }, { 0.4 } };
            Matrix< DDRMat > tDirection = { { -1.0 }, { -1.0 } };
            tDirection                  = tDirection;

            // Setup perturbation matrix for FD
            Matrix< DDRMat > tPerturbation( 2, 1, 0.0 );

            // Check the region of this point
            mtk::Mesh_Region tRegion = tSurfaceMesh.get_region_from_raycast( tTestPoint );
            REQUIRE( tRegion == mtk::Mesh_Region::INSIDE );

            // Get the distance and sensitivities
            bool                tWarning;
            Intersection_Vector tDistance    = tSurfaceMesh.cast_single_ray( tTestPoint, tDirection, tWarning );
            real                tRefDistance = tDistance( 0 ).second;
            REQUIRE( tDistance.size() == 1 );
            CHECK( tRefDistance == Approx( 0.2250 ) );    // distance
            CHECK( tDistance( 0 ).first == 1 );           // intersected facet
            CHECK( not tWarning );                        // no warning

            // Get the sensitivity wrt to the ray origin
            Matrix< DDRMat > tdRdO = tSurfaceMesh.compute_draycast_dorigin( tTestPoint, tDirection, tDistance( 0 ).first );

            // Finite difference raycast sensitivity wrt origin
            // Loop over dimensions
            for ( uint iDim = 0; iDim < 2; iDim++ )
            {
                // Perturb origin positively
                tPerturbation( iDim )                = tEps;
                Matrix< DDRMat >    tPerturbedOrigin = tTestPoint + tPerturbation;
                Intersection_Vector tDistancePlus    = tSurfaceMesh.cast_single_ray( tPerturbedOrigin, tDirection, tWarning );
                CHECK( tDistancePlus.size() == 1 );
                REQUIRE( tDistancePlus( 0 ).first == 1 );    // intersected facet
                REQUIRE( not tWarning );                     // no warning
                real tDistancePlusValue = tDistancePlus( 0 ).second;

                // Perturb origin negatively
                tPerturbation( iDim )              = -tEps;
                tPerturbedOrigin                   = tTestPoint + tPerturbation;
                Intersection_Vector tDistanceMinus = tSurfaceMesh.cast_single_ray( tPerturbedOrigin, tDirection, tWarning );
                CHECK( tDistancePlus.size() == 1 );
                REQUIRE( tDistancePlus( 0 ).first == 1 );    // intersected facet
                REQUIRE( not tWarning );                     // no warning
                real tDistanceMinusValue = tDistanceMinus( 0 ).second;

                // Reset perturbation
                tPerturbation( iDim ) = 0.0;

                // Compute FD sensitivity
                real tRForward  = ( tDistancePlusValue - tRefDistance ) / tEps;
                real tRBackward = ( tRefDistance - tDistanceMinusValue ) / tEps;
                real tRCentral  = ( tDistancePlusValue - tDistanceMinusValue ) / ( 2.0 * tEps );

                // Check sensitivities
                CHECK( tdRdO( iDim ) == Approx( tRForward ) );
                CHECK( tdRdO( iDim ) == Approx( tRBackward ) );
                CHECK( tdRdO( iDim ) == Approx( tRCentral ) );
            }

            // Get the sensitivity wrt to the ray direction
            Matrix< DDRMat > tdRdD = tSurfaceMesh.compute_draycast_ddirection( tTestPoint, tDirection, tDistance( 0 ).first );

            // Finite difference raycast sensitivity wrt direction
            // Loop over dimensions
            for ( uint iDim = 0; iDim < 2; iDim++ )
            {
                // Perturb direction positively
                tPerturbation( iDim )                   = tEps;
                Matrix< DDRMat >    tPerturbedDirection = tDirection + tPerturbation;
                Intersection_Vector tDistancePlus       = tSurfaceMesh.cast_single_ray( tTestPoint, tPerturbedDirection, tWarning );
                real                tDistancePlusValue  = tDistancePlus( 0 ).second;

                // Perturb direction negatively
                tPerturbation( iDim )                      = -tEps;
                Matrix< DDRMat >    tPerturbedDirectionNeg = tDirection + tPerturbation;
                Intersection_Vector tDistanceMinus         = tSurfaceMesh.cast_single_ray( tTestPoint, tPerturbedDirectionNeg, tWarning );
                real                tDistanceMinusValue    = tDistanceMinus( 0 ).second;

                // Reset perturbation
                tPerturbation( iDim ) = 0.0;

                // Compute FD sensitivity
                real tRForward  = ( tDistancePlusValue - tRefDistance ) / tEps;
                real tRBackward = ( tRefDistance - tDistanceMinusValue ) / tEps;
                real tRCentral  = ( tDistancePlusValue - tDistanceMinusValue ) / ( 2.0 * tEps );

                // Check sensitivities
                CHECK( tdRdD( iDim ) == Approx( tRForward ) );
                CHECK( tdRdD( iDim ) == Approx( tRBackward ) );
                CHECK( tdRdD( iDim ) == Approx( tRCentral ) );
            }

            // Get the sensitivity wrt to the vertices of the intersected facet
            Matrix< DDRMat > tdRdV = tSurfaceMesh.compute_draycast_dvertices( tTestPoint, tDirection, tDistance( 0 ).first );

            // Finite difference raycast sensitivity wrt vertices
            Vector< moris_index > tFacetVertices = tSurfaceMesh.get_facets_vertex_indices( tDistance( 0 ).first );
            // Loop over vertices
            for ( uint iVLocal = 0; iVLocal < tFacetVertices.size(); iVLocal++ )
            {
                moris_index iV = tFacetVertices( iVLocal );
                // Loop over dimensions
                for ( uint iDim = 0; iDim < 2; iDim++ )
                {
                    // Perturb vertex positively
                    tPerturbation( iDim ) = tEps;
                    tSurfaceMesh.set_vertex_displacement( iV, tPerturbation );
                    Intersection_Vector tDistancePlus    = tSurfaceMesh.cast_single_ray( tTestPoint, tDirection, tWarning );
                    real                tDistancePlusVal = tDistancePlus( 0 ).second;

                    // Perturb vertex negatively
                    tPerturbation( iDim ) = -tEps;
                    tSurfaceMesh.set_vertex_displacement( iV, tPerturbation );
                    Intersection_Vector tDistanceMinus    = tSurfaceMesh.cast_single_ray( tTestPoint, tDirection, tWarning );
                    real                tDistanceMinusVal = tDistanceMinus( 0 ).second;

                    // Reset perturbation
                    tPerturbation( iDim ) = 0.0;
                    tSurfaceMesh.set_vertex_displacement( iV, tPerturbation );

                    // Compute FD sensitivity
                    real tRForward  = ( tDistancePlusVal - tRefDistance ) / tEps;
                    real tRBackward = ( tRefDistance - tDistanceMinusVal ) / tEps;
                    real tRCentral  = ( tDistancePlusVal - tDistanceMinusVal ) / ( 2.0 * tEps );

                    // Check sensitivities
                    CHECK( tdRdV( iVLocal, iDim ) == Approx( tRForward ) );
                    CHECK( tdRdV( iVLocal, iDim ) == Approx( tRBackward ) );
                    CHECK( tdRdV( iVLocal, iDim ) == Approx( tRCentral ) );
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
                Vector< real > tOffsets  = { 0.0, 0.0, 0.0 };
                Vector< real > tScales   = { 1.0, 1.0, 1.0 };
                Surface_Mesh   tSurfaceMesh( load_vertices_from_object_file( tFilePath, tOffsets, tScales ), load_facets_from_object_file( tFilePath ) );

                // define test point that is inside the object
                Matrix< DDRMat > tTestPoint = { { 0.9, 0.6, 0.7 } };

                // cast a bunch of random rays and ensure they all return the correct result
                mtk::Mesh_Region tPointIsInside = tSurfaceMesh.get_region_from_raycast( tTestPoint );
                REQUIRE( tPointIsInside == mtk::Mesh_Region::INSIDE );

                // repeat test for point that is outside
                tTestPoint     = { { 0.2, 0.6, 0.7 } };
                tPointIsInside = tSurfaceMesh.get_region_from_raycast( tTestPoint );
                REQUIRE( tPointIsInside == mtk::Mesh_Region::OUTSIDE );

                // Repeat for all of them at the same time using batching
                tTestPoint = { { 0.9, 0.2 }, { 0.6, 0.6 }, { 0.7, 0.7 } };

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
                // create object from object file
                std::string    tSurfaceMeshPath = tMorisRoot + "projects/GEN/SDF/test/data/rhombus.obj";
                Vector< real > tOffsets         = { 0.0, 0.0, 0.0 };
                Vector< real > tScales          = { 1.0, 1.0, 1.0 };
                Surface_Mesh   tSurfaceMesh( load_vertices_from_object_file( tSurfaceMeshPath, tOffsets, tScales ), load_facets_from_object_file( tSurfaceMeshPath ) );

                // define test point
                Matrix< DDRMat > tTestPoint = { { -.25 }, { -0.3 } };

                mtk::Mesh_Region tRegion = tSurfaceMesh.get_region_from_raycast( tTestPoint );
                REQUIRE( tRegion == mtk::Mesh_Region::OUTSIDE );

                // repeat for a point inside the surface
                tTestPoint = { { -.25 }, { 0.2 } };

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
