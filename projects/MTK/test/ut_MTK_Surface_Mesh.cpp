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
        Matrix< DDRMat >                tCoordsExpected  = { { 2.25, 0.25, 1.0, 1.25 }, { 1.25, 0.5, -0.25, 0.5 } };
        Vector< Vector< moris_index > > tConnExpected    = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 } };
        Matrix< DDRMat >                tNormalsExpected = { { -0.35112344, -0.70710678, 0.94868330, 0.60000000 }, { 0.93632918, -0.70710678, -0.31622777, -0.80000000 } };
        Matrix< DDRMat >                tCenterExpected  = { { 1.25000000, 0.62500000, 1.12500000, 1.75000000 }, { 0.87500000, 0.12500000, 0.12500000, 0.87500000 } };
        Vector< real >                  tMeasureExpected = { 2.13600094, 1.06066017, 0.79056942, 1.25000000 };

        // load a surface mesh from file
        std::string    tFilePath = tMorisRoot + "/projects/GEN/test/data/triangle_sensitivity_oblique.obj";
        Vector< real > tOffsets  = { 0.0, 0.0 };
        Vector< real > tScales   = { 1.0, 1.0 };
        Surface_Mesh   tSurfaceMesh( load_vertices_from_object_file( tFilePath, tOffsets, tScales ), load_facets_from_object_file( tFilePath ) );

        SECTION( "Basic functionality and sensitivities" )
        {
            real tRaycastDiameterExpected   = 2.220283351;
            real tInscribedDiameterExpected = 0.1986737884;
            real tShortestDiameterExpected  = 0.12377533;

            // Config for shape diameter computation
            uint                              tNumRays            = 6;       // For raycast method
            uint                              tNumSamples         = 1;       // For inscribed circle and shortest distance methods
            real                              tRaycastConeAngle   = 60.0;    // degrees
            real                              tRelativeChord      = 0.25;
            real                              tInscribedConeAngle = 120.0;    // degrees
            real                              tShortestConeAngle  = 40.0;     // degrees
            Agglomeration_Parameters_Constant tAgglom( 4.0, 2.0, 0.0 );       // Same for all methods

            // Check number of vertices and facets
            REQUIRE( tSurfaceMesh.get_number_of_vertices() == tCoordsExpected.n_cols() );
            REQUIRE( tSurfaceMesh.get_number_of_facets() == tConnExpected.size() );

            // Check vertex coordinates
            check_equal( tSurfaceMesh.get_all_vertex_coordinates(), tCoordsExpected );

            // Compute the actual facet measure and centers
            Vector< real >   tFacetMeasures = tSurfaceMesh.compute_facet_measure();
            Matrix< DDRMat > tFacetCenters  = tSurfaceMesh.compute_facet_centroids();

            // Compute the nodal and global shape diameter
            Vector< real > tNodalShapeDiameter = tSurfaceMesh.compute_raycast_shape_diameter( tAgglom, tRaycastConeAngle, tNumRays, 1 );

            // Setup perturbation matrix for FD
            Matrix< DDRMat > tPerturbation( 2, 1, 0.0 );

            // Loop over the surface mesh facets
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

                // Loop over vertices and finite difference the quanitities
                for ( uint iV = 0; iV < tSurfaceMesh.get_number_of_vertices(); iV++ )
                {
                    // Compute analytic sensitivities
                    Matrix< DDRMat > tNormalSens  = tSurfaceMesh.compute_dfacet_normal_dvertex( iF, iV );
                    Matrix< DDRMat > tCenterSens  = tSurfaceMesh.compute_dfacet_centroid_dvertex( iF, iV );
                    Matrix< DDRMat > tMeasureSens = tSurfaceMesh.compute_dfacet_measure_dvertex( iF, iV );

                    // Compute different shape diameters and their sensitivities. Have to recompute the forward solve as it stores intermediate values that are needed for the sensitivity computations. Yes this is expensive, but it's just a test to verify the sensitivities are correct.
                    real             tRaycastGlobalDiameter   = tSurfaceMesh.compute_global_shape_diameter_raycast( tAgglom, tRaycastConeAngle, tNumRays, (uint)1 );
                    Matrix< DDRMat > tRaycastDiameterSens     = tSurfaceMesh.compute_ddiameter_dvertex( iV );
                    real             tInscribedGlobalDiameter = tSurfaceMesh.compute_global_shape_diameter_inscribed_circle( tAgglom, tInscribedConeAngle, tRelativeChord, tNumSamples );
                    Matrix< DDRMat > tInscribedDiameterSens   = tSurfaceMesh.compute_ddiameter_dvertex( iV );
                    real             tShortestGlobalDiameter  = tSurfaceMesh.compute_global_shape_diameter_shortest_distance( tAgglom, tShortestConeAngle, tNumSamples );
                    Matrix< DDRMat > tShortestDiameterSens    = tSurfaceMesh.compute_ddiameter_dvertex( iV );

                    // Check global shape diameter
                    CHECK( tRaycastGlobalDiameter == Approx( tRaycastDiameterExpected ) );
                    CHECK( tInscribedGlobalDiameter == Approx( tInscribedDiameterExpected ) );
                    CHECK( tShortestGlobalDiameter == Approx( tShortestDiameterExpected ) );

                    // Loop over dimensions
                    for ( uint iDim = 0; iDim < 2; iDim++ )
                    {
                        // Perturb positively
                        tPerturbation( iDim ) = tEps;
                        tSurfaceMesh.set_vertex_displacement( iV, tPerturbation );
                        Matrix< DDRMat > tNormalPlus    = tSurfaceMesh.get_facet_normal( iF );
                        Matrix< DDRMat > tCenterPlus    = tSurfaceMesh.compute_facet_centroid( iF );
                        real             tMeasurePlus   = tSurfaceMesh.compute_facet_measure( iF );
                        real             tRaycastPlus   = tSurfaceMesh.compute_global_shape_diameter_raycast( tAgglom, tRaycastConeAngle, tNumRays, 1 );
                        real             tInscribedPlus = tSurfaceMesh.compute_global_shape_diameter_inscribed_circle( tAgglom, tInscribedConeAngle, tRelativeChord, tNumSamples );
                        real             tShortestPlus  = tSurfaceMesh.compute_global_shape_diameter_shortest_distance( tAgglom, tShortestConeAngle, tNumSamples );

                        // Perturb negatively
                        tPerturbation( iDim ) = -tEps;
                        tSurfaceMesh.set_vertex_displacement( iV, tPerturbation );
                        Matrix< DDRMat > tNormalMinus    = tSurfaceMesh.get_facet_normal( iF );
                        Matrix< DDRMat > tCenterMinus    = tSurfaceMesh.compute_facet_centroid( iF );
                        real             tMeasureMinus   = tSurfaceMesh.compute_facet_measure( iF );
                        real             tRaycastMinus   = tSurfaceMesh.compute_global_shape_diameter_raycast( tAgglom, tRaycastConeAngle, tNumRays, 1 );
                        real             tInscribedMinus = tSurfaceMesh.compute_global_shape_diameter_inscribed_circle( tAgglom, tInscribedConeAngle, tRelativeChord, tNumSamples );
                        real             tShortestMinus  = tSurfaceMesh.compute_global_shape_diameter_shortest_distance( tAgglom, tShortestConeAngle, tNumSamples );

                        // Reset perturbation
                        tPerturbation( iDim ) = 0.0;
                        tSurfaceMesh.set_vertex_displacement( iV, tPerturbation );

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

                        // Compute finite difference results - global shape diameter
                        real tSDForward  = ( tRaycastPlus - tRaycastGlobalDiameter ) / tEps;
                        real tSDBackward = ( tRaycastGlobalDiameter - tRaycastMinus ) / tEps;
                        real tSDCentral  = ( tRaycastPlus - tRaycastMinus ) / ( 2.0 * tEps );

                        // Check sensitivities - raycast shape diameter
                        CHECK( tRaycastDiameterSens( iDim ) == Approx( tSDForward ).epsilon( 1e-4 ) );
                        CHECK( tRaycastDiameterSens( iDim ) == Approx( tSDBackward ).epsilon( 1e-4 ) );
                        CHECK( tRaycastDiameterSens( iDim ) == Approx( tSDCentral ).epsilon( 1e-4 ) );

                        // Compute finite difference results - inscribed circle shape diameter
                        tSDForward  = ( tInscribedPlus - tInscribedGlobalDiameter ) / tEps;
                        tSDBackward = ( tInscribedGlobalDiameter - tInscribedMinus ) / tEps;
                        tSDCentral  = ( tInscribedPlus - tInscribedMinus ) / ( 2.0 * tEps );

                        // Check sensitivities - inscribed circle shape diameter
                        CHECK( tInscribedDiameterSens( iDim ) == Approx( tSDForward ).epsilon( 1e-4 ) );
                        CHECK( tInscribedDiameterSens( iDim ) == Approx( tSDBackward ).epsilon( 1e-4 ) );
                        CHECK( tInscribedDiameterSens( iDim ) == Approx( tSDCentral ).epsilon( 1e-4 ) );

                        // Compute finite difference results - shortest distance shape diameter
                        tSDForward  = ( tShortestPlus - tShortestGlobalDiameter ) / tEps;
                        tSDBackward = ( tShortestGlobalDiameter - tShortestMinus ) / tEps;
                        tSDCentral  = ( tShortestPlus - tShortestMinus ) / ( 2.0 * tEps );

                        // Check sensitivities - shortest distance shape diameter
                        CHECK( tShortestDiameterSens( iDim ) == Approx( tSDForward ).epsilon( 1e-6 ) );
                        CHECK( tShortestDiameterSens( iDim ) == Approx( tSDBackward ).epsilon( 1e-6 ) );
                        CHECK( tShortestDiameterSens( iDim ) == Approx( tSDCentral ).epsilon( 1e-6 ) );
                    }
                }
            }

            // Compute vertex normals and FD sensitivities
            Matrix< DDRMat > tVertexNormals = tSurfaceMesh.compute_vertex_normals();

            // Loop over the surface mesh vertices
            for ( uint iV = 0; iV < tSurfaceMesh.get_number_of_vertices(); iV++ )
            {
                // Loop over vertices again
                for ( uint iVC = 0; iVC < tSurfaceMesh.get_number_of_vertices(); iVC++ )
                {
                    // Get the vertex normal sensitivity
                    Matrix< DDRMat > tVertexNormalSens = tSurfaceMesh.compute_dvertex_normal_dvertex( iV, iVC );

                    // Loop over dimensions
                    for ( uint iDim = 0; iDim < 2; iDim++ )
                    {
                        // Perturb positively
                        tPerturbation( iDim ) = tEps;
                        tSurfaceMesh.set_vertex_displacement( iVC, tPerturbation );
                        Matrix< DDRMat > tVertexNormalsPlus = tSurfaceMesh.compute_vertex_normals();

                        // Perturb negatively
                        tPerturbation( iDim ) = -tEps;
                        tSurfaceMesh.set_vertex_displacement( iVC, tPerturbation );
                        Matrix< DDRMat > tVertexNormalsMinus = tSurfaceMesh.compute_vertex_normals();

                        // Reset perturbation
                        tPerturbation( iDim ) = 0.0;
                        tSurfaceMesh.set_vertex_displacement( iVC, tPerturbation );

                        // Compute FD sensitivity
                        Matrix< DDRMat > tVNSForward  = ( tVertexNormalsPlus.get_column( iV ) - tVertexNormals.get_column( iV ) ) / tEps;
                        Matrix< DDRMat > tVNSBackward = ( tVertexNormals.get_column( iV ) - tVertexNormalsMinus.get_column( iV ) ) / tEps;
                        Matrix< DDRMat > tVNSCentral  = ( tVertexNormalsPlus.get_column( iV ) - tVertexNormalsMinus.get_column( iV ) ) / ( 2.0 * tEps );

                        // Check sensitivities for normal
                        for ( uint iComp = 0; iComp < 2; iComp++ )
                        {
                            if ( std::abs( tVNSCentral( iComp ) ) > 1e-8 )
                            {
                                CHECK( tVertexNormalSens( iComp, iDim ) == Approx( tVNSForward( iComp ) ).epsilon( 1e-6 ) );
                                CHECK( tVertexNormalSens( iComp, iDim ) == Approx( tVNSBackward( iComp ) ).epsilon( 1e-6 ) );
                                CHECK( tVertexNormalSens( iComp, iDim ) == Approx( tVNSCentral( iComp ) ).epsilon( 1e-6 ) );
                            }
                        }
                    }
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
        SECTION( "Dot Product Sensitivity" )
        {
            // These can be any arbitrary vectors for testing purposes
            Matrix< DDRMat > tDirection   = { { 0.70710678118 }, { 0.70710678118 } };
            Matrix< DDRMat > tFacetNormal = { { 0.0 }, { 1.0 } };

            real tAngle = tSurfaceMesh.compute_dot_absolute( tDirection, tFacetNormal );
            CHECK( tAngle == Approx( 0.7071067812 ) );

            // Finite difference sensitivity
            Matrix< DDRMat > tPerturbation( 2, 1, 0.0 );
            Matrix< DDRMat > tDAdD = tSurfaceMesh.compute_ddot_absolute( tDirection, tFacetNormal );

            // Loop over dimensions
            for ( uint iDim = 0; iDim < 2; iDim++ )
            {
                // Perturb direction positively
                tPerturbation( iDim )          = tEps;
                Matrix< DDRMat > tPerturbedDir = tDirection + tPerturbation;
                real             tAnglePlus    = tSurfaceMesh.compute_dot_absolute( tPerturbedDir, tFacetNormal );

                // Perturb direction negatively
                tPerturbation( iDim )             = -tEps;
                Matrix< DDRMat > tPerturbedDirNeg = tDirection + tPerturbation;
                real             tAngleMinus      = tSurfaceMesh.compute_dot_absolute( tPerturbedDirNeg, tFacetNormal );

                // Reset perturbation
                tPerturbation( iDim ) = 0.0;

                // Compute FD sensitivity
                real tAForward  = ( tAnglePlus - tAngle ) / tEps;
                real tABackward = ( tAngle - tAngleMinus ) / tEps;
                real tACentral  = ( tAnglePlus - tAngleMinus ) / ( 2.0 * tEps );

                // Check sensitivities
                CHECK( tDAdD( iDim ) == Approx( tAForward ) );
                CHECK( tDAdD( iDim ) == Approx( tABackward ) );
                CHECK( tDAdD( iDim ) == Approx( tACentral ) );
            }
        }
        SECTION( "Agglomeration Function Sensitivity" )
        {
            // Agglomeration function sensitivity test
            Agglomeration_Parameters_Constant tAgglom( 2.0, 0.5, 0.0 );
            uint                              tNumTests = 5;
            Matrix< DDRMat >                  tDummy;

            // Loop over some test values
            for ( uint iTest = 0; iTest < tNumTests; iTest++ )
            {
                real tValue = 0.05 + iTest * 0.05;

                // Compute analytic value and sensitivity
                real tAgglomeratedValue        = tSurfaceMesh.tanh_clip( tValue, tAgglom, tDummy );
                real tDAgglomeratedValueDValue = tSurfaceMesh.dtanh_clip( tValue, tAgglom, tDummy );

                // Finite difference sensitivity
                real tValuePlus             = tValue + tEps;
                real tAgglomeratedValuePlus = tSurfaceMesh.tanh_clip( tValuePlus, tAgglom, tDummy );

                real tValueMinus             = tValue - tEps;
                real tAgglomeratedValueMinus = tSurfaceMesh.tanh_clip( tValueMinus, tAgglom, tDummy );

                real tDForward  = ( tAgglomeratedValuePlus - tAgglomeratedValue ) / tEps;
                real tDBackward = ( tAgglomeratedValue - tAgglomeratedValueMinus ) / tEps;
                real tDCentral  = ( tAgglomeratedValuePlus - tAgglomeratedValueMinus ) / ( 2.0 * tEps );

                // Check sensitivities
                CHECK( tDAgglomeratedValueDValue == Approx( tDForward ) );
                CHECK( tDAgglomeratedValueDValue == Approx( tDBackward ) );
                CHECK( tDAgglomeratedValueDValue == Approx( tDCentral ) );

                // Compute analytic sensitivity
                tAgglomeratedValue        = tSurfaceMesh.max_clip( tValue, tAgglom, tDummy );
                tDAgglomeratedValueDValue = tSurfaceMesh.dmax_clip( tValue, tAgglom, tDummy );

                // Finite difference sensitivity
                tValuePlus             = tValue + tEps;
                tAgglomeratedValuePlus = tSurfaceMesh.max_clip( tValuePlus, tAgglom, tDummy );

                tValueMinus             = tValue - tEps;
                tAgglomeratedValueMinus = tSurfaceMesh.max_clip( tValueMinus, tAgglom, tDummy );

                tDForward  = ( tAgglomeratedValuePlus - tAgglomeratedValue ) / tEps;
                tDBackward = ( tAgglomeratedValue - tAgglomeratedValueMinus ) / tEps;
                tDCentral  = ( tAgglomeratedValuePlus - tAgglomeratedValueMinus ) / ( 2.0 * tEps );

                // Check sensitivities
                CHECK( tDAgglomeratedValueDValue == Approx( tDForward ) );
                CHECK( tDAgglomeratedValueDValue == Approx( tDBackward ) );
                CHECK( tDAgglomeratedValueDValue == Approx( tDCentral ) );
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
    // TEST_CASE( "Sinusoid Shape Diameter Sweep - Cone Angle", "[MTK], [MTK_Surface_Mesh], [Raycast]" )    // brendan delete
    // {
    //     // real tAmp = 0.0;
    //     // real tPeriod  = 17.5;
    //     uint tHgt     = 50;
    //     real tWid     = 10.0;
    //     uint tNPoints = 90;

    //     // // create surface mesh from object file
    //     std::string tPath = "/home/chong/work/AU25/Input_Files/SD_tests/meshes/";
    //     // std::string    tFile            = "sin_" + std::to_string( static_cast< int >( 10 * tAmp ) ) + "a_" + std::to_string( static_cast< int >( tPeriod ) ) + "p_" + std::to_string( static_cast< int >( 10 * tHgt ) ) + "h_" + std::to_string( static_cast< int >( 10 * tWid ) ) + "w_" + std::to_string( static_cast< int >( tNpoints ) ) + "np";
    //     // std::string    tSurfaceMeshPath = tPath + tFile + ".obj";
    //     Vector< real > tOffsets = { 0.0, 0.0, 0.0 };
    //     Vector< real > tScales  = { 1.0, 1.0, 1.0 };
    //     // Surface_Mesh   tSurfaceMesh( load_vertices_from_object_file( tSurfaceMeshPath, tOffsets, tScales ), load_facets_from_object_file( tSurfaceMeshPath ) );

    //     // Config for shape diameter computation
    //     uint                     tNumRaysBash               = 30;
    //     real                     tConeAngleBash             = 60.0;
    //     real                     tAgglomerationExponentBash = 2.0;
    //     real                     tAgglomerationRefBash      = 1.0;
    //     real                     tAgglomerationShiftBash    = 0.0;
    //     Agglomeration_Parameters_Constant tAgglom( tAgglomerationExponentBash, tAgglomerationRefBash, tAgglomerationShiftBash );

    //     // Open file for printing
    //     std::ofstream tSDFile;
    //     tSDFile.open( tPath + "../data/apsweepsfine_maxsmoother.txt" );
    //     tSDFile << "Amp Period Height Width Num_points Agglomeration_Exponent Agglomeration_Reference Agglomeration_Shift Num_rays Cone_Angle Shape_Diameter dSD_dVx dSD_dVy" << std::endl;

    //     // for ( uint tCase = 0; tCase < tHgtValues.size(); tCase++ )
    //     // {
    //     //     real tHgt = tHgtValues( tCase );
    //     //     real tWid = tWidValues( tCase );
    //     // std::cout << "{ ";

    //     // for ( uint tHgt = 2; tHgt <= 50; tHgt++ )
    //     // {
    //     //     for ( real tWid = 0.1; tWid < 10.0; tWid += 0.1 )
    //     //     {
    //     for ( uint tPeriod = 1700; tPeriod < 1800; tPeriod++ )
    //     {
    //         for ( uint tAmp = 300; tAmp <= 480; tAmp++ )
    //         {
    //             // for ( uint tNPoints = 80; tNPoints <= 200; tNPoints += 10 )
    //             // {
    //             //     // create surface mesh from object file
    //             std::string  tFile            = "sin_" + std::to_string( tAmp ) + "a_" + std::to_string( tPeriod ) + "p_" + std::to_string( tHgt ) + "h_" + std::to_string( static_cast< int >( 10 * tWid ) ) + "w_" + std::to_string( static_cast< int >( tNPoints ) ) + "np";
    //             std::string  tSurfaceMeshPath = tPath + tFile + ".obj";
    //             Surface_Mesh tSurfaceMesh( load_vertices_from_object_file( tSurfaceMeshPath, tOffsets, tScales ), load_facets_from_object_file( tSurfaceMeshPath ) );

    //             std::cout << "Processing: " << tFile << std::endl;

    //             // for ( real tAgglomerationExponentBash = 2.0; tAgglomerationExponentBash <= 8.0; tAgglomerationExponentBash += 2.0 )
    //             // {
    //             // for ( real tAgglomerationRefBash = 4.8; tAgglomerationRefBash <= 6.2; tAgglomerationRefBash += 0.1 )
    //             // {
    //             //         for ( real tAgglomerationShiftBash = 0.0; tAgglomerationShiftBash <= 2.0; tAgglomerationShiftBash += 1.0 )
    //             //         {
    //             // for ( uint tNumRaysBash = 10; tNumRaysBash <= 70; tNumRaysBash += 2 )
    //             // {
    //             //     for ( real tConeAngleBash = 20.0; tConeAngleBash <= 150.0; tConeAngleBash += 1.0 )
    //             //     {
    //             // tAgglom.mExp = tAgglomerationExponentBash;
    //             // tAgglom.mRef = tAgglomerationRefBash;
    //             // tAgglom.mShift = tAgglomerationShiftBash;
    //             // real tConeAngleBash = (real)tNumRaysBash * 2.0;

    //             // Compute the nodal and global shape diameter
    //             tSurfaceMesh.compute_global_shape_diameter( tAgglom, tConeAngleBash, tNumRaysBash, 1 );

    //             // Print the results
    //             tSDFile << tAmp << " " << tPeriod << " " << tHgt << " " << tWid << " " << tNPoints << " " << tAgglom.mRef << " " << tAgglom.mRef << " " << tAgglom.mShift << " " << tNumRaysBash << " " << tConeAngleBash << " " << tSurfaceMesh.get_facet_shape_diameter( 0 ) << " " << tSurfaceMesh.compute_ddiameter_dvertex( 0 )( 0 ) << " " << tSurfaceMesh.compute_ddiameter_dvertex( 0 )( 1 ) << std::endl;
    //         }
    //     }
    //     // }

    //     // std::cout << " };" << std::endl;
    //     // }
    //     //         }
    //     //     }
    //     // }
    //     //     }
    //     // }
    //     tSDFile.close();
    // }
    // TEST_CASE( "raycast debug", "[MTK], [MTK_Surface_Mesh], [Raycast]" )    // brendan delete
    // {
    //     // real tAmp     = 0.1;
    //     // real tPeriod  = 5;
    //     // real tHgt     = 10.0;
    //     // real tWid     = 1.5;
    //     // uint tNpoints = 100;
    //     // Processing: sin_1a_81p_21h_11w_90np

    //     // // create surface mesh from object file
    //     // std::string    tPath            = "/home/chong/codes/moris/build_dbg/projects/MTK/test/bin/";
    //     std::string    tPath            = "/home/chong/work/SP26/Input_Files/Bar_SD/new/";
    //     std::string    tFile            = "integ_mesh_iter_241";
    //     std::string    tSurfaceMeshPath = tPath + tFile + ".obj";
    //     Vector< real > tOffsets         = { 0.0, 0.0, 0.0 };
    //     Vector< real > tScales          = { 1.0, 1.0, 1.0 };
    //     Surface_Mesh   tSurfaceMesh( load_vertices_from_object_file( tSurfaceMeshPath, tOffsets, tScales ), load_facets_from_object_file( tSurfaceMeshPath ) );

    //     // // Config for shape diameter computation
    //     uint                     tNumRaysBash   = 30;
    //     real                     tConeAngleBash = 120.0;
    //     Agglomeration_Parameters_Constant tAgglom( 2.0, 1.5, 0.0 );

    //     tSurfaceMesh.compute_global_shape_diameter( tAgglom, tConeAngleBash, tNumRaysBash, 1 );

    //     std::ofstream tSDFile;
    //     tSDFile.open( tPath + "../sd_data.txt" );

    //     for ( uint i = 0; i < tSurfaceMesh.get_number_of_vertices(); i++ )
    //     {
    //         tSDFile << ( i > tSurfaceMesh.get_number_of_facets() ? std::numeric_limits< real >::quiet_NaN() : tSurfaceMesh.get_facet_shape_diameter( i ) ) << " " << tSurfaceMesh.compute_ddiameter_dvertex( i )( 0 ) << " " << tSurfaceMesh.compute_ddiameter_dvertex( i )( 1 ) << std::endl;
    //     }
    // }
}    // namespace moris::mtk
