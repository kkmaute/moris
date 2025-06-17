/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_GEN_Intersection_Surface_Mesh.cpp
 *
 */

#include "catch.hpp"
#include <cmath>
#include "paths.hpp"
#include <Kokkos_Core.hpp>

#include "parameters.hpp"
#include "cl_GEN_Surface_Mesh_Geometry.hpp"

#include "fn_check_equal.hpp"

#if MORIS_HAVE_ARBORX
namespace moris::gen
{
    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Surface Mesh Regularization - Isotropic Laplacian", "[gen], [pdv], [regularization], [isotropic laplacian]," )
    {
        if ( par_size() == 1 )
        {
            // get root from environment
            std::string tMorisRoot = moris::get_base_moris_dir();

            // Create parameter list for surface mesh geometry
            Submodule_Parameter_Lists tFieldParameterLists( "GEOMETRIES" );
            tFieldParameterLists.add_parameter_list( prm::create_surface_mesh_geometry_parameter_list() );
            tFieldParameterLists.set( "file_path", tMorisRoot + "projects/GEN/test/data/triangle_sensitivity_oblique.obj" );
            tFieldParameterLists.set( "intersection_tolerance", 1e-9 );
            tFieldParameterLists.set( "regularization_type", gen::Regularization_Type::ISOTROPIC_LAPLACIAN );

            // Create surface mesh parameters struct from parameter list
            Surface_Mesh_Parameters tParameters( tFieldParameterLists( 0 ) );

            // Dummy mesh
            mtk::Mesh* tMesh = nullptr;

            // Other parameters needed to create surface mesh, all dummy
            Vector< ADV >                 tADVs;
            ADV_Manager                   tADVManager;
            Node_Manager                  tNodeManager( tMesh );
            std::shared_ptr< Library_IO > tLibrary = nullptr;

            // Create surface mesh geometry
            Surface_Mesh_Geometry tGeom( tParameters, tNodeManager, tADVs, tADVManager, tLibrary );

            Matrix< DDRMat > tVertexCoordinatesExpected = { { 0.75, 1.625, 0.75, 1.625 }, { 0.5, 0.5, 0.5, 0.5 } };

            // Regularize the surface mesh
            tGeom.regularize();

            // Get the regularized vertex coordinates
            Matrix< DDRMat > tNewCoords = tGeom.get_all_vertex_coordinates();

            CHECK_EQUAL( tNewCoords, tVertexCoordinatesExpected, );
        }
    }
    // TEST_CASE( "Engine 3D Surface Mesh Intersections", "[gen], [pdv], [intersection], [surface mesh geometry 3d]," )
    // {
    //     // get root from environment
    //     std::string tMorisRoot = moris::get_base_moris_dir();

    //     // HMR parameters
    //     Module_Parameter_Lists tHMRParameters( Module_Type::HMR );

    //     tHMRParameters.set( "number_of_elements_per_dimension", 2, 2, 2 );
    //     tHMRParameters.set( "domain_dimensions", 2.0, 2.0, 2.0 );
    //     tHMRParameters.set( "domain_offset", 0.0, 0.0, 0.0 );
    //     tHMRParameters.set( "lagrange_output_meshes", "0" );

    //     tHMRParameters.set( "lagrange_orders", "1" );
    //     tHMRParameters.set( "lagrange_pattern", "0" );
    //     tHMRParameters.set( "bspline_orders", "1" );
    //     tHMRParameters.set( "bspline_pattern", "0" );

    //     tHMRParameters.set( "lagrange_to_bspline", "0" );

    //     tHMRParameters.set( "truncate_bsplines", true );
    //     tHMRParameters.set( "refinement_buffer", 1 );
    //     tHMRParameters.set( "staircase_buffer", 1 );
    //     tHMRParameters.set( "pattern_initial_refinement", 0 );

    //     tHMRParameters.set( "severity_level", 2 );

    //     // Create HMR
    //     hmr::HMR tHMR( tHMRParameters );

    //     // initial refinement
    //     tHMR.perform_initial_refinement();
    //     tHMR.finalize();

    //     mtk::Interpolation_Mesh* tMesh = tHMR.create_interpolation_mesh( 0 );

    //     // surface mesh
    //     Submodule_Parameter_Lists tSurfaceMeshParameterList( "GEOMETRIES" );
    //     tSurfaceMeshParameterList.add_parameter_list( prm::create_surface_mesh_geometry_parameter_list() );
    //     tSurfaceMeshParameterList.set( "file_path", tMorisRoot + "projects/GEN/test/data/tetra.obj" );
    //     tSurfaceMeshParameterList.set( "intersection_tolerance", 1e-8 );

    //     // Create geometry engine
    //     Geometry_Engine_Parameters tGeometryEngineParameters;
    //     ADV_Manager                tADVManager;
    //     Design_Factory             tDesignFactory( tSurfaceMeshParameterList, tADVManager );
    //     tGeometryEngineParameters.mGeometries = tDesignFactory.get_geometries();
    //     Geometry_Engine tGeometryEngine( tMesh, tGeometryEngineParameters );

    //     // Solution to is_intersected per element, per edge
    //     Vector< Vector< bool > > tIsEdgeIntersected = {
    //         { false, false, false, false, false, false, false, false, false, false, false, false },    // Element 0
    //         { false, false, false, false, false, false, false, false, false, false, false, false },    // Element 1
    //         { false, false, false, false, false, false, false, false, false, false, false, false },    // Element 2
    //         { false, true, true, false, false, false, false, false, false, false, true, false },       // Element 3
    //         { false, false, false, false, false, false, false, false, false, false, false, false },    // Element 4
    //         { false, false, false, false, false, false, false, false, false, false, false, false },    // Element 5
    //         { false, false, false, false, false, false, false, false, false, false, false, false },    // Element 6
    //         { false, false, false, false, false, false, false, false, false, false, false, false },    // Element 7
    //     };

    //     // Intersection local coordinates solutions
    //     Vector< real > tIntersectionLocalCoordinates = { 0.045454545454545454, 0.291666666666667, -0.642857142857142 };

    //     // Intersection global coordinates solutions
    //     Vector< Matrix< DDRMat > > tIntersectionGlobalCoordinates = {
    //         { { 2.0, 1.522727272727272727272727, 0.0 } },
    //         { { 1.354166666666666666, 2.0, 0.0 } },
    //         { { 2.0, 2.0, 0.178571428571429 } }
    //     };

    //     uint tIntersectionCount = 0;
    //     for ( uint iElementIndex = 0; iElementIndex < 8; iElementIndex++ )
    //     {
    //         // Node indices per element
    //         Matrix< IndexMat > tSignedNodeIndices = tMesh->get_nodes_connected_to_element_loc_inds( iElementIndex );
    //         Matrix< DDUMat >   tNodeIndices( 8, 1 );
    //         for ( uint iNode = 0; iNode < 8; iNode++ )
    //         {
    //             tNodeIndices( iNode ) = tSignedNodeIndices( iNode );
    //         }

    //         for ( uint iEdgeNumber = 0; iEdgeNumber < 12; iEdgeNumber++ )
    //         {

    //             // Get the geometry engine result
    //             bool tIntersected = tGeometryEngine.is_intersected_by_active_geometry( { { tSignedNodeIndices( tEdgeOrder( iEdgeNumber )( 0 ) ), tSignedNodeIndices( tEdgeOrder( iEdgeNumber )( 1 ) ) } } );

    //             // check that these are equal
    //             REQUIRE( tIntersected == tIsEdgeIntersected( iElementIndex )( iEdgeNumber ) );

    //             // Check queued intersection
    //             if ( tIntersected )
    //             {
    //                 // Queue intersection
    //                 bool tQueryIntersected = tGeometryEngine.queue_intersection(
    //                         tSignedNodeIndices( tEdgeOrder( iEdgeNumber )( 0 ) ),
    //                         tSignedNodeIndices( tEdgeOrder( iEdgeNumber )( 1 ) ),
    //                         tHexParametricCoordinates( tEdgeOrder( iEdgeNumber )( 0 ) ),
    //                         tHexParametricCoordinates( tEdgeOrder( iEdgeNumber )( 1 ) ),
    //                         tNodeIndices,
    //                         mtk::Geometry_Type::HEX,
    //                         mtk::Interpolation_Order::LINEAR );

    //                 // Check that the query was successful
    //                 CHECK( tQueryIntersected == tIsEdgeIntersected( iElementIndex )( iEdgeNumber ) );

    //                 // Check parents
    //                 bool tFirstParentOnInterface  = false;
    //                 bool tSecondParentOnInterface = false;

    //                 CHECK( tGeometryEngine.queued_intersection_first_parent_on_interface() == tFirstParentOnInterface );
    //                 CHECK( tGeometryEngine.queued_intersection_second_parent_on_interface() == tSecondParentOnInterface );

    //                 // See if local coordinate is a number
    //                 real tLocalCoordinate = tGeometryEngine.get_queued_intersection_local_coordinate();

    //                 // Ensure the local coordinate is a number
    //                 REQUIRE( !std::isnan( tLocalCoordinate ) );

    //                 CHECK( tLocalCoordinate == Approx( tIntersectionLocalCoordinates( tIntersectionCount ) ).margin( 1e-9 ) );

    //                 CHECK_EQUAL( tGeometryEngine.get_queued_intersection_global_coordinates(), tIntersectionGlobalCoordinates( tIntersectionCount ), );

    //                 // Admit intersection
    //                 tGeometryEngine.admit_queued_intersection();
    //                 tIntersectionCount++;
    //             }
    //         }
    //     }
    // }
}    // namespace moris::gen

#endif
