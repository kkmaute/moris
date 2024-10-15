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
#include "fn_eye.hpp"
#include "paths.hpp"
#include <Kokkos_Core.hpp>

#include "cl_GEN_Geometry_Engine_Test.hpp"
#include "cl_GEN_PDV_Host_Manager.hpp"
#include "cl_GEN_Design_Factory.hpp"
#include "fn_GEN_create_simple_mesh.hpp"
#include "cl_GEN_Intersection_Node.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "cl_GEN_Background_Node.hpp"
#include "cl_GEN_Parent_Node.hpp"

#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Mesh_Integration.hpp"
#include "fn_PRM_HMR_Parameters.hpp"

#include "fn_check_equal.hpp"

#ifdef MORIS_HAVE_ARBORX
// initialize Kokkos for the use in the spatial tree library ArborX
std::unique_ptr< Kokkos::ScopeGuard > guard = !Kokkos::is_initialized() && !Kokkos::is_finalized() ? std::make_unique< Kokkos::ScopeGuard >() : nullptr;
namespace moris::gen
{

#ifdef MORIS_HAVE_ARBORX
    // initialize Kokkos for the use in the spatial tree library ArborX
    std::unique_ptr< Kokkos::ScopeGuard > guard = !Kokkos::is_initialized() && !Kokkos::is_finalized() ? std::make_unique< Kokkos::ScopeGuard >() : nullptr;
#endif
    //--------------------------------------------------------------------------------------------------------------

    static Vector< Matrix< DDRMat > > tQuadParametricCoordinates = {
        { { -1.0, -1.0 } },
        { { 1.0, -1.0 } },
        { { 1.0, 1.0 } },
        { { -1.0, 1.0 } }
    };

    static Vector< Matrix< DDRMat > > tHexParametricCoordinates = {
        { { -1.0, -1.0, -1.0 } },
        { { 1.0, -1.0, -1.0 } },
        { { 1.0, 1.0, -1.0 } },
        { { -1.0, 1.0, -1.0 } },
        { { -1.0, -1.0, 1.0 } },
        { { 1.0, -1.0, 1.0 } },
        { { 1.0, 1.0, 1.0 } },
        { { -1.0, 1.0, 1.0 } }
    };

    static Vector< Vector< uint > > tEdgeOrder = {
        { 0, 1 },
        { 1, 2 },
        { 2, 3 },
        { 3, 0 },
        { 4, 5 },
        { 5, 6 },
        { 6, 7 },
        { 7, 4 },
        { 0, 4 },
        { 1, 5 },
        { 2, 6 },
        { 3, 7 },
    };

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Engine 2D Surface Mesh Intersections", "[gen], [pdv], [intersection], [surface mesh geometry 2d]," )
    {
        if ( par_size() == 1 )
        {
            // get root from environment
            std::string tMorisRoot = moris::get_base_moris_dir();

            // Create mesh
            Parameter_List tHMRParameters = prm::create_hmr_parameter_list();

            tHMRParameters.set( "number_of_elements_per_dimension", "2,1" );
            tHMRParameters.set( "domain_dimensions", "2, 1" );
            tHMRParameters.set( "domain_offset", "0.0, 0.0" );
            tHMRParameters.set( "domain_sidesets", "1,2,3,4" );
            tHMRParameters.set( "lagrange_output_meshes", "0" );

            tHMRParameters.set( "lagrange_orders", "1" );
            tHMRParameters.set( "lagrange_pattern", "0" );
            tHMRParameters.set( "bspline_orders", "1" );
            tHMRParameters.set( "bspline_pattern", "0" );
            tHMRParameters.set( "lagrange_to_bspline", "0" );

            tHMRParameters.set( "initial_refinement", "0" );
            tHMRParameters.set( "truncate_bsplines", 1 );
            tHMRParameters.set( "refinement_buffer", 1 );
            tHMRParameters.set( "staircase_buffer", 1 );

            tHMRParameters.set( "severity_level", 2 );

            hmr::HMR tHMR( tHMRParameters );

            // initial refinement
            tHMR.perform_initial_refinement();
            tHMR.finalize();

            mtk::Interpolation_Mesh* tMesh = tHMR.create_interpolation_mesh( 0 );

            // surface mesh
            Submodule_Parameter_Lists tFieldParameterLists( "FIELDS" );
            tFieldParameterLists.add_parameter_list( prm::create_surface_mesh_geometry_parameter_list() );
            tFieldParameterLists.set( "file_path", tMorisRoot + "projects/GEN/SDF/test/data/triangle_sensitivity_oblique.obj" );
            tFieldParameterLists.set( "intersection_tolerance", 1e-9 );


            // Create geometry engine
            Geometry_Engine_Parameters tGeometryEngineParameters;
            ADV_Manager tADVManager;
            Design_Factory tDesignFactory( tFieldParameterLists, tADVManager );
            tGeometryEngineParameters.mGeometries = tDesignFactory.get_geometries();
            Geometry_Engine tGeometryEngine( tMesh, tGeometryEngineParameters );

            // Solution to is_intersected per element, per edge
            Vector< Vector< bool > > tIsEdgeIntersected = {
                { true, true, false, false },    // Element 0
                { true, false, false, true }     // Element 1
            };

            // Intersection local coordinates solutions
            Vector< real > tIntersectionLocalCoordinates = { { 0.5, 0.5625, -5.0 / 6.0, -0.5625 } };

            // Intersection global coordinates solutions
            Vector< Matrix< DDRMat > > tIntersectionGlobalCoordinates = {
                { { 0.75, 0.0 } },
                { { 1.0, 0.78125 } },
                { { 13.0 / 12.0, 0.0 } },
                { { 1.0, 0.78125 } }
            };

            uint tIntersectionCount = 0;
            for ( uint iElementIndex = 0; iElementIndex < 2; iElementIndex++ )
            {
                // Node indices per element
                Matrix< IndexMat > tSignedNodeIndices = tMesh->get_nodes_connected_to_element_loc_inds( iElementIndex );
                Matrix< DDUMat >   tNodeIndices( 4, 1 );
                for ( uint iNode = 0; iNode < 4; iNode++ )
                {
                    tNodeIndices( iNode ) = tSignedNodeIndices( iNode );
                }

                for ( uint iNodeNumber = 0; iNodeNumber < 4; iNodeNumber++ )
                {
                    // Get the geometry engine result
                    bool tIntersected = tGeometryEngine.is_intersected_by_active_geometry( { { tSignedNodeIndices( iNodeNumber ), tSignedNodeIndices( ( iNodeNumber + 1 ) % 4 ) } } );

                    // check that these are equal
                    REQUIRE( tIntersected == tIsEdgeIntersected( iElementIndex )( iNodeNumber ) );

                    // Check queued intersection
                    if ( tIntersected )
                    {
                        // Queue intersection
                        bool tQueryIntersected = tGeometryEngine.queue_intersection(
                                tSignedNodeIndices( iNodeNumber ),
                                tSignedNodeIndices( ( iNodeNumber + 1 ) % 4 ),
                                tQuadParametricCoordinates( iNodeNumber ),
                                tQuadParametricCoordinates( ( iNodeNumber + 1 ) % 4 ),
                                tNodeIndices,
                                mtk::Geometry_Type::QUAD,
                                mtk::Interpolation_Order::LINEAR );

                        // Check that the query was successful
                        CHECK( tQueryIntersected == tIsEdgeIntersected( iElementIndex )( iNodeNumber ) );

                        // Check parents
                        bool tFirstParentOnInterface  = false;
                        bool tSecondParentOnInterface = false;

                        CHECK( tGeometryEngine.queued_intersection_first_parent_on_interface() == tFirstParentOnInterface );
                        CHECK( tGeometryEngine.queued_intersection_second_parent_on_interface() == tSecondParentOnInterface );

                        // See if local coordinate is a number
                        real tLocalCoordinate = tGeometryEngine.get_queued_intersection_local_coordinate();

                        // Ensure the local coordinate is a number
                        REQUIRE( !std::isnan( tLocalCoordinate ) );

                        CHECK( tLocalCoordinate == Approx( tIntersectionLocalCoordinates( tIntersectionCount ) ).margin( 1e-9 ) );

                        // Check global coordinates
                        CHECK_EQUAL( tGeometryEngine.get_queued_intersection_global_coordinates(), tIntersectionGlobalCoordinates( tIntersectionCount ), );

                        // Admit intersection
                        tGeometryEngine.admit_queued_intersection();
                        tIntersectionCount++;
                    }
                }
            }
        }
    }
    TEST_CASE( "Engine 3D Surface Mesh Intersections", "[gen], [pdv], [intersection], [surface mesh geometry 3d]," )
    {
        // get root from environment
        std::string tMorisRoot = moris::get_base_moris_dir();

        // HMR parameters
        Parameter_List tHMRParameters = prm::create_hmr_parameter_list();

        tHMRParameters.set( "number_of_elements_per_dimension", "2,2,2" );
        tHMRParameters.set( "domain_dimensions", "2,2,2" );
        tHMRParameters.set( "domain_sidesets", "1,2,3,4,5,6" );
        tHMRParameters.set( "domain_offset", "0.0,0.0,0.0" );
        tHMRParameters.set( "lagrange_output_meshes", "0" );

        tHMRParameters.set( "lagrange_orders", "1" );
        tHMRParameters.set( "lagrange_pattern", "0" );
        tHMRParameters.set( "bspline_orders", "1" );
        tHMRParameters.set( "bspline_pattern", "0" );

        tHMRParameters.set( "lagrange_to_bspline", "0" );

        tHMRParameters.set( "truncate_bsplines", 1 );
        tHMRParameters.set( "refinement_buffer", 1 );
        tHMRParameters.set( "staircase_buffer", 1 );
        tHMRParameters.set( "initial_refinement", "0" );
        tHMRParameters.set( "initial_refinement_pattern", "0" );

        tHMRParameters.set( "severity_level", 2 );

        // Create HMR
        hmr::HMR tHMR( tHMRParameters );

        // initial refinement
        tHMR.perform_initial_refinement();
        tHMR.finalize();

        mtk::Interpolation_Mesh* tMesh = tHMR.create_interpolation_mesh( 0 );

        // surface mesh
        Parameter_List tSurfaceMeshParameterList = prm::create_surface_mesh_geometry_parameter_list();
        tSurfaceMeshParameterList.set( "file_path", tMorisRoot + "projects/GEN/test/data/tetra.obj" );
        tSurfaceMeshParameterList.set( "intersection_tolerance", 1e-8 );

        // Create geometry engine
        Geometry_Engine_Parameters tGeometryEngineParameters;
        ADV_Manager                tADVManager;
        Design_Factory             tDesignFactory( { tSurfaceMeshParameterList }, tADVManager );
        tGeometryEngineParameters.mGeometries = tDesignFactory.get_geometries();
        Geometry_Engine tGeometryEngine( tMesh, tGeometryEngineParameters );

        // Solution to is_intersected per element, per edge
        Vector< Vector< bool > > tIsEdgeIntersected = {
            { false, false, false, false, false, false, false, false, false, false, false, false },    // Element 0
            { false, false, false, false, false, false, false, false, false, false, false, false },    // Element 1
            { false, false, false, false, false, false, false, false, false, false, false, false },    // Element 2
            { false, true, true, false, false, false, false, false, false, false, true, false },       // Element 3
            { false, false, false, false, false, false, false, false, false, false, false, false },    // Element 4
            { false, false, false, false, false, false, false, false, false, false, false, false },    // Element 5
            { false, false, false, false, false, false, false, false, false, false, false, false },    // Element 6
            { false, false, false, false, false, false, false, false, false, false, false, false },    // Element 7
        };

        // Intersection local coordinates solutions
        Vector< real > tIntersectionLocalCoordinates = { 0.045454545454545454, 0.291666666666667, -0.642857142857142 };

        // Intersection global coordinates solutions
        Vector< Matrix< DDRMat > > tIntersectionGlobalCoordinates = {
            { { 2.0, 1.522727272727272727272727, 0.0 } },
            { { 1.354166666666666666, 2.0, 0.0 } },
            { { 2.0, 2.0, 0.178571428571429 } }
        };

        uint tIntersectionCount = 0;
        for ( uint iElementIndex = 0; iElementIndex < 8; iElementIndex++ )
        {
            // Node indices per element
            Matrix< IndexMat > tSignedNodeIndices = tMesh->get_nodes_connected_to_element_loc_inds( iElementIndex );
            Matrix< DDUMat >   tNodeIndices( 8, 1 );
            for ( uint iNode = 0; iNode < 8; iNode++ )
            {
                tNodeIndices( iNode ) = tSignedNodeIndices( iNode );
            }

            for ( uint iEdgeNumber = 0; iEdgeNumber < 12; iEdgeNumber++ )
            {

                // Get the geometry engine result
                bool tIntersected = tGeometryEngine.is_intersected_by_active_geometry( { { tSignedNodeIndices( tEdgeOrder( iEdgeNumber )( 0 ) ), tSignedNodeIndices( tEdgeOrder( iEdgeNumber )( 1 ) ) } } );

                // check that these are equal
                REQUIRE( tIntersected == tIsEdgeIntersected( iElementIndex )( iEdgeNumber ) );

                // Check queued intersection
                if ( tIntersected )
                {
                    // Queue intersection
                    bool tQueryIntersected = tGeometryEngine.queue_intersection(
                            tSignedNodeIndices( tEdgeOrder( iEdgeNumber )( 0 ) ),
                            tSignedNodeIndices( tEdgeOrder( iEdgeNumber )( 1 ) ),
                            tHexParametricCoordinates( tEdgeOrder( iEdgeNumber )( 0 ) ),
                            tHexParametricCoordinates( tEdgeOrder( iEdgeNumber )( 1 ) ),
                            tNodeIndices,
                            mtk::Geometry_Type::HEX,
                            mtk::Interpolation_Order::LINEAR );

                    // Check that the query was successful
                    CHECK( tQueryIntersected == tIsEdgeIntersected( iElementIndex )( iEdgeNumber ) );

                    // Check parents
                    bool tFirstParentOnInterface  = false;
                    bool tSecondParentOnInterface = false;

                    CHECK( tGeometryEngine.queued_intersection_first_parent_on_interface() == tFirstParentOnInterface );
                    CHECK( tGeometryEngine.queued_intersection_second_parent_on_interface() == tSecondParentOnInterface );

                    // See if local coordinate is a number
                    real tLocalCoordinate = tGeometryEngine.get_queued_intersection_local_coordinate();

                    // Ensure the local coordinate is a number
                    REQUIRE( !std::isnan( tLocalCoordinate ) );

                    CHECK( tLocalCoordinate == Approx( tIntersectionLocalCoordinates( tIntersectionCount ) ).margin( 1e-9 ) );

                    CHECK_EQUAL( tGeometryEngine.get_queued_intersection_global_coordinates(), tIntersectionGlobalCoordinates( tIntersectionCount ), );

                    // Admit intersection
                    tGeometryEngine.admit_queued_intersection();
                    tIntersectionCount++;
                }
            }
        }
    }
}    // namespace moris::gen

#endif
