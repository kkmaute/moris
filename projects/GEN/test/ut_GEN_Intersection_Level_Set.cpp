/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_GEN_Intersection_Level_Set.cpp
 *
 */

#include "catch.hpp"
#include <cmath>
#include "fn_eye.hpp"

#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Geometry_Engine_Test.hpp"
#include "cl_GEN_PDV_Host_Manager.hpp"
#include "cl_GEN_BSpline_Field.hpp"
#include "cl_GEN_Design_Factory.hpp"
#include "fn_check_equal.hpp"
#include "fn_GEN_create_simple_mesh.hpp"
#include "cl_GEN_Mesh_Field.hpp"
#include "cl_MTK_Mesh_Factory.hpp"

#include "fn_PRM_GEN_Parameters.hpp"

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    static Vector< Matrix< DDRMat > > tQuadParametricCoordinates = {
        { { -1.0, -1.0 } },
        { { 1.0, -1.0 } },
        { { 1.0, 1.0 } },
        { { -1.0, 1.0 } }
    };

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    get_quad_local_coordinates( uint tNodeNumber )
    {
        real tRadians = M_PI * ( 5 + 2 * tNodeNumber ) / 4;
        return sqrt( 2 ) * Matrix< DDRMat >( { { cos( tRadians ), sin( tRadians ) } } );
    }

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Linear Intersections", "[gen], [pdv], [intersection], [linear intersection]" )
    {
        if ( par_size() == 1 )
        {
            // Create mesh
            mtk::Interpolation_Mesh* tMesh = create_simple_mesh( 2, 2 );

            // Set up geometry
            Vector< real > tADVs = { 0.25, 0.0, 1.0, 0.0 };

            // Circle
            Parameter_List tCircleParameterList = prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE );
            tCircleParameterList.set( "center_x", -0.25 );
            tCircleParameterList.set( "center_y", 0.0 );
            tCircleParameterList.set( "radius", 0.7499999999 );
            tCircleParameterList.set( "discretization_mesh_index", 0 );

            // Plane 1
            Parameter_List tPlane1ParameterList = prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE );
            tPlane1ParameterList.set( "center_x", 0.25, 0.25, 0.25 );
            tPlane1ParameterList.set( "center_y", 0.0, 0.0, 0.0 );
            tPlane1ParameterList.set( "normal_x", 1.0, 1.0, 1.0 );
            tPlane1ParameterList.set( "normal_y", 0.0, 0.0, 0.0 );

            // Plane 2
            Parameter_List tPlane2ParameterList = prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE );
            tPlane2ParameterList.set( "center_x", 1.0 );
            tPlane2ParameterList.set( "center_y", 0.0 );
            tPlane2ParameterList.set( "normal_x", 1.0 );
            tPlane2ParameterList.set( "normal_y", 0.0 );

            // Create geometry engine
            Geometry_Engine_Parameters tGeometryEngineParameters;
            ADV_Manager tADVManager;
            tADVManager.mADVs = tADVs;
            Design_Factory tDesignFactory( { tCircleParameterList, tPlane1ParameterList, tPlane2ParameterList }, tADVManager );
            tGeometryEngineParameters.mADVManager = tADVManager;
            tGeometryEngineParameters.mGeometries = tDesignFactory.get_geometries();
            Geometry_Engine_Test tGeometryEngine( tMesh, tGeometryEngineParameters );

            // TODO ensure this writes the mesh/fields correctly instead of just relying on no errors being thrown
            tGeometryEngine.output_fields_on_mesh( tMesh, "intersection_test.exo" );

            // Solution for is_intersected() per geometry and per element
            Vector< Vector< bool > > tIsElementIntersected = {
                { true, true, true, true },      // Geometry 0
                { false, true, false, true },    // Geometry 1
                { false, true, false, true }
            };    // Geometry 2

            // Per geometry, per element, per edge
            Vector< Vector< Vector< bool > > > tIsEdgeIntersected = {
                // Geometry 0
                { { false, true, true, false },            // Element 0
                        { false, false, true, true },      // Element 1
                        { true, true, false, false },      // Element 2
                        { true, false, false, true } },    // Element 3
                // Geometry 1
                { { false, false, false, false },
                        { true, false, true, false },
                        { false, false, false, false },
                        { true, false, true, false } },
                // Geometry 2
                { { false, false, false, false },
                        { true, true, true, false },
                        { false, false, false, false },
                        { true, true, true, false } }
            };

            // Intersection coordinates
            real tFrac = 2.0 / ( 3.0 + sqrt( 17.0 ) );

            Matrix< DDRMat > tIntersectionLocalCoordinates = { { -tFrac, 1.0, 0.0, tFrac, -1.0, tFrac, 0.0, -tFrac, -0.5, 0.5, -0.5, 0.5, 0.0000000002, -0.0000000002, 1.0, 0.0, -1.0, 1.0, 0.0, -1.0 } };

            Vector< Matrix< DDRMat > > tIntersectionGlobalCoordinates = {
                { { 0.0, -0.5 - ( tFrac / 2.0 ) } },
                { { -1.0, 0.0 } },
                { { 0.5, 0.0 } },
                { { 0.0, -0.5 - ( tFrac / 2.0 ) } },
                { { -1.0, 0.0 } },
                { { 0.0, 0.5 + ( tFrac / 2.0 ) } },
                { { 0.5, 0.0 } },
                { { 0.0, 0.5 + ( tFrac / 2.0 ) } },
                { { 0.25, -1.0 } },
                { { 0.25, 0.0 } },
                { { 0.25, 0.0 } },
                { { 0.25, 1.0 } },
                { { 0.25, -0.25 - tFrac / 4.0 } },
                { { 0.25, 0.25 + tFrac / 4.0 } },
                { { 1.0, -1.0 } },
                { { 1.0, -1.0 } },
                { { 1.0, 0.0 } },
                { { 1.0, 0.0 } },
                { { 1.0, 0.0 } },
                { { 1.0, 1.0 } }
            };

            // Check element intersections
            uint tIntersectionCount = 0;
            for ( uint tGeometryIndex = 0; tGeometryIndex < 3; tGeometryIndex++ )
            {
                for ( uint tElementIndex = 0; tElementIndex < 4; tElementIndex++ )
                {
                    // Node indices per element
                    Matrix< IndexMat > tSignedNodeIndices = tMesh->get_nodes_connected_to_element_loc_inds( tElementIndex );
                    Matrix< DDUMat >   tNodeIndices( 4, 1 );
                    for ( uint iNode = 0; iNode < 4; iNode++ )
                    {
                        tNodeIndices( iNode ) = tSignedNodeIndices( iNode );
                    }

                    // Check edges for properly queued intersections
                    Matrix< DDRMat > tNodeCoordinates( 4, 2 );
                    for ( uint tNodeNumber = 0; tNodeNumber < 4; tNodeNumber++ )
                    {
                        // Node coordinates
                        Matrix< DDRMat > tFirstNodeCoordinates  = tMesh->get_node_coordinate( tNodeIndices( tNodeNumber ) );
                        Matrix< DDRMat > tSecondNodeCoordinates = tMesh->get_node_coordinate( tNodeIndices( ( tNodeNumber + 1 ) % 4 ) );

                        // Queue intersection
                        bool tIntersectionQueued = tGeometryEngine.queue_intersection(
                                tNodeIndices( tNodeNumber ),
                                tNodeIndices( ( tNodeNumber + 1 ) % 4 ),
                                tQuadParametricCoordinates( tNodeNumber ),
                                tQuadParametricCoordinates( ( tNodeNumber + 1 ) % 4 ),
                                tNodeIndices,
                                mtk::Geometry_Type::QUAD,
                                mtk::Interpolation_Order::LINEAR );
                        REQUIRE( tIntersectionQueued == tIsEdgeIntersected( tGeometryIndex )( tElementIndex )( tNodeNumber ) );

                        // Check queued intersection
                        if ( tIntersectionQueued )
                        {
                            // Check parents
                            bool tFirstParentOnInterface  = false;
                            bool tSecondParentOnInterface = false;

                            // Parent nodes on interface
                            if ( tGeometryIndex == 2 and tFirstNodeCoordinates( 0 ) == Approx( 1.0 ) )
                            {
                                tFirstParentOnInterface = true;
                            }

                            if ( tGeometryIndex == 2 and tSecondNodeCoordinates( 0 ) == Approx( 1.0 ) )
                            {
                                tSecondParentOnInterface = true;
                            }
                            CHECK( tGeometryEngine.queued_intersection_first_parent_on_interface() == tFirstParentOnInterface );
                            CHECK( tGeometryEngine.queued_intersection_second_parent_on_interface() == tSecondParentOnInterface );

                            // See if local coordinate is a number
                            real tLocalCoordinate = tGeometryEngine.get_queued_intersection_local_coordinate();
                            if ( not std::isnan( tLocalCoordinate ) )
                            {
                                // Check local coordinate; note reference values are for radius of exactly 0.75; thus higher margin needed
                                CHECK( tLocalCoordinate == Approx( tIntersectionLocalCoordinates( tIntersectionCount ) ).margin( 1e-9 ) );

                                // Check global coordinates
                                CHECK_EQUAL(
                                        tGeometryEngine.get_queued_intersection_global_coordinates(),
                                        tIntersectionGlobalCoordinates( tIntersectionCount ), );
                            }

                            // Admit intersection
                            tGeometryEngine.admit_queued_intersection();
                            tIntersectionCount++;
                        }

                        // Set node coordinates for element checking
                        tNodeCoordinates.set_row( tNodeNumber, tMesh->get_node_coordinate( tNodeIndices( tNodeNumber ) ) );
                    }

                    // Check with solution
                    bool tIsIntersected = tGeometryEngine.is_intersected( tGeometryIndex, tSignedNodeIndices, tNodeCoordinates );
                    CHECK( tIsIntersected == tIsElementIntersected( tGeometryIndex )( tElementIndex ) );
                }

                // Intersection on intersection
                if ( tGeometryIndex == 1 )
                {
                    // Queue intersection on intersection 1
                    bool tIntersectionQueued = tGeometryEngine.queue_intersection(
                            9,
                            11,
                            { { -1.0, -tFrac } },
                            { { 0.0, 1.0 } },
                            { { 1, 4, 5, 2 } },
                            mtk::Geometry_Type::QUAD,
                            mtk::Interpolation_Order::LINEAR );

                    // Check intersection on intersection 1
                    REQUIRE( tIntersectionQueued == true );
                    CHECK( tGeometryEngine.queued_intersection_first_parent_on_interface() == false );
                    CHECK( tGeometryEngine.queued_intersection_second_parent_on_interface() == false );
                    CHECK( tGeometryEngine.get_queued_intersection_local_coordinate() == Approx( tIntersectionLocalCoordinates( tIntersectionCount ) ) );
                    CHECK_EQUAL(
                            tGeometryEngine.get_queued_intersection_global_coordinates(),
                            tIntersectionGlobalCoordinates( tIntersectionCount ), );

                    // Admit intersection on intersection 1
                    tGeometryEngine.admit_queued_intersection();
                    tIntersectionCount++;

                    // Queue intersection on intersection 2
                    tIntersectionQueued = tGeometryEngine.queue_intersection(
                            11,
                            14,
                            { { 0.0, -1.0 } },
                            { { -1.0, tFrac } },
                            { { 2, 5, 8, 6 } },
                            mtk::Geometry_Type::QUAD,
                            mtk::Interpolation_Order::LINEAR );

                    // Check intersection on intersection 1
                    REQUIRE( tIntersectionQueued == true );
                    CHECK( tGeometryEngine.queued_intersection_first_parent_on_interface() == false );
                    CHECK( tGeometryEngine.queued_intersection_second_parent_on_interface() == false );
                    CHECK( tGeometryEngine.get_queued_intersection_local_coordinate() == Approx( tIntersectionLocalCoordinates( tIntersectionCount ) ) );
                    CHECK_EQUAL(
                            tGeometryEngine.get_queued_intersection_global_coordinates(),
                            tIntersectionGlobalCoordinates( tIntersectionCount ), );

                    // Admit intersection on intersection 1
                    tGeometryEngine.admit_queued_intersection();
                    tIntersectionCount++;
                }

                // Advance geometry index
                if ( tGeometryIndex < 2 )
                {
                    tGeometryEngine.advance_geometry_index();
                }
            }

            // Check total number of intersections
            CHECK( tIntersectionCount == 20 );

            // Test the new child nodes on the level set field (geometry 0)
            CHECK( tGeometryEngine.get_geometric_region( 0, 9, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 10, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 11, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 12, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 13, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 14, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 15, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 16, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 17, { {} } ) == Geometric_Region::POSITIVE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 18, { {} } ) == Geometric_Region::NEGATIVE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 19, { {} } ) == Geometric_Region::NEGATIVE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 20, { {} } ) == Geometric_Region::POSITIVE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 21, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 22, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 23, { {} } ) == Geometric_Region::POSITIVE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 25, { {} } ) == Geometric_Region::POSITIVE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 26, { {} } ) == Geometric_Region::POSITIVE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 28, { {} } ) == Geometric_Region::POSITIVE );

            // Get the PDV host manager
            auto tPDVHostManager = dynamic_cast< PDV_Host_Manager* >( tGeometryEngine.get_design_variable_interface() );

            // Test that the new intersections have been added to the PDV host manager, but ONLY for the circle
            Vector< Matrix< DDRMat > > tPDVValues( 0 );
            Vector< Vector< bool > >     tIsActive( 0 );
            tPDVHostManager->get_ig_pdv_value(
                    { { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28 } },
                    { PDV_Type::X_COORDINATE, PDV_Type::Y_COORDINATE },
                    tPDVValues,
                    tIsActive );

            // Background nodes
            for ( uint tNodeIndex = 0; tNodeIndex < 9; tNodeIndex++ )
            {
                // Check if not active
                CHECK( tIsActive( 0 )( tNodeIndex ) == false );
                CHECK( tIsActive( 1 )( tNodeIndex ) == false );
            }

            // Nodes on the circle and first plane interfaces (depend on ADVs, active)
            for ( uint tNodeIndex = 9; tNodeIndex < 23; tNodeIndex++ )
            {
                // Check if active
                CHECK( tIsActive( 0 )( tNodeIndex ) == true );
                CHECK( tIsActive( 1 )( tNodeIndex ) == true );

                // Check for node coordinates as PDV values
                CHECK( tPDVValues( 0 )( tNodeIndex ) == Approx( tIntersectionGlobalCoordinates( tNodeIndex - 9 )( 0 ) ) );
                CHECK( tPDVValues( 1 )( tNodeIndex ) == Approx( tIntersectionGlobalCoordinates( tNodeIndex - 9 )( 1 ) ) );
            }

            // Nodes on the second plane interface (inactive)
            for ( uint tNodeIndex = 23; tNodeIndex < 27; tNodeIndex++ )
            {
                // Check if not active
                CHECK( tIsActive( 0 )( tNodeIndex ) == false );
                CHECK( tIsActive( 1 )( tNodeIndex ) == false );
            }

            // Check sensitivities
            Vector< Matrix< DDRMat > > tIntersectionSensitivities = {
                { { 0.0, 0.0 }, { ( 9 + sqrt( 17 ) ) / 16, ( 3 * sqrt( 17 ) - 5 ) / 16 } },
                { { 0.0, 2.0 }, { 0.0, 0.0 } },
                { { -0.5, -0.5 }, { 0.0, 0.0 } },
                { { 0.0, 0.0 }, { ( 3 * sqrt( 17 ) - 5 ) / 16, ( 9 + sqrt( 17 ) ) / 16 } },
                { { 2.0, 0.0 }, { 0.0, 0.0 } },
                { { 0.0, 0.0 }, { -( 3 * sqrt( 17 ) - 5 ) / 16, -( 9 + sqrt( 17 ) ) / 16 } },
                { { -0.5, -0.5 }, { 0.0, 0.0 } },
                { { 0.0, 0.0 }, { -( 9 + sqrt( 17 ) ) / 16, -( 3 * sqrt( 17 ) - 5 ) / 16 } },
                { { 0.75, 0.0, 0.1875, 0.75, 0.25, 0.0, -0.1875, 0.25 }, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } },
                { { 0.25, 0.0, -0.1875, 0.0, 0.75, 0.0, 0.1875, 0.0 }, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } },
                { { 0.75, 0.0, 0.1875, 0.0, 0.25, 0.0, -0.1875, 0.0 }, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } },
                { { 0.25, 0.0, -0.1875, -0.25, 0.75, 0.0, 0.1875, -0.75 }, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } },
                { { 0.5, 0.0, 0.125, ( 1 + sqrt( 17 ) ) / 16, 0.5, 0.0, -0.125, 0.0, 0.0, 0.0, 0.0, 0.0 }, { ( 1 + sqrt( 17 ) ) / 8, 0.0, ( 1 + sqrt( 17 ) ) / 32, ( 9 + sqrt( 17 ) ) / 32, ( 1 + sqrt( 17 ) ) / 8, 0.0, -( 1 + sqrt( 17 ) ) / 32, 0.0, ( 9 + sqrt( 17 ) ) / 32, ( 3 * sqrt( 17 ) - 5 ) / 32, ( 1 + sqrt( 17 ) ) / 16, ( 1 + sqrt( 17 ) ) / 16 } },
                { { 0.5, 0.0, -0.125, 0.0, 0.5, 0.0, 0.125, -( 1 + sqrt( 17 ) ) / 16, 0.0, 0.0, 0.0, 0.0 }, { -( 1 + sqrt( 17 ) ) / 8, 0.0, ( 1 + sqrt( 17 ) ) / 32, 0.0, -( 1 + sqrt( 17 ) ) / 8, 0.0, -( 1 + sqrt( 17 ) ) / 32, ( 9 + sqrt( 17 ) ) / 32, -( 1 + sqrt( 17 ) ) / 16, -( 1 + sqrt( 17 ) ) / 16, -( 3 * sqrt( 17 ) - 5 ) / 32, -( 9 + sqrt( 17 ) ) / 32 } }
            };

            Vector< Vector< sint > > tIntersectionIDs = {
                { 6, 9 },
                { 9, 8 },
                { 10, 9 },
                { 9, 6 },
                { 8, 9 },
                { 9, 12 },
                { 9, 10 },
                { 12, 9 },
                { 0, 1, 2, 3, 0, 1, 2, 3 },
                { 0, 1, 2, 3, 0, 1, 2, 3 },
                { 0, 1, 2, 3, 0, 1, 2, 3 },
                { 0, 1, 2, 3, 0, 1, 2, 3 },
                { 0, 1, 2, 3, 0, 1, 2, 3, 6, 9, 10, 9 },
                { 0, 1, 2, 3, 0, 1, 2, 3, 10, 9, 9, 12 }
            };

            Matrix< DDRMat > tHostADVSensitivities;
            Matrix< DDRMat > tI;
            Node_Manager&    tNodeManager = tGeometryEngine.get_node_manager();
            for ( uint iNodeIndex = 9; iNodeIndex < 23; iNodeIndex++ )
            {
                tHostADVSensitivities.set_size( 0.0, 0.0 );
                eye( 2, 2, tI );
                tNodeManager.append_dcoordinate_dadv_from_derived_node( iNodeIndex, tHostADVSensitivities, tI );
                CHECK_EQUAL(
                        tHostADVSensitivities,
                        tIntersectionSensitivities( iNodeIndex - 9 ),
                        1E8, );
                CHECK_EQUAL(
                        tNodeManager.get_coordinate_determining_adv_ids_from_derived_node( iNodeIndex ),
                        tIntersectionIDs( iNodeIndex - 9 ), );
            }

            //------------------------------------------------------------------------------------------------------
            // Start second check
            //------------------------------------------------------------------------------------------------------

            // Set new ADVs, level set field now has no intersections
            mtk::Mesh_Pair tMeshPair( tMesh, create_integration_mesh_from_interpolation_mesh( mtk::MeshType::HMR, tMesh ) );
            tADVs = { { 0.25, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 } };
            tGeometryEngine.set_advs( tADVs );
            tGeometryEngine.reset_mesh_information( tMeshPair.get_interpolation_mesh() );

            // Solution for is_intersected() per geometry and per element
            tIsElementIntersected = {
                { false, false, false, false },    // Geometry 0
                { false, true, false, true },      // Geometry 1
                { false, true, false, true }
            };    // Geometry 2

            // Per geometry, per element, per edge
            tIsEdgeIntersected = {
                // Geometry 0
                { { false, false, false, false },            // Element 0
                        { false, false, false, false },      // Element 1
                        { false, false, false, false },      // Element 2
                        { false, false, false, false } },    // Element 3
                // Geometry 1
                { { false, false, false, false },
                        { true, false, true, false },
                        { false, false, false, false },
                        { true, false, true, false } },
                // Geometry 2
                { { false, false, false, false },
                        { true, true, true, false },
                        { false, false, false, false },
                        { true, true, true, false } }
            };

            // Intersection coordinates
            tIntersectionLocalCoordinates  = { { -0.5, 0.5, -0.5, 0.5, 1.0, 0.0, -1.0, 1.0, 0.0, -1.0 } };
            tIntersectionGlobalCoordinates = {
                { { 0.25, -1.0 } },
                { { 0.25, 0.0 } },
                { { 0.25, 0.0 } },
                { { 0.25, 1.0 } },
                { { 1.0, -1.0 } },
                { { 1.0, -1.0 } },
                { { 1.0, 0.0 } },
                { { 1.0, 0.0 } },
                { { 1.0, 0.0 } },
                { { 1.0, 1.0 } }
            };

            // Check element intersections
            tIntersectionCount = 0;
            for ( uint tGeometryIndex = 0; tGeometryIndex < 3; tGeometryIndex++ )
            {
                for ( uint tElementIndex = 0; tElementIndex < 4; tElementIndex++ )
                {
                    // Node indices per element
                    Matrix< IndexMat > tSignedNodeIndices = tMesh->get_nodes_connected_to_element_loc_inds( tElementIndex );
                    Matrix< DDUMat >   tNodeIndices( 4, 1 );
                    for ( uint iNode = 0; iNode < 4; iNode++ )
                    {
                        tNodeIndices( iNode ) = tSignedNodeIndices( iNode );
                    }

                    // Check edges for properly queued intersections
                    Matrix< DDRMat > tNodeCoordinates( 4, 2 );
                    for ( uint tNodeNumber = 0; tNodeNumber < 4; tNodeNumber++ )
                    {
                        // Node coordinates
                        Matrix< DDRMat > tFirstNodeCoordinates  = tMesh->get_node_coordinate( tNodeIndices( tNodeNumber ) );
                        Matrix< DDRMat > tSecondNodeCoordinates = tMesh->get_node_coordinate( tNodeIndices( ( tNodeNumber + 1 ) % 4 ) );

                        // Queue intersection
                        bool tIntersectionQueued = tGeometryEngine.queue_intersection(
                                tNodeIndices( tNodeNumber ),
                                tNodeIndices( ( tNodeNumber + 1 ) % 4 ),
                                tQuadParametricCoordinates( tNodeNumber ),
                                tQuadParametricCoordinates( ( tNodeNumber + 1 ) % 4 ),
                                tNodeIndices,
                                mtk::Geometry_Type::QUAD,
                                mtk::Interpolation_Order::LINEAR );
                        REQUIRE( tIntersectionQueued == tIsEdgeIntersected( tGeometryIndex )( tElementIndex )( tNodeNumber ) );

                        // Check queued intersection
                        if ( tIntersectionQueued )
                        {
                            // Check parents
                            bool tFirstParentOnInterface  = false;
                            bool tSecondParentOnInterface = false;

                            // Parent nodes on interface
                            if ( tGeometryIndex == 2 and tFirstNodeCoordinates( 0 ) == Approx( 1.0 ) )
                            {
                                tFirstParentOnInterface = true;
                            }
                            if ( tGeometryIndex == 2 and tSecondNodeCoordinates( 0 ) == Approx( 1.0 ) )
                            {
                                tSecondParentOnInterface = true;
                            }
                            CHECK( tGeometryEngine.queued_intersection_first_parent_on_interface() == tFirstParentOnInterface );
                            CHECK( tGeometryEngine.queued_intersection_second_parent_on_interface() == tSecondParentOnInterface );

                            // See if local coordinate is a number
                            real tLocalCoordinate = tGeometryEngine.get_queued_intersection_local_coordinate();
                            if ( not std::isnan( tLocalCoordinate ) )
                            {
                                // Check local coordinate
                                CHECK( tLocalCoordinate == Approx( tIntersectionLocalCoordinates( tIntersectionCount ) ) );

                                // Check global coordinates
                                CHECK_EQUAL(
                                        tGeometryEngine.get_queued_intersection_global_coordinates(),
                                        tIntersectionGlobalCoordinates( tIntersectionCount ), );
                            }

                            // Admit intersection
                            tGeometryEngine.admit_queued_intersection();
                            tIntersectionCount++;
                        }

                        // Set node coordinates for element checking
                        tNodeCoordinates.set_row( tNodeNumber, tMesh->get_node_coordinate( tNodeIndices( tNodeNumber ) ) );
                    }

                    // Check with solution
                    CHECK( tGeometryEngine.is_intersected( tGeometryIndex, tSignedNodeIndices, tNodeCoordinates ) == tIsElementIntersected( tGeometryIndex )( tElementIndex ) );
                }

                // Advance geometry index
                if ( tGeometryIndex < 2 )
                {
                    tGeometryEngine.advance_geometry_index();
                }
            }

            // Check total number of intersections
            CHECK( tIntersectionCount == 10 );

            // Test the new child nodes on the level set field (geometry 0)
            for ( uint iNodeIndex = 9; iNodeIndex < 19; iNodeIndex++ )
            {
                CHECK( tGeometryEngine.get_geometric_region( 0, iNodeIndex, { {} } ) == Geometric_Region::POSITIVE );
            }

            // Clean up
            delete tMesh;
        }
    }


    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Bilinear Intersections", "[gen], [pdv], [intersection], [bilinear intersection]" )
    {
        if ( par_size() == 1 )
        {
            // Create mesh
            mtk::Interpolation_Mesh* tMesh = create_simple_mesh( 2, 2 );

            // Set up circle
            Parameter_List tCircleParameterList = prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE );
            tCircleParameterList.set( "center_x", -0.25 );
            tCircleParameterList.set( "center_y", 0.0 );
            tCircleParameterList.set( "radius", 0.7499999999 );
            tCircleParameterList.set( "discretization_mesh_index", 0 );
            tCircleParameterList.set( "use_multilinear_interpolation", true );

            // Create geometry engine
            Geometry_Engine_Parameters tGeometryEngineParameters;
            ADV_Manager tADVManager;
            Design_Factory tDesignFactory( { tCircleParameterList }, tADVManager );
            tGeometryEngineParameters.mGeometries = tDesignFactory.get_geometries();
            Geometry_Engine_Test tGeometryEngine( tMesh, tGeometryEngineParameters );

            // Solution for is_intersected() per geometry and per element
            Vector< bool > tIsElementIntersected = { true, true, true, true };

            // Per element, per edge
            Vector< Vector< bool > > tIsEdgeIntersected = {
                { false, true, true, false },    // Element 0
                { false, false, true, true },    // Element 1
                { true, true, false, false },    // Element 2
                { true, false, false, true }
            };    // Element 3

            // Intersection coordinates
            real tFrac = 2.0 / ( 3.0 + sqrt( 17.0 ) );

            Matrix< DDRMat > tIntersectionLocalCoordinates = {
                { -tFrac, 1.0, 0.0, tFrac, -1.0, tFrac, 0.0, -tFrac, -0.5, 0.5, -0.5, 0.5 }
            };

            Vector< Matrix< DDRMat > > tIntersectionGlobalCoordinates = {
                { { 0.0, -0.5 - ( tFrac / 2.0 ) } },
                { { -1.0, 0.0 } },
                { { 0.5, 0.0 } },
                { { 0.0, -0.5 - ( tFrac / 2.0 ) } },
                { { -1.0, 0.0 } },
                { { 0.0, 0.5 + ( tFrac / 2.0 ) } },
                { { 0.5, 0.0 } },
                { { 0.0, 0.5 + ( tFrac / 2.0 ) } }
            };

            // Initialize sensitivity variables
            real             tEpsilon = 1E-12;
            Matrix< DDRMat > tHostADVSensitivities;
            Matrix< DDRMat > tI;
            eye( 2, 2, tI );

            // Check element intersections
            uint tIntersectionCount = 0;
            for ( uint tElementIndex = 0; tElementIndex < 4; tElementIndex++ )
            {
                // Get element info
                Matrix< IndexMat > tSignedNodeIndices = tMesh->get_nodes_connected_to_element_loc_inds( tElementIndex );
                Matrix< DDUMat >   tNodeIndices( 4, 1 );
                for ( uint tNodeNumber = 0; tNodeNumber < 4; tNodeNumber++ )
                {
                    tNodeIndices( tNodeNumber ) = tSignedNodeIndices( tNodeNumber );
                }

                // Check edges for properly queued intersections
                for ( uint tNodeNumber = 0; tNodeNumber < 4; tNodeNumber++ )
                {
                    // Queue intersection
                    bool tIntersectionQueued = tGeometryEngine.queue_intersection(
                            tNodeIndices( tNodeNumber ),
                            tNodeIndices( ( tNodeNumber + 1 ) % 4 ),
                            get_quad_local_coordinates( tNodeNumber ),
                            get_quad_local_coordinates( ( tNodeNumber + 1 ) % 4 ),
                            tNodeIndices,
                            mtk::Geometry_Type::QUAD,
                            mtk::Interpolation_Order::LINEAR );
                    REQUIRE( tIntersectionQueued == tIsEdgeIntersected( tElementIndex )( tNodeNumber ) );

                    // Check queued intersection
                    if ( tIntersectionQueued )
                    {
                        // Check parents
                        CHECK( not tGeometryEngine.queued_intersection_first_parent_on_interface() );
                        CHECK( not tGeometryEngine.queued_intersection_second_parent_on_interface() );

                        // Check local coordinates
                        CHECK( tGeometryEngine.get_queued_intersection_local_coordinate() == Approx( tIntersectionLocalCoordinates( tIntersectionCount ) ).margin( 1e-9 ) );

                        // Check global coordinates
                        CHECK_EQUAL( tGeometryEngine.get_queued_intersection_global_coordinates(), tIntersectionGlobalCoordinates( tIntersectionCount ), );

                        // Admit intersection
                        tGeometryEngine.admit_queued_intersection();

                        // Get node manager
                        Node_Manager& tNodeManager = tGeometryEngine.get_node_manager();

                        // Check sensitivities
                        tHostADVSensitivities.set_size( 0.0, 0.0 );
                        tNodeManager.append_dcoordinate_dadv_from_derived_node( 9 + tIntersectionCount, tHostADVSensitivities, tI );
                        Vector< sint > tADVIDs = tNodeManager.get_coordinate_determining_adv_ids_from_derived_node( 9 + tIntersectionCount );

                        // Finite difference sensitivities by queueing dummy nodes
                        Matrix< DDRMat > tFDSensitivities( tHostADVSensitivities.n_rows(), tHostADVSensitivities.n_cols(), 0.0 );
                        Vector< real > tADVs = tGeometryEngine.get_advs();
                        for ( uint iADVIndex = 0; iADVIndex < tADVIDs.size(); iADVIndex++ )
                        {
                            // Get ADV ID
                            sint tADVID = tADVIDs( iADVIndex );

                            // Positive perturbation
                            tADVs( tADVID - 1 ) += tEpsilon;
                            tGeometryEngine.set_advs( tADVs );
                            tGeometryEngine.queue_intersection(
                                    tNodeIndices( tNodeNumber ),
                                    tNodeIndices( ( tNodeNumber + 1 ) % 4 ),
                                    get_quad_local_coordinates( tNodeNumber ),
                                    get_quad_local_coordinates( ( tNodeNumber + 1 ) % 4 ),
                                    tNodeIndices,
                                    mtk::Geometry_Type::QUAD,
                                    mtk::Interpolation_Order::LINEAR );
                            Matrix< DDRMat > tPositiveGlobalCoordinates = tGeometryEngine.get_queued_intersection_global_coordinates();

                            // Negative perturbation
                            tADVs( tADVID - 1 ) -= 2.0 * tEpsilon;
                            tGeometryEngine.set_advs( tADVs );
                            tGeometryEngine.queue_intersection(
                                    tNodeIndices( tNodeNumber ),
                                    tNodeIndices( ( tNodeNumber + 1 ) % 4 ),
                                    get_quad_local_coordinates( tNodeNumber ),
                                    get_quad_local_coordinates( ( tNodeNumber + 1 ) % 4 ),
                                    tNodeIndices,
                                    mtk::Geometry_Type::QUAD,
                                    mtk::Interpolation_Order::LINEAR );
                            Matrix< DDRMat > tNegativeGlobalCoordinates = tGeometryEngine.get_queued_intersection_global_coordinates();

                            // Reset
                            tADVs( tADVID - 1 ) += tEpsilon;
                            tGeometryEngine.set_advs( tADVs );
                            tFDSensitivities( 0, iADVIndex ) = ( tPositiveGlobalCoordinates( 0 ) - tNegativeGlobalCoordinates( 0 ) ) / ( 2.0 * tEpsilon );
                            tFDSensitivities( 1, iADVIndex ) = ( tPositiveGlobalCoordinates( 1 ) - tNegativeGlobalCoordinates( 1 ) ) / ( 2.0 * tEpsilon );
                        }

                        // Check sensitivities, needs a generous error factor
                        CHECK_EQUAL( tHostADVSensitivities, tFDSensitivities, 1E12, );

                        // Increment intersection count
                        tIntersectionCount++;
                    }
                }
            }

            // Check total number of intersections
            CHECK( tIntersectionCount == 8 );

            // Test the new child nodes on the level set field
            CHECK( tGeometryEngine.get_geometric_region( 0, 9, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 10, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 11, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 12, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 13, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 14, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 15, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 16, { {} } ) == Geometric_Region::INTERFACE );

            // Get full element info for element 0
            Matrix< IndexMat >       tSignedNodeIndices = tMesh->get_nodes_connected_to_element_loc_inds( 0 );
            Matrix< DDUMat >         tNodeIndices( 4, 1 );
            Vector< Matrix< DDRMat > > tNodeCoordinates( 4 );
            for ( uint tNodeNumber = 0; tNodeNumber < 4; tNodeNumber++ )
            {
                tNodeIndices( tNodeNumber )     = tSignedNodeIndices( tNodeNumber );
                tNodeCoordinates( tNodeNumber ) = tMesh->get_node_coordinate( tNodeIndices( tNodeNumber ) );
            }

            // Queue custom intersection 1 and check for bilinear intersection
            bool tIntersectionQueued = tGeometryEngine.queue_intersection(
                    0, 2, { { -1.0, -1.0 } }, { { 1.0, 1.0 } }, tNodeIndices, mtk::Geometry_Type::QUAD, mtk::Interpolation_Order::LINEAR );
            REQUIRE( tIntersectionQueued );

            // Queue custom intersection 2 and check for no bilinear intersection
            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    1, 3, { { 1.0, -1.0 } }, { { -1.0, 1.0 } }, tNodeIndices, mtk::Geometry_Type::QUAD, mtk::Interpolation_Order::LINEAR );
            REQUIRE( not tIntersectionQueued );

            // Queue custom intersection 3 and check for bilinear intersection
            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    9, 10, { { 1.0, tFrac } }, { { -1.0, 1.0 } }, tNodeIndices, mtk::Geometry_Type::QUAD, mtk::Interpolation_Order::LINEAR );
            REQUIRE( tIntersectionQueued );

            // Clean up
            delete tMesh;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Bilinear Intersections with Nonzero Threshold", "[gen], [pdv], [intersection], [bilinear intersection], [nonzero_bilinear_threshold]" )
    {
        if ( par_size() == 1 )
        {
            // This tests a pathological case where the curved bi-linear interpolation results in an edge with the same signed parent nodes to have an intersection point
            // My mesh is one element where the geometry is a user defined nodal field (best way to control everything in this problem).
            // The isocontour threshold is away from 0. If this test fails, it is probably because someone modified the logic in cl_GEN_Intersection_Node.cpp which determins
            // if an entity is intersected.

            // Generate mesh
            uint             aNumDim = 2;
            Matrix< DDRMat > aCoords( 6, 2 );
            aCoords( 0, 0 ) = 0.0, aCoords( 0, 1 ) = 0.0;
            aCoords( 1, 0 ) = 0.25, aCoords( 1, 1 ) = 0.0;
            aCoords( 2, 0 ) = 0.25, aCoords( 2, 1 ) = 0.25;
            aCoords( 3, 0 ) = 0.0, aCoords( 3, 1 ) = 1.0;
            aCoords( 4, 0 ) = 0.0, aCoords( 4, 1 ) = 0.0;
            Matrix< IdMat > aElemConn( 1, 4 );

            // 0D to 3D connectivity (node to element)
            aElemConn( 0, 0 ) = 1;
            aElemConn( 0, 1 ) = 2;
            aElemConn( 0, 2 ) = 3;
            aElemConn( 0, 3 ) = 4;

            Matrix< IdMat > aElemLocaltoGlobal = { { 1 } };

            // Create MORIS mesh using MTK database
            moris::mtk::MtkMeshData aMeshData;
            aMeshData.CreateAllEdgesAndFaces    = true;
            aMeshData.SpatialDim                = &aNumDim;
            aMeshData.ElemConn( 0 )             = &aElemConn;
            aMeshData.NodeCoords                = &aCoords;
            aMeshData.LocaltoGlobalElemMap( 0 ) = &aElemLocaltoGlobal;
            moris::mtk::Scalar_Field_Info< DDRMat > tLSF;
            std::string                             tLSFName = "lsf1";
            tLSF.set_field_name( tLSFName );
            tLSF.set_field_entity_rank( mtk::EntityRank::NODE );

            // Add to mesh input field container
            moris::mtk::MtkFieldsInfo tFieldsInfo;
            add_field_for_mesh_input( &tLSF, tFieldsInfo );
            aMeshData.FieldsInfo = &tFieldsInfo;

            moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( mtk::MeshType::STK, aMeshData );

            moris::uint                    tNumNodes = tMeshData->get_num_entities( mtk::EntityRank::NODE );
            moris::Matrix< moris::DDRMat > tLevelsetVal( tNumNodes, 1, -1.3 );

            moris_id tIndexOfNodeId1 = tMeshData->get_loc_entity_ind_from_entity_glb_id( 1, mtk::EntityRank::NODE );
            moris_id tIndexOfNodeId2 = tMeshData->get_loc_entity_ind_from_entity_glb_id( 2, mtk::EntityRank::NODE );
            moris_id tIndexOfNodeId3 = tMeshData->get_loc_entity_ind_from_entity_glb_id( 3, mtk::EntityRank::NODE );
            moris_id tIndexOfNodeId4 = tMeshData->get_loc_entity_ind_from_entity_glb_id( 4, mtk::EntityRank::NODE );

            tLevelsetVal( tIndexOfNodeId1 ) = 5.71398452074828e-01;
            tLevelsetVal( tIndexOfNodeId2 ) = 5.06776630434012e-01;
            tLevelsetVal( tIndexOfNodeId3 ) = 4.25145951405591e-01;
            tLevelsetVal( tIndexOfNodeId4 ) = 5.00394812706282e-01;

            tMeshData->add_mesh_field_real_scalar_data_loc_inds( tLSFName, mtk::EntityRank::NODE, tLevelsetVal );

            Vector< std::shared_ptr< gen::Geometry > > tGeometry( 1 );
            Level_Set_Parameters                     tLevelSetParameters;
            tLevelSetParameters.mUseMultilinearInterpolation = true;
            tLevelSetParameters.mIsocontourThreshold         = 0.5;
            tLevelSetParameters.mIsocontourTolerance         = 1E-13;
            auto tField                                      = std::make_shared< Mesh_Field >( tMeshData, tLSFName );
            tGeometry( 0 )                                   = std::make_shared< Level_Set_Geometry >( tField, tLevelSetParameters );

            Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometry;

            Geometry_Engine tGeometryEngine( tMeshData, tGeometryEngineParameters );

            // get the cell
            moris::mtk::Cell& tCell = tMeshData->get_mtk_cell( 0 );

            // verify that it says that the cell is intersected
            bool tIsIntersected = tGeometryEngine.is_intersected( 0, tCell.get_vertex_inds(), tCell.get_vertex_coords() );

            CHECK( tIsIntersected );

            Vector< Matrix< IndexMat > > tVertexIndices = { tCell.get_vertex_inds() };

            Vector< Matrix< DDRMat > > tLocalCoords( 1 );

            tLocalCoords( 0 ) = { { +0.000000000000000e+00, +0.000000000000000e+00 } };

            tGeometryEngine.create_new_derived_nodes(
                    tVertexIndices,
                    tLocalCoords,
                    mtk::Geometry_Type::QUAD,
                    mtk::Interpolation_Order::LINEAR );

            Matrix< DDRMat >   tVertexCoords  = tCell.get_vertex_coords();
            Matrix< IndexMat > tVertexInds    = tCell.get_vertex_inds();
            Matrix< DDUMat >   tVertexIndsDDU = { { (uint)tVertexInds( 0 ), (uint)tVertexInds( 1 ), (uint)tVertexInds( 2 ), (uint)tVertexInds( 3 ) } };

            Vector< Matrix< DDRMat > > tBGCellCoords( 4 );
            tBGCellCoords( 0 ) = tVertexCoords.get_row( 0 );
            tBGCellCoords( 1 ) = tVertexCoords.get_row( 1 );
            tBGCellCoords( 2 ) = tVertexCoords.get_row( 2 );
            tBGCellCoords( 3 ) = tVertexCoords.get_row( 3 );

            Matrix< DDRMat > tLocalCoordsMat = {
                { -1.000000000000000e+00, -1.000000000000000e+00 },
                { +1.000000000000000e+00, -1.000000000000000e+00 },
                { +1.000000000000000e+00, +1.000000000000000e+00 },
                { -1.000000000000000e+00, +1.000000000000000e+00 },
                { +0.000000000000000e+00, +0.000000000000000e+00 }
            };

            // check that the cell is intersected
            moris_index tNodeIndex1 = 0;
            moris_index tNodeIndex2 = 1;

            bool tIntersectionQueued = tGeometryEngine.queue_intersection(
                    tNodeIndex1,
                    tNodeIndex2,
                    tLocalCoordsMat.get_row( tNodeIndex1 ),
                    tLocalCoordsMat.get_row( tNodeIndex2 ),
                    tVertexIndsDDU,
                    mtk::Geometry_Type::QUAD,
                    mtk::Interpolation_Order::LINEAR );

            CHECK( !tIntersectionQueued );

            // check that the cell is intersected
            tNodeIndex1 = 1;
            tNodeIndex2 = 2;

            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    tNodeIndex1,
                    tNodeIndex2,
                    tLocalCoordsMat.get_row( tNodeIndex1 ),
                    tLocalCoordsMat.get_row( tNodeIndex2 ),
                    tVertexIndsDDU,
                    mtk::Geometry_Type::QUAD,
                    mtk::Interpolation_Order::LINEAR );

            CHECK( tIntersectionQueued );

            tNodeIndex1 = 2;
            tNodeIndex2 = 3;

            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    tNodeIndex1,
                    tNodeIndex2,
                    tLocalCoordsMat.get_row( tNodeIndex1 ),
                    tLocalCoordsMat.get_row( tNodeIndex2 ),
                    tVertexIndsDDU,
                    mtk::Geometry_Type::QUAD,
                    mtk::Interpolation_Order::LINEAR );
            CHECK( tIntersectionQueued );

            tNodeIndex1 = 0;
            tNodeIndex2 = 3;

            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    tNodeIndex1,
                    tNodeIndex2,
                    tLocalCoordsMat.get_row( tNodeIndex1 ),
                    tLocalCoordsMat.get_row( tNodeIndex2 ),
                    tVertexIndsDDU,
                    mtk::Geometry_Type::QUAD,
                    mtk::Interpolation_Order::LINEAR );
            CHECK( !tIntersectionQueued );

            // check that the cell is intersected
            tNodeIndex1 = 0;
            tNodeIndex2 = 4;

            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    tNodeIndex1,
                    tNodeIndex2,
                    tLocalCoordsMat.get_row( tNodeIndex1 ),
                    tLocalCoordsMat.get_row( tNodeIndex2 ),
                    tVertexIndsDDU,
                    mtk::Geometry_Type::QUAD,
                    mtk::Interpolation_Order::LINEAR );

            CHECK( !tIntersectionQueued );

            // check that the cell is intersected
            tNodeIndex1 = 1;
            tNodeIndex2 = 4;

            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    tNodeIndex1,
                    tNodeIndex2,
                    tLocalCoordsMat.get_row( tNodeIndex1 ),
                    tLocalCoordsMat.get_row( tNodeIndex2 ),
                    tVertexIndsDDU,
                    mtk::Geometry_Type::QUAD,
                    mtk::Interpolation_Order::LINEAR );

            CHECK( !tIntersectionQueued );

            // check that the cell is intersected
            tNodeIndex1 = 2;
            tNodeIndex2 = 4;

            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    tNodeIndex1,
                    tNodeIndex2,
                    tLocalCoordsMat.get_row( tNodeIndex1 ),
                    tLocalCoordsMat.get_row( tNodeIndex2 ),
                    tVertexIndsDDU,
                    mtk::Geometry_Type::QUAD,
                    mtk::Interpolation_Order::LINEAR );

            CHECK( tIntersectionQueued );

            // check that the cell is intersected
            tNodeIndex1 = 3;
            tNodeIndex2 = 4;

            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    tNodeIndex1,
                    tNodeIndex2,
                    tLocalCoordsMat.get_row( tNodeIndex1 ),
                    tLocalCoordsMat.get_row( tNodeIndex2 ),
                    tVertexIndsDDU,
                    mtk::Geometry_Type::QUAD,
                    mtk::Interpolation_Order::LINEAR );
            CHECK( !tIntersectionQueued );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::gen
