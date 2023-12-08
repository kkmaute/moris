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
#include "cl_GEN_Pdv_Host_Manager.hpp"
#include "cl_GEN_BSpline_Field.hpp"
#include "cl_GEN_Design_Factory.hpp"
#include "fn_check_equal.hpp"
#include "fn_GEN_create_simple_mesh.hpp"
#include "cl_GEN_Mesh_Field.hpp"
#include "cl_MTK_Mesh_Factory.hpp"

#include "fn_PRM_GEN_Parameters.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    static Cell< Matrix< DDRMat > > tQuadParametricCoordinates = {
        {{ -1.0, -1.0 }},
        {{ 1.0, -1.0 }},
        {{ 1.0, 1.0 }},
        {{ -1.0, 1.0 }}
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
            Matrix< DDRMat > tADVs = { { 0.25, 0.0, 1.0, 0.0 } };

            // Circle
            ParameterList tCircleParameterList = prm::create_level_set_geometry_parameter_list();
            tCircleParameterList.set( "field_type", "circle" );
            tCircleParameterList.set( "constant_parameters", "-0.25, 0.0, 0.7499999999" );
            tCircleParameterList.set( "discretization_mesh_index", 0 );

            // Plane 1
            ParameterList tPlane1ParameterList = prm::create_level_set_geometry_parameter_list();
            tPlane1ParameterList.set( "field_type", "plane" );
            tPlane1ParameterList.set( "field_variable_indices", "all" );
            tPlane1ParameterList.set( "adv_indices", "all" );

            // Plane 2
            ParameterList tPlane2ParameterList = prm::create_level_set_geometry_parameter_list();
            tPlane2ParameterList.set( "field_type", "plane" );
            tPlane2ParameterList.set( "constant_parameters", "1.0, 0.0, 1.0, 0.0" );

            // Create geometry engine
            Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mADVs = tADVs;
            Design_Factory tDesignFactory( { tCircleParameterList, tPlane1ParameterList, tPlane2ParameterList }, tADVs );
            tGeometryEngineParameters.mGeometries = tDesignFactory.get_geometries();
            Geometry_Engine tGeometryEngine( tMesh, tGeometryEngineParameters );

            // TODO ensure this writes the mesh/fields correctly instead of just relying on no errors being thrown
            tGeometryEngine.output_fields_on_mesh( tMesh, "intersection_test.exo" );

            // Solution for is_intersected() per geometry and per element
            Cell< Cell< bool > > tIsElementIntersected = { //
                { true, true, true, true },
                { false, true, false, true },
                { false, true, false, true }
            };

            // Per geometry, per element, per edge
            Cell< Cell< Cell< bool > > > tIsEdgeIntersected = {
                { { false, true, true, false },          // Geometry 0, Element 0
                        { false, false, true, true },    // Geometry 0, Element 1
                        { true, true, false, false },    // Geometry 0, Element 2
                        { true, false, false, true } },
                {                                          // Geometry 0, Element 3
                        { false, false, false, false },    // Geometry 1, Element 0
                        { true, false, true, false },      // Geometry 1, Element 1
                        { false, false, false, false },    // Geometry 1, Element 2
                        { true, false, true, false } },
                {                                          // Geometry 1, Element 3
                        { false, false, false, false },    // Geometry 2, Element 0
                        { true, true, true, false },       // Geometry 2, Element 1
                        { false, false, false, false },    // Geometry 2, Element 2
                        { true, true, true, false } }      // Geometry 2, Element 3
            };

            // Intersection coordinates
            real tFrac = 2.0 / ( 3.0 + sqrt( 17.0 ) );

            Matrix< DDRMat > tIntersectionLocalCoordinates = {
                { -tFrac, 1.0, 0.0, tFrac, -1.0, tFrac, 0.0, -tFrac, -0.5, 0.5, -0.5, 0.5, 0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 0.0, -1.0 }
            };

            Cell< Matrix< DDRMat > > tIntersectionGlobalCoordinates = {
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
                    Matrix< DDUMat > tNodeIndices( 4, 1 );
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
                                tNodeIndices );
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
                            if ( not isnan( tLocalCoordinate ) )
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
                    bool tIsIntersected = tGeometryEngine.is_intersected( tSignedNodeIndices, tNodeCoordinates );
                    CHECK( tIsIntersected == tIsElementIntersected( tGeometryIndex )( tElementIndex ) );
                }

                // Intersection on intersection
                if ( tGeometryIndex == 1 )
                {
                    // Queue intersection on intersection 1
                    bool tIntersectionQueued = tGeometryEngine.queue_intersection(
                            9,
                            11,
                            {{ -1.0, -tFrac }},
                            {{ 0.0, 1.0 }},
                            {{ 1, 4, 5, 2 }} );

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
                            {{ 0.0, -1.0 }},
                            {{ -1.0, tFrac }},
                            {{ 2, 5, 8, 6 }} );

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
            CHECK( tGeometryEngine.get_field_value( 0, 9, { {} } ) == Approx( 0.0 ).margin( MORIS_REAL_EPS ) );
            CHECK( tGeometryEngine.get_field_value( 0, 10, { {} } ) == Approx( 0.0 ).margin( MORIS_REAL_EPS ) );
            CHECK( tGeometryEngine.get_field_value( 0, 11, { {} } ) == Approx( 0.0 ).margin( MORIS_REAL_EPS ) );
            CHECK( tGeometryEngine.get_field_value( 0, 12, { {} } ) == Approx( 0.0 ).margin( MORIS_REAL_EPS ) );
            CHECK( tGeometryEngine.get_field_value( 0, 13, { {} } ) == Approx( 0.0 ).margin( MORIS_REAL_EPS ) );
            CHECK( tGeometryEngine.get_field_value( 0, 14, { {} } ) == Approx( 0.0 ).margin( MORIS_REAL_EPS ) );
            CHECK( tGeometryEngine.get_field_value( 0, 15, { {} } ) == Approx( 0.0 ).margin( MORIS_REAL_EPS ) );
            CHECK( tGeometryEngine.get_field_value( 0, 16, { {} } ) == Approx( 0.0 ).margin( MORIS_REAL_EPS ) );
            CHECK( tGeometryEngine.get_field_value( 0, 17, { {} } ) == Approx( 0.423278 ).margin( MORIS_REAL_EPS ) );
            CHECK( tGeometryEngine.get_field_value( 0, 18, { {} } ) == Approx( -0.25 ).margin( MORIS_REAL_EPS ) );
            CHECK( tGeometryEngine.get_field_value( 0, 19, { {} } ) == Approx( -0.25 ).margin( MORIS_REAL_EPS ) );
            CHECK( tGeometryEngine.get_field_value( 0, 20, { {} } ) == Approx( 0.423278 ).margin( MORIS_REAL_EPS ) );
            CHECK( tGeometryEngine.get_field_value( 0, 21, { {} } ) == Approx( 0.0 ).margin( MORIS_REAL_EPS ) );
            CHECK( tGeometryEngine.get_field_value( 0, 22, { {} } ) == Approx( 0.0 ).margin( MORIS_REAL_EPS ) );
            CHECK( tGeometryEngine.get_field_value( 0, 23, { {} } ) == Approx( ( sqrt( 41.0 ) - 3.0 ) / 4.0 ).margin( MORIS_REAL_EPS ) );
            CHECK( tGeometryEngine.get_field_value( 0, 25, { {} } ) == Approx( 0.5 ).margin( MORIS_REAL_EPS ) );
            CHECK( tGeometryEngine.get_field_value( 0, 26, { {} } ) == Approx( 0.5 ).margin( MORIS_REAL_EPS ) );
            CHECK( tGeometryEngine.get_field_value( 0, 28, { {} } ) == Approx( ( sqrt( 41.0 ) - 3.0 ) / 4.0 ).margin( MORIS_REAL_EPS ) );

            // Get the PDV host manager and set the number of total nodes
            auto tPDVHostManager = dynamic_cast< Pdv_Host_Manager* >( tGeometryEngine.get_design_variable_interface() );

            // Test that the new intersections have been added to the PDV host manager, but ONLY for the circle
            Cell< Matrix< DDRMat > > tPdvValues( 0 );
            Cell< Matrix< DDSMat > > tIsActive( 0 );
            tPDVHostManager->get_ig_pdv_value(
                    { { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28 } },
                    { PDV_Type::X_COORDINATE, PDV_Type::Y_COORDINATE },
                    tPdvValues,
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
                CHECK( tPdvValues( 0 )( tNodeIndex ) == Approx( tIntersectionGlobalCoordinates( tNodeIndex - 9 )( 0 ) ) );
                CHECK( tPdvValues( 1 )( tNodeIndex ) == Approx( tIntersectionGlobalCoordinates( tNodeIndex - 9 )( 1 ) ) );
            }

            // Nodes on the second plane interface (inactive)
            for ( uint tNodeIndex = 23; tNodeIndex < 27; tNodeIndex++ )
            {
                // Check if not active
                CHECK( tIsActive( 0 )( tNodeIndex ) == false );
                CHECK( tIsActive( 1 )( tNodeIndex ) == false );
            }

            // Check sensitivities
            Cell< Matrix< DDRMat > > tIntersectionSensitivities = {
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

            Cell< Matrix< DDSMat > > tIntersectionIDs = {
                { { 6, 9 } },
                { { 9, 8 } },
                { { 10, 9 } },
                { { 9, 6 } },
                { { 8, 9 } },
                { { 9, 12 } },
                { { 9, 10 } },
                { { 12, 9 } },
                { { 0, 1, 2, 3, 0, 1, 2, 3 } },
                { { 0, 1, 2, 3, 0, 1, 2, 3 } },
                { { 0, 1, 2, 3, 0, 1, 2, 3 } },
                { { 0, 1, 2, 3, 0, 1, 2, 3 } },
                { { 0, 1, 2, 3, 0, 1, 2, 3, 6, 9, 10, 9 } },
                { { 0, 1, 2, 3, 0, 1, 2, 3, 10, 9, 9, 12 } }
            };

            Matrix< DDRMat > tHostADVSensitivities;
            Matrix< DDRMat > tI;
            for ( uint tNodeIndex = 9; tNodeIndex < 23; tNodeIndex++ )
            {
                tHostADVSensitivities.set_size( 0.0, 0.0 );
                eye( 2, 2, tI );
                tPDVHostManager->get_intersection_node( tNodeIndex )->append_dcoordinate_dadv( tHostADVSensitivities, tI );
                CHECK_EQUAL(
                        tHostADVSensitivities,
                        tIntersectionSensitivities( tNodeIndex - 9 ),
                        1E8, );
                CHECK_EQUAL(
                        tPDVHostManager->get_intersection_node( tNodeIndex )->get_coordinate_determining_adv_ids(),
                        tIntersectionIDs( tNodeIndex - 9 ), );
            }

            //------------------------------------------------------------------------------------------------------
            // Start second check
            //------------------------------------------------------------------------------------------------------

            // Set new ADVs, level set field now has no intersections
            mtk::Mesh_Pair tMeshPair( tMesh, create_integration_mesh_from_interpolation_mesh( mtk::MeshType::HMR, tMesh ) );
            tGeometryEngine.set_advs( Matrix< DDRMat >( 16, 1, 1.0 ) );
            tGeometryEngine.distribute_advs( tMeshPair );

            // Solution for is_intersected() per geometry and per element
            tIsElementIntersected = { { false, false, false, false }, { false, true, false, true }, { false, true, false, true } };

            // Per geometry, per element, per edge
            tIsEdgeIntersected = {
                { { false, false, false, false },          // Geometry 0, Element 0
                        { false, false, false, false },    // Geometry 0, Element 1
                        { false, false, false, false },    // Geometry 0, Element 2
                        { false, false, false, false } },
                {                                          // Geometry 0, Element 3
                        { false, false, false, false },    // Geometry 1, Element 0
                        { true, false, true, false },      // Geometry 1, Element 1
                        { false, false, false, false },    // Geometry 1, Element 2
                        { true, false, true, false } },
                {                                          // Geometry 1, Element 3
                        { false, false, false, false },    // Geometry 2, Element 0
                        { true, true, true, false },       // Geometry 2, Element 1
                        { false, false, false, false },    // Geometry 2, Element 2
                        { true, true, true, false } }      // Geometry 2, Element 3
            };

            // Intersection coordinates
            tIntersectionLocalCoordinates = {
                { -0.5, 0.5, -0.5, 0.5, 1.0, 0.0, -1.0, 1.0, 0.0, -1.0 }
            };

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
                    Matrix< DDUMat > tNodeIndices( 4, 1 );
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
                                tFirstNodeCoordinates,
                                tSecondNodeCoordinates,
                                tNodeIndices );
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
                            if ( not isnan( tLocalCoordinate ) )
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
                    CHECK( tGeometryEngine.is_intersected( tSignedNodeIndices, tNodeCoordinates ) == tIsElementIntersected( tGeometryIndex )( tElementIndex ) );
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
            CHECK( tGeometryEngine.get_field_value( 0, 9, { {} } ) == Approx( 1.0 ) );
            CHECK( tGeometryEngine.get_field_value( 0, 10, { {} } ) == Approx( 1.0 ) );
            CHECK( tGeometryEngine.get_field_value( 0, 11, { {} } ) == Approx( 1.0 ) );
            CHECK( tGeometryEngine.get_field_value( 0, 12, { {} } ) == Approx( 1.0 ) );

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

            Cell< std::shared_ptr< ge::Level_Set_Geometry > > tGeometry( 1 );
            Level_Set_Parameters tLevelSetParameters;
            tLevelSetParameters.mUseMultilinearInterpolation = true;
            tLevelSetParameters.mIsocontourThreshold = 0.5;
            tLevelSetParameters.mIsocontourTolerance = 1E-13;
            auto tField = std::make_shared< Mesh_Field >( tMeshData, tLSFName );
            tGeometry( 0 ) = std::make_shared< Level_Set_Geometry >( tField, tLevelSetParameters );
            
            Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometry;

            Geometry_Engine tGeometryEngine( tMeshData, tGeometryEngineParameters );

            // get the cell
            moris::mtk::Cell& tCell = tMeshData->get_mtk_cell( 0 );

            // verify that it says that the cell is intersected
            bool tIsIntersected = tGeometryEngine.is_intersected( tCell.get_vertex_inds(), tCell.get_vertex_coords() );

            CHECK( tIsIntersected );

            moris::Cell< Matrix< IndexMat > > tVertexIndices = { tCell.get_vertex_inds() };

            moris::Cell< Matrix< DDRMat > > tLocalCoords( 1 );

            tLocalCoords( 0 ) = { { +0.000000000000000e+00, +0.000000000000000e+00 } };

            Cell< moris_index > tNewNodeIndices = { 4 };

            tGeometryEngine.create_new_derived_nodes( tNewNodeIndices, tVertexIndices, tLocalCoords );

            Matrix< DDRMat >   tVertexCoords  = tCell.get_vertex_coords();
            Matrix< IndexMat > tVertexInds    = tCell.get_vertex_inds();
            Matrix< DDUMat >   tVertexIndsDDU = { { (uint)tVertexInds( 0 ), (uint)tVertexInds( 1 ), (uint)tVertexInds( 2 ), (uint)tVertexInds( 3 ) } };

            Cell< Matrix< DDRMat > > tBGCellCoords( 4 );
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
                    tVertexIndsDDU );

            CHECK( !tIntersectionQueued );

            // check that the cell is intersected
            tNodeIndex1 = 1;
            tNodeIndex2 = 2;

            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    tNodeIndex1,
                    tNodeIndex2,
                    tLocalCoordsMat.get_row( tNodeIndex1 ),
                    tLocalCoordsMat.get_row( tNodeIndex2 ),
                    tVertexIndsDDU );

            CHECK( tIntersectionQueued );

            tNodeIndex1 = 2;
            tNodeIndex2 = 3;

            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    tNodeIndex1,
                    tNodeIndex2,
                    tLocalCoordsMat.get_row( tNodeIndex1 ),
                    tLocalCoordsMat.get_row( tNodeIndex2 ),
                    tVertexIndsDDU );
            CHECK( tIntersectionQueued );

            tNodeIndex1 = 0;
            tNodeIndex2 = 3;

            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    tNodeIndex1,
                    tNodeIndex2,
                    tLocalCoordsMat.get_row( tNodeIndex1 ),
                    tLocalCoordsMat.get_row( tNodeIndex2 ),
                    tVertexIndsDDU );
            CHECK( !tIntersectionQueued );

            // check that the cell is intersected
            tNodeIndex1 = 0;
            tNodeIndex2 = 4;

            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    tNodeIndex1,
                    tNodeIndex2,
                    tLocalCoordsMat.get_row( tNodeIndex1 ),
                    tLocalCoordsMat.get_row( tNodeIndex2 ),
                    tVertexIndsDDU );

            CHECK( !tIntersectionQueued );

            // check that the cell is intersected
            tNodeIndex1 = 1;
            tNodeIndex2 = 4;

            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    tNodeIndex1,
                    tNodeIndex2,
                    tLocalCoordsMat.get_row( tNodeIndex1 ),
                    tLocalCoordsMat.get_row( tNodeIndex2 ),
                    tVertexIndsDDU );

            CHECK( !tIntersectionQueued );

            // check that the cell is intersected
            tNodeIndex1 = 2;
            tNodeIndex2 = 4;

            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    tNodeIndex1,
                    tNodeIndex2,
                    tLocalCoordsMat.get_row( tNodeIndex1 ),
                    tLocalCoordsMat.get_row( tNodeIndex2 ),
                    tVertexIndsDDU );

            CHECK( tIntersectionQueued );

            // check that the cell is intersected
            tNodeIndex1 = 3;
            tNodeIndex2 = 4;

            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    tNodeIndex1,
                    tNodeIndex2,
                    tLocalCoordsMat.get_row( tNodeIndex1 ),
                    tLocalCoordsMat.get_row( tNodeIndex2 ),
                    tVertexIndsDDU );
            CHECK( !tIntersectionQueued );
        }
    }

    TEST_CASE( "Bilinear Intersections", "[gen], [pdv], [intersection], [bilinear intersection]" )
    {
        if ( par_size() == 1 )
        {
            // Create mesh
            mtk::Interpolation_Mesh* tMesh = create_simple_mesh( 2, 2 );

            // Set up circle
            ParameterList tCircleParameterList = prm::create_level_set_geometry_parameter_list();
            tCircleParameterList.set( "field_type", "circle" );
            tCircleParameterList.set( "constant_parameters", "-0.25, 0.0, 0.7499999999" );
            tCircleParameterList.set( "discretization_mesh_index", 0 );
            tCircleParameterList.set( "use_multilinear_interpolation", true );
            Matrix< DDRMat > tADVs( 0, 0 );

            // Create geometry engine
            Geometry_Engine_Parameters tGeometryEngineParameters;
            Design_Factory tDesignFactory( { tCircleParameterList }, tADVs );
            tGeometryEngineParameters.mGeometries = tDesignFactory.get_geometries();
            Geometry_Engine tGeometryEngine( tMesh, tGeometryEngineParameters );

            // Solution for is_intersected() per geometry and per element
            Cell< bool > tIsElementIntersected = { true, true, true, true };

            // Per element, per edge
            Cell< Cell< bool > > tIsEdgeIntersected = {
                { false, true, true, false },    // Element 0
                { false, false, true, true },    // Element 1
                { true, true, false, false },    // Element 2
                { true, false, false, true }     // Element 3
            };

            // Intersection coordinates
            real tFrac = 2.0 / ( 3.0 + sqrt( 17.0 ) );

            Matrix< DDRMat > tIntersectionLocalCoordinates = {
                { -tFrac, 1.0, 0.0, tFrac, -1.0, tFrac, 0.0, -tFrac, -0.5, 0.5, -0.5, 0.5 }
            };

            Cell< Matrix< DDRMat > > tIntersectionGlobalCoordinates = {
                { { 0.0, -0.5 - ( tFrac / 2.0 ) } },
                { { -1.0, 0.0 } },
                { { 0.5, 0.0 } },
                { { 0.0, -0.5 - ( tFrac / 2.0 ) } },
                { { -1.0, 0.0 } },
                { { 0.0, 0.5 + ( tFrac / 2.0 ) } },
                { { 0.5, 0.0 } },
                { { 0.0, 0.5 + ( tFrac / 2.0 ) } }
            };

            // Check element intersections
            uint tIntersectionCount = 0;
            for ( uint tElementIndex = 0; tElementIndex < 4; tElementIndex++ )
            {
                // Get element info
                Matrix< IndexMat > tSignedNodeIndices = tMesh->get_nodes_connected_to_element_loc_inds( tElementIndex );
                Matrix< DDUMat >   tNodeIndices( 4, 1 );

                Cell< Matrix< DDRMat > > tNodeCoordinates( 4 );

                for ( uint tNodeNumber = 0; tNodeNumber < 4; tNodeNumber++ )
                {
                    tNodeIndices( tNodeNumber )     = tSignedNodeIndices( tNodeNumber );
                    tNodeCoordinates( tNodeNumber ) = tMesh->get_node_coordinate( tNodeIndices( tNodeNumber ) );
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
                            tNodeIndices );
                    REQUIRE( tIntersectionQueued == tIsEdgeIntersected( tElementIndex )( tNodeNumber ) );

                    // Check queued intersection
                    if ( tIntersectionQueued )
                    {
                        // Check parents
                        bool tFirstParentOnInterface  = false;
                        bool tSecondParentOnInterface = false;

                        CHECK( tGeometryEngine.queued_intersection_first_parent_on_interface() == tFirstParentOnInterface );
                        CHECK( tGeometryEngine.queued_intersection_second_parent_on_interface() == tSecondParentOnInterface );

                        // Check local coordinates
                        CHECK( tGeometryEngine.get_queued_intersection_local_coordinate() ==    //
                                Approx( tIntersectionLocalCoordinates( tIntersectionCount ) ).margin( 1e-9 ) );

                        // Check global coordinates
                        CHECK( tGeometryEngine.get_queued_intersection_global_coordinates()( 0 ) ==    //
                                Approx( tIntersectionGlobalCoordinates( tIntersectionCount )( 0 ) ).margin( 1e-9 ) );

                        CHECK( tGeometryEngine.get_queued_intersection_global_coordinates()( 1 ) ==    //
                                Approx( tIntersectionGlobalCoordinates( tIntersectionCount )( 1 ) ).margin( 1e-9 ) );

                        // Admit intersection
                        tGeometryEngine.admit_queued_intersection();

                        // Increment intersection count
                        tIntersectionCount++;
                    }
                }
            }

            // Check total number of intersections
            CHECK( tIntersectionCount == 8 );

            // Test the new child nodes on the level set field
            CHECK( tGeometryEngine.get_field_value( 0, 9, { {} } ) == Approx( 0.0 ).margin( MORIS_REAL_EPS ) );
            CHECK( tGeometryEngine.get_field_value( 0, 10, { {} } ) == Approx( 0.0 ).margin( MORIS_REAL_EPS ) );
            CHECK( tGeometryEngine.get_field_value( 0, 11, { {} } ) == Approx( 0.0 ).margin( MORIS_REAL_EPS ) );
            CHECK( tGeometryEngine.get_field_value( 0, 12, { {} } ) == Approx( 0.0 ).margin( MORIS_REAL_EPS ) );
            CHECK( tGeometryEngine.get_field_value( 0, 13, { {} } ) == Approx( 0.0 ).margin( MORIS_REAL_EPS ) );
            CHECK( tGeometryEngine.get_field_value( 0, 14, { {} } ) == Approx( 0.0 ).margin( MORIS_REAL_EPS ) );
            CHECK( tGeometryEngine.get_field_value( 0, 15, { {} } ) == Approx( 0.0 ).margin( MORIS_REAL_EPS ) );
            CHECK( tGeometryEngine.get_field_value( 0, 16, { {} } ) == Approx( 0.0 ).margin( MORIS_REAL_EPS ) );

            // Get full element info for element 0
            Matrix< IndexMat > tSignedNodeIndices = tMesh->get_nodes_connected_to_element_loc_inds( 0 );
            Matrix< DDUMat >   tNodeIndices( 4, 1 );
            Cell< Matrix< DDRMat > > tNodeCoordinates( 4 );
            for ( uint tNodeNumber = 0; tNodeNumber < 4; tNodeNumber++ )
            {
                tNodeIndices( tNodeNumber )     = tSignedNodeIndices( tNodeNumber );
                tNodeCoordinates( tNodeNumber ) = tMesh->get_node_coordinate( tNodeIndices( tNodeNumber ) );
            }

            // Queue custom intersection 1 and check for bilinear intersection
            bool tIntersectionQueued = tGeometryEngine.queue_intersection(
                    0, 2, { { -1.0, -1.0 } }, { { 1.0, 1.0 } }, tNodeIndices );
            REQUIRE( tIntersectionQueued );

            // Queue custom intersection 2 and check for no bilinear intersection
            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    1, 3, { { 1.0, -1.0 } }, { { -1.0, 1.0 } }, tNodeIndices );
            REQUIRE( not tIntersectionQueued );

            // Queue custom intersection 3 and check for bilinear intersection
            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    0, 0, { { 0.75, 0.0 } }, { { -0.75, 0.0 } }, tNodeIndices );
            REQUIRE( tIntersectionQueued );

            // Clean up
            delete tMesh;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}
