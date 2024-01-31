/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_GEN_PDV.cpp
 *
 */

#include "catch.hpp"
#include <cmath>
#include "cl_Matrix.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_MSI_Design_Variable_Interface.hpp"

#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Design_Factory.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "cl_GEN_Background_Node.hpp"
#include "cl_GEN_Parent_Node.hpp"

#define protected public
#define private public
#include "cl_GEN_PDV_Host_Manager.hpp"
#include "cl_GEN_Intersection_Node_Linear.hpp"
#undef protected
#undef private

#include "fn_check_equal.hpp"

namespace moris::ge
{

    class Node_Manager_Test : public Node_Manager
    {
      public:
        Node_Manager_Test()
                : Node_Manager( nullptr )
        {
            mMeshGiven = true;
        }
    };

    //------------------------------------------------------------------------------------------------------------------

    // Dummy IQI sensitivity values so a FEM model doesn't have to be created
    uint             gNumPDVs = 16;
    Matrix< DDRMat > gdIQIdPDV1( 1, gNumPDVs );
    Matrix< DDRMat > gdIQIdPDV2( 1, gNumPDVs );

    class PDV_Host_Manager_Test : public PDV_Host_Manager
    {
      public:
        explicit PDV_Host_Manager_Test( Node_Manager& aNodeManager = Node_Manager::get_trivial_instance() )
                : PDV_Host_Manager( aNodeManager )
        {
        }

        sol::Dist_Vector*
        get_dQIdp() override
        {
            // Factory
            sol::Matrix_Vector_Factory tDistributedFactory;

            // IQI/PDV sensitivities
            Matrix< DDSMat > tFullPDVIds( gNumPDVs, 1 );
            if ( par_rank() == 0 )
            {
                for ( uint tPDVIndex = 0; tPDVIndex < gNumPDVs; tPDVIndex++ )
                {
                    tFullPDVIds( tPDVIndex ) = tPDVIndex;
                    gdIQIdPDV1( tPDVIndex )  = 1.0;
                    gdIQIdPDV2( tPDVIndex )  = (real)tPDVIndex;
                }
            }

            // PDV IDs
            Matrix< DDSMat > tOwnedPDVIds = this->get_my_local_global_map();

            // IQI sensitivity vector
            sol::Dist_Map*    tPDVMap   = tDistributedFactory.create_map( tOwnedPDVIds );
            sol::Dist_Vector* tdIQIdPDV = tDistributedFactory.create_vector( tPDVMap, 2, false, true );

            // Fill values
            if ( par_rank() == 0 )
            {
                tdIQIdPDV->replace_global_values( tFullPDVIds, gdIQIdPDV1, 0 );
                tdIQIdPDV->replace_global_values( tFullPDVIds, gdIQIdPDV2, 1 );
            }
            tdIQIdPDV->vector_global_assembly();

            return tdIQIdPDV;
        }
    };

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Interpolation PDV Creation", "[gen], [pdv], [interpolation pdv], [interpolation pdv serial]" )
    {
        if ( par_size() == 1 )
        {
            // Create PDV_Type host manager
            PDV_Host_Manager_Test tPDVHostManager;

            // PDV type list
            tPDVHostManager.mPDVTypeList = { PDV_Type::DENSITY, PDV_Type::TEMPERATURE, PDV_Type::ELASTIC_MODULUS };
            tPDVHostManager.mPDVTypeMap.set_size( 10, 1, -1 );
            tPDVHostManager.mPDVTypeMap( 3 ) = 0;
            tPDVHostManager.mPDVTypeMap( 4 ) = 1;
            tPDVHostManager.mPDVTypeMap( 5 ) = 2;

            // Node indices per set
            Cell< Cell< uint > > tNodeIndicesPerSet = { { 0, 1, 2, 3 }, { 2, 3, 4, 5 } };

            // Node indices per set to request
            Cell< Matrix< IndexMat > > tRequestNodeIndicesPerSet( 2 );
            tRequestNodeIndicesPerSet( 0 ) = { { 0, 1, 2, 3 } };
            tRequestNodeIndicesPerSet( 1 ) = { { 2, 3, 4, 5 } };

            // Node IDs per set
            Cell< Cell< sint > > tNodeIdsPerSet = { { 0, 1, 2, 3 }, { 2, 3, 4, 5 } };

            // Node ownership per set
            Cell< Cell< uint > > tNodeOwnersPerSet = { { 0, 0, 0, 0 }, { 0, 0, 0, 0 } };

            // Node coordinates per set
            Cell< Matrix< DDRMat > > tNodeCoordinatesPerSet( 2 );
            tNodeCoordinatesPerSet( 0 ).set_size( 4, 3, 0.0 );
            tNodeCoordinatesPerSet( 1 ).set_size( 4, 3, 0.0 );

            // PDV_Type types per set
            Cell< Cell< Cell< PDV_Type > > > tIpPDVTypes( 2 );
            tIpPDVTypes( 0 ) = { { PDV_Type::DENSITY }, { PDV_Type::TEMPERATURE } };
            tIpPDVTypes( 1 ) = { { PDV_Type::TEMPERATURE }, { PDV_Type::ELASTIC_MODULUS } };

            // Create PDV_Type hosts
            tPDVHostManager.set_interpolation_pdv_types( tIpPDVTypes );
            tPDVHostManager.create_interpolation_pdv_hosts(
                    tNodeIndicesPerSet,
                    tNodeIdsPerSet,
                    tNodeOwnersPerSet,
                    tNodeCoordinatesPerSet );

            // Set PDVs
            for ( uint tMeshSetIndex = 0; tMeshSetIndex < 2; tMeshSetIndex++ )
            {
                for ( uint tNodeIndex = 0; tNodeIndex < 4; tNodeIndex++ )
                {
                    for ( uint tPDVIndex = 0; tPDVIndex < 2; tPDVIndex++ )
                    {
                        tPDVHostManager.create_interpolation_pdv(
                                (uint)tNodeIndicesPerSet( tMeshSetIndex )( tNodeIndex ),
                                tIpPDVTypes( tMeshSetIndex )( tPDVIndex )( 0 ),
                                (real)tNodeIndicesPerSet( tMeshSetIndex )( tNodeIndex ) + ( tIpPDVTypes( tMeshSetIndex )( tPDVIndex )( 0 ) == PDV_Type::TEMPERATURE ) );
                    }
                }
            }

            // Create PDV IDs
            tPDVHostManager.create_pdv_ids();

            // Check PDVs
            Cell< Matrix< DDRMat > > tPDVValues;
            for ( uint tMeshSetIndex = 0; tMeshSetIndex < 2; tMeshSetIndex++ )
            {
                for ( uint tPDVIndex = 0; tPDVIndex < 2; tPDVIndex++ )
                {
                    tPDVValues.clear();
                    tPDVHostManager.get_ip_pdv_value( tRequestNodeIndicesPerSet( tMeshSetIndex ), tIpPDVTypes( tMeshSetIndex )( tPDVIndex ), tPDVValues );

                    for ( uint tNodeIndex = 0; tNodeIndex < 4; tNodeIndex++ )
                    {
                        CHECK( tPDVValues( 0 )( tNodeIndex ) == Approx( (real)tNodeIndicesPerSet( tMeshSetIndex )( tNodeIndex ) + ( tIpPDVTypes( tMeshSetIndex )( tPDVIndex )( 0 ) == PDV_Type::TEMPERATURE ) ) );
                    }
                }
            }

            // ------------------- Check global map ----------------------- //
            const Matrix< DDSMat >& tLocalGlobalMap = tPDVHostManager.get_my_local_global_map();

            REQUIRE( tLocalGlobalMap.length() == 14 );
            for ( sint tGlobalPDVIndex = 0; tGlobalPDVIndex < 14; tGlobalPDVIndex++ )
            {
                CHECK( tLocalGlobalMap( tGlobalPDVIndex ) == tGlobalPDVIndex );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Parallel Interpolation PDV Creation", "[gen], [pdv], [interpolation pdv], [interpolation pdv parallel]" )
    {
        if ( par_size() == 2 )
        {
            // Create PDV_Type host manager
            PDV_Host_Manager_Test tPDVHostManager;

            tPDVHostManager.mPDVTypeList = { PDV_Type::DENSITY, PDV_Type::TEMPERATURE };
            tPDVHostManager.mPDVTypeMap.set_size( 10, 1, -1 );
            tPDVHostManager.mPDVTypeMap( 3 ) = 0;
            tPDVHostManager.mPDVTypeMap( 4 ) = 1;

            // ----------------- Interpolation PDVs ---------------------- //
            Cell< Cell< uint > > tNodeIndicesPerSet( 2 );
            Cell< Cell< sint > > tNodeIdsPerSet( 2 );
            Cell< Cell< uint > > tNodeOwnersPerSet( 2 );
            Cell< Matrix< DDRMat > > tNodeCoordinatesPerSet( 2 );
            Cell< Cell< Cell< PDV_Type > > > tIpPDVTypes( 2 );

            // Get my rank and other rank
            uint tMyRank = par_rank();
            uint tOtherRank = ( tMyRank + 1 ) % 2;

            // Node indices per set
            tNodeIndicesPerSet( tMyRank ) = { 0, 1, 2, 3 };
            tNodeIndicesPerSet( tOtherRank ) = { tOtherRank * 2, tOtherRank * 2 + 1 };

            // Node coordinates per set
            tNodeCoordinatesPerSet( tMyRank ).set_size( 4, 3, 0.0 );
            tNodeCoordinatesPerSet( tOtherRank ).set_size( 2, 3, 0.0 );

            // Communication table
            tPDVHostManager.mCommTable.set_size( 2, 1, 1 );
            tPDVHostManager.mCommTable( tMyRank ) = 0;

            // PDV types on set 0
            tIpPDVTypes( 0 ).resize( 1 );
            tIpPDVTypes( 0 )( 0 ).resize( 1 );
            tIpPDVTypes( 0 )( 0 )( 0 ) = PDV_Type::DENSITY;

            // PDV types on set 1
            tIpPDVTypes( 1 ).resize( 2 );
            tIpPDVTypes( 1 )( 0 ).resize( 1 );
            tIpPDVTypes( 1 )( 1 ).resize( 1 );
            tIpPDVTypes( 1 )( 0 )( 0 ) = PDV_Type::TEMPERATURE;
            tIpPDVTypes( 1 )( 1 )( 0 ) = PDV_Type::DENSITY;

            if ( par_rank() == 0 )
            {
                // Node indices, IDs and owners
                tNodeIdsPerSet = { { 0, 1, 2, 3 }, { 2, 3 } };
                tNodeOwnersPerSet = { { 0, 0, 0, 1 }, { 0, 1 } };

                // Vertex ID to index map
                tPDVHostManager.mIPVertexIdtoIndMap[ 2 ] = 2;
            }
            else if ( par_rank() == 1 )
            {
                // Node IDs and owners
                tNodeIdsPerSet = { { 2, 3 }, { 2, 3, 4, 5 } };
                tNodeOwnersPerSet = { { 0, 1 }, { 0, 1, 1, 1 } };

                // Vertex ID to index map
                tPDVHostManager.mIPVertexIdtoIndMap[ 3 ] = 1;
            }

            // Create PDV_Type hosts
            tPDVHostManager.set_interpolation_pdv_types( tIpPDVTypes );
            tPDVHostManager.create_interpolation_pdv_hosts(
                    tNodeIndicesPerSet,
                    tNodeIdsPerSet,
                    tNodeOwnersPerSet,
                    tNodeCoordinatesPerSet );

            // Set PDVs
            for ( uint tMeshSetIndex = 0; tMeshSetIndex < 2; tMeshSetIndex++ )
            {
                for ( uint tNodeIndex = 0; tNodeIndex < tNodeIndicesPerSet( tMeshSetIndex ).size(); tNodeIndex++ )
                {
                    for ( uint tPDVGroupIndex = 0; tPDVGroupIndex < tIpPDVTypes( tMeshSetIndex ).size(); tPDVGroupIndex++ )
                    {
                        tPDVHostManager.create_interpolation_pdv(
                                (uint)tNodeIndicesPerSet( tMeshSetIndex )( tNodeIndex ),
                                tIpPDVTypes( tMeshSetIndex )( tPDVGroupIndex )( 0 ),
                                (real)tMeshSetIndex );
                    }
                }
            }

            // Create PDV IDs
            tPDVHostManager.create_pdv_ids();

            // ------------------- Check global map ----------------------- //
            const Matrix< DDSMat >& tLocalGlobalMap   = tPDVHostManager.get_my_local_global_map();
            const Matrix< DDSMat >& tLocalGlobalOverlappingMap = tPDVHostManager.get_my_local_global_overlapping_map();

            if ( par_rank() == 0 )
            {
                REQUIRE( tLocalGlobalMap.length() == 4 );
                CHECK( tLocalGlobalMap( 0 ) == 0 );
                CHECK( tLocalGlobalMap( 1 ) == 1 );
                CHECK( tLocalGlobalMap( 2 ) == 2 );
                CHECK( tLocalGlobalMap( 3 ) == 3 );

                REQUIRE( tLocalGlobalOverlappingMap.length() == 6 );
                CHECK( tLocalGlobalOverlappingMap( 0 ) == 0 );
                CHECK( tLocalGlobalOverlappingMap( 1 ) == 1 );
                CHECK( tLocalGlobalOverlappingMap( 2 ) == 2 );
                CHECK( tLocalGlobalOverlappingMap( 3 ) == 4 );
                CHECK( tLocalGlobalOverlappingMap( 4 ) == 3 );
                CHECK( tLocalGlobalOverlappingMap( 5 ) == 7 );
            }
            if ( par_rank() == 1 )
            {
                REQUIRE( tLocalGlobalMap.length() == 6 );
                CHECK( tLocalGlobalMap( 0 ) == 4 );
                CHECK( tLocalGlobalMap( 1 ) == 5 );
                CHECK( tLocalGlobalMap( 2 ) == 6 );
                CHECK( tLocalGlobalMap( 3 ) == 7 );
                CHECK( tLocalGlobalMap( 4 ) == 8 );
                CHECK( tLocalGlobalMap( 5 ) == 9 );

                REQUIRE( tLocalGlobalOverlappingMap.length() == 8 );
                CHECK( tLocalGlobalOverlappingMap( 0 ) == 2 );
                CHECK( tLocalGlobalOverlappingMap( 1 ) == 4 );
                CHECK( tLocalGlobalOverlappingMap( 2 ) == 5 );
                CHECK( tLocalGlobalOverlappingMap( 3 ) == 6 );
                CHECK( tLocalGlobalOverlappingMap( 4 ) == 3 );
                CHECK( tLocalGlobalOverlappingMap( 5 ) == 7 );
                CHECK( tLocalGlobalOverlappingMap( 6 ) == 8 );
                CHECK( tLocalGlobalOverlappingMap( 7 ) == 9 );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Intersection PDV Creation", "[gen], [pdv], [intersection pdv]" )
    {
        if ( par_size() == 2 )
        {
            // Create circle
            Matrix< DDRMat >            tADVs   = { { 0.0, 0.0, 1.0 } };
            auto tCircleField = std::make_shared< Circle >(
                    tADVs,
                    Matrix< DDUMat >( { { 0, 1, 2 } } ),
                    Matrix< DDUMat >( { { 0, 1, 2 } } ),
                    Matrix< DDRMat >( { {} } ) );

            auto tCircleGeometry = std::make_shared< Level_Set_Geometry >( tCircleField );

            // Create node manager with no background nodes
            Node_Manager_Test tNodeManager;

            // Create PDV_Type host manager (here: test version defined above)
            PDV_Host_Manager_Test tPDVHostManager( tNodeManager );

            // Node IDs/owners per set
            uint             tNumOwnedNodes = 0;
            Matrix< IdMat >  tIpNodeIdsPerSet( 4, 1 );
            Matrix< DDSMat > tIpNodeOwnersPerSet( 4, 1 );

            if ( par_rank() == 0 )
            {
                tIpNodeIdsPerSet = { { 0 }, { 1 }, { 2 }, { 3 } };

                tIpNodeOwnersPerSet = { { 0 }, { 0 }, { 0 }, { 1 } };
                tNumOwnedNodes      = 3;

                tPDVHostManager.mCommTable.set_size( 2, 1, 0 );
                tPDVHostManager.mCommTable( 1, 0 ) = 1;

                tPDVHostManager.mIGVertexIdtoIndMap[ 2 ] = 2;
            }
            else if ( par_rank() == 1 )
            {
                tIpNodeIdsPerSet = { { 2 }, { 3 }, { 4 }, { 5 }, { 6 }, { 7 } };
                tIpNodeOwnersPerSet = { { 0 }, { 1 }, { 1 }, { 1 }, { 1 }, { 1 } };
                tNumOwnedNodes      = 5;

                tPDVHostManager.mCommTable.set_size( 2, 1, 1 );
                tPDVHostManager.mCommTable( 1, 0 ) = 0;

                tPDVHostManager.mIGVertexIdtoIndMap[ 3 ] = 1;
            }

            // Loop over all node indices
            Cell< Background_Node* > tTemporaryBackgroundNodes;
            for ( uint tNodeIndex = 0; tNodeIndex < tIpNodeIdsPerSet.length(); tNodeIndex++ )
            {
                // Go around a circle to create parent coordinates
                real tRadians = tIpNodeIdsPerSet( tNodeIndex ) * M_PI / 2.0;
                Matrix< DDRMat > tFirstParentCoordinates  = { { 0.5 * cos( tRadians ), 0.5 * sin( tRadians ) } };
                Matrix< DDRMat > tSecondParentCoordinates = { { 1.5 * cos( tRadians ), 1.5 * sin( tRadians ) } };

                // Create parent nodes
                auto tFirstNode = new Background_Node( 0, tFirstParentCoordinates );
                Parent_Node tFirstParentNode( *tFirstNode, {{ -1.0 }} );
                auto tSecondNode = new Background_Node( 0, tSecondParentCoordinates );
                Parent_Node tSecondParentNode( *tSecondNode, {{ 1.0 }} );

                // Background nodes need to be deleted later
                tTemporaryBackgroundNodes.push_back( tFirstNode );
                tTemporaryBackgroundNodes.push_back( tSecondNode );

                // Assign as base nodes
                Cell< Background_Node* > tBackgroundNodes( { tFirstNode, tSecondNode } );

                // Create intersection node
                auto tIntersectionNode = new Intersection_Node_Linear(
                        tNodeIndex,
                        tBackgroundNodes,
                        tFirstParentNode,
                        tSecondParentNode,
                        mtk::Geometry_Type::LINE,
                        mtk::Interpolation_Order::LINEAR,
                        *tCircleGeometry );

                // Add to node manager
                tNodeManager.add_derived_node( tIntersectionNode );
                tNodeManager.update_derived_node(
                        tNodeIndex,
                        tIpNodeIdsPerSet( tNodeIndex ),
                        tIpNodeOwnersPerSet( tNodeIndex ) );
            }

            // Create PDV IDs
            tPDVHostManager.create_pdv_ids();

            // ------------------- Check global map ----------------------- //
            const Matrix< DDSMat >& tLocalGlobalMap   = tPDVHostManager.get_my_local_global_map();
            const Matrix< DDSMat >& tLocalGlobalOverlappingMap = tPDVHostManager.get_my_local_global_overlapping_map();

            REQUIRE( tLocalGlobalMap.length() == tNumOwnedNodes * 2 );
            REQUIRE( tLocalGlobalOverlappingMap.length() == ( tNumOwnedNodes + 1 ) * 2 );

            if ( par_rank() == 0 )
            {
                CHECK( tLocalGlobalMap( 0 ) == 0 );
                CHECK( tLocalGlobalMap( 1 ) == 1 );
                CHECK( tLocalGlobalMap( 2 ) == 2 );
                CHECK( tLocalGlobalMap( 3 ) == 3 );
                CHECK( tLocalGlobalMap( 4 ) == 4 );
                CHECK( tLocalGlobalMap( 5 ) == 5 );

                CHECK( tLocalGlobalOverlappingMap( 0 ) == 0 );
                CHECK( tLocalGlobalOverlappingMap( 1 ) == 1 );
                CHECK( tLocalGlobalOverlappingMap( 2 ) == 2 );
                CHECK( tLocalGlobalOverlappingMap( 3 ) == 3 );
                CHECK( tLocalGlobalOverlappingMap( 4 ) == 4 );
                CHECK( tLocalGlobalOverlappingMap( 5 ) == 5 );
                CHECK( tLocalGlobalOverlappingMap( 6 ) == 6 );
                CHECK( tLocalGlobalOverlappingMap( 7 ) == 7 );
            }

            if ( par_rank() == 1 )
            {
                CHECK( tLocalGlobalMap( 0 ) == 6 );
                CHECK( tLocalGlobalMap( 1 ) == 7 );
                CHECK( tLocalGlobalMap( 2 ) == 8 );
                CHECK( tLocalGlobalMap( 3 ) == 9 );
                CHECK( tLocalGlobalMap( 4 ) == 10 );
                CHECK( tLocalGlobalMap( 5 ) == 11 );

                CHECK( tLocalGlobalOverlappingMap( 0 ) == 4 );
                CHECK( tLocalGlobalOverlappingMap( 1 ) == 5 );
                CHECK( tLocalGlobalOverlappingMap( 2 ) == 6 );
                CHECK( tLocalGlobalOverlappingMap( 3 ) == 7 );
                CHECK( tLocalGlobalOverlappingMap( 4 ) == 8 );
                CHECK( tLocalGlobalOverlappingMap( 5 ) == 9 );
                CHECK( tLocalGlobalOverlappingMap( 6 ) == 10 );
                CHECK( tLocalGlobalOverlappingMap( 7 ) == 11 );
            }

            // Set owned ADV IDs
            Matrix< DDSMat > tOwnedADVIds( 0, 0 );
            if ( par_rank() == 0 )
            {
                tOwnedADVIds = { { 0 }, { 1 }, { 2 } };
            }
            tPDVHostManager.set_owned_adv_ids( tOwnedADVIds );

            // Get sensitivities
            Matrix< DDSMat > tFullADVIds( 0, 0 );
            if ( par_rank() == 0 )
            {
                tFullADVIds = tOwnedADVIds;
            }
            Matrix< DDRMat > tdIQIdADV = tPDVHostManager.compute_diqi_dadv( tFullADVIds );

            // Check sensitivities
            if ( par_rank() == 0 )
            {
                CHECK_EQUAL( tdIQIdADV, Matrix< DDRMat>({ { 4.0, 4.0, 0.0 }, { 24.0, 36.0, -16.0 } }), 1.0E6, );
            }

            // Clean up temporary background nodes
            for ( auto iBackgroundNode : tTemporaryBackgroundNodes )
            {
                delete iBackgroundNode;
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "PDV Sensitivities", "[gen], [pdv], [sensitivity], [pdv sensitivity]" )
    {
        // Create PDV_Type host manager (here: test version defined above)
        PDV_Host_Manager_Test tPDVHostManager;

        tPDVHostManager.mPDVTypeList = { PDV_Type::DENSITY };
        tPDVHostManager.mPDVTypeMap.set_size( 10, 1, -1 );
        tPDVHostManager.mPDVTypeMap( 3 ) = 0;

        // Constant property parameter list
        ParameterList tParameterList = moris::prm::create_gen_property_parameter_list();
        tParameterList.set( "field_type", "constant" );
        tParameterList.set( "field_variable_indices", "0" );
        tParameterList.set( "pdv_type", "DENSITY" );

        // Create ADVs
        uint             tNumADVs = 2 * par_size();
        Matrix< DDRMat > tADVs( tNumADVs, 1 );

        // Create constant properties
        Cell< std::shared_ptr< Property > > tProperties( tNumADVs );
        for ( uint tPropertyIndex = 0; tPropertyIndex < tNumADVs; tPropertyIndex++ )
        {
            tParameterList.set( "adv_indices", std::to_string( tPropertyIndex ) );
            Design_Factory tDesignFactory( { tParameterList }, tADVs );
            tProperties( tPropertyIndex ) = tDesignFactory.get_properties()( 0 );
        }

        // Node indices, IDs, ownership and coordinates per set
        Cell< Cell< uint > > tNodeIndicesPerSet( 1 );
        Cell< Cell< sint > > tNodeIdsPerSet( 1 );
        Cell< Cell< uint > > tNodeOwnersPerSet( 1 );
        Cell< Matrix< DDRMat > > tNodeCoordinatesPerSet( 1 );

        tNodeIndicesPerSet( 0 ).resize( gNumPDVs, 1 );
        tNodeIdsPerSet( 0 ).resize( gNumPDVs, 1 );
        tNodeOwnersPerSet( 0 ).resize( gNumPDVs, par_rank() );
        tNodeCoordinatesPerSet( 0 ).resize( gNumPDVs, 3 );

        for ( uint tNodeIndex = 0; tNodeIndex < gNumPDVs; tNodeIndex++ )
        {
            tNodeIndicesPerSet( 0 )( tNodeIndex ) = tNodeIndex;
            tNodeIdsPerSet( 0 )( tNodeIndex )     = tNodeIndex * ( par_rank() + 1 );
        }

        // fill coordinates for all nodes with zeros
        tNodeCoordinatesPerSet( 0 ).fill( 0.0 );

        // PDV_Type types per set
        Cell< Cell< Cell< PDV_Type > > > tIpPDVTypes( 1 );
        tIpPDVTypes( 0 ).resize( 1 );
        tIpPDVTypes( 0 )( 0 ).resize( 1 );
        tIpPDVTypes( 0 )( 0 )( 0 ) = PDV_Type::DENSITY;

        // Communication table
        Matrix< DDSMat > tCommunicationTable( par_size(), 1, 0 );
        for ( uint tProcessorIndex = 1; tProcessorIndex < (uint)par_size(); tProcessorIndex++ )
        {
            tCommunicationTable( tProcessorIndex ) = tProcessorIndex;
        }
        tPDVHostManager.set_communication_table( tCommunicationTable );

        // Create PDV_Type hosts
        tPDVHostManager.set_interpolation_pdv_types( tIpPDVTypes );
        tPDVHostManager.create_interpolation_pdv_hosts(
                tNodeIndicesPerSet,
                tNodeIdsPerSet,
                tNodeOwnersPerSet,
                tNodeCoordinatesPerSet );

        // Set PDVs
        for ( uint tNodeIndex = 0; tNodeIndex < gNumPDVs; tNodeIndex++ )
        {
            tPDVHostManager.create_interpolation_pdv(
                    (uint)tNodeIndicesPerSet( 0 )( tNodeIndex ),
                    tIpPDVTypes( 0 )( 0 )( 0 ),
                    tProperties( tNodeIndex % tNumADVs ) );
        }
        tPDVHostManager.create_pdv_ids();

        // Owned ADV IDs
        tPDVHostManager.set_owned_adv_ids( { { 2 * par_rank() }, { 2 * par_rank() + 1 } } );

        // Full ADV IDs
        Matrix< DDSMat > tFullADVIds( 0, 0 );
        if ( par_rank() == 0 )
        {
            tFullADVIds.resize( tNumADVs, 1 );
            for ( uint tADVIndex = 0; tADVIndex < tNumADVs; tADVIndex++ )
            {
                tFullADVIds( tADVIndex ) = tADVIndex;
            }
        }

        // Compute sensitivities
        Matrix< DDRMat > tdIQIdADV = tPDVHostManager.compute_diqi_dadv( tFullADVIds );

        // Check sensitivities
        REQUIRE( tdIQIdADV.n_rows() == 2 );
        if ( par_rank() == 0 )
        {
            REQUIRE( tdIQIdADV.n_cols() == tNumADVs );
            for ( uint tADVIndex = 0; tADVIndex < tNumADVs; tADVIndex++ )
            {
                // First IQI
                CHECK( tdIQIdADV( 0, tADVIndex ) == Approx( ( gNumPDVs / tNumADVs ) + ( tADVIndex < gNumPDVs % tNumADVs ) ) );

                // Second IQI
                real tdIQI2dADV = 0.0;
                uint tPDV       = tADVIndex;
                while ( tPDV < gNumPDVs )
                {
                    tdIQI2dADV += tPDV;
                    tPDV += tNumADVs;
                }
                CHECK( tdIQIdADV( 1, tADVIndex ) == Approx( tdIQI2dADV ) );
            }
        }
        else
        {
            REQUIRE( tdIQIdADV.n_cols() == 0 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}
