/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_GEN_Pdv.cpp
 *
 */

#include "catch.hpp"
#include "math.h"
#include "cl_Matrix.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_MSI_Design_Variable_Interface.hpp"

#include "cl_GEN_Circle.hpp"
#include "fn_GEN_create_properties.hpp"
#include "fn_PRM_GEN_Parameters.hpp"

#define protected public
#define private public
#include "cl_GEN_Pdv_Host_Manager.hpp"
#include "cl_GEN_Intersection_Node_Linear.hpp"
#undef protected
#undef private

#include "fn_GEN_check_equal.hpp"

namespace moris
{

    //------------------------------------------------------------------------------------------------------------------

    // Dummy IQI sensitivity values so a FEM model doesn't have to be created
    uint             gNumPDVs = 16;
    Matrix< DDRMat > gdIQIdPDV1( 1, gNumPDVs );
    Matrix< DDRMat > gdIQIdPDV2( 1, gNumPDVs );

    namespace ge
    {
        class Pdv_Host_Manager_Test : public ge::Pdv_Host_Manager
        {
            sol::Dist_Vector*
            get_dQIdp()
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
    }    // namespace ge

    //--------------------------------------------------------------------------------------------------------------

    namespace ge
    {
        TEST_CASE( "Interpolation PDV Creation", "[gen], [pdv], [interpolation pdv], [interpolation pdv serial]" )
        {
            if ( par_size() == 1 )
            {
                // Create PDV_Type host manager
                Pdv_Host_Manager tPDVHostManager;

                // PDV type list
                tPDVHostManager.mPdvTypeList = { PDV_Type::DENSITY, PDV_Type::TEMPERATURE, PDV_Type::ELASTIC_MODULUS };
                tPDVHostManager.mPdvTypeMap.set_size( 10, 1, -1 );
                tPDVHostManager.mPdvTypeMap( 3 ) = 0;
                tPDVHostManager.mPdvTypeMap( 4 ) = 1;
                tPDVHostManager.mPdvTypeMap( 5 ) = 2;

                // Node indices per set
                Cell< Matrix< DDSMat > > tIpNodeIndicesPerSet( 2 );
                tIpNodeIndicesPerSet( 0 ).resize( 4, 1 );
                tIpNodeIndicesPerSet( 1 ).resize( 4, 1 );
                tIpNodeIndicesPerSet( 0 ) = { { 0, 1, 2, 3 } };
                tIpNodeIndicesPerSet( 1 ) = { { 2, 3, 4, 5 } };

                // Node ownership per set
                Cell< Matrix< DDSMat > > tIpNodeOwnersPerSet( 2 );
                tIpNodeOwnersPerSet( 0 ).resize( 4, 1 );
                tIpNodeOwnersPerSet( 1 ).resize( 4, 1 );
                tIpNodeOwnersPerSet( 0 ) = { { 0, 0, 0, 0 } };
                tIpNodeOwnersPerSet( 1 ) = { { 0, 0, 0, 0 } };

                // Node coordinates per set
                Cell< Matrix< DDRMat > > tNodeCoordinatesPerSet( 2 );
                tNodeCoordinatesPerSet( 0 ).resize( 4, 3 );
                tNodeCoordinatesPerSet( 1 ).resize( 4, 3 );
                tNodeCoordinatesPerSet( 0 ).fill( 0.0 );
                tNodeCoordinatesPerSet( 1 ).fill( 0.0 );

                // PDV_Type types per set
                Cell< Cell< Cell< PDV_Type > > > tIpPdvTypes( 2 );
                tIpPdvTypes( 0 ).resize( 2 );
                tIpPdvTypes( 1 ).resize( 2 );
                tIpPdvTypes( 0 )( 0 ).resize( 1 );
                tIpPdvTypes( 0 )( 1 ).resize( 1 );
                tIpPdvTypes( 1 )( 0 ).resize( 1 );
                tIpPdvTypes( 1 )( 1 ).resize( 1 );
                tIpPdvTypes( 0 )( 0 )( 0 ) = PDV_Type::DENSITY;
                tIpPdvTypes( 0 )( 1 )( 0 ) = PDV_Type::TEMPERATURE;
                tIpPdvTypes( 1 )( 0 )( 0 ) = PDV_Type::TEMPERATURE;
                tIpPdvTypes( 1 )( 1 )( 0 ) = PDV_Type::ELASTIC_MODULUS;

                // Create PDV_Type hosts
                tPDVHostManager.create_interpolation_pdv_hosts(
                        tIpNodeIndicesPerSet,
                        tIpNodeIndicesPerSet,
                        tIpNodeOwnersPerSet,
                        tNodeCoordinatesPerSet,
                        tIpPdvTypes );

                // Set PDVs
                for ( uint tMeshSetIndex = 0; tMeshSetIndex < 2; tMeshSetIndex++ )
                {
                    for ( uint tNodeIndex = 0; tNodeIndex < 4; tNodeIndex++ )
                    {
                        for ( uint tPdvIndex = 0; tPdvIndex < 2; tPdvIndex++ )
                        {
                            tPDVHostManager.create_interpolation_pdv(
                                    (uint)tIpNodeIndicesPerSet( tMeshSetIndex )( tNodeIndex ),
                                    tIpPdvTypes( tMeshSetIndex )( tPdvIndex )( 0 ),
                                    (real)tIpNodeIndicesPerSet( tMeshSetIndex )( tNodeIndex ) + ( tIpPdvTypes( tMeshSetIndex )( tPdvIndex )( 0 ) == PDV_Type::TEMPERATURE ) );
                        }
                    }
                }

                // Create PDV IDs
                tPDVHostManager.create_pdv_ids();

                // Check PDVs
                Cell< Matrix< DDRMat > > tPdvValues;
                for ( uint tMeshSetIndex = 0; tMeshSetIndex < 2; tMeshSetIndex++ )
                {
                    for ( uint tPdvIndex = 0; tPdvIndex < 2; tPdvIndex++ )
                    {
                        tPdvValues.clear();
                        tPDVHostManager.get_ip_pdv_value( tIpNodeIndicesPerSet( tMeshSetIndex ), tIpPdvTypes( tMeshSetIndex )( tPdvIndex ), tPdvValues );

                        for ( uint tNodeIndex = 0; tNodeIndex < 4; tNodeIndex++ )
                        {
                            CHECK( tPdvValues( 0 )( tNodeIndex ) == Approx( (real)tIpNodeIndicesPerSet( tMeshSetIndex )( tNodeIndex ) + ( tIpPdvTypes( tMeshSetIndex )( tPdvIndex )( 0 ) == PDV_Type::TEMPERATURE ) ) );
                        }
                    }
                }

                // ------------------- Check global map ----------------------- //
                const Matrix< DDSMat >& tLocalGlobalMap = tPDVHostManager.get_my_local_global_map();

                REQUIRE( tLocalGlobalMap.length() == 14 );
                for ( sint tGlobalPdvIndex = 0; tGlobalPdvIndex < 14; tGlobalPdvIndex++ )
                {
                    CHECK( tLocalGlobalMap( tGlobalPdvIndex ) == tGlobalPdvIndex );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        TEST_CASE( "Parallel Interpolation PDV Creation", "[gen], [pdv], [interpolation pdv], [interpolation pdv parallel]" )
        {
            if ( par_size() == 2 )
            {
                // Create PDV_Type host manager
                Pdv_Host_Manager tPDVHostManager;

                tPDVHostManager.mPdvTypeList = { PDV_Type::DENSITY, PDV_Type::TEMPERATURE };
                tPDVHostManager.mPdvTypeMap.set_size( 10, 1, -1 );
                tPDVHostManager.mPdvTypeMap( 3 ) = 0;
                tPDVHostManager.mPdvTypeMap( 4 ) = 1;

                // ----------------- Interpolation PDVs ---------------------- //

                Cell< Matrix< DDSMat > > tIpNodeIndicesPerSet( 1 );
                tIpNodeIndicesPerSet( 0 ).resize( 4, 1 );

                Cell< Matrix< DDSMat > > tIpNodeIdsPerSet( 1 );
                tIpNodeIdsPerSet( 0 ).resize( 4, 1 );

                Cell< Matrix< DDSMat > > tIpNodeOwnersPerSet( 1 );
                tIpNodeOwnersPerSet( 0 ).resize( 4, 1 );

                // define and fill coordinates for all nodes with zeros
                Cell< Matrix< DDRMat > > tNodeCoordinatesPerSet( 1 );
                tNodeCoordinatesPerSet( 0 ).resize( 4, 3 );
                tNodeCoordinatesPerSet( 0 ).fill( 0.0 );

                // PDV_Type types per set
                Cell< Cell< Cell< PDV_Type > > > tIpPdvTypes( 1 );
                tIpPdvTypes( 0 ).resize( 2 );
                tIpPdvTypes( 0 )( 0 ).resize( 1 );
                tIpPdvTypes( 0 )( 1 ).resize( 1 );
                tIpPdvTypes( 0 )( 0 )( 0 ) = PDV_Type::DENSITY;
                tIpPdvTypes( 0 )( 1 )( 0 ) = PDV_Type::TEMPERATURE;

                if ( par_rank() == 0 )
                {
                    // Node indices per set
                    tIpNodeIndicesPerSet( 0 ) = { { 0, 1, 2, 3 } };

                    tIpNodeIdsPerSet( 0 ) = { { 0, 1, 2, 3 } };

                    tIpNodeOwnersPerSet( 0 ) = { { 0, 0, 0, 1 } };

                    tPDVHostManager.mCommTable.set_size( 2, 1, 0 );
                    tPDVHostManager.mCommTable( 1, 0 ) = 1;

                    tPDVHostManager.mIPVertexIdtoIndMap[ 2 ] = 2;
                }
                else if ( par_rank() == 1 )
                {
                    // Node indices per set
                    tIpNodeIndicesPerSet( 0 ) = { { 0, 1, 2, 3 } };

                    tIpNodeIdsPerSet( 0 ) = { { 2, 3, 4, 5 } };

                    tIpNodeOwnersPerSet( 0 ) = { { 0, 1, 1, 1 } };

                    tPDVHostManager.mCommTable.set_size( 2, 1, 1 );
                    tPDVHostManager.mCommTable( 1, 0 ) = 0;

                    tPDVHostManager.mIPVertexIdtoIndMap[ 3 ] = 1;
                }

                // Create PDV_Type hosts
                tPDVHostManager.create_interpolation_pdv_hosts(
                        tIpNodeIndicesPerSet,
                        tIpNodeIdsPerSet,
                        tIpNodeOwnersPerSet,
                        tNodeCoordinatesPerSet,
                        tIpPdvTypes );

                // Set PDVs
                for ( uint tMeshSetIndex = 0; tMeshSetIndex < 1; tMeshSetIndex++ )
                {
                    for ( uint tNodeIndex = 0; tNodeIndex < 4; tNodeIndex++ )
                    {
                        for ( uint tPdvIndex = 0; tPdvIndex < 2; tPdvIndex++ )
                        {
                            tPDVHostManager.create_interpolation_pdv(
                                    (uint)tIpNodeIndicesPerSet( tMeshSetIndex )( tNodeIndex ),
                                    tIpPdvTypes( tMeshSetIndex )( tPdvIndex )( 0 ),
                                    (real)tMeshSetIndex );
                        }
                    }
                }

                // Create PDV IDs
                tPDVHostManager.create_pdv_ids();

                // ------------------- Check global map ----------------------- //
                const Matrix< DDSMat >& tLocalGlobalMap   = tPDVHostManager.get_my_local_global_map();
                const Matrix< DDSMat >& tLocalGlobalOSMap = tPDVHostManager.get_my_local_global_overlapping_map();

                REQUIRE( tLocalGlobalMap.length() == 6 );
                REQUIRE( tLocalGlobalOSMap.length() == 8 );

                if ( par_rank() == 0 )
                {
                    CHECK( tLocalGlobalMap( 0 ) == 0 );
                    CHECK( tLocalGlobalMap( 1 ) == 1 );
                    CHECK( tLocalGlobalMap( 2 ) == 2 );
                    CHECK( tLocalGlobalMap( 3 ) == 3 );
                    CHECK( tLocalGlobalMap( 4 ) == 4 );
                    CHECK( tLocalGlobalMap( 5 ) == 5 );

                    CHECK( tLocalGlobalOSMap( 0 ) == 0 );
                    CHECK( tLocalGlobalOSMap( 1 ) == 1 );
                    CHECK( tLocalGlobalOSMap( 2 ) == 2 );
                    CHECK( tLocalGlobalOSMap( 3 ) == 6 );
                    CHECK( tLocalGlobalOSMap( 4 ) == 3 );
                    CHECK( tLocalGlobalOSMap( 5 ) == 4 );
                    CHECK( tLocalGlobalOSMap( 6 ) == 5 );
                    CHECK( tLocalGlobalOSMap( 7 ) == 9 );
                }
                if ( par_rank() == 1 )
                {
                    CHECK( tLocalGlobalMap( 0 ) == 6 );
                    CHECK( tLocalGlobalMap( 1 ) == 7 );
                    CHECK( tLocalGlobalMap( 2 ) == 8 );
                    CHECK( tLocalGlobalMap( 3 ) == 9 );
                    CHECK( tLocalGlobalMap( 4 ) == 10 );
                    CHECK( tLocalGlobalMap( 5 ) == 11 );

                    CHECK( tLocalGlobalOSMap( 0 ) == 2 );
                    CHECK( tLocalGlobalOSMap( 1 ) == 6 );
                    CHECK( tLocalGlobalOSMap( 2 ) == 7 );
                    CHECK( tLocalGlobalOSMap( 3 ) == 8 );
                    CHECK( tLocalGlobalOSMap( 4 ) == 5 );
                    CHECK( tLocalGlobalOSMap( 5 ) == 9 );
                    CHECK( tLocalGlobalOSMap( 6 ) == 10 );
                    CHECK( tLocalGlobalOSMap( 7 ) == 11 );
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
                std::shared_ptr< Geometry > tCircle = std::make_shared< Circle >(
                        tADVs,
                        Matrix< DDUMat >( { { 0, 1, 2 } } ),
                        Matrix< DDUMat >( { { 0, 1, 2 } } ),
                        Matrix< DDRMat >( { {} } ) );

                // Create PDV_Type host manager (here: test version defined above)
                Pdv_Host_Manager_Test tPDVHostManager;

                tPDVHostManager.set_num_background_nodes( 0 );

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
                for ( uint tNodeIndex = 0; tNodeIndex < tIpNodeIdsPerSet.length(); tNodeIndex++ )
                {
                    // Go around a circle to create parent coordinates
                    real tRadians = tIpNodeIdsPerSet( tNodeIndex ) * M_PI / 2.0;

                    Matrix< DDRMat > tFirstParentCoordinates  = { { 0.5 * cos( tRadians ), 0.5 * sin( tRadians ) } };
                    Matrix< DDRMat > tSecondParentCoordinates = { { 1.5 * cos( tRadians ), 1.5 * sin( tRadians ) } };

                    // Create intersection node
                    std::shared_ptr< Intersection_Node > tIntersectionNode =    //
                            std::make_shared< Intersection_Node_Linear >(
                                    nullptr,
                                    nullptr,
                                    0,
                                    0,
                                    tFirstParentCoordinates,
                                    tSecondParentCoordinates,
                                    tCircle,
                                    0.0,
                                    0.0 );

                    // Add intersection node to PDV host manager
                    tPDVHostManager.set_intersection_node( tNodeIndex, tIntersectionNode );
                    tPDVHostManager.update_intersection_node(
                            tNodeIndex,
                            tIpNodeIdsPerSet( tNodeIndex ),
                            tIpNodeOwnersPerSet( tNodeIndex ) );
                }

                // Create PDV IDs
                tPDVHostManager.create_pdv_ids();

                // ------------------- Check global map ----------------------- //
                const Matrix< DDSMat >& tLocalGlobalMap   = tPDVHostManager.get_my_local_global_map();
                const Matrix< DDSMat >& tLocalGlobalOSMap = tPDVHostManager.get_my_local_global_overlapping_map();

                REQUIRE( tLocalGlobalMap.length() == tNumOwnedNodes * 2 );
                REQUIRE( tLocalGlobalOSMap.length() == ( tNumOwnedNodes + 1 ) * 2 );

                if ( par_rank() == 0 )
                {
                    CHECK( tLocalGlobalMap( 0 ) == 0 );
                    CHECK( tLocalGlobalMap( 1 ) == 1 );
                    CHECK( tLocalGlobalMap( 2 ) == 2 );
                    CHECK( tLocalGlobalMap( 3 ) == 3 );
                    CHECK( tLocalGlobalMap( 4 ) == 4 );
                    CHECK( tLocalGlobalMap( 5 ) == 5 );

                    CHECK( tLocalGlobalOSMap( 0 ) == 0 );
                    CHECK( tLocalGlobalOSMap( 1 ) == 1 );
                    CHECK( tLocalGlobalOSMap( 2 ) == 2 );
                    CHECK( tLocalGlobalOSMap( 3 ) == 3 );
                    CHECK( tLocalGlobalOSMap( 4 ) == 4 );
                    CHECK( tLocalGlobalOSMap( 5 ) == 5 );
                    CHECK( tLocalGlobalOSMap( 6 ) == 6 );
                    CHECK( tLocalGlobalOSMap( 7 ) == 7 );
                }

                if ( par_rank() == 1 )
                {
                    CHECK( tLocalGlobalMap( 0 ) == 6 );
                    CHECK( tLocalGlobalMap( 1 ) == 7 );
                    CHECK( tLocalGlobalMap( 2 ) == 8 );
                    CHECK( tLocalGlobalMap( 3 ) == 9 );
                    CHECK( tLocalGlobalMap( 4 ) == 10 );
                    CHECK( tLocalGlobalMap( 5 ) == 11 );

                    CHECK( tLocalGlobalOSMap( 0 ) == 4 );
                    CHECK( tLocalGlobalOSMap( 1 ) == 5 );
                    CHECK( tLocalGlobalOSMap( 2 ) == 6 );
                    CHECK( tLocalGlobalOSMap( 3 ) == 7 );
                    CHECK( tLocalGlobalOSMap( 4 ) == 8 );
                    CHECK( tLocalGlobalOSMap( 5 ) == 9 );
                    CHECK( tLocalGlobalOSMap( 6 ) == 10 );
                    CHECK( tLocalGlobalOSMap( 7 ) == 11 );
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
                    check_equal( tdIQIdADV, { { 4.0, 4.0, 0.0 }, { 24.0, 36.0, -16.0 } }, 1e-12 );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        TEST_CASE( "PDV Sensitivities", "[gen], [pdv], [sensitivity], [pdv sensitivity]" )
        {
            // Create PDV_Type host manager (here: test version defined above)
            Pdv_Host_Manager_Test tPDVHostManager;

            tPDVHostManager.mPdvTypeList = { PDV_Type::DENSITY };
            tPDVHostManager.mPdvTypeMap.set_size( 10, 1, -1 );
            tPDVHostManager.mPdvTypeMap( 3 ) = 0;

            // Constant property parameter list
            ParameterList tParameterList = moris::prm::create_gen_property_parameter_list();
            ;
            tParameterList.set( "type", "constant" );
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
                tProperties( tPropertyIndex ) = create_property( tParameterList, tADVs );
            }

            // Node indices, IDs, ownership and coordinates per set
            Cell< Matrix< DDSMat > > tIpNodeIndicesPerSet( 1 );
            Cell< Matrix< DDSMat > > tIpNodeIdsPerSet( 1 );
            Cell< Matrix< DDSMat > > tIpNodeOWnersPerSet( 1 );
            Cell< Matrix< DDRMat > > tNodeCoordinatesPerSet( 1 );

            tIpNodeIndicesPerSet( 0 ).set_size( gNumPDVs, 1 );
            tIpNodeIdsPerSet( 0 ).set_size( gNumPDVs, 1 );
            tIpNodeOWnersPerSet( 0 ).set_size( gNumPDVs, 1, par_rank() );
            tNodeCoordinatesPerSet( 0 ).resize( gNumPDVs, 3 );

            for ( uint tNodeIndex = 0; tNodeIndex < gNumPDVs; tNodeIndex++ )
            {
                tIpNodeIndicesPerSet( 0 )( tNodeIndex ) = tNodeIndex;
                tIpNodeIdsPerSet( 0 )( tNodeIndex )     = tNodeIndex * ( par_rank() + 1 );
            }

            // fill coordinates for all nodes with zeros
            tNodeCoordinatesPerSet( 0 ).fill( 0.0 );

            // PDV_Type types per set
            Cell< Cell< Cell< PDV_Type > > > tIpPdvTypes( 1 );
            tIpPdvTypes( 0 ).resize( 1 );
            tIpPdvTypes( 0 )( 0 ).resize( 1 );
            tIpPdvTypes( 0 )( 0 )( 0 ) = PDV_Type::DENSITY;

            // Communication table
            Matrix< DDSMat > tCommunicationTable( par_size(), 1, 0 );
            for ( uint tProcessorIndex = 1; tProcessorIndex < (uint)par_size(); tProcessorIndex++ )
            {
                tCommunicationTable( tProcessorIndex ) = tProcessorIndex;
            }
            tPDVHostManager.set_communication_table( tCommunicationTable );

            // Create PDV_Type hosts
            tPDVHostManager.create_interpolation_pdv_hosts(
                    tIpNodeIndicesPerSet,
                    tIpNodeIdsPerSet,
                    tIpNodeOWnersPerSet,
                    tNodeCoordinatesPerSet,
                    tIpPdvTypes );

            // Set PDVs
            for ( uint tNodeIndex = 0; tNodeIndex < gNumPDVs; tNodeIndex++ )
            {
                tPDVHostManager.create_interpolation_pdv(
                        (uint)tIpNodeIndicesPerSet( 0 )( tNodeIndex ),
                        tIpPdvTypes( 0 )( 0 )( 0 ),
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

    }    // namespace ge
}    // namespace moris
