/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Mesh_Data.cpp
 *
 */

#include <memory>
#include <mpi.h>
#include "catch.hpp"

// Logger include
#include "cl_Logger.hpp"
#include "paths.hpp"

// XTKL: General Includes
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "mesh/cl_Mesh_Data.hpp"

#include "xtk_typedefs.hpp"
#include "cl_XTK_Matrix_Base_Utilities.hpp"    // For print and equal to for matrices
#include "cl_Mesh_Enums.hpp"
#include "mesh/cl_Mesh_Builder.hpp"
#include "mesh/cl_Mesh_Data_Stk.hpp"
#include "mesh/cl_Mesh_Builder_Stk.hpp"

// Topology includes
#include "topology/cl_XTK_Hexahedron_8_Topology.hpp"

using namespace xtk;

/*
 * Test Cases
 *  - Mesh from file
 *  - Mesh from data
 */
TEST_CASE( "STK Mesh Test Serial", "[MESH][STK]" )
{
    int tProcRank = 0;
    int tProcSize = 0;
    MPI_Comm_rank( MPI_COMM_WORLD, &tProcRank );
    MPI_Comm_size( MPI_COMM_WORLD, &tProcSize );

    if ( tProcSize == 1 )
    {

        // Intialize STK Mesh Builder
        mesh::Mesh_Builder_Stk< xtk::real, xtk::size_t, xtk::moris::DDRMat, xtk::moris::DDSTMat > tMeshBuilder;

        SECTION( "Mesh from file" )
        {
            // Specify Mesh Inputs
            std::string                tMeshFileName = "generated:1x1x2";
            moris::Cell< std::string > tScalarFields( 0 );

            // Generate mesh from file
            std::shared_ptr< mesh::Mesh_Data< xtk::real, xtk::size_t, xtk::moris::DDRMat, xtk::moris::DDSTMat > > tMeshData = tMeshBuilder.build_mesh_from_string( tMeshFileName, tScalarFields, true );

            // Get information from Mesh
            xtk::size_t tNumNodes    = tMeshData->get_num_entities( EntityRank::NODE );
            xtk::size_t tNumEdges    = tMeshData->get_num_entities( EntityRank::EDGE );
            xtk::size_t tNumFaces    = tMeshData->get_num_entities( EntityRank::FACE );
            xtk::size_t tNumElements = tMeshData->get_num_entities( EntityRank::ELEMENT );

            REQUIRE( tNumNodes == 12 );
            REQUIRE( tNumEdges == 20 );
            REQUIRE( tNumFaces == 11 );
            REQUIRE( tNumElements == 2 );

            // Get information about element 1
            moris::Matrix< moris::DDSTMat > tElement1Nodes =
                    tMeshData->get_entity_connected_to_entity_loc_inds( (xtk::size_t)0, EntityRank::ELEMENT, EntityRank::NODE );

            moris::Matrix< moris::DDSTMat > tElement1Faces =
                    tMeshData->get_entity_connected_to_entity_loc_inds( (xtk::size_t)0, EntityRank::ELEMENT, EntityRank::FACE );

            moris::Matrix< moris::DDSTMat > tElement1Edges =
                    tMeshData->get_entity_connected_to_entity_loc_inds( (xtk::size_t)0, EntityRank::ELEMENT, EntityRank::EDGE );

            // Define expected values
            moris::Matrix< moris::DDSTMat > tExpectedElement1Nodes( { { 0, 1, 3, 2, 4, 5, 7, 6 } } );

            moris::Matrix< moris::DDSTMat > tExpectedElement1Edges( { { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 } } );

            moris::Matrix< moris::DDSTMat > tExpectedElement1Faces( { { 0, 1, 2, 3, 4, 5 } } );

            CHECK( xtk::equal_to( tElement1Nodes, tExpectedElement1Nodes ) );
            CHECK( xtk::equal_to( tElement1Edges, tExpectedElement1Edges ) );
            CHECK( xtk::equal_to( tElement1Faces, tExpectedElement1Faces ) );

            moris::Matrix< moris::DDSTMat > tElement2Nodes = tMeshData->get_entity_connected_to_entity_loc_inds( (xtk::size_t)1, EntityRank::ELEMENT, EntityRank::NODE );
            moris::Matrix< moris::DDSTMat > tElement2Edges = tMeshData->get_entity_connected_to_entity_loc_inds( (xtk::size_t)1, EntityRank::ELEMENT, EntityRank::EDGE );
            moris::Matrix< moris::DDSTMat > tElement2Faces = tMeshData->get_entity_connected_to_entity_loc_inds( (xtk::size_t)1, EntityRank::ELEMENT, EntityRank::FACE );
            moris::Matrix< moris::DDSTMat > tExpectedElement2Nodes( { { 4, 5, 7, 6, 8, 9, 11, 10 } } );
            moris::Matrix< moris::DDSTMat > tExpectedElement2Edges( { { 4, 5, 6, 7, 12, 13, 14, 15, 16, 17, 18, 19 } } );
            moris::Matrix< moris::DDSTMat > tExpectedElement2Faces( { { 6, 7, 8, 9, 5, 10 } } );

            CHECK( xtk::equal_to( tElement2Nodes, tExpectedElement2Nodes ) );
            CHECK( xtk::equal_to( tElement2Edges, tExpectedElement2Edges ) );
            CHECK( xtk::equal_to( tElement2Faces, tExpectedElement2Faces ) );

            // Check node coordinates
            moris::Matrix< moris::DDSTMat > tNodeIndex( { { 0, 1, 2, 4 } } );
            moris::Matrix< moris::DDRMat >  tNodeCoordinates = tMeshData->get_selected_node_coordinates_loc_inds( tNodeIndex );
            moris::Matrix< moris::DDRMat >  tExpectedNodeCoordinates( { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } } );

            CHECK( xtk::equal_to( tNodeCoordinates, tExpectedNodeCoordinates ) );
        }
    }
}

TEST_CASE( "Batch Create New Nodes Functions", "[MESH][BATCH_CREATE]" )
{

    /*
     * Tests the following functions:
     * - allocate entity ids
     * - get local to global map
     * - batch create new nodes
     *
     * TODO: The following tests
     * - batch create new nodes with fields (NOT DONE)
     * - parallel versions of the above
     */

    /*
     * Get information about MPI
     */
    int tProcRank = 0;
    int tProcSize = 0;
    MPI_Comm_rank( MPI_COMM_WORLD, &tProcRank );
    MPI_Comm_size( MPI_COMM_WORLD, &tProcSize );

    /*
     * Construct the Mesh
     */
    moris::Cell< std::string >                                                                            tScalarFields( 0 );
    std::string                                                                                           tMeshFileName = "generated:1x1x2";
    mesh::Mesh_Builder_Stk< xtk::real, xtk::size_t, xtk::moris::DDRMat, xtk::moris::DDSTMat >             tMeshBuilder;
    std::shared_ptr< mesh::Mesh_Data< xtk::real, xtk::size_t, xtk::moris::DDRMat, xtk::moris::DDSTMat > > tMeshData = tMeshBuilder.build_mesh_from_string( tMeshFileName, tScalarFields, true );

    /*
     * Test Allocation of Entity Ids
     */

    xtk::size_t tFirstId;
    if ( tProcRank == 0 )
    {
        tFirstId = tMeshData->allocate_entity_ids( 2, EntityRank::NODE );
        CHECK( tFirstId == 13 );
    }
    else
    {
        tFirstId = tMeshData->allocate_entity_ids( 2, EntityRank::NODE );
        CHECK( tFirstId == 15 );
    }

    //    std::cout<<"Proc Rank: "<< tProcRank<< " First Node Id:"<< tFirstId<<std::endl;

    /*
     * Test Batch Creation Of Nodes
     */

    /*
     * Setup pending node data structure
     */
    moris::Cell< xtk::Pending_Node< xtk::real, xtk::size_t, xtk::moris::DDRMat, xtk::moris::DDSTMat > > tPendingNodes( 3 );

    /*
     * Setup Parent Element Topology with Nodes Corresponding to Element 1
     */
    moris::Matrix< moris::DDSTMat >                                                               tElementNodesForTopology = tMeshData->get_entity_connected_to_entity_loc_inds( 0, EntityRank::ELEMENT, EntityRank::NODE );
    xtk::Hexahedron_8_Topology< xtk::real, xtk::size_t, xtk::moris::DDRMat, xtk::moris::DDSTMat > tDummyTopology( tElementNodesForTopology );

    /*
     * Pending Node 0
     */
    xtk::size_t                    tNodeIndex1 = 12;
    xtk::size_t                    tNodeId1    = 51;
    moris::Matrix< moris::DDRMat > tNodeCoords1( { { 1.24, 1.3, 1.5 } } );
    moris::Matrix< moris::DDRMat > tLocalCoords1( { { 0.0, 0.0, 0.0 } } );

    tPendingNodes( 0 ).set_pending_node_info( &tNodeIndex1, &tNodeId1, tNodeCoords1, tDummyTopology, tLocalCoords1 );

    /*
     * Pending Node 1
     */
    xtk::size_t                    tNodeIndex2 = 14;
    xtk::size_t                    tNodeId2    = 94;
    moris::Matrix< moris::DDRMat > tNodeCoords2( { { -3.24, -0.3, -2.5 } } );
    moris::Matrix< moris::DDRMat > tLocalCoords2( { { 0.0, 0.0, 0.0 } } );

    tPendingNodes( 1 ).set_pending_node_info( &tNodeIndex2, &tNodeId2, tNodeCoords2, tDummyTopology, tLocalCoords2 );

    /*
     * Pending Node 2
     */
    xtk::size_t                    tNodeIndex3 = 13;
    xtk::size_t                    tNodeId3    = 200;
    moris::Matrix< moris::DDRMat > tNodeCoords3( { { 1.9, -2.3, 5.5 } } );
    moris::Matrix< moris::DDRMat > tLocalCoords3( { { 0.0, 0.0, 0.0 } } );

    tPendingNodes( 2 ).set_pending_node_info( &tNodeIndex3, &tNodeId3, tNodeCoords3, tDummyTopology, tLocalCoords3 );

    /*
     * Check the map prior to modifying the mesh
     */

    moris::Matrix< moris::DDSTMat > tExpectedMap( { { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 } } );

    if ( tProcSize == 1 )
    {
        CHECK( xtk::equal_to( tExpectedMap, tMeshData->get_local_to_global_map( EntityRank::NODE ) ) );
    }
    //       tMeshData->batch_create_new_nodes(tPendingNodes);

    tExpectedMap = moris::Matrix< moris::DDSTMat >( { { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 51, 200, 94 } } );

    if ( tProcSize == 1 )
    {
        CHECK( xtk::equal_to( tExpectedMap, tMeshData->get_local_to_global_map( EntityRank::NODE ) ) );
    }
    /*
     * Check to see that the coordinates are correct
     *
     */
    moris::Matrix< moris::DDSTMat > tNodeIndices( { { tNodeIndex1, tNodeIndex2, tNodeIndex3 } } );

    moris::Matrix< moris::DDRMat > tExpectedNodeCoordinates( { { 1.24, 1.3, 1.5 }, { 1.9, -2.3, 5.5 }, { -3.24, -0.3, -2.5 } } );

    if ( tProcSize == 1 )
    {
        CHECK( xtk::equal_to( tExpectedNodeCoordinates, ( tMeshData->get_selected_node_coordinates_loc_inds( tNodeIndices ) ) ) );
    }
    /*
     * Do a second round of batch creation with a different number of nodes
     */

    tPendingNodes = moris::Cell< xtk::Pending_Node< xtk::real, xtk::size_t, xtk::moris::DDRMat, xtk::moris::DDSTMat > >( 2 );

    /*
     * Pending Node 0
     */
    tNodeIndex1   = 15;
    tNodeId1      = 100;
    tNodeCoords1  = moris::Matrix< moris::DDRMat >( { { 15, 15, 15 } } );
    tLocalCoords1 = moris::Matrix< moris::DDRMat >( { { 0.0, 0.0, 0.0 } } );

    tPendingNodes( 0 ).set_pending_node_info( &tNodeIndex1, &tNodeId1, tNodeCoords1, tDummyTopology, tLocalCoords1 );

    /*
     * Pending Node 1
     */
    tNodeIndex2 = 16;
    tNodeId2    = 200;
    moris::Matrix< moris::DDRMat > tNodeCoords( { { 16, 16, 16 } } );
    tLocalCoords2 = moris::Matrix< moris::DDRMat >( { { 0.0, 0.0, 0.0 } } );

    tPendingNodes( 1 ).set_pending_node_info( &tNodeIndex2, &tNodeId2, tNodeCoords2, tDummyTopology, tLocalCoords2 );

    //         tMeshData->batch_create_new_nodes(tPendingNodes);
}

TEST_CASE( "Part Ordinals", "[MESH][PARTS][ORDINALS]" )
{

    /*
     * Tests the following functions:
     *  - get_num_buckets
     *  - get_entities_in_bucket_loc_index
     *  - get_entity_part_membership_ordinals
     *  - get_part_name_from_part_ordinals
     *  - get_all_parts by entity rank
     *
     *  Currently, this is only testing parts for elements but will be extended to faces and sides
     */

    /*
     * Load Mesh which has 3 block sets. These blocks are named:
     *  - top_bread
     *  - bacon_lettuce_tomato
     *  - bottom_bread
     *
     * Side Sets will eventually be named
     *  - top_crust
     *  - bottom_crust
     *
     * Node Sets will be named
     *  - sesame_seeds
     */
    std::string                                                                                           tPrefix       = moris::get_base_moris_dir();
    std::string                                                                                           tMeshFileName = tPrefix + "/TestExoFiles/sandwich.e";
    mesh::Mesh_Builder_Stk< xtk::real, xtk::size_t, xtk::moris::DDRMat, xtk::moris::DDSTMat >             tMeshBuilder;
    std::shared_ptr< mesh::Mesh_Data< xtk::real, xtk::size_t, xtk::moris::DDRMat, xtk::moris::DDSTMat > > tMeshData = tMeshBuilder.build_mesh_from_string( tMeshFileName, {}, true );

    /*
     * Get number of buckets in the mesh
     */

    xtk::size_t tNumBuckets = tMeshData->get_num_buckets( EntityRank::ELEMENT );

    /*
     * Iterate over buckets
     */
    moris::Cell< std::string > tPartNames;
    moris::Cell< xtk::size_t > tPartOrdinals;

    for ( xtk::size_t i = 0; i < tNumBuckets; i++ )
    {
        moris::Matrix< moris::DDSTMat > tEntitiesInBucket = tMeshData->get_entities_in_bucket_loc_index( i, EntityRank::ELEMENT );
        tMeshData->get_entity_part_membership_ordinals( (tEntitiesInBucket)( 0, 1 ), EntityRank::ELEMENT, tPartOrdinals );
        tMeshData->get_part_name_from_part_ordinals( tPartOrdinals, tPartNames );
    }

    /*
     * Ask for part ordinal of entity
     */
    size_t     tElementIndex = 14;
    EntityRank tEntityRank   = EntityRank::ELEMENT;

    tMeshData->get_entity_part_membership_ordinals( tElementIndex, tEntityRank, tPartOrdinals );
    tMeshData->get_part_name_from_part_ordinals( tPartOrdinals, tPartNames );

    CHECK( tPartOrdinals.size() == 1 );
    CHECK( tPartNames( 0 ).compare( "meat" ) == 0 );

    /*
     * Do it for a different element
     */
    tElementIndex = 1500;
    tEntityRank   = EntityRank::ELEMENT;

    tMeshData->get_entity_part_membership_ordinals( tElementIndex, tEntityRank, tPartOrdinals );
    tMeshData->get_part_name_from_part_ordinals( tPartOrdinals, tPartNames );

    CHECK( tPartOrdinals.size() == 1 );
    CHECK( tPartNames( 0 ).compare( "bottom_bread" ) == 0 );

    /*
     * And a third element
     */

    tElementIndex = 3000;
    tEntityRank   = EntityRank::ELEMENT;

    tMeshData->get_entity_part_membership_ordinals( tElementIndex, tEntityRank, tPartOrdinals );
    tMeshData->get_part_name_from_part_ordinals( tPartOrdinals, tPartNames );

    CHECK( tPartOrdinals.size() == 1 );
    CHECK( tPartNames( 0 ).compare( "top_bread" ) == 0 );

    tMeshData->get_all_part_names( EntityRank::ELEMENT, tPartNames );

    CHECK( tPartNames.size() == 3 );
    CHECK( tPartNames( 0 ).compare( "meat" ) == 0 );
    CHECK( tPartNames( 1 ).compare( "bottom_bread" ) == 0 );
    CHECK( tPartNames( 2 ).compare( "top_bread" ) == 0 );
}

/*
 * This tests the interface with side sets in a mesh
 */
TEST_CASE( "STK Mesh with Side Set", "[STK][SIDE_SET]" )
{
    /*
     * Load Mesh which is a unit cube with 2 faces belonging to a side set
     */
    std::string                                                                                           tPrefix       = moris::get_base_moris_dir();
    std::string                                                                                           tMeshFileName = tPrefix + "/TestExoFiles/cube_1x1x1_with_side_set.e";
    mesh::Mesh_Builder_Stk< xtk::real, xtk::size_t, xtk::moris::DDRMat, xtk::moris::DDSTMat >             tMeshBuilder;
    std::shared_ptr< mesh::Mesh_Data< xtk::real, xtk::size_t, xtk::moris::DDRMat, xtk::moris::DDSTMat > > tMeshData = tMeshBuilder.build_mesh_from_string( tMeshFileName, {}, true );

    /*
     * Iterate over buckets
     */
    xtk::size_t                tNumBuckets = tMeshData->get_num_buckets( EntityRank::ELEMENT );
    moris::Cell< std::string > tPartNames;
    moris::Cell< xtk::size_t > tPartOrdinals;
    //    xtk::size_t tNumParts;

    for ( xtk::size_t i = 0; i < tNumBuckets; i++ )
    {

        moris::Matrix< moris::DDSTMat > tEntitiesInBucket = tMeshData->get_entities_in_bucket_loc_index( i, EntityRank::ELEMENT );

        if ( tEntitiesInBucket.n_cols() != 0 )
        {
            tMeshData->get_entity_part_membership_ordinals( (tEntitiesInBucket)( 0, 0 ), EntityRank::ELEMENT, tPartOrdinals );
            tMeshData->get_part_name_from_part_ordinals( tPartOrdinals, tPartNames );
        }
    }

    tMeshData->get_all_part_names( EntityRank::FACE, tPartNames );

    //    CHECK(tPartNames.size());
    //    CHECK(tPartNames(0).compare("meat") == 0);
    //    CHECK(tPartNames(1).compare("bottom_bread") == 0);
    //    CHECK(tPartNames(2).compare("top_bread") == 0);
    //    tNumParts = tPartNames.size();
    //    for(size_t iPart = 0; iPart<tNumParts; iPart++)
    //    {
    //        std::cout<<iPart<<": "<< tPartNames(iPart)<< " "<<std::endl;
    //    }

    tPrefix = std::getenv( "XTKOUTPUT" );

    MORIS_ERROR( tPrefix.size() > 0,
            "Environment variable XTKOUTPUT not set." );

    std::string tMeshOutputFile = tPrefix + "/sideset_test.e";
    tMeshData->write_output_mesh( tMeshOutputFile );
}

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

//#include <percept/xfer/STKMeshTransferSetup.hpp>

namespace xtk
{

    TEST_CASE( "Discretizing level set field onto a fine mesh and transferring to a coarse mesh via Percept", "[MESH_FIELDS]" )
    {
        //    // Specify Mesh inputs and outputs as well as fields
        //    const std::string tFineMeshInput        = "../TestExoFiles/fine10x10x10cube.e";
        //    const std::string tCoarseMeshInput      = "generated:10x10x10";
        //    const std::string tMeshNameOutputFine   = "../TestExoFiles/Outputs/fine10x10x10cube_with_lsf.e";
        //    const std::string tMeshNameOutputCoarse = "../TestExoFiles/Outputs/coarse10x10x10cube_with_lsf.e";
        //    const std::string tFieldName1           = "LEVELSET_FIELD_01"; // NOTE: Each mesh declares this field but the level set field is only discretized on the fine mesh and then transfered to the coarse mesh
        //
        //    // Initialize Matrix Factory---------------------------------------------
        //    Matrix_Factory<real, size_t> tMatrixFactory;
        //
        //    // Create Mesh (w/o edge and face data)----------------------------------
        //    Cell<std::string> tScalarFieldNames = {tFieldName1};
        //    mesh::Mesh_Builder_Stk<real, size_t, xtk::moris::DDRMat, xtk::moris::DDSTMat> tMeshBuilder;
        //
        //    // Setup geometry and discretize onto a levelset mesh--------------------
        //    real tRadius =  5;
        //    real tXCenter = 10.0;
        //    real tYCenter = 10.0;
        //    real tZCenter = 10.0;
        //    Analytic_Level_Set_Sphere<real, size_t, xtk::moris::DDRMat, xtk::moris::DDSTMat> tLevelsetSphere1(tRadius, tXCenter, tYCenter, tZCenter);
        //    Cell<Geometry<real, size_t, xtk::moris::DDRMat, xtk::moris::DDSTMat>*> tLevelSetFunctions = {&tLevelsetSphere1};
        //    Mesh_Field_Geometry<real,size_t, xtk::moris::DDRMat, xtk::moris::DDSTMat> tLevelSetMeshManager(tMatrixFactory,tLevelSetFunctions,tFineMeshInput,tScalarFieldNames,tMeshBuilder);
        //
        //    // Build a coarse mesh which the levelset mesh will be transfered to-----
        //    std::shared_ptr<mesh::Mesh_Data<real, size_t,xtk::moris::DDRMat, xtk::moris::DDSTMat>> tCoarseMesh = tMeshBuilder.build_mesh_from_string(tMatrixFactory,tCoarseMeshInput,tScalarFieldNames,false);
        //
        //    // Get Reference to the fine mesh
        //    mesh::Mesh_Data<real, size_t,xtk::moris::DDRMat, xtk::moris::DDSTMat> & tFineMesh = tLevelSetMeshManager.get_level_set_mesh();
        //
        //    // Transfer level set field-----------------------------------------------
        //    // Setup STK Transfer Function (boost shared pointer because this it what Percept uses)
        //    boost::shared_ptr<percept::STKMeshTransfer> tTransfer = percept::buildSTKMeshTransfer<percept::STKMeshTransfer>(
        //                                                                                   tFineMesh.mesh_bulk_data(),
        //                                                                                   tFineMesh.get_coordinate_field(),
        //                                                                                   tFineMesh.get_field(EntityRank::NODE,tFieldName1),
        //                                                                                   tCoarseMesh->mesh_bulk_data(),
        //                                                                                   tCoarseMesh->get_coordinate_field(),
        //                                                                                   tCoarseMesh->get_field(EntityRank::NODE,tFieldName1),
        //                                                                                   "transfer_lsf");
        //
        //    initializeSTKMeshTransfer(&(tTransfer));
        //
        //    tTransfer->apply();
        //
        //    // Write mesh to outputs to exodus files
        //    size_t tTime  = 1;
        //    tFineMesh.write_output_mesh(tMeshNameOutputFine,tScalarFieldNames,tTime);
        //    tCoarseMesh->write_output_mesh(tMeshNameOutputCoarse,tLevelSetMeshManager.get_level_set_field_name(),tTime);
    }

#include <stk_mesh/base/MetaData.hpp>     // for MetaData
#include <stk_mesh/base/BulkData.hpp>     // for BulkData
#include <stk_mesh/base/FieldBase.hpp>    // for FieldBase
#include <stk_mesh/base/Types.hpp>        // for EntityRank, PropertyBase, etc
    TEST_CASE( "MESH FIELDS TESTING", "[MESH_FIELDS]" )
    {
        //    Matrix_Factory<real, size_t> tMatrixFactory;
        //    std::string tMeshInputFile = "/TestExoFiles/mesh_test_fields.e";
        //    mesh::Mesh_Builder_Stk<real, size_t, xtk::moris::DDRMat, xtk::moris::DDSTMat> tMeshBuilder;
        //    std::shared_ptr<mesh::Mesh_Data<xtk::real, xtk::size_t,xtk::moris::DDRMat, xtk::moris::DDSTMat>> tMeshData = tMeshBuilder.build_mesh_from_string(tMatrixFactory, tMeshInputFile, {}, true);
        //    stk::mesh::MetaData & tMeta = tMeshData->mesh_meta_data();
        //    stk::mesh::FieldVector const & tFields = tMeta.get_fields();

        //    CHECK(tFields.size() == 2);

        /*
         * Check to see the field from the exodus file is included and the data has been loaded
         * also see if the interface function accesses the value correctly
         */
        //    moris::Cell<std::string> tFieldNames = {"levelset_field_01"};
        //    CHECK(tMeshData->get_entity_field_value(0,tFieldNames(0), EntityRank::NODE)==Approx(275));
    }

#include <stk_mesh/base/GetEntities.hpp>    // for count_entities

    TEST_CASE( "Mesh Side Set and Block Set Functions",
            "[MESH_SETS]" )
    {
        //    std::string tMeshInputFile = "/TestExoFiles/mesh_test_nodeset_sideset_1x1x1.e";
        //    std::string tMeshOutputFile = "/mesh_test_nodeset_sideset_1x1x1_output.e";
        //    Matrix_Factory<> tMatrixFactory;
        //    mesh::Mesh_Builder_Stk<real, size_t, xtk::moris::DDRMat, xtk::moris::DDSTMat> tMeshBuilder;
        //    std::shared_ptr<mesh::Mesh_Data<xtk::real, xtk::size_t,xtk::moris::DDRMat, xtk::moris::DDSTMat>> tMeshData = tMeshBuilder.build_mesh_from_string(tMatrixFactory, tMeshInputFile, {}, true);
        //
        //    stk::mesh::MetaData & tMeta = tMeshData->mesh_meta_data();
        //    stk::mesh::BulkData & tBulk = tMeshData->mesh_bulk_data();
        //
        //    const stk::mesh::PartVector tParts = tMeta.get_mesh_parts();
        //
        ////    std::cout<<"Number of Parts: "<<tParts.size()<<std::endl;
        ////    std::cout<<"Part 0 Name: "<< tParts[0]->name()<<std::endl;
        ////    std::cout<<"Part 1 Name: "<< tParts[1]->name()<<std::endl;
        ////    std::cout<<"Part 2 Name: "<< tParts[2]->name()<<std::endl;
        ////    std::cout<<"Part 3 Name: "<< tParts[3]->name()<<std::endl;
        //
        //    stk::mesh::Part * tPartSideSetName = tMeta.get_part("surface_1");
        //    unsigned tPartOrdinal = tPartSideSetName->mesh_meta_data_ordinal();
        //    stk::mesh::Part & tPartSideSetOrdinal = tMeta.get_part(tPartOrdinal);
        //
        //    stk::mesh::Selector tSideSetSelectorName(tPartSideSetName);
        //    stk::mesh::Selector tSideSetSelectorOrdinal(tPartSideSetOrdinal);
        //
        ////    xtk::size_t tNumNodes = stk::mesh::count_selected_entities(tSideSetSelectorName, tBulk.buckets(stk::topology::NODE_RANK));
        ////    xtk::size_t tNumEdges = stk::mesh::count_selected_entities(tSideSetSelectorName, tBulk.buckets(stk::topology::EDGE_RANK));
        ////    xtk::size_t tNumFaces = stk::mesh::count_selected_entities(tSideSetSelectorName, tBulk.buckets(stk::topology::FACE_RANK));
        ////    xtk::size_t tNumElems = stk::mesh::count_selected_entities(tSideSetSelectorName, tBulk.buckets(stk::topology::ELEM_RANK));
        //////
        ////    std::cout<<"\nNumber of nodes: "<<tNumNodes<<std::endl;
        ////    std::cout<<"Number of edges: "<<tNumEdges<<std::endl;
        ////    std::cout<<"Number of faces: "<<tNumFaces<<std::endl;
        ////    std::cout<<"Number of elems: "<<tNumElems<<std::endl;
        //
        ////    tNumNodes = stk::mesh::count_selected_entities(tSideSetSelectorOrdinal, tBulk.buckets(stk::topology::NODE_RANK));
        ////    tNumEdges = stk::mesh::count_selected_entities(tSideSetSelectorOrdinal, tBulk.buckets(stk::topology::EDGE_RANK));
        ////    tNumFaces = stk::mesh::count_selected_entities(tSideSetSelectorOrdinal, tBulk.buckets(stk::topology::FACE_RANK));
        ////    tNumElems = stk::mesh::count_selected_entities(tSideSetSelectorOrdinal, tBulk.buckets(stk::topology::ELEM_RANK));
        //
        ////    std::cout<<"\nNumber of nodes: "<<tNumNodes<<std::endl;
        ////    std::cout<<"Number of edges: "<<tNumEdges<<std::endl;
        ////    std::cout<<"Number of faces: "<<tNumFaces<<std::endl;
        ////    std::cout<<"Number of elems: "<<tNumElems<<std::endl;
        //
        //    stk::mesh::Part * tPartNodeSet = tMeta.get_part("nodelist_1");
        //
        //    stk::mesh::Selector tNodeSetSelector(tPartNodeSet);
        ////
        ////    tNumNodes = stk::mesh::count_selected_entities(tNodeSetSelector, tBulk.buckets(stk::topology::NODE_RANK));
        ////    tNumEdges = stk::mesh::count_selected_entities(tNodeSetSelector, tBulk.buckets(stk::topology::EDGE_RANK));
        ////    tNumFaces = stk::mesh::count_selected_entities(tNodeSetSelector, tBulk.buckets(stk::topology::FACE_RANK));
        ////    tNumElems = stk::mesh::count_selected_entities(tNodeSetSelector, tBulk.buckets(stk::topology::ELEM_RANK));
        //
        ////    std::cout<<"\nNumber of nodes: "<<tNumNodes<<std::endl;
        ////    std::cout<<"Number of edges: "<<tNumEdges<<std::endl;
        ////    std::cout<<"Number of faces: "<<tNumFaces<<std::endl;
        ////    std::cout<<"Number of elems: "<<tNumElems<<std::endl;
        //
        //
        //    const stk::mesh::BucketVector & tNodeBuckets = tBulk.buckets(stk::topology::NODE_RANK);
        //    const stk::mesh::BucketVector & tEdgeBuckets = tBulk.buckets(stk::topology::EDGE_RANK);
        //    const stk::mesh::BucketVector & tFaceBuckets = tBulk.buckets(stk::topology::FACE_RANK);
        //    const stk::mesh::BucketVector & tElemBuckets = tBulk.buckets(stk::topology::ELEM_RANK);
        //
        ////    std::cout<<"\nNumber of Node Buckets:" << tNodeBuckets.size() << std::endl;
        ////    std::cout<<"Number of Edge Buckets:" << tEdgeBuckets.size() << std::endl;
        ////    std::cout<<"Number of Face Buckets:" << tFaceBuckets.size() << std::endl;
        ////    std::cout<<"Number of Elem Buckets:" << tElemBuckets.size() << std::endl;
        //
        //    tMeshData->write_output_mesh(tMeshOutputFile);
    }

}    // namespace xtk
