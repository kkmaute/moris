/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Mesh_From_Data.cpp
 *
 */

#include <catch.hpp>
#include <iostream>

// MORIS project header files.
#include "algorithms.hpp"
#include "cl_MTK_Mesh.hpp" // MTK/src
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"

#include "cl_Communication_Tools.hpp" // COM/src

// Third-party header files.
#include <stk_io/StkMeshIoBroker.hpp>     // for StkMeshIoBroker
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities
#include <stk_mesh/base/MetaData.hpp>     // for MetaData
#include <stk_mesh/base/BulkData.hpp>     // for BulkData
#include <stk_mesh/base/Selector.hpp>     // for Selector
#include <stk_mesh/base/FEMHelpers.hpp>   // for Selector
#include "stk_io/DatabasePurpose.hpp"     // for DatabasePurpose::READ_MESH
#include "stk_mesh/base/CoordinateSystems.hpp" // for Cartesian
#include "stk_mesh/base/CreateFaces.hpp"  // for handling faces
#include "stk_mesh/base/CreateEdges.hpp"  // for handling faces
#include "stk_mesh/base/Bucket.hpp"       // for buckets
#include "stk_mesh/base/Field.hpp"    // for coordinates
#include "stk_mesh/base/GetEntities.hpp"    // for coordinates
#include "stk_mesh/base/FieldParallel.hpp"  // for handling parallel fields

// ----------------------------------------------------------------------------

namespace moris
{
namespace mtk
{

TEST_CASE( "Creating a 2D mesh from data in serial", "[Mesh_from_data_1]" )
            {
	if (par_size() == 1)
	{
    // Parallel
    uint p_rank = moris::par_rank();
    uint p_size = moris::par_size();

    if(p_rank == 0 && p_size == 1 ) // specify it is a serial test only
    {
        // Generate data for test
        uint aNumDim = 2;
        Matrix< DDRMat >  aCoords(6,2);
        aCoords(0,0) = 0.0, aCoords(0,1) = 0.0;
        aCoords(1,0) = 1.0, aCoords(1,1) = 0.0;
        aCoords(2,0) = 1.0, aCoords(2,1) = 1.0;
        aCoords(3,0) = 0.0, aCoords(3,1) = 1.0;
        aCoords(4,0) = 2.0, aCoords(4,1) = 0.0;
        aCoords(5,0) = 2.0, aCoords(5,1) = 1.0;
        Matrix< IdMat >     aElemConn( 2, 4 );

        SECTION( "using consecutive node and element ids" )
        {
            // 0D to 3D connectivity (node to element)
            aElemConn( 0, 0 ) = 1; aElemConn( 0, 1 ) = 2; aElemConn( 0, 2 ) = 3; aElemConn( 0, 3 ) = 4;
            aElemConn( 1, 0 ) = 2; aElemConn( 1, 1 ) = 5; aElemConn( 1, 2 ) = 6; aElemConn( 1, 3 ) = 3;

            Matrix< IdMat >  aElemLocaltoGlobal = {{1},{2}};

            // No need of an element map since elements in connectivity table are assumed to be contiguous

            // Create MORIS mesh using MTK database
            MtkMeshData aMeshData;
            aMeshData.CreateAllEdgesAndFaces = false;
            aMeshData.SpatialDim = &aNumDim;
            aMeshData.ElemConn(0)= &aElemConn;
            aMeshData.NodeCoords = &aCoords;
            aMeshData.LocaltoGlobalElemMap(0) = &aElemLocaltoGlobal;

            Mesh* tMesh2D_QUADs  = create_interpolation_mesh( MeshType::STK, aMeshData );

            // =============================
            // Testing basic functionalities
            // =============================

            uint tNumElements = tMesh2D_QUADs->get_num_elems();
            REQUIRE(moris::equal_to(tNumElements,2));

            uint tNumNodes    = tMesh2D_QUADs->get_num_nodes();
            REQUIRE(moris::equal_to(tNumNodes,6));

            Matrix< IdMat >  tAvailableNodeIDs = tMesh2D_QUADs->generate_unique_node_ids(4);
            REQUIRE(moris::equal_to(tAvailableNodeIDs(0,0),7));
            REQUIRE(moris::equal_to(tAvailableNodeIDs(1,0),8));
            REQUIRE(moris::equal_to(tAvailableNodeIDs(2,0),9));
            REQUIRE(moris::equal_to(tAvailableNodeIDs(3,0),10));

            delete tMesh2D_QUADs;
        }
    }
	}
 }
TEST_CASE( "Creating a 3D 2 element mesh from data in serial using non-consecutive node and element id maps")
{
    // Parallel
    uint p_rank = moris::par_rank();
    uint p_size = moris::par_size();

    if(p_rank == 0 && p_size == 1 ) // specify it is a serial test only
    {
                uint aNumDim = 3;
                Matrix< IdMat >  aElemConn = {{1000,2,4,38,543,6,8,77},{543,6,8,77,93,10,12,111}};
                Matrix< IdMat >  aNodeLocaltoGlobalNC = {{1000},{2},{38},{4},{543},{6},{77},{8},{93},{10},{111},{12}};
                Matrix< DDRMat >  aCoords   = {{0.0, 0.0, 0.0},
                                               {1.0, 0.0, 0.0},
                                               {0.0, 1.0, 0.0},
                                               {1.0, 1.0, 0.0},
                                               {0.0, 0.0, 1.0},
                                               {1.0, 0.0, 1.0},
                                               {0.0, 1.0, 1.0},
                                               {1.0, 1.0, 1.0},
                                               {0.0, 0.0, 2.0},
                                               {1.0, 0.0, 2.0},
                                               {0.0, 1.0, 2.0},
                                               {1.0, 1.0, 2.0}};
                Matrix< IdMat >  aNodeLocaltoGlobalnc = {{1000},{2},{38},{4},{543},{6},{77},{8},{93},{10},{111},{12}};
                Matrix< IdMat >  aElemLocaltoGlobalNC = {{51},{32}};

                // Create MORIS mesh using MTK database
                moris::mtk::MtkMeshData aMeshData;
                aMeshData.CreateAllEdgesAndFaces  = false;
                aMeshData.SpatialDim              = &aNumDim;
                aMeshData.ElemConn(0)             = &aElemConn;
                aMeshData.NodeCoords              = &aCoords;
                aMeshData.LocaltoGlobalElemMap(0) = &aElemLocaltoGlobalNC;
                aMeshData.LocaltoGlobalNodeMap    = &aNodeLocaltoGlobalNC;

                moris::mtk::Mesh* tMesh3D_HEXs = create_interpolation_mesh( MeshType::STK, aMeshData );

                // =============================
                // Testing basic functionalities
                // =============================
                uint tNumElems = tMesh3D_HEXs->get_num_elems();

                //TODO: ADD more tests here (previous implementation does not test an existing function)
                //
                REQUIRE(moris::equal_to(tNumElems,2));
                //                    REQUIRE(moris::equal_to(tNodeIDs(0),2));
                //                    REQUIRE(moris::equal_to(tNodeIDs(1),4));
                //                    REQUIRE(moris::equal_to(tNodeIDs(2),6));
                //                    REQUIRE(moris::equal_to(tNodeIDs(3),8));
                //                    REQUIRE(moris::equal_to(tNodeIDs(4),10));
                //                    REQUIRE(moris::equal_to(tNodeIDs(5),12));
                //                    REQUIRE(moris::equal_to(tNodeIDs(6),38));
                //                    REQUIRE(moris::equal_to(tNodeIDs(7),77));
                //                    REQUIRE(moris::equal_to(tNodeIDs(8),93));
                //                    REQUIRE(moris::equal_to(tNodeIDs(9),111));
                //                    REQUIRE(moris::equal_to(tNodeIDs(10),543));
                //                    REQUIRE(moris::equal_to(tNodeIDs(11),1000));
                //
                //                    REQUIRE(moris::equal_to(tElemIDs(0),32));
                //                    REQUIRE(moris::equal_to(tElemIDs(1),51));

                delete tMesh3D_HEXs;
    }
}

TEST_CASE( "with 2 block sets, 1 node set, and 1 side set","[Mesh_with_blocks]" )
{
    if( par_size() == 1 )
    {
    // Parallel
    uint p_rank = moris::par_rank();
    uint p_size = moris::par_size();

    if(p_rank == 0 && p_size == 1 ) // specify it is a serial test only
    {
        // Generate data for test
        uint aNumDim = 3;
        Matrix< IdMat >  aElemConn = {{1,2,4,3,5,6,8,7},{5,6,8,7,9,10,12,11}};
        Matrix< DDRMat >  aCoords   = {{0.0, 0.0, 0.0},
                                       {1.0, 0.0, 0.0},
                                       {0.0, 1.0, 0.0},
                                       {1.0, 1.0, 0.0},
                                       {0.0, 0.0, 1.0},
                                       {1.0, 0.0, 1.0},
                                       {0.0, 1.0, 1.0},
                                       {1.0, 1.0, 1.0},
                                       {0.0, 0.0, 2.0},
                                       {1.0, 0.0, 2.0},
                                       {0.0, 1.0, 2.0},
                                       {1.0, 1.0, 2.0}};
        Matrix< IdMat >  aNodeLocaltoGlobal = {{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12}};
        Matrix< IdMat >  aElemLocaltoGlobal = {{1},{2}};

        // Create 2  block sets (one over each element) and a node set that contains only 4 nodes
        // NOTE: A side set requires a two column matrix. The first column contains the ids of the elements to which
        // the faces are associated (for reference), and the second column represents the ordinal of the face,
        // which should go from 0 to n-1, being n the number of faces. numbering goes counterclockwise.

        // Initialize Sets information structure
        MtkSetsInfo tMtkMeshSets;

        // Declare block sets
        //////////////////////
        Matrix< IdMat >tCellIdsBS1({{1}});

        MtkBlockSetInfo tBlockSet1;
        tBlockSet1.mCellIdsInSet = &tCellIdsBS1;
        tBlockSet1.mBlockSetName = "blockset_1";
        tBlockSet1.mBlockSetTopo = CellTopology::HEX8;

        Matrix< IdMat >tCellIdsBS2({{2}});
        MtkBlockSetInfo tBlockSet2;
        tBlockSet2.mCellIdsInSet = &tCellIdsBS2;
        tBlockSet2.mBlockSetName = "blockset_2";
        tBlockSet2.mBlockSetTopo = CellTopology::HEX8;

        // Add block sets to mtk mesh sets
        tMtkMeshSets.add_block_set(&tBlockSet1);
        tMtkMeshSets.add_block_set(&tBlockSet2);

        // Declare side sets
        /////////////////////
        Matrix< IdMat > tElemIdsAndSideOrdsSS1  = { {1, 3},
                                                    {1, 4},
                                                    {1, 5},
                                                    {2, 1},
                                                    {2, 2} };

        MtkSideSetInfo tSideSetStruc;
        tSideSetStruc.mElemIdsAndSideOrds = &tElemIdsAndSideOrdsSS1;
        tSideSetStruc.mSideSetName        = "Sideset_1" ;

        // Add side side set to mesh sets
        tMtkMeshSets.add_side_set(&tSideSetStruc);

        // Declare node sets
        /////////////////////
        Matrix< IdMat > tNodeIdsNS1  = { {1}, {3}, {5}, {6} };

        MtkNodeSetInfo tNodeSet1;
        tNodeSet1.mNodeIds     = &tNodeIdsNS1;
        tNodeSet1.mNodeSetName = "Nodeset_1" ;

        // Add node set to Mtk mesh sets
        tMtkMeshSets.add_node_set(&tNodeSet1);

        MtkMeshData aMeshData;
        aMeshData.ElemConn = moris::Cell<Matrix < IdMat >*>(1);
        aMeshData.LocaltoGlobalElemMap = moris::Cell<Matrix < IdMat >*>(1);

        aMeshData.SpatialDim              = &aNumDim;
        aMeshData.ElemConn(0)             = &aElemConn;
        aMeshData.NodeCoords              = &aCoords;
        aMeshData.SetsInfo                = &tMtkMeshSets;
        aMeshData.LocaltoGlobalElemMap(0) = &aElemLocaltoGlobal;
        aMeshData.LocaltoGlobalNodeMap    = &aNodeLocaltoGlobal;

        moris::mtk::Mesh* tMesh = create_interpolation_mesh( MeshType::STK, aMeshData );

        // ========================
        // Testing sets information
        // ========================

        Matrix< IndexMat >  tBlockIndices1   = tMesh->get_set_entity_loc_inds( EntityRank::ELEMENT, "blockset_1" );
        Matrix< IndexMat >  tBlockIndices2   = tMesh->get_set_entity_loc_inds( EntityRank::ELEMENT, "blockset_2" );
        Matrix< IndexMat >  tNodeSetIndices1 = tMesh->get_set_entity_loc_inds( EntityRank::NODE,    "Nodeset_1" );
        Matrix< IndexMat >  tSideSetIndices1 = tMesh->get_set_entity_loc_inds( EntityRank::FACE,    "Sideset_1" );

        //        Matrix< DDUMat >  tNodesInBlockSet1 = tMesh.get_nodes_in_block_set( 1 );
//        Matrix< DDUMat >  tNodesInBlockSet2 = tMesh.get_nodes_in_block_set( 2 );
//        Matrix< DDUMat >  tNodesInNodeSet1  = tMesh.get_nodes_in_node_set( 1 );
//        Matrix< DDUMat >  tFacesInFaceSet1  = tMesh.get_faces_in_side_set( 1 );
//
//        Matrix< DDUMat >  tFacesInFaceSet1AndBlockSet1 = tMesh.get_intersected_entities_field_set(EntityRank::FACE, "Sideset_1", "blockset_1");
//        Matrix< DDUMat >  tFacesInFaceSet1AndBlockSet2 = tMesh.get_intersected_entities_field_set(EntityRank::FACE, "Sideset_1", "blockset_2");
//        Matrix< DDUMat >  tNodesInFaceSet1AndBlockSet2 = tMesh.get_intersected_entities_field_set(EntityRank::NODE, "Sideset_1", "blockset_2");
//
//        uint tNumFaces = tMesh.get_num_faces();
//        uint tNumEdges = tMesh.get_num_edges();
//
        REQUIRE(moris::equal_to(tBlockIndices1(0,0),0));
        REQUIRE(moris::equal_to(tBlockIndices2(0,0),1));

        Matrix< IdMat > tNodeSetIds1 = moris::mtk::convert_entity_indices_to_ids(tNodeSetIndices1,EntityRank::NODE,tMesh);
        REQUIRE(moris::equal_to(tNodeSetIds1(0),tNodeIdsNS1(0)));
        REQUIRE(moris::equal_to(tNodeSetIds1(1),tNodeIdsNS1(1)));
        REQUIRE(moris::equal_to(tNodeSetIds1(2),tNodeIdsNS1(2)));
        REQUIRE(moris::equal_to(tNodeSetIds1(3),tNodeIdsNS1(3)));
//
//        REQUIRE(moris::equal_to(tNodesInNodeSet1(0,0),1));
//        REQUIRE(moris::equal_to(tNodesInNodeSet1(1,0),3));
//        REQUIRE(moris::equal_to(tNodesInNodeSet1(2,0),5));
//        REQUIRE(moris::equal_to(tNodesInNodeSet1(3,0),6));
//
//        REQUIRE(moris::equal_to(tNodesInBlockSet1(0,0),1));
//        REQUIRE(moris::equal_to(tNodesInBlockSet1(1,0),2));
//        REQUIRE(moris::equal_to(tNodesInBlockSet1(2,0),3));
//        REQUIRE(moris::equal_to(tNodesInBlockSet1(3,0),4));
//        REQUIRE(moris::equal_to(tNodesInBlockSet1(4,0),5));
//        REQUIRE(moris::equal_to(tNodesInBlockSet1(5,0),6));
//        REQUIRE(moris::equal_to(tNodesInBlockSet1(6,0),7));
//        REQUIRE(moris::equal_to(tNodesInBlockSet1(7,0),8));
//
//        REQUIRE(moris::equal_to(tNodesInBlockSet2(0,0),5));
//        REQUIRE(moris::equal_to(tNodesInBlockSet2(1,0),6));
//        REQUIRE(moris::equal_to(tNodesInBlockSet2(2,0),7));
//        REQUIRE(moris::equal_to(tNodesInBlockSet2(3,0),8));
//        REQUIRE(moris::equal_to(tNodesInBlockSet2(4,0),9));
//        REQUIRE(moris::equal_to(tNodesInBlockSet2(5,0),10));
//        REQUIRE(moris::equal_to(tNodesInBlockSet2(6,0),11));
//        REQUIRE(moris::equal_to(tNodesInBlockSet2(7,0),12));
//
//        REQUIRE(moris::equal_to(tSideSetIndices1(0),14));
//        REQUIRE(moris::equal_to(tSideSetIndices1(1),15));
//        REQUIRE(moris::equal_to(tFacesInFaceSet1(2),16));
//        REQUIRE(moris::equal_to(tFacesInFaceSet1(3),22));
//        REQUIRE(moris::equal_to(tFacesInFaceSet1(4),23));
//
//        REQUIRE(moris::equal_to(tNumFaces,5));
//        REQUIRE(moris::equal_to(tNumEdges,0));
//
//        REQUIRE(moris::equal_to(tFacesInFaceSet1AndBlockSet1(0),14));
//        REQUIRE(moris::equal_to(tFacesInFaceSet1AndBlockSet1(1),15));
//        REQUIRE(moris::equal_to(tFacesInFaceSet1AndBlockSet1(2),16));
//
//        REQUIRE(moris::equal_to(tFacesInFaceSet1AndBlockSet2(0),22));
//        REQUIRE(moris::equal_to(tFacesInFaceSet1AndBlockSet2(1),23));
//
//        REQUIRE(moris::equal_to(tNodesInFaceSet1AndBlockSet2(0),5));
//        REQUIRE(moris::equal_to(tNodesInFaceSet1AndBlockSet2(2),7));
//        REQUIRE(moris::equal_to(tNodesInFaceSet1AndBlockSet2(3),8));
//        REQUIRE(moris::equal_to(tNodesInFaceSet1AndBlockSet2(5),11));

        delete tMesh;
    }
    }
}

TEST_CASE("parallel 4 element mesh","[PAR_MTK_FROM_DATA]")
{
    if(par_size() == 4)
    {
        std::string tPrefix = std::getenv("MORISOUTPUT");
        std::string tMeshOutputFile = tPrefix + "/mtk_par_mtk_from_data.e";

        bool tAura = false;
        bool tCreateEdgesAndFaces = true;
        if(par_rank() == 0)
        {
            uint aNumDim = 3;
            Matrix< IdMat >  aElemConn = {{1,10,11,2,4,13,14,5}};
            Matrix< IdMat >  aNodeLocaltoGlobalNC = {{1,2,5,4,10,11,14,13}};
            Matrix< DDRMat >  aCoords   = {{0.0,0.0,0.0},
                                           {1.0,0.0,0.0},
                                           {1.0,1.0,0.0},
                                           {0.0,1.0,0.0},
                                           {0.0,0.0,1.0},
                                           {1.0,0.0,1.0},
                                           {1.0,1.0,1.0},
                                           {0.0,1.0,1.0}};
            Matrix< IdMat >  aElemLocaltoGlobalNC = {{1}};
            Matrix< IdMat > aNodeSharedProcs(8,3,MORIS_ID_MAX);

            // Node 1 is not shared
            aNodeSharedProcs(1,0) = 1;
            aNodeSharedProcs(2,0) = 1; aNodeSharedProcs(2,1) = 2; aNodeSharedProcs(2,2) = 3;
            aNodeSharedProcs(3,0) = 2;
            // Node 10 is not shared
            aNodeSharedProcs(5,0) = 1;
            aNodeSharedProcs(6,0) = 1; aNodeSharedProcs(6,1) = 2; aNodeSharedProcs(6,2) = 3;
            aNodeSharedProcs(7,0) = 2;

            // Create MORIS mesh using MTK database
            moris::mtk::MtkMeshData aMeshData;
            aMeshData.CreateAllEdgesAndFaces  = tCreateEdgesAndFaces;
            aMeshData.AutoAuraOptionInSTK     = tAura;
            aMeshData.SpatialDim              = &aNumDim;
            aMeshData.ElemConn(0)             = &aElemConn;
            aMeshData.NodeCoords              = &aCoords;
            aMeshData.LocaltoGlobalElemMap(0) = &aElemLocaltoGlobalNC;
            aMeshData.LocaltoGlobalNodeMap    = &aNodeLocaltoGlobalNC;
            aMeshData.NodeProcsShared         = &aNodeSharedProcs;

            moris::mtk::Mesh* tParMesh = create_interpolation_mesh( MeshType::STK, aMeshData );
            tParMesh->create_output_mesh(tMeshOutputFile);

            delete tParMesh;

        }
        else if(par_rank() == 1)
        {
            uint aNumDim = 3;
            Matrix< IdMat >  aElemConn = {{2,11,12,3,5,14,15,6}};
            Matrix< IdMat >  aNodeLocaltoGlobalNC = {{2,5,11,14,3,6,12,15}};
            Matrix< DDRMat >  aCoords   = {{1.0,0.0,0.0},
                                           {1.0,1.0,0.0},
                                           {1.0,0.0,1.0},
                                           {1.0,1.0,1.0},
                                           {2.0,0.0,0.0},
                                           {2.0,1.0,0.0},
                                           {2.0,0.0,1.0},
                                           {2.0,1.0,1.0}};
            Matrix< IdMat >  aElemLocaltoGlobalNC = {{2}};
            Matrix< IdMat > aNodeSharedProcs(8,3,MORIS_ID_MAX);

            aNodeSharedProcs(0,0) = 0;
            aNodeSharedProcs(1,0) = 0; aNodeSharedProcs(1,1) = 2; aNodeSharedProcs(1,2) = 3;
            aNodeSharedProcs(2,0) = 0;
            aNodeSharedProcs(3,0) = 0; aNodeSharedProcs(3,1) = 2; aNodeSharedProcs(3,2) = 3;
            aNodeSharedProcs(4,0) = MORIS_ID_MAX;
            aNodeSharedProcs(5,0) = 3;
            aNodeSharedProcs(6,0) = MORIS_ID_MAX;
            aNodeSharedProcs(7,0) = 3;

            // Create MORIS mesh using MTK database
            moris::mtk::MtkMeshData aMeshData;
            aMeshData.CreateAllEdgesAndFaces  = tCreateEdgesAndFaces;
            aMeshData.AutoAuraOptionInSTK     = tAura;
            aMeshData.SpatialDim              = &aNumDim;
            aMeshData.ElemConn(0)             = &aElemConn;
            aMeshData.NodeCoords              = &aCoords;
            aMeshData.LocaltoGlobalElemMap(0) = &aElemLocaltoGlobalNC;
            aMeshData.LocaltoGlobalNodeMap    = &aNodeLocaltoGlobalNC;
            aMeshData.NodeProcsShared         = &aNodeSharedProcs;

            moris::mtk::Mesh* tParMesh = create_interpolation_mesh( MeshType::STK, aMeshData );
            tParMesh->create_output_mesh(tMeshOutputFile);
            delete tParMesh;
        }
        else if(par_rank() == 2)
        {
            uint aNumDim = 3;
            Matrix< IdMat >  aElemConn = {{4,13,14,5,7,16,17,8}};
            Matrix< IdMat >  aNodeLocaltoGlobalNC = {{5,4,14,13,8,7,17,16}};
            Matrix< DDRMat >  aCoords   = {{1.0,1.0,0.0},
                                           {0.0,1.0,0.0},
                                           {1.0,1.0,1.0},
                                           {0.0,1.0,1.0},
                                           {1.0,2.0,0.0},
                                           {0.0,2.0,0.0},
                                           {1.0,2.0,1.0},
                                           {0.0,2.0,1.0}};
            Matrix< IdMat >  aElemLocaltoGlobalNC = {{4}};
            Matrix< IdMat > aNodeSharedProcs(8,3,MORIS_ID_MAX);

            aNodeSharedProcs(0,0) = 0; aNodeSharedProcs(0,1) = 1; aNodeSharedProcs(0,2) = 3;
            aNodeSharedProcs(1,0) = 0;
            aNodeSharedProcs(2,0) = 0; aNodeSharedProcs(2,1) = 1; aNodeSharedProcs(2,2) = 3;
            aNodeSharedProcs(3,0) = 0;
            aNodeSharedProcs(4,0) = 3;
            aNodeSharedProcs(5,0) = MORIS_ID_MAX;
            aNodeSharedProcs(6,0) = 3;
            aNodeSharedProcs(7,0) = MORIS_ID_MAX;

            // Create MORIS mesh using MTK database
            moris::mtk::MtkMeshData aMeshData;
            aMeshData.CreateAllEdgesAndFaces  = tCreateEdgesAndFaces;
            aMeshData.AutoAuraOptionInSTK     = tAura;
            aMeshData.SpatialDim              = &aNumDim;
            aMeshData.ElemConn(0)             = &aElemConn;
            aMeshData.NodeCoords              = &aCoords;
            aMeshData.LocaltoGlobalElemMap(0) = &aElemLocaltoGlobalNC;
            aMeshData.LocaltoGlobalNodeMap    = &aNodeLocaltoGlobalNC;
            aMeshData.NodeProcsShared         = &aNodeSharedProcs;

            moris::mtk::Mesh* tParMesh = create_interpolation_mesh( MeshType::STK, aMeshData );
            tParMesh->create_output_mesh(tMeshOutputFile);
            delete tParMesh;
        }

        else if(par_rank() == 3)
        {
            uint aNumDim = 3;
            Matrix< IdMat >  aElemConn = {{5,14,15,6,8,17,18,9}};
            Matrix< IdMat >  aNodeLocaltoGlobalNC = {{5,6,9,8,14,15,18,17}};
            Matrix< DDRMat >  aCoords   = {{1.0,1.0,0.0},
                                           {2.0,1.0,0.0},
                                           {2.0,2.0,0.0},
                                           {1.0,2.0,0.0},
                                           {1.0,1.0,1.0},
                                           {2.0,1.0,1.0},
                                           {2.0,2.0,1.0},
                                           {1.0,2.0,1.0}};
            Matrix< IdMat >  aElemLocaltoGlobalNC = {{3}};
            Matrix< IdMat > aNodeSharedProcs(8,3,MORIS_ID_MAX);

            aNodeSharedProcs(0,0) = 0; aNodeSharedProcs(0,1) = 1; aNodeSharedProcs(0,2) = 2;
            aNodeSharedProcs(1,0) = 1;
            aNodeSharedProcs(2,0) = MORIS_ID_MAX;
            aNodeSharedProcs(3,0) = 2;
            aNodeSharedProcs(4,0) = 0; aNodeSharedProcs(4,1) = 1; aNodeSharedProcs(4,2) = 2;
            aNodeSharedProcs(5,0) = 1;
            aNodeSharedProcs(6,0) = MORIS_ID_MAX;
            aNodeSharedProcs(7,0) = 2;

            // Create MORIS mesh using MTK database
            moris::mtk::MtkMeshData aMeshData;
            aMeshData.CreateAllEdgesAndFaces  = tCreateEdgesAndFaces;
            aMeshData.AutoAuraOptionInSTK     = tAura;
            aMeshData.SpatialDim              = &aNumDim;
            aMeshData.ElemConn(0)             = &aElemConn;
            aMeshData.NodeCoords              = &aCoords;
            aMeshData.LocaltoGlobalElemMap(0) = &aElemLocaltoGlobalNC;
            aMeshData.LocaltoGlobalNodeMap    = &aNodeLocaltoGlobalNC;
            aMeshData.NodeProcsShared         = &aNodeSharedProcs;

            moris::mtk::Mesh* tParMesh = create_interpolation_mesh( MeshType::STK, aMeshData );
            tParMesh->create_output_mesh(tMeshOutputFile);
            delete tParMesh;
        }
    }
}

}
}

