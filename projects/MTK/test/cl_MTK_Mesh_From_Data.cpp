/*
 * cl_MTK_Mesh_From_Data.cpp
 *
 *  Created on: Sep 28, 2018
 *      Author: barrera/doble
 */
// Third-party header files.
#include <catch.hpp>
#include <iostream>

// MORIS project header files.
#include "algorithms.hpp"
#include "cl_MTK_Mesh.hpp" // MTK/src
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_Mesh_Factory.hpp"

#include "cl_Communication_Tools.hpp" // COM/src


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

            Mesh* tMesh2D_QUADs  = create_mesh( MeshType::STK, aMeshData );

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

                moris::mtk::Mesh* tMesh3D_HEXs = create_mesh( MeshType::STK, aMeshData );

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
    }
}
TEST_CASE( "Creating a 3D 2 element mesh from data in serial ")
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

        SECTION( "with 2 block sets, 1 node set, and 1 side set" )
        {
            // Create 2  block sets (one over each element) and a node set that contains only 4 nodes
            // NOTE: A side set requires a two column matrix. The first column contains the ids of the elements to which
            // the faces are associated (for reference), and the second column represents the ordinal of the face,
            // which should go from 0 to n-1, being n the number of faces. numbering goes counterclockwise.

            // Declare block sets
            //////////////////////
//            Matrix<IndexMat> tBlockSetsPartOwners = { {0},{1} };
//
//            MtkBlockSetsInfo tBlockSetStruc;
//            tBlockSetStruc.BSetInds = &tBlockSetsPartOwners;
//            tBlockSetStruc.BSetNames   = { "blockset_1", "blockset_2" };
//
//            // Declare side sets
//            /////////////////////
//            Matrix< IdMat > tSideset_1  = { {1, 3}, {1, 4}, {1, 5}, {2, 1}, {2, 2} };
//            moris::Cell< Matrix< IdMat >  > tSideSetsInfo = { tSideset_1 };
//
//            MtkSideSetsInfo tSideSetStruc;
//            tSideSetStruc.ElemIdsAndSideOrds = &tSideSetsInfo;
//            tSideSetStruc.SSetNames   = { "Sideset_1" };
//
//            // Declare node sets
//            /////////////////////
//            Matrix<IdMat> tNodeSet_1  = { {1}, {3}, {5}, {6} };
//            moris::Cell< Matrix< IdMat >  > tNodeSetsEntIds = { tNodeSet_1 };
//
//            MtkNodeSetsInfo tNodeSetStruc;
//            tNodeSetStruc.EntIds = &tNodeSetsEntIds;
//            tNodeSetStruc.NSetNames = { "Nodeset_1" };
//
//            // Create MORIS mesh using MTK database
//            ///////////////////////////////////////
//            MtkSetsInfo aMeshSets;
//            aMeshSets.NodeSetsInfo = &tNodeSetStruc;
//            aMeshSets.SideSetsInfo   = &tSideSetStruc;
//            aMeshSets.BlockSetsInfo   = &tBlockSetStruc;
//
//            MtkMeshData aMeshData;
//            aMeshData.SpatialDim              = & aNumDim;
//            aMeshData.ElemConn(0)             = & aElemConn;
//            aMeshData.NodeCoords              = & aCoords;
//            aMeshData.SetsInfo                = & aMeshSets;
//            aMeshData.LocaltoGlobalElemMap(0) = & aElemLocaltoGlobal;
//            aMeshData.LocaltoGlobalNodeMap    = & aNodeLocaltoGlobal;
//
//            std::cout<<" create mesh w/ block "<<std::endl;
//            moris::mtk::Mesh* tMesh = create_mesh( MeshType::STK, aMeshData );

            // ========================
            // Testing sets information
            // ========================
            // TODO: ADD BLOCK INFORMATION ACCESS TO STK
//            Matrix< DDUMat >  tBlock1 = tMesh.get_set_entity_ids( EntityRank::ELEMENT, "blockset_1" );
//            Matrix< DDUMat >  tBlock2 = tMesh.get_set_entity_ids( EntityRank::ELEMENT, "blockset_2" );
//            Matrix< DDUMat >  tNodeSet1 = tMesh.get_set_entity_ids( EntityRank::NODE, "Nodeset_1" );
//            Matrix< DDUMat >  tSideSet1 = tMesh.get_set_entity_ids( EntityRank::FACE, "Sideset_1" );
//
//            Matrix< DDUMat >  tNodesInBlockSet1 = tMesh.get_nodes_in_block_set( 1 );
//            Matrix< DDUMat >  tNodesInBlockSet2 = tMesh.get_nodes_in_block_set( 2 );
//            Matrix< DDUMat >  tNodesInNodeSet1 = tMesh.get_nodes_in_node_set( 1 );
//            Matrix< DDUMat >  tFacesInFaceSet1 = tMesh.get_faces_in_side_set( 1 );
//
//            uint tNumFaces = tMesh.get_num_faces();
//            uint tNumEdges = tMesh.get_num_edges();
//
//            REQUIRE(moris::equal_to(tBlock1(0,0),1));
//            REQUIRE(moris::equal_to(tBlock2(0,0),2));
//
//            REQUIRE(moris::equal_to(tNodeSet1(0,0),1));
//            REQUIRE(moris::equal_to(tNodeSet1(1,0),3));
//            REQUIRE(moris::equal_to(tNodeSet1(2,0),5));
//            REQUIRE(moris::equal_to(tNodeSet1(3,0),6));
//
//            REQUIRE(moris::equal_to(tNodesInNodeSet1(0,0),1));
//            REQUIRE(moris::equal_to(tNodesInNodeSet1(1,0),3));
//            REQUIRE(moris::equal_to(tNodesInNodeSet1(2,0),5));
//            REQUIRE(moris::equal_to(tNodesInNodeSet1(3,0),6));
//
//            REQUIRE(moris::equal_to(tNodesInBlockSet1(0,0),1));
//            REQUIRE(moris::equal_to(tNodesInBlockSet1(1,0),2));
//            REQUIRE(moris::equal_to(tNodesInBlockSet1(2,0),3));
//            REQUIRE(moris::equal_to(tNodesInBlockSet1(3,0),4));
//            REQUIRE(moris::equal_to(tNodesInBlockSet1(4,0),5));
//            REQUIRE(moris::equal_to(tNodesInBlockSet1(5,0),6));
//            REQUIRE(moris::equal_to(tNodesInBlockSet1(6,0),7));
//            REQUIRE(moris::equal_to(tNodesInBlockSet1(7,0),8));
//
//            REQUIRE(moris::equal_to(tNodesInBlockSet2(0,0),5));
//            REQUIRE(moris::equal_to(tNodesInBlockSet2(1,0),6));
//            REQUIRE(moris::equal_to(tNodesInBlockSet2(2,0),7));
//            REQUIRE(moris::equal_to(tNodesInBlockSet2(3,0),8));
//            REQUIRE(moris::equal_to(tNodesInBlockSet2(4,0),9));
//            REQUIRE(moris::equal_to(tNodesInBlockSet2(5,0),10));
//            REQUIRE(moris::equal_to(tNodesInBlockSet2(6,0),11));
//            REQUIRE(moris::equal_to(tNodesInBlockSet2(7,0),12));
//
//            REQUIRE(moris::equal_to(tSideSet1(0),14));
//            REQUIRE(moris::equal_to(tSideSet1(1),15));
//            REQUIRE(moris::equal_to(tFacesInFaceSet1(2),16));
//            REQUIRE(moris::equal_to(tFacesInFaceSet1(3),22));
//            REQUIRE(moris::equal_to(tFacesInFaceSet1(4),23));
//
//            REQUIRE(moris::equal_to(tNumFaces,5));
//            REQUIRE(moris::equal_to(tNumEdges,0));
//
//            REQUIRE(moris::equal_to(tFacesInFaceSet1AndBlockSet1(0),14));
//            REQUIRE(moris::equal_to(tFacesInFaceSet1AndBlockSet1(1),15));
//            REQUIRE(moris::equal_to(tFacesInFaceSet1AndBlockSet1(2),16));
//
//            REQUIRE(moris::equal_to(tFacesInFaceSet1AndBlockSet2(0),22));
//            REQUIRE(moris::equal_to(tFacesInFaceSet1AndBlockSet2(1),23));
//
//            REQUIRE(moris::equal_to(tNodesInFaceSet1AndBlockSet2(0),5));
//            REQUIRE(moris::equal_to(tNodesInFaceSet1AndBlockSet2(2),7));
//            REQUIRE(moris::equal_to(tNodesInFaceSet1AndBlockSet2(3),8));
//            REQUIRE(moris::equal_to(tNodesInFaceSet1AndBlockSet2(5),11));
        }



    } // if end (par check)

} // test end

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

        moris::mtk::Mesh* tMesh = create_mesh( MeshType::STK, aMeshData );

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

        REQUIRE(moris::equal_to(tNodeSetIndices1(0,0),0));
        REQUIRE(moris::equal_to(tNodeSetIndices1(1,0),2));
        REQUIRE(moris::equal_to(tNodeSetIndices1(2,0),4));
        REQUIRE(moris::equal_to(tNodeSetIndices1(3,0),5));
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


TEST_CASE("parallel test 4 element cluster","[PAR_MTK_FROM_DATA]")
{
    if(par_size() == 4)
    {
        std::string tPrefix = std::getenv("MORISOUTPUT");
        std::string tMeshOutputFile = tPrefix + "/mtk_par_mtk_from_data.e";
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
            aMeshData.CreateAllEdgesAndFaces  = false;
            aMeshData.AutoAuraOptionInSTK     = true;
            aMeshData.SpatialDim              = &aNumDim;
            aMeshData.ElemConn(0)             = &aElemConn;
            aMeshData.NodeCoords              = &aCoords;
            aMeshData.LocaltoGlobalElemMap(0) = &aElemLocaltoGlobalNC;
            aMeshData.LocaltoGlobalNodeMap    = &aNodeLocaltoGlobalNC;
            aMeshData.NodeProcsShared         = &aNodeSharedProcs;

            moris::mtk::Mesh* tParMesh = create_mesh( MeshType::STK, aMeshData );
            tParMesh->create_output_mesh(tMeshOutputFile);
            std::cout<<"Mesh outputted to: "<<tMeshOutputFile<<std::endl;
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
            aMeshData.CreateAllEdgesAndFaces  = false;
            aMeshData.SpatialDim              = &aNumDim;
            aMeshData.ElemConn(0)             = &aElemConn;
            aMeshData.NodeCoords              = &aCoords;
            aMeshData.LocaltoGlobalElemMap(0) = &aElemLocaltoGlobalNC;
            aMeshData.LocaltoGlobalNodeMap    = &aNodeLocaltoGlobalNC;
            aMeshData.NodeProcsShared         = &aNodeSharedProcs;

            moris::mtk::Mesh* tParMesh = create_mesh( MeshType::STK, aMeshData );
            tParMesh->create_output_mesh(tMeshOutputFile);
            std::cout<<"Num elements = "<<tParMesh->get_num_entities(EntityRank::ELEMENT)<<std::endl;
            std::cout<<"Num nodes = "<<tParMesh->get_num_entities(EntityRank::NODE)<<std::endl;

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
            aMeshData.CreateAllEdgesAndFaces  = false;
            aMeshData.SpatialDim              = &aNumDim;
            aMeshData.ElemConn(0)             = &aElemConn;
            aMeshData.NodeCoords              = &aCoords;
            aMeshData.LocaltoGlobalElemMap(0) = &aElemLocaltoGlobalNC;
            aMeshData.LocaltoGlobalNodeMap    = &aNodeLocaltoGlobalNC;
            aMeshData.NodeProcsShared         = &aNodeSharedProcs;

            moris::mtk::Mesh* tParMesh = create_mesh( MeshType::STK, aMeshData );
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
            aMeshData.CreateAllEdgesAndFaces  = false;
            aMeshData.SpatialDim              = &aNumDim;
            aMeshData.ElemConn(0)             = &aElemConn;
            aMeshData.NodeCoords              = &aCoords;
            aMeshData.LocaltoGlobalElemMap(0) = &aElemLocaltoGlobalNC;
            aMeshData.LocaltoGlobalNodeMap    = &aNodeLocaltoGlobalNC;
            aMeshData.NodeProcsShared         = &aNodeSharedProcs;

            moris::mtk::Mesh* tParMesh = create_mesh( MeshType::STK, aMeshData );
            tParMesh->create_output_mesh(tMeshOutputFile);
            delete tParMesh;
        }

    }
}


//    }
//

//

//
//            SECTION( "with multiple sets and fields allocated to sets instead of the entire mesh")
//            {
//                // Create 2  block sets (one over each element) and a node set that contains only 4 nodes
//                // NOTE: A side set requires a two column matrix. The first column contains the ids of the elements to which
//                // the faces are associated (for reference), and the second column represents the ordinal of the face,
//                // which should go from 0 to n-1, being n the number of faces. numbering goes counterclockwise.
//
//                // Declare side sets
//                /////////////////////
//                Mat<uint> tSideset_1  = { {1, 5}, {2, 0} };
//                Cell< Matrix< DDUMat >  > tSideSetsInfo = { tSideset_1 };
//
//                MtkSideSetsInfo tSideSetStruc;
//                tSideSetStruc.ElemIdsAndSideOrds = &tSideSetsInfo;
//                tSideSetStruc.SSetNames   = { "Sideset_1" };
//
//                // Declare node sets
//                /////////////////////
//                Mat<uint> tNodeSet_1  = { {4}, {5}, {7} };
//                Cell< Matrix< DDUMat >  > tNodeSetsEntIds = { tNodeSet_1 };
//
//                MtkNodeSetsInfo tNodeSetStruc;
//                tNodeSetStruc.EntIds = &tNodeSetsEntIds;
//                tNodeSetStruc.NSetNames = { "Nodeset_1" };
//
//                // Declare fields
//                /////////////////
//
//                // Create nodal and elemental fields of reals
//                Matrix< DDRMat >  tElemFieldData_1      = {{123.45},{678.9}};
//                Matrix< DDRMat >  tNodalFieldData_1     = {{-1.2},{-2.3},{-3.4}};
//                Matrix< DDRMat >  tSideFieldNodalData_1 = {{555.555}};
//                Matrix< DDRMat >  tSideFieldFaceData_1  = {{111.11},{222.22}};
//
//                Cell < Matrix< DDRMat >  > aFieldData      = { tElemFieldData_1, tNodalFieldData_1, tSideFieldNodalData_1, tSideFieldFaceData_1 };
//                // Field set owner. If not provided, it will use the universal part.
//                Cell < std::string > aSetOwnerName   = { "", "NodeSet_1", "Sideset_1", "Sideset_1" };
//                Cell < std::string > aFieldName      = { "ElementField_1", "NodesetField_1", "SidesetField_1", "SidesetField_2" };
//                Cell < enum EntityRank > aFieldRanks = { EntityRank::ELEMENT, EntityRank::NODE, EntityRank::NODE, EntityRank::FACE };
//
//                // Create MORIS mesh using MTK database
//                ///////////////////////////////////////
//                MtkSetsInfo aMeshSets;
//                aMeshSets.NodeSetsInfo = &tNodeSetStruc;
//                aMeshSets.SideSetsInfo   = &tSideSetStruc;
//
//                moris::MtkFieldsInfo aFieldsInfo;
//                aFieldsInfo.FieldsData = &aFieldData;
//                aFieldsInfo.FieldsName = aFieldName;
//                aFieldsInfo.FieldsRank = aFieldRanks;
//                aFieldsInfo.SetsOwner  = &aSetOwnerName;
//
//                moris::MtkMeshData aMeshData;
//                aMeshData.SpatialDim = &aNumDim;
//                aMeshData.ElemConn   = &aElemConn;
//                aMeshData.NodeCoords = &aCoords;
//                aMeshData.SetsInfo   = &aMeshSets;
//                aMeshData.FieldsInfo = &aFieldsInfo;
//                aMeshData.LocaltoGlobalElemMap = &aElemLocaltoGlobal;
//                aMeshData.LocaltoGlobalNodeMap = &aNodeLocaltoGlobal;
//
//                moris::mesh tMesh( MeshType::MTK, aMeshData );
//
//                // ========================
//                // Testing outputting field
//                // ========================
//
//                Matrix< DDRMat >  tField_1  = tMesh.get_field_values( EntityRank::ELEMENT, "ElementField_1" );
//                Matrix< DDRMat >  tField_2  = tMesh.get_field_values( EntityRank::NODE, "NodesetField_1" );
//                Matrix< DDRMat >  tField_3  = tMesh.get_field_values( EntityRank::NODE, "SidesetField_1" );
//                Matrix< DDRMat >  tField_4  = tMesh.get_field_values( EntityRank::FACE, "SidesetField_2" );
//
//                REQUIRE(moris::equal_to(tField_1(0,0),123.45));
//                REQUIRE(moris::equal_to(tField_1(1,0),678.9));
//
//                REQUIRE(moris::equal_to(tField_2(0,0),-1.2));
//                REQUIRE(moris::equal_to(tField_2(1,0),-2.3));
//                REQUIRE(moris::equal_to(tField_2(2,0),-3.4));
//
//                REQUIRE(moris::equal_to(tField_3(0,0),555.555));
//                REQUIRE(moris::equal_to(tField_3(1,0),555.555));
//                REQUIRE(moris::equal_to(tField_3(2,0),555.555));
//                REQUIRE(moris::equal_to(tField_3(3,0),555.555));
//                REQUIRE(moris::equal_to(tField_3(4,0),555.555));
//
//                REQUIRE(moris::equal_to(tField_4(0,0),111.11));
//                REQUIRE(moris::equal_to(tField_4(1,0),222.22));
//            }
//        }
//    }
//
//    SECTION( "Creating a 3D 3 element mesh in parallel using element processor owner list")
//    {
//        if( p_size == 2 ) // specify it is a 2 processor test
//        {
//            uint aNumDim = 3;
//            Mat<uint >  aElemProcs;
//            Mat<uint >  aElemConn;
//            Mat<real >  aNodeCoords;
//            Mat<uint>    aElemLocaltoGlobal;
//            Mat<uint>    aNodeLocaltoGlobal;
//
//            if ( p_rank == 0 )
//            {
//                // Generate data for test
//                Mat< uint>  tElemProcsDummy  = {{0}, {0}, {1}};
//                Matrix< DDUMat >  tElemConnDummy   = {{1000,2,4,3000,5,6,8,7000},
//                        {5,6,8,7000,9000,10,12,11000},
//                        {9000,10,12,11000,13,14,16,15000}};
//                Matrix< DDRMat >  tNodeCoordsDummy = {{0.0, 0.0, 0.0},
//                        {1.0, 0.0, 0.0},
//                        {0.0, 1.0, 0.0},
//                        {1.0, 1.0, 0.0},
//                        {0.0, 0.0, 1.0},
//                        {1.0, 0.0, 1.0},
//                        {0.0, 1.0, 1.0},
//                        {1.0, 1.0, 1.0},
//                        {0.0, 0.0, 2.0},
//                        {1.0, 0.0, 2.0},
//                        {0.0, 1.0, 2.0},
//                        {1.0, 1.0, 2.0},
//                        {0.0, 0.0, 3.0},
//                        {1.0, 0.0, 3.0},
//                        {0.0, 1.0, 3.0},
//                        {1.0, 1.0, 3.0}};
//                Matrix< DDUMat >  tElemLocaltoGlobalDummy = {{15},{7},{31}};
//                Matrix< DDUMat >  tNodeLocaltoGlobalDummy = {{1000},{2},{3000},{4},{5},
//                        {6},{7000},{8},{9000},{10},
//                        {11000},{12},{13},{14},{15000},{16}};
//                aElemProcs         = tElemProcsDummy;
//                aElemConn          = tElemConnDummy;
//                aNodeCoords        = tNodeCoordsDummy;
//                aElemLocaltoGlobal = tElemLocaltoGlobalDummy;
//                aNodeLocaltoGlobal = tNodeLocaltoGlobalDummy;
//
//            }
//            else
//            {
//
//                // Generate data for test
//                moris::Mat< moris::uint>  tElemProcsDummy  = {{0}, {1}};
//                moris::Mat< moris::uint > tElemConnDummy = {{5,6,8,7000,9000,10,12,11000},
//                        {9000,10,12,11000,13,14,16,15000}};
//                moris::Mat< moris::real > tNodeCoordsDummy = {{0.0, 0.0, 1.0},
//                        {1.0, 0.0, 1.0},
//                        {0.0, 1.0, 1.0},
//                        {1.0, 1.0, 1.0},
//                        {0.0, 0.0, 2.0},
//                        {1.0, 0.0, 2.0},
//                        {0.0, 1.0, 2.0},
//                        {1.0, 1.0, 2.0},
//                        {0.0, 0.0, 3.0},
//                        {1.0, 0.0, 3.0},
//                        {0.0, 1.0, 3.0},
//                        {1.0, 1.0, 3.0}};
//                moris::Mat< moris::uint > tElemLocaltoGlobalDummy = {{7},{31}};
//                moris::Mat< moris::uint > tNodeLocaltoGlobalDummy = {{5},{6},{7000},{8},{9000},{10},
//                        {11000},{12},{13},{14},{15000},{16}};
//
//                aElemProcs         = tElemProcsDummy;
//                aElemConn          = tElemConnDummy;
//                aNodeCoords        = tNodeCoordsDummy;
//                aElemLocaltoGlobal = tElemLocaltoGlobalDummy;
//                aNodeLocaltoGlobal = tNodeLocaltoGlobalDummy;
//            }
//
//            SECTION( "mesh without sets or fields")
//            {
//                moris::MtkMeshData aMeshData;
//                aMeshData.SpatialDim   = &aNumDim;
//                aMeshData.ElemConn     = &aElemConn;
//                aMeshData.NodeCoords   = &aNodeCoords;
//                aMeshData.EntProcOwner = &aElemProcs;
//                aMeshData.CreateAllEdgesAndFaces = true;
//                aMeshData.LocaltoGlobalElemMap   = &aElemLocaltoGlobal;
//                aMeshData.LocaltoGlobalNodeMap   = &aNodeLocaltoGlobal;
//
//                moris::mesh tParallelMesh( MeshType::MTK, aMeshData );
//
//                // ========================================
//                // Testing get_entities_owned_current_proc
//                // ========================================
//
//                Matrix< DDUMat >  tElemsOwned  = tParallelMesh.get_entities_owned_current_proc(EntityRank::ELEMENT);
//                Matrix< DDUMat >  tNodesOwned  = tParallelMesh.get_entities_owned_current_proc(EntityRank::NODE);
//                Matrix< DDUMat >  tProcsSharedByEntityId  = tParallelMesh.get_procs_sharing_entity_by_id(9000,EntityRank::NODE);
//
//                // Index = 5-> location of node with id 9000 for processor 1.
//                Matrix< DDUMat >  tProcsSharedByEntityInd = tParallelMesh.get_procs_sharing_entity_by_index(5,EntityRank::NODE);
//
//                uint tEntityOwner = tParallelMesh.parallel_owner_rank_by_entity_index(5,EntityRank::NODE);
//
//                // Access maps (all of the local to global maps for all entity types).
//                Mat < uint > lcl2GblNodeMap = tParallelMesh.get_nodal_local_map();
//                Mat < uint > lcl2GblElemMap = tParallelMesh.get_elemental_local_map();
//                Mat < uint > lcl2GblEdgeMap = tParallelMesh.get_edge_local_map();
//                Mat < uint > lcl2GblFaceMap = tParallelMesh.get_face_local_map();
//
//                // Access maps (all of the local to global maps for all entity types).
//                Mat < uint > lcl2GblNodeOwnerProc = tParallelMesh.get_nodal_owner_proc_map();
//                Mat < uint > lcl2GblElemOwnerProc = tParallelMesh.get_elemental_owner_proc_map();
//                Mat < uint > lcl2GblEdgeOwnerProc = tParallelMesh.get_edge_owner_proc_map();
//                Mat < uint > lcl2GblFaceOwnerProc = tParallelMesh.get_face_owner_proc_map();
//
//                Cell < Cell < uint > > tNodesSharedPerProc = tParallelMesh.get_nodes_shared_processors();
//                Cell < Cell < uint > > tElemsSharedPerProc = tParallelMesh.get_elements_shared_processors();
//                Cell < Cell < uint > > tEdgesSharedPerProc = tParallelMesh.get_edges_shared_processors();
//                Cell < Cell < uint > > tFacesSharedPerProc = tParallelMesh.get_faces_shared_processors();
//
//                //Check entity owner
//                REQUIRE(moris::equal_to(tEntityOwner,0));
//
//                if( p_rank == 0 )
//                {
//                    // Check elements locally owned
//                    REQUIRE(moris::equal_to(tElemsOwned(0,0),7));
//                    REQUIRE(moris::equal_to(tElemsOwned(1,0),15));
//
//                    // Check nodes locally owned
//                    REQUIRE(moris::equal_to(tNodesOwned(0,0),2));
//                    REQUIRE(moris::equal_to(tNodesOwned(1,0),4));
//                    REQUIRE(moris::equal_to(tNodesOwned(2,0),5));
//                    REQUIRE(moris::equal_to(tNodesOwned(3,0),6));
//                    REQUIRE(moris::equal_to(tNodesOwned(4,0),8));
//                    REQUIRE(moris::equal_to(tNodesOwned(5,0),10));
//                    REQUIRE(moris::equal_to(tNodesOwned(6,0),12));
//                    REQUIRE(moris::equal_to(tNodesOwned(7,0),1000));
//                    REQUIRE(moris::equal_to(tNodesOwned(8,0),3000));
//                    REQUIRE(moris::equal_to(tNodesOwned(9,0),7000));
//                    REQUIRE(moris::equal_to(tNodesOwned(10,0),9000));
//                    REQUIRE(moris::equal_to(tNodesOwned(11,0),11000));
//
//                    //Check processors sharing inquired node
//                    REQUIRE(moris::equal_to(tProcsSharedByEntityId(0,0),1));
//                    REQUIRE(moris::equal_to(tProcsSharedByEntityInd(0,0),UINT_MAX)); // Node for local index 5 not shared
//
//                    // Nodal maps
//                    REQUIRE(moris::equal_to(lcl2GblNodeMap(0),1000));
//                    REQUIRE(moris::equal_to(lcl2GblNodeMap(1),2));
//                    REQUIRE(moris::equal_to(lcl2GblNodeOwnerProc(7),0));
//                    REQUIRE(moris::equal_to(lcl2GblNodeOwnerProc(14),1));
//
//                    // Elemental maps
//                    REQUIRE(moris::equal_to(lcl2GblElemMap(0),15));
//                    REQUIRE(moris::equal_to(lcl2GblElemOwnerProc(0),0));
//                    REQUIRE(moris::equal_to(lcl2GblElemOwnerProc(2),1));
//
//                    // Edge maps
//                    REQUIRE(moris::equal_to(lcl2GblEdgeMap(19),20));
//                    REQUIRE(moris::equal_to(lcl2GblEdgeMap(24),45));
//                    REQUIRE(moris::equal_to(lcl2GblEdgeOwnerProc(10),0));
//                    REQUIRE(moris::equal_to(lcl2GblEdgeOwnerProc(25),1));
//
//                    // Face maps
//                    REQUIRE(moris::equal_to(lcl2GblFaceMap(11),19));
//                    REQUIRE(moris::equal_to(lcl2GblFaceMap(15),24));
//                    REQUIRE(moris::equal_to(lcl2GblFaceOwnerProc(3),0));
//                    REQUIRE(moris::equal_to(lcl2GblFaceOwnerProc(14),1));
//
//                    // Check shared structure (same as proc 1 since it is shared info)
//                    REQUIRE(moris::equal_to(tNodesSharedPerProc.size(),1));
//                    REQUIRE(moris::equal_to(tNodesSharedPerProc(0).size(),4));
//                    REQUIRE(moris::equal_to(tNodesSharedPerProc(0)(0),10));
//                    REQUIRE(moris::equal_to(tNodesSharedPerProc(0)(1),12));
//                    REQUIRE(moris::equal_to(tNodesSharedPerProc(0)(2),9000));
//                    REQUIRE(moris::equal_to(tNodesSharedPerProc(0)(3),11000));
//                    REQUIRE(moris::equal_to(tElemsSharedPerProc(0).size(),0)); // no aura, no shared elements
//                    REQUIRE(moris::equal_to(tEdgesSharedPerProc.size(),1));
//                    REQUIRE(moris::equal_to(tEdgesSharedPerProc(0).size(),4));
//                    REQUIRE(moris::equal_to(tEdgesSharedPerProc(0)(0),5));
//                    REQUIRE(moris::equal_to(tEdgesSharedPerProc(0)(1),6));
//                    REQUIRE(moris::equal_to(tEdgesSharedPerProc(0)(2),7));
//                    REQUIRE(moris::equal_to(tEdgesSharedPerProc(0)(3),8)); //
//                    REQUIRE(moris::equal_to(tFacesSharedPerProc.size(),1));
//                    REQUIRE(moris::equal_to(tFacesSharedPerProc(0).size(),1));
//                    REQUIRE(moris::equal_to(tFacesSharedPerProc(0)(0),6));
//
//                    // ===================================================
//                    // Testing entities connected to element with ID = 7
//                    // ===================================================
//                    uint elementID = 7;
//
//                    // Initialize and fill cells to store IDs of faces, edges and nodes connected to current element (elementID = 1)
//
//                    Matrix< DDUMat >  elemsConnectedToElement = tParallelMesh.get_elements_connected_to_element(elementID);
//
//                    // Get number of elements, faces and edges connected to current node
//                    uint NumberOfElemsConnectedToElement = elemsConnectedToElement.length();
//
//                    // Check the number of elements and its IDs connected to current element
//                    REQUIRE( moris::equal_to(NumberOfElemsConnectedToElement, 2) );
//                    REQUIRE( moris::equal_to(elemsConnectedToElement(0), 15) );
//                    REQUIRE( moris::equal_to(elemsConnectedToElement(1), 31) );
//                }
//                else
//                {
//                    // Check elements locally owned
//                    REQUIRE(moris::equal_to(tElemsOwned(0,0),31));
//
//                    // Check nodes locally owned
//                    REQUIRE(moris::equal_to(tNodesOwned(0,0),13));
//                    REQUIRE(moris::equal_to(tNodesOwned(1,0),14));
//                    REQUIRE(moris::equal_to(tNodesOwned(2,0),16));
//                    REQUIRE(moris::equal_to(tNodesOwned(3,0),15000));
//
//                    //Check processors sharing inquired node
//                    REQUIRE(moris::equal_to(tProcsSharedByEntityId(0,0),0));
//                    REQUIRE(moris::equal_to(tProcsSharedByEntityInd(0,0),0));
//
//                    // Nodal map
//                    REQUIRE(moris::equal_to(lcl2GblNodeMap(2),7000));
//                    REQUIRE(moris::equal_to(lcl2GblNodeMap(5),10));
//                    REQUIRE(moris::equal_to(lcl2GblNodeOwnerProc(7),0));
//                    REQUIRE(moris::equal_to(lcl2GblNodeOwnerProc(9),1));
//
//                    // Elemental map
//                    REQUIRE(moris::equal_to(lcl2GblElemMap(0),7));
//                    REQUIRE(moris::equal_to(lcl2GblElemOwnerProc(0),0));
//                    REQUIRE(moris::equal_to(lcl2GblElemOwnerProc(1),1));
//
//                    // Edge map
//                    REQUIRE(moris::equal_to(lcl2GblEdgeMap(12),41));
//                    REQUIRE(moris::equal_to(lcl2GblEdgeMap(17),46));
//                    REQUIRE(moris::equal_to(lcl2GblEdgeOwnerProc(10),0));
//                    REQUIRE(moris::equal_to(lcl2GblEdgeOwnerProc(15),1));
//
//                    // Face map
//                    REQUIRE(moris::equal_to(lcl2GblFaceMap(6),19));
//                    REQUIRE(moris::equal_to(lcl2GblFaceMap(2),3));
//                    REQUIRE(moris::equal_to(lcl2GblFaceOwnerProc(3),0));
//                    REQUIRE(moris::equal_to(lcl2GblFaceOwnerProc(9),1));
//
//                    // Check shared structure (same as proc 0 since it is shared info)
//                    REQUIRE(moris::equal_to(tNodesSharedPerProc.size(),1));
//                    REQUIRE(moris::equal_to(tNodesSharedPerProc(0).size(),4));
//                    REQUIRE(moris::equal_to(tNodesSharedPerProc(0)(0),10));
//                    REQUIRE(moris::equal_to(tNodesSharedPerProc(0)(1),12));
//                    REQUIRE(moris::equal_to(tNodesSharedPerProc(0)(2),9000));
//                    REQUIRE(moris::equal_to(tNodesSharedPerProc(0)(3),11000));
//                    REQUIRE(moris::equal_to(tElemsSharedPerProc(0).size(),0)); // no aura, no shared elements
//                    REQUIRE(moris::equal_to(tEdgesSharedPerProc.size(),1));
//                    REQUIRE(moris::equal_to(tEdgesSharedPerProc(0).size(),4));
//                    REQUIRE(moris::equal_to(tEdgesSharedPerProc(0)(0),5));
//                    REQUIRE(moris::equal_to(tEdgesSharedPerProc(0)(1),6));
//                    REQUIRE(moris::equal_to(tEdgesSharedPerProc(0)(2),7));
//                    REQUIRE(moris::equal_to(tEdgesSharedPerProc(0)(3),8)); //
//                    REQUIRE(moris::equal_to(tFacesSharedPerProc.size(),1));
//                    REQUIRE(moris::equal_to(tFacesSharedPerProc(0).size(),1));
//                    REQUIRE(moris::equal_to(tFacesSharedPerProc(0)(0),6));
//                }
//            }
//        }
//    }
//
//    SECTION( "Creating a 3D 3 element mesh in parallel using nodal processor owner list in parallel")
//    {
//        if( p_size == 2 ) // specify it is a 2 processor test
//        {
//            uint aNumDim = 3;
//            Mat< moris::uint >  aNodeProcs;
//            Mat< moris::uint >  aElemConn;
//            Mat< moris::real >  aNodeCoords;
//            Mat<moris::uint>    aElemLocaltoGlobal;
//            Mat<moris::uint>    aNodeLocaltoGlobal;
//
//            Cell < Matrix< DDRMat >  > aFieldData;
//            Cell < std::string > aFieldName;
//            Cell < enum EntityRank > aFieldRanks ;
//
//            if ( p_rank == 0 )
//            {
//                // Mesh information
//                // Generate data for test
//                Mat< moris::uint>  tNodeProcsDummy  = {{0}, {0}, {0}, {0},
//                        {0}, {0}, {0}, {0},
//                        {0}, {0}, {0}, {0}};
//                Mat< moris::uint > tElemConnDummy   = {{1000,2,4,3000,5,6,8,7000},
//                        {5,6,8,7000,9000,10,12,11000}};
//                Mat< moris::real > tNodeCoordsDummy = {{0.0, 0.0, 0.0},
//                        {1.0, 0.0, 0.0},
//                        {0.0, 1.0, 0.0},
//                        {1.0, 1.0, 0.0},
//                        {0.0, 0.0, 1.0},
//                        {1.0, 0.0, 1.0},
//                        {0.0, 1.0, 1.0},
//                        {1.0, 1.0, 1.0},
//                        {0.0, 0.0, 2.0},
//                        {1.0, 0.0, 2.0},
//                        {0.0, 1.0, 2.0},
//                        {1.0, 1.0, 2.0}};
//                Mat< moris::uint > tElemLocaltoGlobalDummy = {{15},{7}};
//                Mat< moris::uint > tNodeLocaltoGlobalDummy = {{1000},{2},{3000},{4},{5},{6},{7000},{8},{9000},
//                        {10},{11000},{12}};
//                aNodeProcs         = tNodeProcsDummy;
//                aElemConn          = tElemConnDummy;
//                aNodeCoords        = tNodeCoordsDummy;
//                aElemLocaltoGlobal = tElemLocaltoGlobalDummy;
//                aNodeLocaltoGlobal = tNodeLocaltoGlobalDummy;
//
//                // Field information
//                // Create nodal and elemental fields of reals
//                Matrix< DDRMat >  tNodalFieldData_1 = {{10.0},{10.0},{10.0},{10.0},
//                        {50.0},{50.0},{35.0},{35.0},
//                        {40.0},{45.0},{45.0},{40.0}};
//                Matrix< DDRMat >  aElemFieldData_1 = {{111.11, 222.22, 333.33},{0.253, -0.253, 0.976}};
//
//                aFieldData  = { tNodalFieldData_1, aElemFieldData_1};
//                aFieldName  = { "nodalField_1", "elementField_1"};
//                aFieldRanks = { EntityRank::NODE, EntityRank::ELEMENT };
//
//            } else {
//                // Mesh information
//                // Generate data for test
//                Mat< moris::uint>  tNodeProcsDummy  = {{0}, {0}, {0}, {0},
//                        {1}, {1}, {1}, {1}};
//                Mat< moris::uint > tElemConnDummy   = {{9000,10,12,11000,13,14,16,15000}};
//                Mat< moris::real > tNodeCoordsDummy = {{0.0, 0.0, 2.0},
//                        {1.0, 0.0, 2.0},
//                        {0.0, 1.0, 2.0},
//                        {1.0, 1.0, 2.0},
//                        {0.0, 0.0, 3.0},
//                        {1.0, 0.0, 3.0},
//                        {0.0, 1.0, 3.0},
//                        {1.0, 1.0, 3.0}};
//                Mat< moris::uint > tElemLocaltoGlobalDummy = {{31}};
//                Mat< moris::uint > tNodeLocaltoGlobalDummy = {{9000},{10},{11000},{12},{13},{14},{15000},{16}};
//
//                aNodeProcs         = tNodeProcsDummy;
//                aElemConn          = tElemConnDummy;
//                aNodeCoords        = tNodeCoordsDummy;
//                aElemLocaltoGlobal = tElemLocaltoGlobalDummy;
//                aNodeLocaltoGlobal = tNodeLocaltoGlobalDummy;
//
//                // Field information
//                // Create nodal and elemental fields of reals
//                Matrix< DDRMat >  tNodalFieldData_1 = {{40.0},{45.0},{45.0},{40.0},
//                        {60.0},{65.0},{65.0},{60.0}};
//                Matrix< DDRMat >  aElemFieldData_1 = {{-1.0, -3.0, -4.0}};
//
//                aFieldData  = { tNodalFieldData_1, aElemFieldData_1};
//                aFieldName  = { "nodalField_1", "elementField_1"};
//                aFieldRanks = { EntityRank::NODE, EntityRank::ELEMENT };
//            }
//
//            SECTION( "with multiple fields")
//            {
//                // Create MORIS mesh using MTK database
//                moris::MtkFieldsInfo aFieldsInfo;
//                aFieldsInfo.FieldsData = &aFieldData;
//                aFieldsInfo.FieldsName = aFieldName;
//                aFieldsInfo.FieldsRank = aFieldRanks;
//
//                moris::MtkMeshData aMeshData;
//                aMeshData.SpatialDim    = &aNumDim;
//                aMeshData.ElemConn      = &aElemConn;
//                aMeshData.NodeCoords    = &aNodeCoords;
//                aMeshData.EntProcOwner  = &aNodeProcs;
//                aMeshData.FieldsInfo     = &aFieldsInfo;
//                aMeshData.LocaltoGlobalElemMap = &aElemLocaltoGlobal;
//                aMeshData.LocaltoGlobalNodeMap = &aNodeLocaltoGlobal;
//
//                moris::mesh tParallelMesh( MeshType::MTK, aMeshData );
//
//                // ========================================
//                // Testing get_entities_owned_current_proc
//                // ========================================
//
//                Matrix< DDUMat >  tElemsOwned  = tParallelMesh.get_entities_owned_current_proc(EntityRank::ELEMENT);
//                Matrix< DDUMat >  tNodesOwned  = tParallelMesh.get_entities_owned_current_proc(EntityRank::NODE);
//                Matrix< DDRMat >  tNodalField_1 = tParallelMesh.get_field_values( EntityRank::NODE, "nodalField_1" );
//                Matrix< DDRMat >  tElemField_1  = tParallelMesh.get_field_values( EntityRank::ELEMENT, "elementField_1" );
//
//                if( p_rank == 0 )
//                {
//                    // Check some entities locally owned
//                    REQUIRE(moris::equal_to(tNodesOwned(0,0),2));
//                    REQUIRE(moris::equal_to(tNodesOwned(2,0),5));
//                    REQUIRE(moris::equal_to(tNodesOwned(5,0),10));
//                    REQUIRE(moris::equal_to(tElemsOwned(0,0),7));
//                    REQUIRE(moris::equal_to(tElemsOwned(1,0),15));
//
//                    // Check some field values
//                    REQUIRE(moris::equal_to(tNodalField_1(0,0),10.0));
//                    REQUIRE(moris::equal_to(tNodalField_1(2,0),50.0));
//                    REQUIRE(moris::equal_to(tNodalField_1(4,0),35.0));
//                    REQUIRE(moris::equal_to(tNodalField_1(10,0),40.0));
//                    REQUIRE(moris::equal_to(tElemField_1(0,0),0.253));
//                    REQUIRE(moris::equal_to(tElemField_1(1,1),222.22));
//                }
//                else
//                {
//                    // Check some entities locally owned
//                    REQUIRE(moris::equal_to(tNodesOwned(0,0),13));
//                    REQUIRE(moris::equal_to(tNodesOwned(3,0),15000));
//                    REQUIRE(moris::equal_to(tElemsOwned(0,0),31));
//
//                    // Check some field values
//                    REQUIRE(moris::equal_to(tNodalField_1(0,0),45.0));
//                    REQUIRE(moris::equal_to(tNodalField_1(2,0),60.0));
//                    REQUIRE(moris::equal_to(tNodalField_1(4,0),60.0));
//                    REQUIRE(moris::equal_to(tElemField_1(0,1),-3.0));
//                    REQUIRE(moris::equal_to(tElemField_1(0,2),-4.0));
//                }
//            }
//        }
//    }
//
//    SECTION( "Creating a 3D 3 element mesh in parallel")
//    {
//        if( p_size == 2 ) // specify it is a serial test only
//        {
//            uint aNumDim = 3;
//            Matrix< DDUMat >  aNodeProcs;
//            Matrix< DDUMat >  aElemConn;
//            Matrix< DDRMat >  aNodeCoords;
//            Matrix< DDUMat >  aElemLocaltoGlobal;
//            Matrix< DDUMat >  aNodeLocaltoGlobal;
//
//            if ( p_rank == 0 )
//            {
//                // Mesh information
//                // Generate data for test
//                Mat< moris::uint>  tNodeProcsDummy  = {{0}, {0}, {0}, {0},
//                        {0}, {0}, {0}, {0},
//                        {0}, {0}, {0}, {0}};
//                Mat< moris::uint > tElemConnDummy   = {{1000,2,4,3000,5,6,8,7000},
//                        {5,6,8,7000,9000,10,12,11000}};
//                Mat< moris::real > tNodeCoordsDummy = {{0.0, 0.0, 0.0},
//                        {1.0, 0.0, 0.0},
//                        {0.0, 1.0, 0.0},
//                        {1.0, 1.0, 0.0},
//                        {0.0, 0.0, 1.0},
//                        {1.0, 0.0, 1.0},
//                        {0.0, 1.0, 1.0},
//                        {1.0, 1.0, 1.0},
//                        {0.0, 0.0, 2.0},
//                        {1.0, 0.0, 2.0},
//                        {0.0, 1.0, 2.0},
//                        {1.0, 1.0, 2.0}};
//                Mat< moris::uint > tElemLocaltoGlobalDummy = {{15},{7}};
//                Mat< moris::uint > tNodeLocaltoGlobalDummy = {{1000},{2},{3000},{4},{5},{6},{7000},{8},{9000},
//                        {10},{11000},{12}};
//                aNodeProcs         = tNodeProcsDummy;
//                aElemConn          = tElemConnDummy;
//                aNodeCoords        = tNodeCoordsDummy;
//                aElemLocaltoGlobal = tElemLocaltoGlobalDummy;
//                aNodeLocaltoGlobal = tNodeLocaltoGlobalDummy;
//            }
//            else
//            {
//                // Mesh information
//                // Generate data for test
//                Mat< moris::uint>  tNodeProcsDummy  = {{0}, {0}, {0}, {0},
//                        {1}, {1}, {1}, {1}};
//                Mat< moris::uint > tElemConnDummy   = {{9000,10,12,11000,13,14,16,15000}};
//                Mat< moris::real > tNodeCoordsDummy = {{0.0, 0.0, 2.0},
//                        {1.0, 0.0, 2.0},
//                        {0.0, 1.0, 2.0},
//                        {1.0, 1.0, 2.0},
//                        {0.0, 0.0, 3.0},
//                        {1.0, 0.0, 3.0},
//                        {0.0, 1.0, 3.0},
//                        {1.0, 1.0, 3.0}};
//                Mat< moris::uint > tElemLocaltoGlobalDummy = {{31}};
//                Mat< moris::uint > tNodeLocaltoGlobalDummy = {{9000},{10},{11000},{12},{13},{14},{15000},{16}};
//
//                aNodeProcs         = tNodeProcsDummy;
//                aElemConn          = tElemConnDummy;
//                aNodeCoords        = tNodeCoordsDummy;
//                aElemLocaltoGlobal = tElemLocaltoGlobalDummy;
//                aNodeLocaltoGlobal = tNodeLocaltoGlobalDummy;
//            }
//
//            SECTION( "without fields or sets for parallel data check")
//            {
//                moris::MtkMeshData aMeshData;
//                aMeshData.SpatialDim    = &aNumDim;
//                aMeshData.ElemConn      = &aElemConn;
//                aMeshData.NodeCoords    = &aNodeCoords;
//                aMeshData.EntProcOwner  = &aNodeProcs;
//                aMeshData.LocaltoGlobalElemMap = &aElemLocaltoGlobal;
//                aMeshData.LocaltoGlobalNodeMap = &aNodeLocaltoGlobal;
//
//                moris::mesh tParallelMesh( MeshType::MTK, aMeshData );
//
//                // ========================================
//                // Testing get_entities_owned_current_proc
//                // ========================================
//
//                Matrix< DDUMat >  tNodeOwned  = tParallelMesh.get_entities_owned_current_proc(EntityRank::NODE);
//                Matrix< DDUMat >  tNodeShared = tParallelMesh.get_entities_glb_shared_current_proc(EntityRank::NODE);
//                Matrix< DDUMat >  tNodeAura   = tParallelMesh.get_entities_in_aura(EntityRank::NODE);
//                Matrix< DDUMat >  tNodeUniv   = tParallelMesh.get_entities_universal(EntityRank::NODE);
//                Matrix< DDUMat >  tNodeOwnedAndShared  = tParallelMesh.get_entities_owned_and_shared_by_current_proc(EntityRank::NODE);
//
//                Matrix< DDUMat >  tElemOwned = tParallelMesh.get_entities_owned_current_proc(EntityRank::ELEMENT);
//                Matrix< DDUMat >  tElemAura  = tParallelMesh.get_entities_in_aura(EntityRank::ELEMENT);
//                Matrix< DDUMat >  tElemUniv  = tParallelMesh.get_entities_universal(EntityRank::ELEMENT);
//                Matrix< DDUMat >  tElemOwnedAndShared  = tParallelMesh.get_entities_owned_and_shared_by_current_proc(EntityRank::ELEMENT);
//
//                if( p_rank == 0 )
//                {
//                    // Check some of the requested data
//                    // For this processor we don't have shared nodes
//                    REQUIRE(moris::equal_to(tNodeOwned(3),6));
//                    REQUIRE(moris::equal_to(tNodeOwnedAndShared(10),11000));
//                    REQUIRE(moris::equal_to(tNodeAura(0),13));
//                    REQUIRE(moris::equal_to(tNodeUniv(15),15000));
//
//                    // Elements don't have shared entities
//                    REQUIRE(moris::equal_to(tElemOwned(1),15));
//                    REQUIRE(moris::equal_to(tElemAura(0),31));
//                    REQUIRE(moris::equal_to(tElemUniv(2),31));
//                }
//                else
//                {
//                    // Check some of the requested data
//                    // For this processor we don't have shared nodes
//                    REQUIRE(moris::equal_to(tNodeOwned(1),14));
//                    REQUIRE(moris::equal_to(tNodeShared(2),9000));
//                    REQUIRE(moris::equal_to(tNodeOwnedAndShared(7),16));
//                    REQUIRE(moris::equal_to(tNodeUniv(7),15000));
//
//                    // Elements don't have shared entities
//                    REQUIRE(moris::equal_to(tElemOwned(0),31));
//                    REQUIRE(moris::equal_to(tElemUniv(0),31));
//                }
//            }
//
//            SECTION( "with multiple fields on the entire domain")
//            {
//                // Test for showing how to provide information of multiple fields and output them into
//                // a mesh file. In this example 1 nodal and 1 element fields of reals are provided for
//                // a three element (HEX8) mesh generated from data.
//
//                Cell < Matrix< DDRMat >  > aFieldData;
//                Cell < std::string > aFieldName;
//                Cell < enum EntityRank > aFieldRanks ;
//
//                if ( p_rank == 0 )
//                {
//                    // Field information
//                    // Create nodal and elemental fields of reals
//                    Matrix< DDRMat >  tNodalFieldData_1 = {{10.0},{10.0},{10.0},{10.0},
//                            {50.0},{50.0},{35.0},{35.0},
//                            {40.0},{45.0},{45.0},{40.0}};
//                    Matrix< DDRMat >  aElemFieldData_1 = {{111.11, 222.22, 333.33},{0.253, -0.253, 0.976}};
//
//                    aFieldData  = { tNodalFieldData_1, aElemFieldData_1};
//                    aFieldName  = { "nodalField_1", "elementField_1"};
//                    aFieldRanks = { EntityRank::NODE, EntityRank::ELEMENT };
//                }
//                else
//                {
//                    // Field information
//                    // Create nodal and elemental fields of reals
//                    Matrix< DDRMat >  tNodalFieldData_1 = {{40.0},{45.0},{45.0},{40.0},
//                            {60.0},{65.0},{65.0},{60.0}};
//                    Matrix< DDRMat >  aElemFieldData_1 = {{-1.0, -3.0, -4.0}};
//
//                    aFieldData  = { tNodalFieldData_1, aElemFieldData_1};
//                    aFieldName  = { "nodalField_1", "elementField_1"};
//                    aFieldRanks = { EntityRank::NODE, EntityRank::ELEMENT };
//                }
//
//                // Create MORIS mesh using MTK database
//                moris::MtkFieldsInfo aFieldsInfo;
//                aFieldsInfo.FieldsData = &aFieldData;
//                aFieldsInfo.FieldsName = aFieldName;
//                aFieldsInfo.FieldsRank = aFieldRanks;
//
//                moris::MtkMeshData aMeshData;
//                aMeshData.SpatialDim    = &aNumDim;
//                aMeshData.ElemConn      = &aElemConn;
//                aMeshData.NodeCoords    = &aNodeCoords;
//                aMeshData.EntProcOwner  = &aNodeProcs;
//                aMeshData.FieldsInfo     = &aFieldsInfo;
//                aMeshData.LocaltoGlobalElemMap = &aElemLocaltoGlobal;
//                aMeshData.LocaltoGlobalNodeMap = &aNodeLocaltoGlobal;
//
//                moris::mesh tParallelMesh( MeshType::MTK, aMeshData );
//
//                // ========================================
//                // Testing get_entities_owned_current_proc
//                // ========================================
//
//                Matrix< DDUMat >  tElemsOwned   = tParallelMesh.get_entities_owned_current_proc(EntityRank::ELEMENT);
//                Matrix< DDUMat >  tNodesOwned   = tParallelMesh.get_entities_owned_current_proc(EntityRank::NODE);
//                Matrix< DDRMat >  tNodalField_1 = tParallelMesh.get_field_values( EntityRank::NODE, "nodalField_1" );
//                Matrix< DDRMat >  tElemField_1  = tParallelMesh.get_field_values( EntityRank::ELEMENT, "elementField_1" );
//                Matrix< DDUMat >  tFieldEnts_1  = tParallelMesh.get_field_entities( EntityRank::NODE, "nodalField_1" );
//                Matrix< DDUMat >  tFieldEnts_2  = tParallelMesh.get_field_entities( EntityRank::ELEMENT, "elementField_1" );
//
//                if( p_rank == 0 )
//                {
//                    // Check some entities owned
//                    REQUIRE(moris::equal_to(tFieldEnts_1(0),2));
//                    REQUIRE(moris::equal_to(tFieldEnts_1(2),5));
//                    REQUIRE(moris::equal_to(tFieldEnts_1(4),8));
//                    REQUIRE(moris::equal_to(tFieldEnts_1(8),3000));
//                    REQUIRE(moris::equal_to(tFieldEnts_1(10),9000));
//                    REQUIRE(moris::equal_to(tFieldEnts_2(0),7));
//                    REQUIRE(moris::equal_to(tFieldEnts_2(1),15));
//
//                    // Check some field values
//                    REQUIRE(moris::equal_to(tNodalField_1(0,0),10.0));
//                    REQUIRE(moris::equal_to(tNodalField_1(4,0),35.0));
//                    REQUIRE(moris::equal_to(tNodalField_1(6,0),40.0));
//                    REQUIRE(moris::equal_to(tNodalField_1(8,0),10.0));
//                    REQUIRE(moris::equal_to(tElemField_1(1,0),111.11));
//                    REQUIRE(moris::equal_to(tElemField_1(0,1),-0.253));
//
//                }
//                else
//                {
//                    // Check some entities owned
//                    REQUIRE(moris::equal_to(tFieldEnts_1(0),10));
//                    REQUIRE(moris::equal_to(tFieldEnts_1(2),13));
//                    REQUIRE(moris::equal_to(tFieldEnts_1(7),15000));
//                    REQUIRE(moris::equal_to(tFieldEnts_2(0),31));
//
//                    // Check some field values
//                    REQUIRE(moris::equal_to(tNodalField_1(0,0),45.0));
//                    REQUIRE(moris::equal_to(tNodalField_1(2,0),60.0));
//                    REQUIRE(moris::equal_to(tNodalField_1(4,0),60.0));
//                    REQUIRE(moris::equal_to(tElemField_1(0,2),-4.0));
//                }
//            }
//
//            SECTION( "with sets allocated in separate processors")
//            {
//                Matrix< DDUMat >  tBlockSetsPartOwners;
//                Cell< Matrix< DDUMat >  > tNodeSetsEntIds;
//                Cell< Matrix< DDUMat >  > tSideSetsInfo;
//
//                if ( p_rank == 0 )
//                {
//                    // Populate block sets
//                    Mat<uint> tBlockSetsInds = { {1}, {0} };
//                    tBlockSetsPartOwners = tBlockSetsInds;
//
//                    // Populate side sets
//                    Matrix< DDUMat >  tSideset_1 = { {15, 5}, {7, 0} };
//                    tSideSetsInfo = { tSideset_1 };
//
//                    // Populate node sets
//                    Matrix< DDUMat >  tNodeSet_1 = { {1000}, {4}, {6} };
//                    tNodeSetsEntIds = { tNodeSet_1 };
//                }
//                else
//                {
//                    // Populate block sets
//                    Mat<uint> tBlockSetsInds = { {0} };
//                    tBlockSetsPartOwners = tBlockSetsInds;
//
//                    // Populate side sets
//                    Matrix< DDUMat >  tSideset_1 = { {31, 3} };
//                    tSideSetsInfo = { tSideset_1 };
//
//                    // Populate node sets
//                    Matrix< DDUMat >  tNodeSet_1 = { };
//                    tNodeSetsEntIds = { tNodeSet_1 };
//                }
//
//                // Declare block sets
//                /////////////////////
//                MtkBlockSetsInfo tBlockSetStruc;
//                tBlockSetStruc.BSetInds = &tBlockSetsPartOwners;
//                tBlockSetStruc.BSetNames   = { "blockset_001", "blockset_002" };
//
//                // Declare side sets
//                /////////////////////
//                MtkSideSetsInfo tSideSetStruc;
//                tSideSetStruc.ElemIdsAndSideOrds = &tSideSetsInfo;
//                tSideSetStruc.SSetNames   = { "Sideset_001" };
//
//                // Declare node sets
//                /////////////////////
//                MtkNodeSetsInfo tNodeSetStruc;
//                tNodeSetStruc.EntIds = &tNodeSetsEntIds;
//                tNodeSetStruc.NSetNames = { "Nodeset_001" };
//
//                // Create MORIS mesh using MTK database
//                ///////////////////////////////////////
//                MtkSetsInfo aMeshSets;
//                aMeshSets.NodeSetsInfo = &tNodeSetStruc;
//                aMeshSets.SideSetsInfo   = &tSideSetStruc;
//                aMeshSets.BlockSetsInfo   = &tBlockSetStruc;
//
//                moris::MtkMeshData aMeshData;
//                aMeshData.SpatialDim   = &aNumDim;
//                aMeshData.ElemConn     = &aElemConn;
//                aMeshData.NodeCoords   = &aNodeCoords;
//                aMeshData.EntProcOwner = &aNodeProcs;
//                aMeshData.SetsInfo     = &aMeshSets;
//                aMeshData.LocaltoGlobalElemMap = &aElemLocaltoGlobal;
//                aMeshData.LocaltoGlobalNodeMap = &aNodeLocaltoGlobal;
//
//                moris::mesh tParallelMesh( MeshType::MTK, aMeshData );
//
//                // ========================================
//                // Testing get_entities_owned_current_proc
//                // ========================================
//
//                Matrix< DDUMat >  tElemsOwned = tParallelMesh.get_entities_owned_current_proc(EntityRank::ELEMENT);
//                Matrix< DDUMat >  tNodesOwned = tParallelMesh.get_entities_owned_current_proc(EntityRank::NODE);
//
//                Matrix< DDUMat >  tNodeSet = tParallelMesh.get_set_entity_ids( EntityRank::NODE, "NodeSet_001" );
//                Matrix< DDUMat >  tSideSet = tParallelMesh.get_set_entity_ids( EntityRank::FACE, "Sideset_001" );
//                Matrix< DDUMat >  tBlockSet1 = tParallelMesh.get_set_entity_ids( EntityRank::ELEMENT, "Blockset_001" );
//                Matrix< DDUMat >  tBlockSet2 = tParallelMesh.get_set_entity_ids( EntityRank::ELEMENT, "Blockset_002" );
//
//                Matrix< DDUMat >  tNodesInNodeSet = tParallelMesh.get_nodes_in_node_set( 1 );
//                Matrix< DDUMat >  tNodesInSideSet = tParallelMesh.get_nodes_in_side_set( 1 );
//
//                if( p_rank == 0 )
//                {
//                    // Check some entities owned
//                    REQUIRE(moris::equal_to(tNodeSet(0),4));
//                    REQUIRE(moris::equal_to(tNodeSet(1),6));
//                    REQUIRE(moris::equal_to(tNodeSet(2),1000));
//                    REQUIRE(moris::equal_to(tSideSet(0),71));
//                    REQUIRE(moris::equal_to(tSideSet(1),156));
//                    REQUIRE(moris::equal_to(tBlockSet2(0),15));
//                    REQUIRE(moris::equal_to(tBlockSet1(0),7));
//                }
//                else
//                {
//                    // zero nodes in node set for this processor
//                    // Check entity owned
//                    REQUIRE(moris::equal_to(tSideSet(0),314));
//                    REQUIRE(moris::equal_to(tBlockSet1(0),31));
//                }
//
//                // Output mesh
//                std::string OutputFileName1 = "meshWMultipleSetsPara.exo";  //output name
//                tParallelMesh.create_output_mesh(OutputFileName1);
//            }
//        }
//    }
//
//    SECTION( "Creating a 1D mesh from data in parallel" )
//    {
//        if( p_size ==2 ) // specify it is a serial test only
//        {
//            // Generate data for test
//            uint aNumDim = 1;
//            Mat<uint >  aNodeProcs;
//            Mat<uint >  aElemConn;
//            Mat<real >  aNodeCoords;
//            Mat<uint>    aElemLocaltoGlobal;
//            Mat<uint>    aNodeLocaltoGlobal;
//
//            if ( p_rank == 0 )
//            {
//                // Generate data for test
//                Mat< uint>  tNodeProcsDummy  = {{0}, {0}, {0}, {0}};
//                Matrix< DDUMat >  tElemConnDummy   = {{1,2},{2,5},{5,6}};
//                Matrix< DDRMat >  tNodeCoordsDummy = {{0.0,0.0},{1.0,0.0},{2.0,0.0},{3.0,0.0}};
//                Matrix< DDUMat >  tElemLocaltoGlobalDummy = {{1},{3},{10}};
//                Matrix< DDUMat >  tNodeLocaltoGlobalDummy = {{1},{2},{5},{6}};
//                aNodeProcs         = tNodeProcsDummy;
//                aElemConn          = tElemConnDummy;
//                aNodeCoords        = tNodeCoordsDummy;
//                aElemLocaltoGlobal = tElemLocaltoGlobalDummy;
//                aNodeLocaltoGlobal = tNodeLocaltoGlobalDummy;
//
//            } else {
//
//                // Generate data for test
//                moris::Mat< uint>  tNodeProcsDummy  = {{0}, {1}};
//                Matrix< DDUMat >  tElemConnDummy   = {{6,3}};
//                Matrix< DDRMat >  tNodeCoordsDummy = {{3.0,0.0},{4.0,0.0}};
//                Matrix< DDUMat >  tElemLocaltoGlobalDummy = {{2}};
//                Matrix< DDUMat >  tNodeLocaltoGlobalDummy = {{6},{3}};
//
//                aNodeProcs         = tNodeProcsDummy;
//                aElemConn          = tElemConnDummy;
//                aNodeCoords        = tNodeCoordsDummy;
//                aElemLocaltoGlobal = tElemLocaltoGlobalDummy;
//                aNodeLocaltoGlobal = tNodeLocaltoGlobalDummy;
//            }
//
//            moris::MtkMeshData aMeshData;
//            aMeshData.SpatialDim   = &aNumDim;
//            aMeshData.ElemConn     = &aElemConn;
//            aMeshData.NodeCoords   = &aNodeCoords;
//            aMeshData.EntProcOwner = &aNodeProcs;
//            aMeshData.LocaltoGlobalElemMap = &aElemLocaltoGlobal;
//            aMeshData.LocaltoGlobalNodeMap = &aNodeLocaltoGlobal;
//
//            moris::mesh tMesh1DParallel( MeshType::MTK, aMeshData );
//
//            // =============================
//            // Testing basic functionalities
//            // =============================
//
//            Matrix< DDUMat >  tElemsOwned  = tMesh1DParallel.get_entities_owned_current_proc(EntityRank::ELEMENT);
//            Matrix< DDUMat >  tNodesOwned  = tMesh1DParallel.get_entities_owned_current_proc(EntityRank::NODE);
//
//            if ( p_rank == 0 )
//            {
//                REQUIRE(moris::equal_to(tNodesOwned(0),1));
//                REQUIRE(moris::equal_to(tNodesOwned(1),2));
//                REQUIRE(moris::equal_to(tNodesOwned(2),5));
//                REQUIRE(moris::equal_to(tNodesOwned(3),6));
//
//                REQUIRE(moris::equal_to(tElemsOwned(0),1));
//                REQUIRE(moris::equal_to(tElemsOwned(1),3));
//                REQUIRE(moris::equal_to(tElemsOwned(2),10));
//            }
//            else
//            {
//                REQUIRE(moris::equal_to(tNodesOwned(0),3));
//
//                REQUIRE(moris::equal_to(tElemsOwned(0),2));
//            }
//
//
//            // Output mesh
//            std::string OutputFileName1 = "mesh1DProblemParallel.exo";  //output name
//            tMesh1DParallel.create_output_mesh(OutputFileName1);
//        }
//    }
//
//
//    SECTION( "Creating a 2D mesh with 'hanging nodes' in serial")
//    {
//        if( p_size == 1 ) // serial test
//        {
//            // Verification of mesh behavior when giving hanging nodes in mesh
//
//            //   *-----*-----*  nodes: 1,2,3
//            //   |     |     |  elements: 1,2
//            //   |     |     |
//            //   *-----*-----*  nodes: 4,5,6
//            //   |           |  element: 3
//            //   |           |
//            //   *-----------*  nodes 7,8
//
//            moris::uint aNumDim = 2;
//            moris::Mat< moris::uint >  aElemConn   = { {1,2,5,4}, {2,3,6,5}, {4,6,8,7} };
//            moris::Cell< std::string > aPartNames  = { "block_1" };
//            moris::Mat< moris::real >  aNodeCoords = {{ 0.0, 2.0, 0.0 },
//                                                      { 1.0, 2.0, 0.0 },
//                                                      { 2.0, 2.0, 0.0 },
//                                                      { 0.0, 1.0, 0.0 },
//                                                      { 1.0, 1.0, 0.0 },
//                                                      { 2.0, 1.0, 0.0 },
//                                                      { 0.0, 0.0, 0.0 },
//                                                      { 2.0, 0.0, 0.0 }};
//
//            moris::MtkMeshData aMeshData;
//            aMeshData.SpatialDim             = &aNumDim;
//            aMeshData.ElemConn               = &aElemConn;
//            aMeshData.NodeCoords             = &aNodeCoords;
//            aMeshData.CreateAllEdgesAndFaces = true;
//
//            // Creating MTK mesh
//            moris::mesh tHangingNodesMesh( MeshType::MTK, aMeshData );
//
//            // Get node and element entities own by current processor.
//            Matrix< DDUMat >  nodesInElem1 = tHangingNodesMesh.get_nodes_connected_to_element( 1 );
//            Matrix< DDUMat >  nodesInElem2 = tHangingNodesMesh.get_nodes_connected_to_element( 2 );
//            Matrix< DDUMat >  nodesInElem3 = tHangingNodesMesh.get_nodes_connected_to_element( 3 );
//
//            Matrix< DDUMat >  edgesInElem1 = tHangingNodesMesh.get_edges_connected_to_element( 1 );
//            Matrix< DDUMat >  edgesInElem2 = tHangingNodesMesh.get_edges_connected_to_element( 2 );
//            Matrix< DDUMat >  edgesInElem3 = tHangingNodesMesh.get_edges_connected_to_element( 3 );
//
//            Matrix< DDUMat >  elemsConnToHangingNodes = tHangingNodesMesh.get_elements_connected_to_node( 5 );
//
//
//            Matrix< DDUMat >  elemsInEdge1 = tHangingNodesMesh.get_elements_connected_to_edge( 8 );
//            Matrix< DDUMat >  elemsInEdge2 = tHangingNodesMesh.get_elements_connected_to_edge( 9 );
//            Matrix< DDUMat >  elemsInEdge3 = tHangingNodesMesh.get_elements_connected_to_edge( 10 );
//            Matrix< DDUMat >  elemsInEdge4 = tHangingNodesMesh.get_elements_connected_to_edge( 11 );
//
//            uint numEdgesInMesh = tHangingNodesMesh.get_num_edges( );
//
//            // ========================================
//            // Testing output related to special mesh
//            // ========================================
//
//            REQUIRE(moris::equal_to( nodesInElem1( 0 ), 1 ) );
//            REQUIRE(moris::equal_to( nodesInElem1( 1 ), 2 ) );
//            REQUIRE(moris::equal_to( nodesInElem1( 2 ), 5 ) );
//            REQUIRE(moris::equal_to( nodesInElem1( 3 ), 4 ) );
//
//            REQUIRE(moris::equal_to( nodesInElem2( 0 ), 2 ) );
//            REQUIRE(moris::equal_to( nodesInElem2( 1 ), 3 ) );
//            REQUIRE(moris::equal_to( nodesInElem2( 2 ), 6 ) );
//            REQUIRE(moris::equal_to( nodesInElem2( 3 ), 5 ) );
//
//            REQUIRE(moris::equal_to( nodesInElem3( 0 ), 4 ) );
//            REQUIRE(moris::equal_to( nodesInElem3( 1 ), 6 ) );
//            REQUIRE(moris::equal_to( nodesInElem3( 2 ), 8 ) );
//            REQUIRE(moris::equal_to( nodesInElem3( 3 ), 7 ) );
//
//            REQUIRE(moris::equal_to( edgesInElem1( 0 ), 1 ) );
//            REQUIRE(moris::equal_to( edgesInElem1( 1 ), 2 ) );
//            REQUIRE(moris::equal_to( edgesInElem1( 2 ), 3 ) );
//            REQUIRE(moris::equal_to( edgesInElem1( 3 ), 4 ) );
//
//            REQUIRE(moris::equal_to( edgesInElem2( 0 ), 5 ) );
//            REQUIRE(moris::equal_to( edgesInElem2( 1 ), 6 ) );
//            REQUIRE(moris::equal_to( edgesInElem2( 2 ), 7 ) );
//            REQUIRE(moris::equal_to( edgesInElem2( 3 ), 2 ) );
//
//            REQUIRE(moris::equal_to( edgesInElem3( 0 ), 8  ) );
//            REQUIRE(moris::equal_to( edgesInElem3( 1 ), 9  ) );
//            REQUIRE(moris::equal_to( edgesInElem3( 2 ), 10 ) );
//            REQUIRE(moris::equal_to( edgesInElem3( 3 ), 11 ) );
//
//            REQUIRE(moris::equal_to( elemsInEdge1( 0 ), 3 ) );
//            REQUIRE(moris::equal_to( elemsInEdge2( 0 ), 3 ) );
//            REQUIRE(moris::equal_to( elemsInEdge3( 0 ), 3 ) );
//            REQUIRE(moris::equal_to( elemsInEdge4( 0 ), 3 ) );
//
//            REQUIRE(moris::equal_to( elemsConnToHangingNodes( 0 ), 1 ) );
//            REQUIRE(moris::equal_to( elemsConnToHangingNodes( 1 ), 2 ) );
//
//            REQUIRE(moris::equal_to( numEdgesInMesh, 11 ) );
//        }
//    }
//}

}
}
