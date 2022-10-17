/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Mesh_New.cpp
 *
 */

#include "catch.hpp"

#include "cl_MTK_Mesh.hpp"    // MTK/src
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Data_STK.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"
#include "cl_MTK_Interpolation_Mesh_STK.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Integration_Mesh_STK.hpp"

#include "cl_Communication_Tools.hpp"
#include "cl_Matrix.hpp"
#include "op_equal_equal.hpp"
#include "fn_all_true.hpp"
#include "fn_equal_to.hpp"

namespace moris
{
    namespace mtk
    {

        TEST_CASE( "Surrogate XTK Mesh", "[MTK_Surrogate]" )
        {
            uint p_size = moris::par_size();

            if ( p_size == 1 )    // specify it is a serial test only
            {
                // Define the Integration Mesh (from data from xtk)
                Matrix< DDRMat > tNodeCoordinates = { { 0, 0, 0 },
                    { 1, 0, 0 },
                    { 0, 1, 0 },
                    { 1, 1, 0 },
                    { 0, 0, 1 },
                    { 1, 0, 1 },
                    { 0, 1, 1 },
                    { 1, 1, 1 },
                    { 0, 0, 2 },
                    { 1, 0, 2 },
                    { 0, 1, 2 },
                    { 1, 1, 2 },
                    { 0, 0, 3 },
                    { 1, 0, 3 },
                    { 0, 1, 3 },
                    { 1, 1, 3 },
                    { 0, 0, 4 },
                    { 1, 0, 4 },
                    { 0, 1, 4 },
                    { 1, 1, 4 },
                    { 0.5, 0, 3.5 },
                    { 1, 0.5, 3.5 },
                    { 0.5, 1, 3.5 },
                    { 0, 0.5, 3.5 },
                    { 0.5, 0.5, 3 },
                    { 0.5, 0.5, 4 },
                    { 0.5, 0.5, 3.5 },
                    { 0.1, 0, 3.1 },
                    { 0.9, 0, 3.1 },
                    { 0.1, 0.1, 3.1 },
                    { 0.9, 0.1, 3.1 },
                    { 1, 0, 3.1 },
                    { 0, 0, 3.1 },
                    { 1, 0.1, 3.1 },
                    { 1, 0.9, 3.1 },
                    { 0.9, 0.9, 3.1 },
                    { 1, 1, 3.1 },
                    { 0.9, 1, 3.1 },
                    { 0.1, 1, 3.1 },
                    { 0.1, 0.9, 3.1 },
                    { 0, 1, 3.1 },
                    { 0, 0.9, 3.1 },
                    { 0, 0.1, 3.1 },
                    { 0.5, 0.5, 3.1 } };

                Matrix< IndexMat > tLocalToGlobalNodeMap = { { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44 } };

                Matrix< IndexMat > tInterpElemsAsIntegCellIds     = { { 1, 2, 3, 4 } };
                Matrix< IndexMat > tInterpElemsAsIntegCellToNodes = { { 1, 2, 4, 3, 5, 6, 8, 7 },
                    { 5, 6, 8, 7, 9, 10, 12, 11 },
                    { 9, 10, 12, 11, 13, 14, 16, 15 },
                    { 13, 14, 16, 15, 17, 18, 20, 19 } };

                // Tetrathedral cells in material phase 1
                Matrix< IndexMat > tCellIdsPhase0    = { { 6, 8, 10, 12, 14, 16, 17, 18, 20, 31, 32, 33, 42, 43, 44, 53, 54, 55, 62, 63, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84 } };
                Matrix< IndexMat > tCellToNodePhase0 = { { 29, 14, 31, 32 },
                    { 13, 28, 30, 33 },
                    { 35, 16, 36, 37 },
                    { 14, 31, 32, 34 },
                    { 16, 36, 37, 38 },
                    { 39, 15, 40, 41 },
                    { 13, 30, 40, 43 },
                    { 15, 40, 41, 42 },
                    { 30, 13, 33, 43 },
                    { 13, 14, 31, 29 },
                    { 28, 13, 31, 29 },
                    { 13, 28, 31, 30 },
                    { 16, 14, 35, 36 },
                    { 14, 31, 34, 36 },
                    { 14, 34, 35, 36 },
                    { 15, 16, 39, 40 },
                    { 16, 36, 38, 40 },
                    { 16, 38, 39, 40 },
                    { 15, 13, 40, 43 },
                    { 15, 40, 42, 43 },
                    { 13, 30, 31, 44 },
                    { 14, 13, 31, 44 },
                    { 25, 13, 14, 44 },
                    { 14, 31, 36, 44 },
                    { 16, 14, 36, 44 },
                    { 25, 14, 16, 44 },
                    { 16, 36, 40, 44 },
                    { 15, 16, 40, 44 },
                    { 25, 16, 15, 44 },
                    { 30, 13, 40, 44 },
                    { 13, 15, 40, 44 },
                    { 13, 25, 15, 44 } };

                // Tetrathedral cells in material phase 1
                Matrix< IndexMat > tCellToNodeGhost0 = { { 21, 27, 31, 30 },
                    { 17, 18, 21, 27 },
                    { 31, 27, 34, 36 },
                    { 18, 20, 22, 27 },
                    { 36, 27, 38, 40 },
                    { 20, 19, 23, 27 },
                    { 17, 24, 19, 27 },
                    { 30, 27, 31, 44 },
                    { 31, 27, 36, 44 },
                    { 36, 27, 40, 44 },
                    { 27, 30, 40, 44 },
                    { 17, 26, 18, 27 },
                    { 18, 26, 20, 27 },
                    { 20, 26, 19, 27 },
                    { 17, 19, 26, 27 },
                    { 21, 28, 30, 31 },
                    { 21, 28, 31, 29 },
                    { 21, 29, 31, 32 },
                    { 27, 21, 31, 32 },
                    { 18, 21, 27, 32 },
                    { 28, 21, 30, 33 },
                    { 21, 27, 30, 33 },
                    { 21, 17, 27, 33 },
                    { 27, 22, 34, 36 },
                    { 34, 22, 35, 36 },
                    { 22, 35, 36, 37 },
                    { 27, 22, 36, 37 },
                    { 20, 22, 27, 37 },
                    { 31, 27, 32, 34 },
                    { 27, 18, 32, 34 },
                    { 27, 22, 18, 34 },
                    { 27, 23, 38, 40 },
                    { 38, 23, 39, 40 },
                    { 36, 27, 37, 38 },
                    { 27, 20, 37, 38 },
                    { 27, 23, 20, 38 },
                    { 23, 39, 40, 41 },
                    { 27, 23, 40, 41 },
                    { 19, 23, 27, 41 },
                    { 27, 24, 42, 43 },
                    { 30, 27, 40, 43 },
                    { 40, 27, 42, 43 },
                    { 40, 27, 41, 42 },
                    { 27, 19, 41, 42 },
                    { 27, 24, 19, 42 },
                    { 27, 30, 33, 43 },
                    { 17, 27, 33, 43 },
                    { 24, 27, 17, 43 } };

                Matrix< IndexMat > tCellIdsGhost0 = { { 5, 7, 9, 11, 13, 15, 19, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 34, 35, 36, 37, 38, 39, 40, 41, 45, 46, 47, 48, 49, 50, 51, 52, 56, 57, 58, 59, 60, 61, 64, 65, 66, 67, 68, 69, 70, 71, 72 } };

                moris::mtk::MtkSetsInfo tMtkMeshSets;
                // Define side sets on the integration mesh (i.e. fixed bc, interface and ghost)

                // interface side set
                moris::mtk::MtkSideSetInfo tInterfaceSideSet;
                Matrix< IndexMat >         tInterfaceElemIdandSideOrd = { { 6, 2 },
                    { 8, 1 },
                    { 9, 2 },
                    { 10, 2 },
                    { 12, 1 },
                    { 13, 2 },
                    { 14, 1 },
                    { 16, 2 },
                    { 17, 1 },
                    { 18, 1 },
                    { 20, 2 },
                    { 21, 2 },
                    { 22, 2 },
                    { 23, 2 },
                    { 24, 1 },
                    { 29, 1 },
                    { 30, 1 },
                    { 32, 2 },
                    { 33, 1 },
                    { 34, 1 },
                    { 37, 2 },
                    { 41, 2 },
                    { 43, 1 },
                    { 44, 1 },
                    { 45, 1 },
                    { 48, 2 },
                    { 52, 2 },
                    { 54, 1 },
                    { 55, 1 },
                    { 56, 2 },
                    { 59, 1 },
                    { 63, 1 },
                    { 65, 2 },
                    { 66, 2 },
                    { 67, 2 },
                    { 70, 1 },
                    { 73, 1 },
                    { 76, 1 },
                    { 79, 1 },
                    { 82, 2 } };
                tInterfaceSideSet.mElemIdsAndSideOrds                 = &tInterfaceElemIdandSideOrd;
                tInterfaceSideSet.mSideSetName                        = "iside";
                tMtkMeshSets.add_side_set( &tInterfaceSideSet );

                // Fixed bc
                moris::mtk::MtkSideSetInfo tFixed;
                Matrix< IndexMat >         tFixedElementsAndOrds = { { 1, 4 } };
                tFixed.mElemIdsAndSideOrds                       = &tFixedElementsAndOrds;
                tFixed.mSideSetName                              = "fixed";
                tMtkMeshSets.add_side_set( &tFixed );

                // Fixed bc
                moris::mtk::MtkSideSetInfo tGhost;
                Matrix< IndexMat >         tGhostCellAndOrds = { { 3, 5 }, { 4, 4 } };
                tGhost.mElemIdsAndSideOrds                   = &tGhostCellAndOrds;
                tGhost.mSideSetName                          = "ghost_facets";
                tMtkMeshSets.add_side_set( &tGhost );

                // add block sets
                // Tet Cells in Omega 0
                moris::mtk::MtkBlockSetInfo tOmega0BlockSetTet;
                tOmega0BlockSetTet.mCellIdsInSet = &tCellIdsPhase0;
                tOmega0BlockSetTet.mBlockSetName = "Omega_0_tets";
                tOmega0BlockSetTet.mBlockSetTopo = CellTopology::TET4;
                tMtkMeshSets.add_block_set( &tOmega0BlockSetTet );

                // Hex Cells in Omega 0
                Matrix< IdMat >             tOmega0HexCellIds = { { 1, 2 } };
                moris::mtk::MtkBlockSetInfo tOmega0BlockSetHex;
                tOmega0BlockSetHex.mCellIdsInSet = &tOmega0HexCellIds;
                tOmega0BlockSetHex.mBlockSetName = "Omega_0_hex";
                tOmega0BlockSetHex.mBlockSetTopo = CellTopology::HEX8;
                tMtkMeshSets.add_block_set( &tOmega0BlockSetHex );

                // Cells in the ghost domain of omega 1
                moris::mtk::MtkBlockSetInfo tOmega0GhostBlockSetTet;
                tOmega0GhostBlockSetTet.mCellIdsInSet = &tCellIdsGhost0;
                tOmega0GhostBlockSetTet.mBlockSetName = "Omega_0_Ghost";
                tOmega0GhostBlockSetTet.mBlockSetTopo = CellTopology::TET4;
                tMtkMeshSets.add_block_set( &tOmega0GhostBlockSetTet );

                // Integration Cells for Ghost penalization only
                Matrix< IdMat >             tGhostCellIds = { { 3, 4 } };
                moris::mtk::MtkBlockSetInfo tCellsForGhost;
                tCellsForGhost.mCellIdsInSet = &tGhostCellIds;
                tCellsForGhost.mBlockSetName = "Ghost_Cells_0";
                tCellsForGhost.mBlockSetTopo = CellTopology::HEX8;
                tMtkMeshSets.add_block_set( &tCellsForGhost );

                // Mesh data input structure
                moris::mtk::MtkMeshData tMeshDataInput( 3 );

                moris::uint     tSpatialDim = 3;
                Matrix< IdMat > tNodeOwner( 1, tNodeCoordinates.n_rows(), moris::par_rank() );
                tMeshDataInput.ElemConn( 0 )             = &tInterpElemsAsIntegCellToNodes;
                tMeshDataInput.ElemConn( 1 )             = &tCellToNodePhase0;
                tMeshDataInput.ElemConn( 2 )             = &tCellToNodeGhost0;
                tMeshDataInput.LocaltoGlobalElemMap( 0 ) = ( &tInterpElemsAsIntegCellIds );
                tMeshDataInput.LocaltoGlobalElemMap( 1 ) = ( &tCellIdsPhase0 );
                tMeshDataInput.LocaltoGlobalElemMap( 2 ) = ( &tCellIdsGhost0 );

                tMeshDataInput.CreateAllEdgesAndFaces = false;
                tMeshDataInput.Verbose                = false;
                tMeshDataInput.SpatialDim             = &tSpatialDim;
                tMeshDataInput.NodeCoords             = &tNodeCoordinates;
                tMeshDataInput.NodeProcOwner          = &tNodeOwner;
                tMeshDataInput.LocaltoGlobalNodeMap   = &tLocalToGlobalNodeMap;

                tMeshDataInput.SetsInfo         = &tMtkMeshSets;
                tMeshDataInput.MarkNoBlockForIO = false;

                Integration_Mesh* tIntegMeshData = moris::mtk::create_integration_mesh( MeshType::STK, tMeshDataInput );

                // write to a file
                std::string tIntegFileOutput = "./integration_mesh.exo";
                tIntegMeshData->create_output_mesh( tIntegFileOutput, "tIntegFileOutput" );

                // Define the relationship between integration mesh and interpolation mesh
                Matrix< IndexMat > tVertexIDs = { { 13, 14, 16, 15, 17, 18, 20, 19, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44 } };

                // element level parameter
                Matrix< DDRMat > tLocalCoordinatesWrtElem4 = { { -1, -1, -1 },
                    { 1, -1, -1 },
                    { 1, 1, -1 },
                    { -1, 1, -1 },
                    { -1, -1, 1 },
                    { 1, -1, 1 },
                    { 1, 1, 1 },
                    { -1, 1, 1 },
                    { 0, -1, 0 },
                    { 1, 0, 0 },
                    { 0, 1, 0 },
                    { -1, 0, 0 },
                    { 0, 0, -1 },
                    { 0, 0, 1 },
                    { 0, 0, 0 },
                    { -0.8, -1, -0.8 },
                    { 0.8, -1, -0.8 },
                    { -0.8, -0.8, -0.8 },
                    { 0.8, -0.8, -0.8 },
                    { 1, -1, -0.8 },
                    { -1, -1, -0.8 },
                    { 1, -0.8, -0.8 },
                    { 1, 0.8, -0.8 },
                    { 0.8, 0.8, -0.8 },
                    { 1, 1, -0.8 },
                    { 0.8, 1, -0.8 },
                    { -0.8, 1, -0.8 },
                    { -0.8, 0.8, -0.8 },
                    { -1, 1, -0.8 },
                    { -1, 0.8, -0.8 },
                    { -1, -0.8, -0.8 },
                    { 0, 0, -0.8 } };

                delete tIntegMeshData;
            }
        }

        TEST_CASE( "Interpolation are the same Integration Mesh", "[NEW_MTK]" )
        {

            uint p_size = moris::par_size();

            if ( p_size == 1 )    // specify it is a serial test only
            {
                // Define the Interpolation Mesh
                std::string tInterpString = "generated:1x1x4";

                // construct the mesh data
                Interpolation_Mesh* tInterpMesh1 = create_interpolation_mesh( MeshType::STK, tInterpString, NULL );
                Integration_Mesh*   tIntegMesh1  = create_integration_mesh_from_interpolation_mesh( MeshType::STK, tInterpMesh1 );

                CHECK( moris::equal_to( tIntegMesh1->get_num_elems(), tInterpMesh1->get_num_elems() ) );
                CHECK( moris::equal_to( tIntegMesh1->get_num_nodes(), tInterpMesh1->get_num_nodes() ) );
                CHECK( moris::equal_to( tIntegMesh1->get_num_edges(), tInterpMesh1->get_num_edges() ) );
                CHECK( moris::equal_to( tIntegMesh1->get_num_faces(), tInterpMesh1->get_num_faces() ) );

                // Check vertex functions
                moris_id tNodeId          = 1;
                moris_id tNodeIndexInterp = tInterpMesh1->get_loc_entity_ind_from_entity_glb_id( tNodeId, EntityRank::NODE );

                mtk::Vertex const & tVertexInterp     = tInterpMesh1->get_mtk_vertex( tNodeIndexInterp );
                Matrix< DDRMat >    tNodeCoordsInterp = tInterpMesh1->get_node_coordinate( tNodeIndexInterp );
                Matrix< DDRMat >    tVertCoordsInterp = tVertexInterp.get_coords();
                CHECK( all_true( tNodeCoordsInterp == tVertCoordsInterp ) );

                // Check vertex functions
                moris_id tNodeIndexInteg = tInterpMesh1->get_loc_entity_ind_from_entity_glb_id( tNodeId, EntityRank::NODE );

                mtk::Vertex const & tVertexInteg     = tIntegMesh1->get_mtk_vertex( tNodeIndexInteg );
                Matrix< DDRMat >    tNodeCoordsInteg = tIntegMesh1->get_node_coordinate( tNodeIndexInteg );
                Matrix< DDRMat >    tVertCoordsInteg = tVertexInteg.get_coords();
                CHECK( all_true( tNodeCoordsInteg == tVertCoordsInteg ) );

                CHECK( all_true( tVertCoordsInteg == tVertCoordsInterp ) );

                // verify addresses match
                CHECK( &tVertexInteg == &tVertexInterp );

                delete tInterpMesh1;
                delete tIntegMesh1;
            }
        }

        TEST_CASE( "Integration Mesh from File", "[Integration Mesh]" )
        {    // Parallel
            uint p_size = moris::par_size();

            if ( p_size == 1 )    // specify it is a serial test only
            {
                const std::string fileName2 = "generated:8x8x8";
                // const std::string fileName2 = "generated:8x8x8|sideset:xXyYzZ";    // 512 elements, 729 nodes, 1944 edges, 1728 faces

                // create an integration mesh
                Integration_Mesh* tMesh3D_HEXs = create_integration_mesh( MeshType::STK, fileName2, NULL, true );

                uint NumElements2 = tMesh3D_HEXs->get_num_elems();
                uint NumNodes2    = tMesh3D_HEXs->get_num_nodes();
                uint NumEdges2    = tMesh3D_HEXs->get_num_edges();
                uint NumFaces2    = tMesh3D_HEXs->get_num_faces();

                CHECK( moris::equal_to( NumElements2, 512 ) );
                CHECK( moris::equal_to( NumNodes2, 729 ) );
                CHECK( moris::equal_to( NumEdges2, 1944 ) );
                CHECK( moris::equal_to( NumFaces2, 1728 ) );

                // ================================================
                // Testing entities connected to node with ID = 1
                // ================================================

                // Initialize and fill cells to store IDs of elements, faces and edges connected to current node (nodeID = 1)
                moris_id        nodeID                  = 1;
                moris_id        nodeIndex               = tMesh3D_HEXs->get_loc_entity_ind_from_entity_glb_id( nodeID, EntityRank::NODE );
                Matrix< IdMat > elementsConnectedToNode = tMesh3D_HEXs->get_elements_connected_to_node_glob_ids( nodeID );
                Matrix< IdMat > facesConnectedToNode    = tMesh3D_HEXs->get_faces_connected_to_node_glob_ids( nodeID );
                Matrix< IdMat > edgesConnectedToNode    = tMesh3D_HEXs->get_edges_connected_to_node_glob_ids( nodeID );

                Matrix< IndexMat > tElementsConnectedToNodeInds = tMesh3D_HEXs->get_elements_connected_to_node_loc_inds( nodeIndex );
                Matrix< IndexMat > tFacesConnectedToNodeInds    = tMesh3D_HEXs->get_faces_connected_to_node_loc_inds( nodeIndex );
                Matrix< IndexMat > tEdgesConnectedToNodeInds    = tMesh3D_HEXs->get_edges_connected_to_node_loc_inds( nodeIndex );

                // make sure the ids match when the indices are converted to ids
                bool tElementIdsMatch = all_true( elementsConnectedToNode == convert_entity_indices_to_ids( tElementsConnectedToNodeInds, EntityRank::ELEMENT, tMesh3D_HEXs ) );
                bool tFaceIdsMatch    = all_true( facesConnectedToNode == convert_entity_indices_to_ids( tFacesConnectedToNodeInds, EntityRank::FACE, tMesh3D_HEXs ) );
                bool tEdgeIdsMatch    = all_true( edgesConnectedToNode == convert_entity_indices_to_ids( tEdgesConnectedToNodeInds, EntityRank::EDGE, tMesh3D_HEXs ) );
                CHECK( tElementIdsMatch );
                CHECK( tFaceIdsMatch );
                CHECK( tEdgeIdsMatch );

                // Get number of elements, faces and edges connected to current node
                uint NumberOfElemsConnectedToNode = elementsConnectedToNode.numel();
                uint NumberOfFacesConnectedToNode = facesConnectedToNode.numel();
                uint NumberOfEdgesConnectedToNode = edgesConnectedToNode.numel();

                // Check the number of elements and its IDs connected to current node
                CHECK( moris::equal_to( NumberOfElemsConnectedToNode, 1 ) );
                CHECK( moris::equal_to( elementsConnectedToNode( 0 ), 1 ) );

                // Check the number of faces and its IDs connected to current node
                CHECK( moris::equal_to( NumberOfFacesConnectedToNode, 3 ) );
                //        CHECK( moris::equal_to(facesConnectedToNode(0), 1) );
                //        CHECK( moris::equal_to(facesConnectedToNode(1), 4) );
                //        CHECK( moris::equal_to(facesConnectedToNode(2), 5) );

                // Check the number of edges and its IDs connected to current node
                CHECK( moris::equal_to( NumberOfEdgesConnectedToNode, 3 ) );
                CHECK( moris::equal_to( edgesConnectedToNode( 0 ), 1 ) );
                CHECK( moris::equal_to( edgesConnectedToNode( 1 ), 4 ) );
                CHECK( moris::equal_to( edgesConnectedToNode( 2 ), 9 ) );

                // ================================================
                // Testing entities connected to edge with ID = 25
                // ================================================
                moris_id edgeID  = 25;
                moris_id edgeIND = tMesh3D_HEXs->get_loc_entity_ind_from_entity_glb_id( edgeID, EntityRank::EDGE );

                // Initialize and fill cells to store IDs of elements and faces connected to current edge (edgeID = 25)
                Matrix< IdMat > elementsConnectedToEdge = tMesh3D_HEXs->get_elements_connected_to_edge_glob_ids( edgeID );
                Matrix< IdMat > facesConnectedToEdge    = tMesh3D_HEXs->get_faces_connected_to_edge_glob_ids( edgeID );

                // Initialize and fill cells to store IDs of elements and faces connected to current edge (edgeIND = 24)
                Matrix< IndexMat > tElementsConnectedToEdgeInd = tMesh3D_HEXs->get_elements_connected_to_edge_loc_inds( edgeIND );
                Matrix< IndexMat > tFacesConnectedToEdgeInd    = tMesh3D_HEXs->get_faces_connected_to_edge_loc_inds( edgeIND );

                tElementIdsMatch = all_true( elementsConnectedToEdge == convert_entity_indices_to_ids( tElementsConnectedToEdgeInd, EntityRank::ELEMENT, tMesh3D_HEXs ) );
                tFaceIdsMatch    = all_true( facesConnectedToEdge == convert_entity_indices_to_ids( tFacesConnectedToEdgeInd, EntityRank::FACE, tMesh3D_HEXs ) );
                CHECK( tElementIdsMatch );
                CHECK( tFaceIdsMatch );

                // Get number of elements, faces and edges connected to current edge
                uint NumberOfElemsConnectedToEdge = elementsConnectedToEdge.numel();
                uint NumberOfFacesConnectedToEdge = facesConnectedToEdge.numel();

                // Check the number of elements and its IDs connected to current edge
                CHECK( moris::equal_to( NumberOfElemsConnectedToEdge, 4 ) );
                CHECK( moris::equal_to( elementsConnectedToEdge( 0 ), 67 ) );
                CHECK( moris::equal_to( elementsConnectedToEdge( 1 ), 68 ) );
                CHECK( moris::equal_to( elementsConnectedToEdge( 2 ), 3 ) );
                CHECK( moris::equal_to( elementsConnectedToEdge( 3 ), 4 ) );

                // Check the number of faces and its IDs connected to current edge
                CHECK( moris::equal_to( NumberOfFacesConnectedToEdge, 4 ) );
                //        CHECK( moris::equal_to(facesConnectedToEdge(0), 283) );
                //        CHECK( moris::equal_to(facesConnectedToEdge(1), 16) );
                //        CHECK( moris::equal_to(facesConnectedToEdge(2), 13) );
                //        CHECK( moris::equal_to(facesConnectedToEdge(3), 21) );

                // ================================================
                // Testing entities connected to face with ID = 25
                // ================================================

                // Get number of elements, faces and edges connected to current face
                moris_id    faceID    = 50;
                moris_index faceIndex = tMesh3D_HEXs->get_loc_entity_ind_from_entity_glb_id( faceID, EntityRank::FACE );

                Matrix< IdMat >    elementsConnectedToFace     = tMesh3D_HEXs->get_elements_connected_to_face_glob_ids( faceID );
                Matrix< IndexMat > tElementsConnectedToFaceInd = tMesh3D_HEXs->get_elements_connected_to_face_loc_inds( faceIndex );

                // verify indices and ids are consistent
                tElementIdsMatch = all_true( elementsConnectedToFace == convert_entity_indices_to_ids( tElementsConnectedToFaceInd.get_row( 0 ), EntityRank::ELEMENT, tMesh3D_HEXs ) );
                CHECK( tElementIdsMatch );

                uint NumberOfElemsConnectedToFace = elementsConnectedToFace.numel();

                // Check the number of elements and its IDs connected to current face
                CHECK( moris::equal_to( NumberOfElemsConnectedToFace, 2 ) );
                //        CHECK( moris::equal_to(elementsConnectedToFace(0), 74) );
                //        CHECK( moris::equal_to(elementsConnectedToFace(1), 10) );

                // ===================================================
                // Testing entities connected to element with ID = 100
                // ===================================================
                moris_id    elementID  = 100;
                moris_index elementInd = tMesh3D_HEXs->get_loc_entity_ind_from_entity_glb_id( elementID, EntityRank::ELEMENT );

                // Initialize and fill cells to store IDs of faces, edges and nodes connected to current element (elementID = 1)
                Matrix< IdMat >    elemsConnectedToElement     = tMesh3D_HEXs->get_element_connected_to_element_glob_ids( elementID );
                Matrix< IndexMat > tElemsConnectedToElementInd = tMesh3D_HEXs->get_elements_connected_to_element_and_face_ord_loc_inds( elementInd );

                // Check consistency of element ids
                tElementIdsMatch = all_true( elemsConnectedToElement.get_row( 0 ) == convert_entity_indices_to_ids( tElemsConnectedToElementInd.get_row( 0 ), EntityRank::ELEMENT, tMesh3D_HEXs ).get_row( 0 ) );
                CHECK( tElementIdsMatch );

                // Check face ordinals
                bool tOrdsMatch = all_true( elemsConnectedToElement.get_row( 1 ) == tElemsConnectedToElementInd.get_row( 1 ) );
                CHECK( tOrdsMatch );

                Matrix< IdMat > facesConnectedToElement = tMesh3D_HEXs->get_faces_connected_to_element_glob_ids( elementID );
                Matrix< IdMat > edgesConnectedToElement = tMesh3D_HEXs->get_edges_connected_to_element_glob_ids( elementID );
                Matrix< IdMat > nodesConnectedToElement = tMesh3D_HEXs->get_nodes_connected_to_element_glob_ids( elementID );

                Matrix< IndexMat > tFacesConnectedToElementInds = tMesh3D_HEXs->get_faces_connected_to_element_loc_inds( elementInd );
                Matrix< IndexMat > tEdgesConnectedToElementInds = tMesh3D_HEXs->get_edges_connected_to_element_loc_inds( elementInd );
                Matrix< IndexMat > tNodesConnectedToElementInds = tMesh3D_HEXs->get_nodes_connected_to_element_loc_inds( elementInd );

                tFaceIdsMatch      = all_true( facesConnectedToElement == convert_entity_indices_to_ids( tFacesConnectedToElementInds, EntityRank::FACE, tMesh3D_HEXs ) );
                tEdgeIdsMatch      = all_true( edgesConnectedToElement == convert_entity_indices_to_ids( tEdgesConnectedToElementInds, EntityRank::EDGE, tMesh3D_HEXs ) );
                bool tNodeIdsMatch = all_true( nodesConnectedToElement == convert_entity_indices_to_ids( tNodesConnectedToElementInds, EntityRank::NODE, tMesh3D_HEXs ) );

                CHECK( tFaceIdsMatch );
                CHECK( tEdgeIdsMatch );
                CHECK( tNodeIdsMatch );

                // Get number of elements, faces and edges connected to current node
                uint NumberOfElemsConnectedToElement = elemsConnectedToElement.n_cols();
                uint NumberOfFacesConnectedToElement = facesConnectedToElement.numel();
                uint NumberOfEdgesConnectedToElement = edgesConnectedToElement.numel();
                uint NumberOfNodesConnectedToElement = nodesConnectedToElement.numel();

                // Check the number of elements and its IDs connected to current element
                CHECK( moris::equal_to( NumberOfElemsConnectedToElement, 6 ) );
                CHECK( moris::equal_to( elemsConnectedToElement( 0, 0 ), 92 ) );
                CHECK( moris::equal_to( elemsConnectedToElement( 0, 1 ), 101 ) );
                CHECK( moris::equal_to( elemsConnectedToElement( 0, 2 ), 108 ) );
                CHECK( moris::equal_to( elemsConnectedToElement( 0, 3 ), 99 ) );
                CHECK( moris::equal_to( elemsConnectedToElement( 0, 4 ), 36 ) );
                CHECK( moris::equal_to( elemsConnectedToElement( 0, 5 ), 164 ) );

                // Check the number of faces and its IDs connected to current element
                CHECK( moris::equal_to( NumberOfFacesConnectedToElement, 6 ) );
                //        CHECK( moris::equal_to(facesConnectedToElement(0), 367) );
                //        CHECK( moris::equal_to(facesConnectedToElement(1), 391) );
                //        CHECK( moris::equal_to(facesConnectedToElement(2), 392) );
                //        CHECK( moris::equal_to(facesConnectedToElement(3), 388) );
                //        CHECK( moris::equal_to(facesConnectedToElement(4), 157) );
                //        CHECK( moris::equal_to(facesConnectedToElement(5), 393) );

                // Check the number of edges and its IDs connected to current element
                CHECK( moris::equal_to( NumberOfEdgesConnectedToElement, 12 ) );
                CHECK( moris::equal_to( edgesConnectedToElement( 0 ), 176 ) );
                CHECK( moris::equal_to( edgesConnectedToElement( 1 ), 218 ) );
                CHECK( moris::equal_to( edgesConnectedToElement( 2 ), 219 ) );
                CHECK( moris::equal_to( edgesConnectedToElement( 3 ), 213 ) );
                CHECK( moris::equal_to( edgesConnectedToElement( 4 ), 477 ) );
                CHECK( moris::equal_to( edgesConnectedToElement( 5 ), 502 ) );
                CHECK( moris::equal_to( edgesConnectedToElement( 6 ), 503 ) );
                CHECK( moris::equal_to( edgesConnectedToElement( 7 ), 499 ) );
                CHECK( moris::equal_to( edgesConnectedToElement( 8 ), 475 ) );
                CHECK( moris::equal_to( edgesConnectedToElement( 9 ), 478 ) );
                CHECK( moris::equal_to( edgesConnectedToElement( 10 ), 504 ) );
                CHECK( moris::equal_to( edgesConnectedToElement( 11 ), 501 ) );

                // Check the number of nodes and its IDs connected to current element
                CHECK( moris::equal_to( NumberOfNodesConnectedToElement, 8 ) );
                CHECK( moris::equal_to( nodesConnectedToElement( 0 ), 121 ) );
                CHECK( moris::equal_to( nodesConnectedToElement( 1 ), 122 ) );
                CHECK( moris::equal_to( nodesConnectedToElement( 2 ), 131 ) );
                CHECK( moris::equal_to( nodesConnectedToElement( 3 ), 130 ) );
                CHECK( moris::equal_to( nodesConnectedToElement( 4 ), 202 ) );
                CHECK( moris::equal_to( nodesConnectedToElement( 5 ), 203 ) );
                CHECK( moris::equal_to( nodesConnectedToElement( 6 ), 212 ) );
                CHECK( moris::equal_to( nodesConnectedToElement( 7 ), 211 ) );

                // ===================================================
                // Testing generate_unique_node_ids for 5 nodes
                // ===================================================

                // Ask for 5 unique node IDS
                Matrix< IdMat > tAvailableNodeIDs = tMesh3D_HEXs->generate_unique_node_ids( 5 );

                // Make sure it returns the correct IDs
                CHECK( moris::equal_to( tAvailableNodeIDs( 0, 0 ), 730 ) );
                CHECK( moris::equal_to( tAvailableNodeIDs( 1, 0 ), 731 ) );
                CHECK( moris::equal_to( tAvailableNodeIDs( 2, 0 ), 732 ) );
                CHECK( moris::equal_to( tAvailableNodeIDs( 3, 0 ), 733 ) );
                CHECK( moris::equal_to( tAvailableNodeIDs( 4, 0 ), 734 ) );

                // ===================================================
                // Test Cell and Vertex API
                // ===================================================
                mtk::Cell const & tCell = tMesh3D_HEXs->get_mtk_cell( elementInd );

                Matrix< IdMat > tCellVertexIds  = tCell.get_vertex_ids();
                Matrix< IdMat > tCellVertexInds = tCell.get_vertex_inds();

                CHECK( all_true( tCellVertexIds == nodesConnectedToElement ) );
                CHECK( all_true( tCellVertexInds == tNodesConnectedToElementInds ) );

                // Check vertex functions
                mtk::Vertex const & tVertex     = tMesh3D_HEXs->get_mtk_vertex( nodeIndex );
                Matrix< DDRMat >    tNodeCoords = tMesh3D_HEXs->get_node_coordinate( nodeIndex );
                Matrix< DDRMat >    tVertCoords = tVertex.get_coords();
                CHECK( all_true( tNodeCoords == tVertCoords ) );

                CHECK( equal_to( tVertex.get_id(), nodeID ) );
                CHECK( equal_to( tVertex.get_index(), nodeIndex ) );

                // ===================================================
                // Dump to file
                // ===================================================
                std::string tFileOutput = "./mtk_generated_ut.exo";
                tMesh3D_HEXs->create_output_mesh( tFileOutput );

                delete tMesh3D_HEXs;
            }
        }

    }    // namespace mtk
}    // namespace moris
