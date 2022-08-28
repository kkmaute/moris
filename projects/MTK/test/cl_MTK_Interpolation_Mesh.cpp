/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Interpolation_Mesh.cpp
 *
 */

#include "catch.hpp"

#include "paths.hpp"

#include "cl_MTK_Mesh.hpp" // MTK/src
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Data_STK.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"
#include "cl_MTK_Interpolation_Mesh_STK.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Integration_Mesh_STK.hpp"
#include "cl_MTK_Mesh_Tools.hpp"

#include "cl_Communication_Tools.hpp"
#include "cl_Matrix.hpp"
#include "op_equal_equal.hpp"
#include "fn_all_true.hpp"
#include "fn_equal_to.hpp"

namespace moris
{
namespace mtk
{
TEST_CASE( "Interpolation Mesh from File STK","[Interpolation Mesh]")
{    // Parallel
    uint p_size = moris::par_size();

    // File prefix
    std::string tMORISROOT = moris::get_base_moris_dir();
    if( p_size == 1 ) // specify it is a serial test only
    {
        const std::string fileName2 = "generated:8x8x8";
        //const std::string fileName2 = "generated:8x8x8|sideset:xXyYzZ";    // 512 elements, 729 nodes, 1944 edges, 1728 faces

        // create an interpolation mesh
        Interpolation_Mesh* tMesh3D_HEXs = create_interpolation_mesh( MeshType::STK, fileName2, NULL );

        uint NumElements2      = tMesh3D_HEXs->get_num_elems();
        uint NumNodes2         = tMesh3D_HEXs->get_num_nodes();
        uint NumEdges2         = tMesh3D_HEXs->get_num_edges( );
        uint NumFaces2         = tMesh3D_HEXs->get_num_faces( );

        CHECK( moris::equal_to(NumElements2, 512) );
        CHECK( moris::equal_to(NumNodes2, 729) );
        CHECK( moris::equal_to(NumEdges2, 1944) );
        CHECK( moris::equal_to(NumFaces2, 1728) );

        // ================================================
        // Testing entities connected to node with ID = 1
        // ================================================

        // Initialize and fill cells to store IDs of elements, faces and edges connected to current node (nodeID = 1)
        moris_id nodeID = 1;
        moris_id nodeIndex = tMesh3D_HEXs->get_loc_entity_ind_from_entity_glb_id(nodeID,EntityRank::NODE);
        Matrix< IdMat > elementsConnectedToNode = tMesh3D_HEXs->get_elements_connected_to_node_glob_ids(nodeID);
        Matrix< IdMat > facesConnectedToNode    = tMesh3D_HEXs->get_faces_connected_to_node_glob_ids(nodeID);
        Matrix< IdMat > edgesConnectedToNode    = tMesh3D_HEXs->get_edges_connected_to_node_glob_ids(nodeID);

        Matrix< IndexMat > tElementsConnectedToNodeInds = tMesh3D_HEXs->get_elements_connected_to_node_loc_inds(nodeIndex);
        Matrix< IndexMat > tFacesConnectedToNodeInds    = tMesh3D_HEXs->get_faces_connected_to_node_loc_inds(nodeIndex);
        Matrix< IndexMat > tEdgesConnectedToNodeInds    = tMesh3D_HEXs->get_edges_connected_to_node_loc_inds(nodeIndex);

        // make sure the ids match when the indices are converted to ids
        bool tElementIdsMatch = all_true(elementsConnectedToNode == convert_entity_indices_to_ids(tElementsConnectedToNodeInds,EntityRank::ELEMENT,tMesh3D_HEXs));
        bool tFaceIdsMatch    = all_true(facesConnectedToNode == convert_entity_indices_to_ids(tFacesConnectedToNodeInds,EntityRank::FACE,tMesh3D_HEXs));
        bool tEdgeIdsMatch    = all_true(edgesConnectedToNode == convert_entity_indices_to_ids(tEdgesConnectedToNodeInds,EntityRank::EDGE,tMesh3D_HEXs));
        CHECK(tElementIdsMatch);
        CHECK(tFaceIdsMatch);
        CHECK(tEdgeIdsMatch);

        // Get number of elements, faces and edges connected to current node
        uint NumberOfElemsConnectedToNode = elementsConnectedToNode.numel();
        uint NumberOfFacesConnectedToNode = facesConnectedToNode.numel();
        uint NumberOfEdgesConnectedToNode = edgesConnectedToNode.numel();

        // Check the number of elements and its IDs connected to current node
        CHECK( moris::equal_to(NumberOfElemsConnectedToNode, 1) );
        CHECK( moris::equal_to(elementsConnectedToNode(0), 1) );

        // Check the number of faces and its IDs connected to current node
        CHECK( moris::equal_to(NumberOfFacesConnectedToNode, 3) );
//        CHECK( moris::equal_to(facesConnectedToNode(0), 1) );
//        CHECK( moris::equal_to(facesConnectedToNode(1), 4) );
//        CHECK( moris::equal_to(facesConnectedToNode(2), 5) );

        // Check the number of edges and its IDs connected to current node
        CHECK( moris::equal_to(NumberOfEdgesConnectedToNode, 3) );
        CHECK( moris::equal_to(edgesConnectedToNode(0), 1) );
        CHECK( moris::equal_to(edgesConnectedToNode(1), 4) );
        CHECK( moris::equal_to(edgesConnectedToNode(2), 9) );

        // ================================================
        // Testing entities connected to edge with ID = 25
        // ================================================
        moris_id edgeID = 25;
        moris_id edgeIND = tMesh3D_HEXs->get_loc_entity_ind_from_entity_glb_id(edgeID,EntityRank::EDGE);

        // Initialize and fill cells to store IDs of elements and faces connected to current edge (edgeID = 25)
        Matrix< IdMat > elementsConnectedToEdge = tMesh3D_HEXs->get_elements_connected_to_edge_glob_ids(edgeID);
        Matrix< IdMat > facesConnectedToEdge    = tMesh3D_HEXs->get_faces_connected_to_edge_glob_ids(edgeID);

        // Initialize and fill cells to store IDs of elements and faces connected to current edge (edgeIND = 24)
        Matrix< IndexMat > tElementsConnectedToEdgeInd = tMesh3D_HEXs->get_elements_connected_to_edge_loc_inds(edgeIND);
        Matrix< IndexMat > tFacesConnectedToEdgeInd    = tMesh3D_HEXs->get_faces_connected_to_edge_loc_inds(edgeIND);

        tElementIdsMatch = all_true(elementsConnectedToEdge == convert_entity_indices_to_ids(tElementsConnectedToEdgeInd,EntityRank::ELEMENT,tMesh3D_HEXs));
        tFaceIdsMatch    = all_true(facesConnectedToEdge == convert_entity_indices_to_ids(tFacesConnectedToEdgeInd,EntityRank::FACE,tMesh3D_HEXs));
        CHECK(tElementIdsMatch);
        CHECK(tFaceIdsMatch);

        // Get number of elements, faces and edges connected to current edge
        uint NumberOfElemsConnectedToEdge = elementsConnectedToEdge.numel();
        uint NumberOfFacesConnectedToEdge = facesConnectedToEdge.numel();

        // Check the number of elements and its IDs connected to current edge
        CHECK( moris::equal_to(NumberOfElemsConnectedToEdge, 4) );
        CHECK( moris::equal_to(elementsConnectedToEdge(0), 67) );
        CHECK( moris::equal_to(elementsConnectedToEdge(1), 68) );
        CHECK( moris::equal_to(elementsConnectedToEdge(2), 3) );
        CHECK( moris::equal_to(elementsConnectedToEdge(3), 4) );

        // Check the number of faces and its IDs connected to current edge
        CHECK( moris::equal_to(NumberOfFacesConnectedToEdge, 4) );
//        CHECK( moris::equal_to(facesConnectedToEdge(0), 283) );
//        CHECK( moris::equal_to(facesConnectedToEdge(1), 16) );
//        CHECK( moris::equal_to(facesConnectedToEdge(2), 13) );
//        CHECK( moris::equal_to(facesConnectedToEdge(3), 21) );

        // ================================================
        // Testing entities connected to face with ID = 25
        // ================================================

        // Get number of elements, faces and edges connected to current face
        moris_id    faceID    = 50;
        moris_index faceIndex = tMesh3D_HEXs->get_loc_entity_ind_from_entity_glb_id(faceID,EntityRank::FACE);

        Matrix< IdMat > elementsConnectedToFace        = tMesh3D_HEXs->get_elements_connected_to_face_glob_ids(faceID);
        Matrix< IndexMat > tElementsConnectedToFaceInd = tMesh3D_HEXs->get_elements_connected_to_face_loc_inds(faceIndex);

        // verify indices and ids are consistent
        tElementIdsMatch = all_true(elementsConnectedToFace == convert_entity_indices_to_ids(tElementsConnectedToFaceInd,EntityRank::ELEMENT,tMesh3D_HEXs));
        CHECK(tElementIdsMatch);

        uint NumberOfElemsConnectedToFace = elementsConnectedToFace.numel();

        // Check the number of elements and its IDs connected to current face
        CHECK( moris::equal_to(NumberOfElemsConnectedToFace, 2) );
//        CHECK( moris::equal_to(elementsConnectedToFace(0), 74) );
//        CHECK( moris::equal_to(elementsConnectedToFace(1), 10) );

        // ===================================================
        // Testing entities connected to element with ID = 100
        // ===================================================
        moris_id elementID = 100;
        moris_index elementInd = tMesh3D_HEXs->get_loc_entity_ind_from_entity_glb_id(elementID,EntityRank::ELEMENT);

        // Initialize and fill cells to store IDs of faces, edges and nodes connected to current element (elementID = 1)
        Matrix< IdMat > elemsConnectedToElement = tMesh3D_HEXs->get_element_connected_to_element_glob_ids(elementID);
        Matrix< IndexMat > tElemsConnectedToElementInd = tMesh3D_HEXs->get_elements_connected_to_element_and_face_ord_loc_inds(elementInd);

        // Check consistency of element ids
        tElementIdsMatch = all_true(elemsConnectedToElement.get_row(0) == convert_entity_indices_to_ids( tElemsConnectedToElementInd.get_row(0), EntityRank::ELEMENT, tMesh3D_HEXs).get_row(0));
        CHECK(tElementIdsMatch);

        // Check face ordinals
        bool tOrdsMatch = all_true(elemsConnectedToElement.get_row(1) == tElemsConnectedToElementInd.get_row(1));
        CHECK(tOrdsMatch);

        Matrix< IdMat > facesConnectedToElement = tMesh3D_HEXs->get_faces_connected_to_element_glob_ids(elementID);
        Matrix< IdMat > edgesConnectedToElement = tMesh3D_HEXs->get_edges_connected_to_element_glob_ids(elementID);
        Matrix< IdMat > nodesConnectedToElement = tMesh3D_HEXs->get_nodes_connected_to_element_glob_ids(elementID);

        Matrix< IndexMat > tFacesConnectedToElementInds = tMesh3D_HEXs->get_faces_connected_to_element_loc_inds(elementInd);
        Matrix< IndexMat > tEdgesConnectedToElementInds = tMesh3D_HEXs->get_edges_connected_to_element_loc_inds(elementInd);
        Matrix< IndexMat > tNodesConnectedToElementInds = tMesh3D_HEXs->get_nodes_connected_to_element_loc_inds(elementInd);

        tFaceIdsMatch = all_true(facesConnectedToElement == convert_entity_indices_to_ids(tFacesConnectedToElementInds,EntityRank::FACE,tMesh3D_HEXs));
        tEdgeIdsMatch = all_true(edgesConnectedToElement == convert_entity_indices_to_ids(tEdgesConnectedToElementInds,EntityRank::EDGE,tMesh3D_HEXs));
        bool tNodeIdsMatch = all_true(nodesConnectedToElement == convert_entity_indices_to_ids(tNodesConnectedToElementInds,EntityRank::NODE,tMesh3D_HEXs));

        CHECK(tFaceIdsMatch);
        CHECK(tEdgeIdsMatch);
        CHECK(tNodeIdsMatch);

        // Get number of elements, faces and edges connected to current node
        uint NumberOfElemsConnectedToElement = elemsConnectedToElement.n_cols();
        uint NumberOfFacesConnectedToElement = facesConnectedToElement.numel();
        uint NumberOfEdgesConnectedToElement = edgesConnectedToElement.numel();
        uint NumberOfNodesConnectedToElement = nodesConnectedToElement.numel();

        // Check the number of elements and its IDs connected to current element
        CHECK( moris::equal_to(NumberOfElemsConnectedToElement, 6) );
        CHECK( moris::equal_to(elemsConnectedToElement(0,0), 92)  );
        CHECK( moris::equal_to(elemsConnectedToElement(0,1), 101) );
        CHECK( moris::equal_to(elemsConnectedToElement(0,2), 108) );
        CHECK( moris::equal_to(elemsConnectedToElement(0,3), 99)  );
        CHECK( moris::equal_to(elemsConnectedToElement(0,4), 36)  );
        CHECK( moris::equal_to(elemsConnectedToElement(0,5), 164) );

        // Check the number of faces and its IDs connected to current element
        CHECK( moris::equal_to(NumberOfFacesConnectedToElement, 6) );
//        CHECK( moris::equal_to(facesConnectedToElement(0), 367) );
//        CHECK( moris::equal_to(facesConnectedToElement(1), 391) );
//        CHECK( moris::equal_to(facesConnectedToElement(2), 392) );
//        CHECK( moris::equal_to(facesConnectedToElement(3), 388) );
//        CHECK( moris::equal_to(facesConnectedToElement(4), 157) );
//        CHECK( moris::equal_to(facesConnectedToElement(5), 393) );

        // Check the number of edges and its IDs connected to current element
        CHECK( moris::equal_to(NumberOfEdgesConnectedToElement, 12) );
        CHECK( moris::equal_to(edgesConnectedToElement(0) , 176) );
        CHECK( moris::equal_to(edgesConnectedToElement(1) , 218) );
        CHECK( moris::equal_to(edgesConnectedToElement(2) , 219) );
        CHECK( moris::equal_to(edgesConnectedToElement(3) , 213) );
        CHECK( moris::equal_to(edgesConnectedToElement(4) , 477) );
        CHECK( moris::equal_to(edgesConnectedToElement(5) , 502) );
        CHECK( moris::equal_to(edgesConnectedToElement(6) , 503) );
        CHECK( moris::equal_to(edgesConnectedToElement(7) , 499) );
        CHECK( moris::equal_to(edgesConnectedToElement(8) , 475) );
        CHECK( moris::equal_to(edgesConnectedToElement(9) , 478) );
        CHECK( moris::equal_to(edgesConnectedToElement(10), 504) );
        CHECK( moris::equal_to(edgesConnectedToElement(11), 501) );

        // Check the number of nodes and its IDs connected to current element
        CHECK( moris::equal_to(NumberOfNodesConnectedToElement, 8) );
        CHECK( moris::equal_to(nodesConnectedToElement(0), 121) );
        CHECK( moris::equal_to(nodesConnectedToElement(1), 122) );
        CHECK( moris::equal_to(nodesConnectedToElement(2), 131) );
        CHECK( moris::equal_to(nodesConnectedToElement(3), 130) );
        CHECK( moris::equal_to(nodesConnectedToElement(4), 202) );
        CHECK( moris::equal_to(nodesConnectedToElement(5), 203) );
        CHECK( moris::equal_to(nodesConnectedToElement(6), 212) );
        CHECK( moris::equal_to(nodesConnectedToElement(7), 211) );

        // ===================================================
        // Testing generate_unique_node_ids for 5 nodes
        // ===================================================

        // Ask for 5 unique node IDS
        Matrix< IdMat > tAvailableNodeIDs = tMesh3D_HEXs->generate_unique_node_ids(5);

        // Make sure it returns the correct IDs
        CHECK(moris::equal_to(tAvailableNodeIDs(0,0),730));
        CHECK(moris::equal_to(tAvailableNodeIDs(1,0),731));
        CHECK(moris::equal_to(tAvailableNodeIDs(2,0),732));
        CHECK(moris::equal_to(tAvailableNodeIDs(3,0),733));
        CHECK(moris::equal_to(tAvailableNodeIDs(4,0),734));

        // ===================================================
        // Test Cell and Vertex API
        // ===================================================
        mtk::Cell const & tCell = tMesh3D_HEXs->get_mtk_cell(elementInd);

        Matrix< IdMat > tCellVertexIds = tCell.get_vertex_ids();
        Matrix< IdMat > tCellVertexInds = tCell.get_vertex_inds();

        CHECK(all_true(tCellVertexIds == nodesConnectedToElement));
        CHECK(all_true(tCellVertexInds == tNodesConnectedToElementInds));

        // Check vertex functions
        mtk::Vertex const & tVertex  = tMesh3D_HEXs->get_mtk_vertex(nodeIndex);
        Matrix< DDRMat > tNodeCoords = tMesh3D_HEXs->get_node_coordinate(nodeIndex);
        Matrix< DDRMat > tVertCoords = tVertex.get_coords();
        CHECK(all_true(tNodeCoords == tVertCoords));

        CHECK(equal_to(tVertex.get_id(),nodeID));
        CHECK(equal_to(tVertex.get_index(),nodeIndex));

        // ===================================================
        // Dump to file
        // ===================================================
        std::string tFileOutput = "./mtk_generated_ut.exo";
        tMesh3D_HEXs->create_output_mesh(tFileOutput);

        delete tMesh3D_HEXs;

    }
}
}
}

