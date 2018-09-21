/*
 * cl_Mesh.cpp
 *
 *  Created on: Sep 17, 2018
 *      Author: doble
 */


#include "catch.hpp"


#include "fn_equal_to.hpp" // ALG/src
#include "cl_Matrix.hpp"
#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_Communication_Tools.hpp"
#include "fn_print.hpp"
#include "op_equal_equal.hpp"
#include "fn_all_true.hpp"
namespace moris
{
namespace mtk
{
TEST_CASE("New Mesh", "[moris],[mesh],[cl_Mesh],[Mesh]")
{

    // Parallel
    uint p_rank = moris::par_rank();
    uint p_size = moris::par_size();

    // File prefix
    std::string tMORISROOT = std::getenv("MORISROOT");

    std::string tPrefix =  tMORISROOT +"projects/STK/test/MeshFiles/";


    SECTION( "Reading 3D mesh from ExodusII file")
    {
        if( p_rank == 0 && p_size == 1 ) // specify it is a serial test only
        {
            // NOTE: Define the path always relative to $MORISROOT
            const std::string fileName =tPrefix + "Cube8Elems.g";    // 8 elements, 27 nodes

            // Create MORIS mesh using MTK database
            Mesh_Temp* Mesh1 = create_mesh( MeshType::STK, fileName, NULL );

            uint NumElements = Mesh1->get_num_entities(EntityRank::ELEMENT);
            uint NumNodes    = Mesh1->get_num_entities(EntityRank::NODE);

            REQUIRE( moris::equal_to(NumElements, 8) );
            REQUIRE( moris::equal_to(NumNodes, 27) );

            // ===================================================
            // Testing generate_unique_node_ids for 4 nodes
            // ===================================================

            // Ask for 5 unique node IDS
            Matrix< IdMat > tAvailableNodeIDs = Mesh1->generate_unique_entity_ids(4,EntityRank::NODE);

            REQUIRE(moris::equal_to(tAvailableNodeIDs(0),28));
            REQUIRE(moris::equal_to(tAvailableNodeIDs(1),29));
            REQUIRE(moris::equal_to(tAvailableNodeIDs(2),30));
            REQUIRE(moris::equal_to(tAvailableNodeIDs(3),31));

            // ===================================================
            // Testing elements connected to element in boundaries
            // ===================================================
            moris_id    tElemId = 5;
            moris_index tElemIndex = 4;
            Matrix< IdMat > tElemsConnectedToElemId = Mesh1->get_element_connected_to_element_glob_ids(tElemId);

            // Element Ids
            REQUIRE(moris::equal_to(tElemsConnectedToElemId(0,0),1));
            REQUIRE(moris::equal_to(tElemsConnectedToElemId(0,1),6));
            REQUIRE(moris::equal_to(tElemsConnectedToElemId(0,2),7));
            // Element face Ordinals
            REQUIRE(moris::equal_to(tElemsConnectedToElemId(1,0),1));
            REQUIRE(moris::equal_to(tElemsConnectedToElemId(1,1),2));
            REQUIRE(moris::equal_to(tElemsConnectedToElemId(1,2),5));

            Matrix< IndexMat > tElemsConnectedToElemInd = Mesh1->get_element_connected_to_element_loc_inds(tElemIndex);

            // convert local using map to global and check if they are the same
            for(uint i =0; i<tElemsConnectedToElemInd.n_cols(); i++)
            {
                REQUIRE(tElemsConnectedToElemId(0,i) == Mesh1->get_glb_entity_id_from_entity_loc_index(tElemsConnectedToElemInd(0,i),EntityRank::ELEMENT));
                REQUIRE(tElemsConnectedToElemId(1,i) == tElemsConnectedToElemInd(1,i));
            }

        }
    }

    SECTION( "Creating 8x8x8 3D mesh generated from a string")
    {
        if(p_rank == 0 && p_size == 1 ) // specify it is a serial test only
        {
            const std::string fileName2 = "generated:8x8x8";    // 512 elements, 729 nodes, 1944 edges, 1728 faces

            // Create MORIS mesh using MTK database
            Mesh_Temp* tMesh3D_HEXs = create_mesh( MeshType::STK, fileName2, NULL );

            uint NumElements2      = tMesh3D_HEXs->get_num_elems();
            uint NumNodes2         = tMesh3D_HEXs->get_num_nodes();
            uint NumEdges2         = tMesh3D_HEXs->get_num_edges( );
            uint NumFaces2         = tMesh3D_HEXs->get_num_faces( );

            REQUIRE( moris::equal_to(NumElements2, 512) );
            REQUIRE( moris::equal_to(NumNodes2, 729) );
            REQUIRE( moris::equal_to(NumEdges2, 1944) );
            REQUIRE( moris::equal_to(NumFaces2, 1728) );

            // ================================================
            // Testing entities connected to node with ID = 1
            // ================================================

            // Initialize and fill cells to store IDs of elements, faces and edges connected to current node (nodeID = 1)
            moris_id nodeID = 1;
            moris_id nodeIndex = 0;
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
            REQUIRE( moris::equal_to(NumberOfElemsConnectedToNode, 1) );
            REQUIRE( moris::equal_to(elementsConnectedToNode(0), 1) );

            // Check the number of faces and its IDs connected to current node
            REQUIRE( moris::equal_to(NumberOfFacesConnectedToNode, 3) );
            REQUIRE( moris::equal_to(facesConnectedToNode(0), 1) );
            REQUIRE( moris::equal_to(facesConnectedToNode(1), 4) );
            REQUIRE( moris::equal_to(facesConnectedToNode(2), 5) );

            // Check the number of edges and its IDs connected to current node
            REQUIRE( moris::equal_to(NumberOfEdgesConnectedToNode, 3) );
            REQUIRE( moris::equal_to(edgesConnectedToNode(0), 1) );
            REQUIRE( moris::equal_to(edgesConnectedToNode(1), 4) );
            REQUIRE( moris::equal_to(edgesConnectedToNode(2), 9) );

            // ================================================
            // Testing entities connected to edge with ID = 25
            // ================================================
            moris_id    edgeID = 25;
            moris_index edgeIND = 24;

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
            REQUIRE( moris::equal_to(NumberOfElemsConnectedToEdge, 4) );
            REQUIRE( moris::equal_to(elementsConnectedToEdge(0), 67) );
            REQUIRE( moris::equal_to(elementsConnectedToEdge(1), 68) );
            REQUIRE( moris::equal_to(elementsConnectedToEdge(2), 3) );
            REQUIRE( moris::equal_to(elementsConnectedToEdge(3), 4) );

            // Check the number of faces and its IDs connected to current edge
            REQUIRE( moris::equal_to(NumberOfFacesConnectedToEdge, 4) );
            REQUIRE( moris::equal_to(facesConnectedToEdge(0), 283) );
            REQUIRE( moris::equal_to(facesConnectedToEdge(1), 16) );
            REQUIRE( moris::equal_to(facesConnectedToEdge(2), 13) );
            REQUIRE( moris::equal_to(facesConnectedToEdge(3), 21) );

            // ================================================
            // Testing entities connected to face with ID = 25
            // ================================================

            // Get number of elements, faces and edges connected to current face
            moris_id    faceID    = 50;
            moris_index faceIndex = 49;
            Matrix< IdMat > elementsConnectedToFace        = tMesh3D_HEXs->get_elements_connected_to_face_glob_ids(faceID);
            Matrix< IndexMat > tElementsConnectedToFaceInd = tMesh3D_HEXs->get_elements_connected_to_face_loc_inds(faceIndex);

            // verify indices and ids are consistent
            tElementIdsMatch = all_true(elementsConnectedToFace == convert_entity_indices_to_ids(tElementsConnectedToFaceInd,EntityRank::ELEMENT,tMesh3D_HEXs));
            CHECK(tElementIdsMatch);

            uint NumberOfElemsConnectedToFace = elementsConnectedToFace.numel();

            // Check the number of elements and its IDs connected to current face
            REQUIRE( moris::equal_to(NumberOfElemsConnectedToFace, 2) );
            REQUIRE( moris::equal_to(elementsConnectedToFace(0), 74) );
            REQUIRE( moris::equal_to(elementsConnectedToFace(1), 10) );

            // ===================================================
            // Testing entities connected to element with ID = 100
            // ===================================================
            moris_id elementID = 100;
            moris_index elementInd = 99;

            // Initialize and fill cells to store IDs of faces, edges and nodes connected to current element (elementID = 1)
            Matrix< IdMat > elemsConnectedToElement = tMesh3D_HEXs->get_element_connected_to_element_glob_ids(elementID);
            Matrix< IndexMat > tElemsConnectedToElementInd = tMesh3D_HEXs->get_element_connected_to_element_loc_inds(elementInd);

            // Check consistency of element ids
            tElementIdsMatch = all_true(elemsConnectedToElement.get_row(0) == convert_entity_indices_to_ids( tElemsConnectedToElementInd, EntityRank::ELEMENT, tMesh3D_HEXs).get_row(0));
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
            REQUIRE( moris::equal_to(NumberOfElemsConnectedToElement, 6) );
            REQUIRE( moris::equal_to(elemsConnectedToElement(0,0), 92)  );
            REQUIRE( moris::equal_to(elemsConnectedToElement(0,1), 101) );
            REQUIRE( moris::equal_to(elemsConnectedToElement(0,2), 108) );
            REQUIRE( moris::equal_to(elemsConnectedToElement(0,3), 99)  );
            REQUIRE( moris::equal_to(elemsConnectedToElement(0,4), 36)  );
            REQUIRE( moris::equal_to(elemsConnectedToElement(0,5), 164) );

            // Check the number of faces and its IDs connected to current element
            REQUIRE( moris::equal_to(NumberOfFacesConnectedToElement, 6) );
            REQUIRE( moris::equal_to(facesConnectedToElement(0), 367) );
            REQUIRE( moris::equal_to(facesConnectedToElement(1), 391) );
            REQUIRE( moris::equal_to(facesConnectedToElement(2), 392) );
            REQUIRE( moris::equal_to(facesConnectedToElement(3), 388) );
            REQUIRE( moris::equal_to(facesConnectedToElement(4), 157) );
            REQUIRE( moris::equal_to(facesConnectedToElement(5), 393) );

            // Check the number of edges and its IDs connected to current element
            REQUIRE( moris::equal_to(NumberOfEdgesConnectedToElement, 12) );
            REQUIRE( moris::equal_to(edgesConnectedToElement(0) , 176) );
            REQUIRE( moris::equal_to(edgesConnectedToElement(1) , 218) );
            REQUIRE( moris::equal_to(edgesConnectedToElement(2) , 219) );
            REQUIRE( moris::equal_to(edgesConnectedToElement(3) , 213) );
            REQUIRE( moris::equal_to(edgesConnectedToElement(4) , 477) );
            REQUIRE( moris::equal_to(edgesConnectedToElement(5) , 502) );
            REQUIRE( moris::equal_to(edgesConnectedToElement(6) , 503) );
            REQUIRE( moris::equal_to(edgesConnectedToElement(7) , 499) );
            REQUIRE( moris::equal_to(edgesConnectedToElement(8) , 475) );
            REQUIRE( moris::equal_to(edgesConnectedToElement(9) , 478) );
            REQUIRE( moris::equal_to(edgesConnectedToElement(10), 504) );
            REQUIRE( moris::equal_to(edgesConnectedToElement(11), 501) );

            // Check the number of nodes and its IDs connected to current element
            REQUIRE( moris::equal_to(NumberOfNodesConnectedToElement, 8) );
            REQUIRE( moris::equal_to(nodesConnectedToElement(0), 121) );
            REQUIRE( moris::equal_to(nodesConnectedToElement(1), 122) );
            REQUIRE( moris::equal_to(nodesConnectedToElement(2), 131) );
            REQUIRE( moris::equal_to(nodesConnectedToElement(3), 130) );
            REQUIRE( moris::equal_to(nodesConnectedToElement(4), 202) );
            REQUIRE( moris::equal_to(nodesConnectedToElement(5), 203) );
            REQUIRE( moris::equal_to(nodesConnectedToElement(6), 212) );
            REQUIRE( moris::equal_to(nodesConnectedToElement(7), 211) );

            // ===================================================
            // Testing generate_unique_node_ids for 5 nodes
            // ===================================================

            // Ask for 5 unique node IDS
            Matrix< IdMat > tAvailableNodeIDs = tMesh3D_HEXs->generate_unique_node_ids(5);

            // Make sure it returns the correct IDs
            REQUIRE(moris::equal_to(tAvailableNodeIDs(0,0),730));
            REQUIRE(moris::equal_to(tAvailableNodeIDs(1,0),731));
            REQUIRE(moris::equal_to(tAvailableNodeIDs(2,0),732));
            REQUIRE(moris::equal_to(tAvailableNodeIDs(3,0),733));
            REQUIRE(moris::equal_to(tAvailableNodeIDs(4,0),734));


            // ===================================================
            // Test Cell and Vertex API
            // ===================================================
            mtk::Cell const & tCell = tMesh3D_HEXs->get_mtk_cell(elementInd);

            Matrix< IdMat > tCellVertexIds = tCell.get_vertex_ids();
            Matrix< IdMat > tCellVertexInds = tCell.get_vertex_inds();

            REQUIRE(all_true(tCellVertexIds == nodesConnectedToElement));
            REQUIRE(all_true(tCellVertexInds == tNodesConnectedToElementInds));


            // Check vertex functions
            mtk::Vertex const & tVertex  = tMesh3D_HEXs->get_mtk_vertex(nodeIndex);
            Matrix< DDRMat > tNodeCoords = tMesh3D_HEXs->get_node_coordinate(nodeIndex);
            Matrix< DDRMat > tVertCoords = tVertex.get_coords();
            REQUIRE(all_true(tNodeCoords == tVertCoords));

            REQUIRE(equal_to(tVertex.get_id(),nodeID));
            REQUIRE(equal_to(tVertex.get_index(),nodeIndex));


        }
    }
}





}
}
