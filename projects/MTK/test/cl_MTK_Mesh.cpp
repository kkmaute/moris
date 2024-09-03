/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Mesh.cpp
 *
 */

#include "catch.hpp"

#include "paths.hpp"

// MTK includes
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"

// Other MORIS includes
#include "cl_Communication_Tools.hpp"
#include "cl_Matrix.hpp"
#include "op_equal_equal.hpp"
#include "fn_all_true.hpp"
#include "fn_equal_to.hpp"

namespace moris::mtk
{
    TEST_CASE( "Reading 3D mesh from ExodusII file", "[moris],[mesh],[cl_Mesh],[Mesh]" )
    {

        // Parallel
        uint p_rank = moris::par_rank();
        uint p_size = moris::par_size();

        // File prefix
        std::string tMORISROOT = moris::get_base_moris_dir();

        std::string tPrefix = tMORISROOT + "projects/MTK/test/Test_Files/";

        if ( p_rank == 0 && p_size == 1 )    // specify it is a serial test only
        {
            // NOTE: Define the path always relative to $MORISROOT
            const std::string fileName = tPrefix + "Cube8Elems.g";    // 8 elements, 27 nodes

            // Create MORIS mesh using MTK database
            Mesh* Mesh1 = create_interpolation_mesh( MeshType::STK, fileName, nullptr );

            uint NumElements = Mesh1->get_num_entities( EntityRank::ELEMENT );
            uint NumNodes    = Mesh1->get_num_entities( EntityRank::NODE );

            CHECK( moris::equal_to( NumElements, (uint)8 ) );
            CHECK( moris::equal_to( NumNodes, 27 ) );

            // ===================================================
            // Testing elements connected to element in boundaries
            // ===================================================
            moris_id    tElemId    = 5;
            moris_index tElemIndex = 4;

            Matrix< IdMat > tElemsConnectedToElemId = Mesh1->get_element_connected_to_element_glob_ids( tElemId );

            // Element Ids
            CHECK( moris::equal_to( tElemsConnectedToElemId( 0, 0 ), 1 ) );
            CHECK( moris::equal_to( tElemsConnectedToElemId( 0, 1 ), 6 ) );
            CHECK( moris::equal_to( tElemsConnectedToElemId( 0, 2 ), 7 ) );
            // Element face Ordinals
            CHECK( moris::equal_to( tElemsConnectedToElemId( 1, 0 ), 1 ) );
            CHECK( moris::equal_to( tElemsConnectedToElemId( 1, 1 ), 2 ) );
            CHECK( moris::equal_to( tElemsConnectedToElemId( 1, 2 ), 5 ) );

            Matrix< IndexMat > tElemsConnectedToElemInd = Mesh1->get_elements_connected_to_element_and_face_ord_loc_inds( tElemIndex );

            // convert local using map to global and check if they are the same
            for ( uint i = 0; i < tElemsConnectedToElemInd.n_cols(); i++ )
            {
                CHECK( tElemsConnectedToElemId( 0, i ) == Mesh1->get_glb_entity_id_from_entity_loc_index( tElemsConnectedToElemInd( 0, i ), EntityRank::ELEMENT ) );
                CHECK( tElemsConnectedToElemId( 1, i ) == tElemsConnectedToElemInd( 1, i ) );
            }

            // ===================================================
            // Testing generate_unique_node_ids for 4 nodes
            // and adding those nodes to the mesh
            // ===================================================

            // Ask for 5 unique node IDS
            Matrix< IdMat > tAvailableNodeIDs = Mesh1->generate_unique_entity_ids( 4, EntityRank::NODE );
            CHECK( moris::equal_to( tAvailableNodeIDs( 0 ), 28 ) );
            CHECK( moris::equal_to( tAvailableNodeIDs( 1 ), 29 ) );
            CHECK( moris::equal_to( tAvailableNodeIDs( 2 ), 30 ) );
            CHECK( moris::equal_to( tAvailableNodeIDs( 3 ), 31 ) );

            delete Mesh1;
        }
    }
    TEST_CASE( "Creating 8x8x8 3D mesh generated from a string", "[MTK_MESH_1]" )
    {    // Parallel
        uint p_size = moris::par_size();

        // File prefix
        std::string tMORISROOT = moris::get_base_moris_dir();
        if ( p_size == 1 )    // specify it is a serial test only
        {
            const std::string fileName2 = "generated:8x8x8";
            // const std::string fileName2 = "generated:8x8x8|sideset:xXyYzZ";    // 512 elements, 729 nodes, 1944 edges, 1728 faces

            // Create MORIS mesh using MTK database
            Mesh* tMesh3D_HEXs = create_interpolation_mesh( MeshType::STK, fileName2, nullptr, true );

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
            tElementIdsMatch = all_true( elementsConnectedToFace == convert_entity_indices_to_ids( tElementsConnectedToFaceInd, EntityRank::ELEMENT, tMesh3D_HEXs ) );
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

    TEST_CASE( "Testing a side set on an 8x8x8 generated mesh", "[MTK_MESH_1_SIDE_SET]" )
    {    // Parallel
        uint p_size = moris::par_size();

        // File prefix
        std::string tMORISROOT = moris::get_base_moris_dir();
        if ( p_size == 1 )    // specify it is a serial test only
        {
            const std::string fileName2 = "generated:8x8x8|sideset:xXyYzZ";    // 512 elements, 729 nodes, 1944 edges, 1728 faces

            // Create MORIS mesh using MTK database
            Mesh* tMeshWithSideSets     = create_interpolation_mesh( MeshType::STK, fileName2, nullptr );
            tMeshWithSideSets->mVerbose = false;

            Vector< std::string > tGoldSideSetNames = { { "surface_1" },
                { "surface_2" },
                { "surface_3" },
                { "surface_4" },
                { "surface_5" },
                { "surface_6" } };

            // ===================================================
            // CHECK SIDE SET FUNCTIONS
            // ===================================================

            // verify names
            Vector< std::string > tSetNames = tMeshWithSideSets->get_set_names( EntityRank::FACE );

            CHECK( tSetNames.size() / 2 == tGoldSideSetNames.size() );

            for ( moris::uint i = 0; i < tSetNames.size() / 2; i++ )
            {

                // verify the expected names show up
                auto tIt = std::find( tGoldSideSetNames.begin(), tGoldSideSetNames.end(), tSetNames( i ) );

                CHECK( tIt != tGoldSideSetNames.end() );

                // get the cells in this side set
                Matrix< IndexMat > tCellIndsInSideSet1;
                Matrix< IndexMat > tCellSideOrdsInSideSet1;
                tMeshWithSideSets->get_sideset_elems_loc_inds_and_ords( tSetNames( i ), tCellIndsInSideSet1, tCellSideOrdsInSideSet1 );

                // verify there are the expected number of elements in the side set
                CHECK( tCellIndsInSideSet1.numel() == 64 );
                CHECK( tCellSideOrdsInSideSet1.numel() == 64 );

                // get the cell ptrs
                Vector< mtk::Cell const * > tCellsInSideSet;
                Matrix< IndexMat >          tCellSideOrdsInSideSet2;
                tMeshWithSideSets->get_sideset_cells_and_ords( tSetNames( i ), tCellsInSideSet, tCellSideOrdsInSideSet2 );

                CHECK( tCellsInSideSet.size() == 64 );
                CHECK( tCellSideOrdsInSideSet2.numel() == 64 );

                CHECK( all_true( tCellSideOrdsInSideSet2 == tCellSideOrdsInSideSet1 ) );

                // retrieve cell inds
                Matrix< IndexMat > tCellIndsInSideSet2( 1, tCellsInSideSet.size() );

                for ( moris::uint iC = 0; iC < tCellsInSideSet.size(); iC++ )
                {
                    tCellIndsInSideSet2( iC ) = tCellsInSideSet( iC )->get_index();
                }

                CHECK( all_true( tCellIndsInSideSet2 == tCellIndsInSideSet1 ) );
            }

            // ===================================================
            // Dump to file
            // ===================================================
            std::string tFileOutput = "./mtk_sideset_ut.exo";
            tMeshWithSideSets->create_output_mesh( tFileOutput );

            delete tMeshWithSideSets;
        }
    }

    TEST_CASE( "Parallel Generated Mesh", "[MTK_2PROC]" )
    {
        if ( par_size() == 2 )
        {
            std::string fileName2   = "generated:1x1x2|sideset:XYZ|nodeset:X";
            std::string tFileOutput = "./mtk_2_proc_test.exo";
            // Create MORIS mesh using MTK database
            Mesh* tParMesh = create_interpolation_mesh( MeshType::STK, fileName2 );

            // Each processor has the full mesh because of the aura
            if ( par_rank() == 0 )
            {

                CHECK( tParMesh->get_num_entities( EntityRank::ELEMENT ) == 2 );
                CHECK( tParMesh->get_num_entities( EntityRank::NODE ) == 12 );
            }

            if ( par_rank() == 1 )
            {
                CHECK( tParMesh->get_num_entities( EntityRank::ELEMENT ) == 2 );
                CHECK( tParMesh->get_num_entities( EntityRank::NODE ) == 12 );
            }
            tParMesh->create_output_mesh( tFileOutput );

            delete tParMesh;
        }
    }

    TEST_CASE( "MTK Mesh from file via STK, with a fields not on the file declared", "[Mesh_File_Field]" )
    {
        // NOTE: Define the path always relative to $MORISROOT
        const std::string fileName = "generated:2x2x4";

        // Declare scalar node field
        moris::mtk::Scalar_Field_Info< DDRMat > tNodeField1;
        std::string                             tFieldName1 = "node_field_1";
        tNodeField1.set_field_name( tFieldName1 );
        tNodeField1.set_field_entity_rank( EntityRank::NODE );

        // Declare scalar node field
        moris::mtk::Scalar_Field_Info< DDRMat > tNodeField2;
        std::string                             tFieldName2 = "node_field_2";
        tNodeField2.set_field_name( tFieldName2 );
        tNodeField2.set_field_entity_rank( EntityRank::NODE );

        // Declare Element field
        moris::mtk::Scalar_Field_Info< DDRMat > tElementField1;
        std::string                             tFieldName3 = "elem_field_1";
        tElementField1.set_field_name( tFieldName3 );
        tElementField1.set_field_entity_rank( EntityRank::ELEMENT );

        // Initialize field information container
        moris::mtk::MtkFieldsInfo tFieldsInfo;

        // Place the node field into the field info container
        add_field_for_mesh_input( &tNodeField1, tFieldsInfo );
        add_field_for_mesh_input( &tNodeField2, tFieldsInfo );
        add_field_for_mesh_input( &tElementField1, tFieldsInfo );

        // Declare some supplementary fields
        MtkMeshData tMeshData;
        tMeshData.FieldsInfo = &tFieldsInfo;

        // Create MORIS mesh using MTK database
        Mesh* Mesh1 = create_interpolation_mesh( MeshType::STK, fileName, &tMeshData );

        // add the field data for node field 1
        Matrix< DDRMat > tFieldData1( Mesh1->get_num_entities( EntityRank::NODE ), 1 );
        tFieldData1.fill( 10.0 );

        Mesh1->add_mesh_field_real_scalar_data_loc_inds( tFieldName1, EntityRank::NODE, tFieldData1 );

        // add the field data for mesh field 2
        Matrix< DDRMat > tFieldData2( Mesh1->get_num_entities( EntityRank::NODE ), 1 );
        tFieldData2.fill( -10.0 );

        Mesh1->add_mesh_field_real_scalar_data_loc_inds( tFieldName2, EntityRank::NODE, tFieldData2 );

        // add the field data for element field 1
        Matrix< DDRMat > tFieldData3( Mesh1->get_num_entities( EntityRank::ELEMENT ), 1 );
        tFieldData3.fill( -11.0 );

        Mesh1->add_mesh_field_real_scalar_data_loc_inds( tFieldName3, EntityRank::ELEMENT, tFieldData3 );

        CHECK( Mesh1->get_entity_field_value_real_scalar( { { 0 } }, tFieldName1, EntityRank::NODE )( 0 ) == 10.0 );
        CHECK( Mesh1->get_entity_field_value_real_scalar( { { 0 } }, tFieldName2, EntityRank::NODE )( 0 ) == -10.0 );
        CHECK( Mesh1->get_entity_field_value_real_scalar( { { 0 } }, tFieldName3, EntityRank::ELEMENT )( 0 ) == -11.0 );

        // Verify Field numbers
        CHECK( Mesh1->get_num_fields( EntityRank::NODE ) == 3 ); /*1 for coordinate field*/
        CHECK( Mesh1->get_num_fields( EntityRank::ELEMENT ) == 1 );

        // verify ordinals
        CHECK( Mesh1->get_field_ind( tFieldName1, EntityRank::NODE ) == 1 );
        CHECK( Mesh1->get_field_ind( tFieldName2, EntityRank::NODE ) == 2 );
        CHECK( Mesh1->get_field_ind( tFieldName3, EntityRank::ELEMENT ) == 0 );

        // output mesh
        std::string tMeshOutputFile = "./MTK_Mesh_File_Data.e";
        Mesh1->create_output_mesh( tMeshOutputFile );

        delete Mesh1;
    }

}    // namespace moris::mtk
