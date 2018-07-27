// Third-party header files.
#include <catch.hpp>
#include <iostream>

// MORIS project header files.
#include "algorithms.hpp"
#include "cl_Hierarchical_Mesh.hpp" // STK/src/Hierarchical
#include "cl_Mesh.hpp" // MTK/src
#include "cl_Database.hpp" // MTK/src
#include "cl_STK_Implementation.hpp" // STK/src

#include "cl_Communication_Tools.hpp" // COM/src

using namespace moris;
// ----------------------------------------------------------------------------

TEST_CASE("MtkMeshFromString", "[moris],[mesh],[cl_Mesh],[Mesh]")
{
    // Parallel
    uint p_rank = 0;
    uint p_size = 1;
#ifdef MORIS_HAVE_PARALLEL
    p_rank = moris::par_rank();
    p_size = moris::par_size();
#endif

    SECTION( "Reading 3D mesh from ExodusII file")
    {
        if( p_rank == 0 && p_size == 1 ) // specify it is a serial test only
        {
            // NOTE: Define the path always relative to $MORISROOT
            const std::string fileName = "projects/STK/test/MeshFiles/Cube8Elems.g";    // 8 elements, 27 nodes

            // Create MORIS mesh using MTK database
            moris::mesh Mesh1( MeshType::MTK, fileName );

            uint NumElements      = Mesh1.get_num_elems();
            uint NumNodes         = Mesh1.get_num_nodes();

            REQUIRE( moris::equal_to(NumElements, 8) );
            REQUIRE( moris::equal_to(NumNodes, 27) );

            // ===================================================
            // Testing generate_unique_node_ids for 4 nodes
            // ===================================================

            // Ask for 5 unique node IDS
            moris::Mat<moris::uint> tAvailableNodeIDs = Mesh1.generate_unique_node_ids(4);

            REQUIRE(moris::equal_to(tAvailableNodeIDs(0,0),28));
            REQUIRE(moris::equal_to(tAvailableNodeIDs(1,0),29));
            REQUIRE(moris::equal_to(tAvailableNodeIDs(2,0),30));
            REQUIRE(moris::equal_to(tAvailableNodeIDs(3,0),31));

            // ===================================================
            // Testing elements connected to element in boundaries
            // ===================================================
            moris::Mat< uint > tElemsConnectedToElem = Mesh1.get_elements_connected_to_element(5);

            REQUIRE(moris::equal_to(tElemsConnectedToElem(0),1));
            REQUIRE(moris::equal_to(tElemsConnectedToElem(1),6));
            REQUIRE(moris::equal_to(tElemsConnectedToElem(2),7));

            tElemsConnectedToElem = Mesh1.get_elements_connected_to_element(3);

            REQUIRE(moris::equal_to(tElemsConnectedToElem(0),2));
            REQUIRE(moris::equal_to(tElemsConnectedToElem(1),8));
            REQUIRE(moris::equal_to(tElemsConnectedToElem(2),4));
        }
    }

    SECTION( "Creating 8x8x8 3D mesh generated from a string")
    {
        if(p_rank == 0 && p_size == 1 ) // specify it is a serial test only
        {
            const std::string fileName2 = "generated:8x8x8";    // 512 elements, 729 nodes, 1944 edges, 1728 faces

            // Create MORIS mesh using MTK database
            moris::mesh tMesh3D_HEXs( MeshType::MTK, fileName2 );

            uint NumElements2      = tMesh3D_HEXs.get_num_elems();
            uint NumNodes2         = tMesh3D_HEXs.get_num_nodes();
            uint NumEdges2         = tMesh3D_HEXs.get_num_edges( );
            uint NumFaces2         = tMesh3D_HEXs.get_num_faces( );

            REQUIRE( moris::equal_to(NumElements2, 512) );
            REQUIRE( moris::equal_to(NumNodes2, 729) );
            REQUIRE( moris::equal_to(NumEdges2, 1944) );
            REQUIRE( moris::equal_to(NumFaces2, 1728) );

            // ================================================
            // Testing entities connected to node with ID = 1
            // ================================================

            // Initialize and fill cells to store IDs of elements, faces and edges connected to current node (nodeID = 1)
            uint nodeID = 1;
            Mat< uint > elementsConnectedToNode = tMesh3D_HEXs.get_elements_connected_to_node(nodeID);
            Mat< uint > facesConnectedToNode    = tMesh3D_HEXs.get_faces_connected_to_node(nodeID);
            Mat< uint > edgesConnectedToNode    = tMesh3D_HEXs.get_edges_connected_to_node(nodeID);

            // Get number of elements, faces and edges connected to current node
            uint NumberOfElemsConnectedToNode = elementsConnectedToNode.n_rows();
            uint NumberOfFacesConnectedToNode = facesConnectedToNode.n_rows();
            uint NumberOfEdgesConnectedToNode = edgesConnectedToNode.n_rows();

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
            uint edgeID = 25;

            // Initialize and fill cells to store IDs of elements and faces connected to current edge (edgeID = 25)
            Mat< uint > elementsConnectedToEdge = tMesh3D_HEXs.get_elements_connected_to_edge(edgeID);
            Mat< uint > facesConnectedToEdge    = tMesh3D_HEXs.get_faces_connected_to_edge(edgeID);

            // Get number of elements, faces and edges connected to current edge
            uint NumberOfElemsConnectedToEdge = elementsConnectedToEdge.n_rows();
            uint NumberOfFacesConnectedToEdge = facesConnectedToEdge.n_rows();

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
            uint faceID = 50;
            Mat< uint > elementsConnectedToFace = tMesh3D_HEXs.get_elements_connected_to_face(faceID);
            uint NumberOfElemsConnectedToFace = elementsConnectedToFace.n_rows();

            // Check the number of elements and its IDs connected to current face
            REQUIRE( moris::equal_to(NumberOfElemsConnectedToFace, 2) );
            REQUIRE( moris::equal_to(elementsConnectedToFace(0), 74) );
            REQUIRE( moris::equal_to(elementsConnectedToFace(1), 10) );

            // ===================================================
            // Testing entities connected to element with ID = 100
            // ===================================================
            uint elementID = 100;

            // Initialize and fill cells to store IDs of faces, edges and nodes connected to current element (elementID = 1)

            Mat< uint > elemsConnectedToElement = tMesh3D_HEXs.get_elements_connected_to_element(elementID);
            Mat< uint > facesConnectedToElement = tMesh3D_HEXs.get_faces_connected_to_element(elementID);
            Mat< uint > edgesConnectedToElement = tMesh3D_HEXs.get_edges_connected_to_element(elementID);
            Mat< uint > nodesConnectedToElement = tMesh3D_HEXs.get_nodes_connected_to_element(elementID);

            // Get number of elements, faces and edges connected to current node
            uint NumberOfElemsConnectedToElement = elemsConnectedToElement.length();
            uint NumberOfFacesConnectedToElement = facesConnectedToElement.length();
            uint NumberOfEdgesConnectedToElement = edgesConnectedToElement.length();
            uint NumberOfNodesConnectedToElement = nodesConnectedToElement.length();

            // Check the number of elements and its IDs connected to current element
            REQUIRE( moris::equal_to(NumberOfElemsConnectedToElement, 6) );
            REQUIRE( moris::equal_to(elemsConnectedToElement(0), 92)  );
            REQUIRE( moris::equal_to(elemsConnectedToElement(1), 101) );
            REQUIRE( moris::equal_to(elemsConnectedToElement(2), 108) );
            REQUIRE( moris::equal_to(elemsConnectedToElement(3), 99)  );
            REQUIRE( moris::equal_to(elemsConnectedToElement(4), 36)  );
            REQUIRE( moris::equal_to(elemsConnectedToElement(5), 164) );

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
            Mat< uint > tAvailableNodeIDs = tMesh3D_HEXs.generate_unique_node_ids(5);

            // Make sure it returns the correct IDs
            REQUIRE(moris::equal_to(tAvailableNodeIDs(0,0),730));
            REQUIRE(moris::equal_to(tAvailableNodeIDs(1,0),731));
            REQUIRE(moris::equal_to(tAvailableNodeIDs(2,0),732));
            REQUIRE(moris::equal_to(tAvailableNodeIDs(3,0),733));
            REQUIRE(moris::equal_to(tAvailableNodeIDs(4,0),734));

        }
    }

    SECTION( "Creating a 3D mesh from string in parallel and checking parallel data in mesh")
    {
        if(  p_size == 2 ) // specify it is a serial test only
        {
            // Test to verify a mesh generated from data gets expected basic shared information
            // and creates the aura and other parts created by default in STK.

            const std::string tFileName = "generated:1x1x2";    // 512 elements, 729 nodes, 1944 edges, 1728 faces
            moris::mesh tParallelMesh( MeshType::MTK, tFileName );

            // ========================================
            // Testing get_entities_owned_current_proc
            // ========================================

            Mat< uint > tNodeOwned  = tParallelMesh.get_entities_owned_current_proc(EntityRank::NODE);
            Mat< uint > tNodeShared  = tParallelMesh.get_entities_glb_shared_current_proc(EntityRank::NODE);
            Mat< uint > tNodeOwnedAndShared  = tParallelMesh.get_entities_owned_and_shared_by_current_proc(EntityRank::NODE);
            Mat< uint > tNodeAura = tParallelMesh.get_entities_in_aura(EntityRank::NODE);
            Mat< uint > tNodeUniv = tParallelMesh.get_entities_universal(EntityRank::NODE);

            Mat< uint > tElemOwned  = tParallelMesh.get_entities_owned_current_proc(EntityRank::ELEMENT);
            Mat< uint > tElemOwnedAndShared  = tParallelMesh.get_entities_owned_and_shared_by_current_proc(EntityRank::ELEMENT);
            Mat< uint > tElemAura = tParallelMesh.get_entities_in_aura(EntityRank::ELEMENT);
            Mat< uint > tElemUniv = tParallelMesh.get_entities_universal(EntityRank::ELEMENT);

            if( p_rank == 0 )
            {
                // Check some of the requested data
                REQUIRE(moris::equal_to(tNodeOwned(3,0),4));
                REQUIRE(moris::equal_to(tNodeShared(1,0),6));
                REQUIRE(moris::equal_to(tNodeOwnedAndShared(7,0),8));
                REQUIRE(moris::equal_to(tNodeAura(3,0),12));
                REQUIRE(moris::equal_to(tNodeUniv(11,0),12));

                // Elements don't have shared entities
                REQUIRE(moris::equal_to(tElemOwned(0,0),1));
                REQUIRE(moris::equal_to(tElemAura(0,0),2));
                REQUIRE(moris::equal_to(tElemUniv(1,0),2));
            } else {

                // Check some of the requested data
                REQUIRE(moris::equal_to(tNodeOwned(0,0),9));
                REQUIRE(moris::equal_to(tNodeShared(0,0),5));
                REQUIRE(moris::equal_to(tNodeOwnedAndShared(7,0),8));
                REQUIRE(moris::equal_to(tNodeAura(3,0),4));
                REQUIRE(moris::equal_to(tNodeUniv(11,0),12));

                // Elements don't have shared entities
                REQUIRE(moris::equal_to(tElemOwned(0,0),2));
                REQUIRE(moris::equal_to(tElemAura(0,0),1));
                REQUIRE(moris::equal_to(tElemUniv(1,0),2));
            }

        }
    }

    SECTION( "Writing Mesh to an Exodus File")
    {
        //        const std::string fileName4 = "generated:2x2x2";    // 512 elements, 729 nodes, 1944 edges, 1728 faces
        //        std::string OutputFileName1 = "test/src/mesh/TestOuputMesh_.e";  //output name
        //
        //
        //        // Create MORIS mesh using MTK database
        //        moris::mesh* Mesh4 = MorisFactory.create_mesh( MeshType::MTK, fileName4 );
        //
        //        Mesh4.create_output_mesh(OutputFileName1);
        //        delete Mesh4;
    }

    SECTION( "String Generated 1x1x2 Mesh With AutoAura Option")
    {
        if(p_size == 2)
        {
            const std::string AuraFile = "generated:1x1x2";    // 512 elements, 729 nodes, 1944 edges, 1728 faces

            // Create MORIS mesh using MTK database
            moris::mesh AuraMesh( MeshType::MTK, AuraFile );

            // Check number of entities in aura
            uint tNumNodesInAura = AuraMesh.get_num_entities_aura(EntityRank::NODE);
            uint tNumEdgesInAura = AuraMesh.get_num_entities_aura(EntityRank::EDGE);
            uint tNumFacesInAura = AuraMesh.get_num_entities_aura(EntityRank::FACE);
            uint tNumElemsInAura = AuraMesh.get_num_entities_aura(EntityRank::ELEMENT);
            REQUIRE(moris::equal_to(tNumNodesInAura,4));
            REQUIRE(moris::equal_to(tNumEdgesInAura,8));
            REQUIRE(moris::equal_to(tNumFacesInAura,5));
            REQUIRE(moris::equal_to(tNumElemsInAura,1));

            // Check Coordinates of nodes in aura
            Mat< real > tNodeCoords = AuraMesh.get_all_nodes_coords_aura();

            // Check Ownership
            // For proc 0 - this is a locally owned part.
            Mat< uint > tProcsForSharedEntity = AuraMesh.get_procs_sharing_entity_by_id(6,EntityRank::NODE);
            Mat< uint > tProcsForNotSharedEntity = AuraMesh.get_procs_sharing_entity_by_id(2,EntityRank::NODE);

            if(p_rank == 0)
            {
                REQUIRE(moris::equal_to(tProcsForSharedEntity(0),1));
                REQUIRE(moris::equal_to(tProcsForNotSharedEntity(0),UINT_MAX));
            }
            else if(p_rank ==1)
            {
                REQUIRE(moris::equal_to(tProcsForSharedEntity(0),0));
                REQUIRE(moris::equal_to(tProcsForNotSharedEntity(0),UINT_MAX));
            }
        }
    }

    SECTION( "Reading 3D Mesh from exodusII file")
    {
        if( p_rank == 0 && p_size == 1 ) // specify it is a serial test only
        {
            // This example shows how to access the nodes for a subset of the mesh (nodesets and sidesets).
            //
            //   *-----*-----*  nodes: 2,1,3
            //   |     |     |  elements: 2,1
            //   |     |     |
            //   *-----*-----*  nodes: 5,9,8
            //   |     |     |  elements: 4,3
            //   |     |     |
            //   *-----*-----*  nodes 4,7,6
            //
            // NodeSet 1 = nodes(1,2,3)
            // NodeSet 2 = nodes(4)
            // SideSet 1 = nodes(3,7,9)
            // SideSet 2 = nodes(6,8)

            // NOTE: Define the path always relative to $MORISROOT
            const std::string fileName = "/test/src/mesh/MeshFiles/SidesAndNodeSets2DMesh.g";    // 4 elements, 9 nodes

            // Create MORIS mesh using MTK database
            moris::mesh MeshSetsTests( MeshType::MTK, fileName );

            // Get basic mesh information
            uint tNumElements = MeshSetsTests.get_num_elems();
            uint tNumNodes    = MeshSetsTests.get_num_nodes();

            // Get nodes contained in node sets
            Mat< uint > tNodesInNodeSet1 = MeshSetsTests.get_nodes_in_node_set(1);
            Mat< uint > tNodesInNodeSet2 = MeshSetsTests.get_nodes_in_node_set(2);

            // Get nodes contained in side sets
            Mat< uint > tNodesInSideSet1 = MeshSetsTests.get_nodes_in_side_set(1);
            Mat< uint > tNodesInSideSet2 = MeshSetsTests.get_nodes_in_side_set(2);

            // Get nodes contained in block sets
           Mat< uint > tNodesInBlockSet1 = MeshSetsTests.get_nodes_in_block_set(1);

            // Get edges contained in side set
           Mat< uint > tEdgesInSideSet1 = MeshSetsTests.get_edges_in_side_set(1);

            // Verify basic information read is correct
            REQUIRE(moris::equal_to(tNumElements,4));
            REQUIRE(moris::equal_to(tNumNodes,9));

            // Verify that node sets contained the expected nodes
            REQUIRE(moris::equal_to(tNodesInNodeSet1(0,0),1));
            REQUIRE(moris::equal_to(tNodesInNodeSet1(1,0),2));
            REQUIRE(moris::equal_to(tNodesInNodeSet1(2,0),3));
            REQUIRE(moris::equal_to(tNodesInNodeSet2(0,0),4));

            // Verify that side sets contained the expected nodes
            REQUIRE(moris::equal_to(tNodesInSideSet1(0,0),3));
            REQUIRE(moris::equal_to(tNodesInSideSet1(1,0),7));
            REQUIRE(moris::equal_to(tNodesInSideSet1(2,0),9));
            REQUIRE(moris::equal_to(tNodesInSideSet2(0,0),6));
            REQUIRE(moris::equal_to(tNodesInSideSet2(1,0),8));

            // Verify that block set contained the expected nodes
            REQUIRE(moris::equal_to(tNodesInBlockSet1(0,0),1));
            REQUIRE(moris::equal_to(tNodesInBlockSet1(1,0),2));
            REQUIRE(moris::equal_to(tNodesInBlockSet1(2,0),3));
            REQUIRE(moris::equal_to(tNodesInBlockSet1(3,0),4));
            REQUIRE(moris::equal_to(tNodesInBlockSet1(4,0),5));
            REQUIRE(moris::equal_to(tNodesInBlockSet1(5,0),6));
            REQUIRE(moris::equal_to(tNodesInBlockSet1(6,0),7));
            REQUIRE(moris::equal_to(tNodesInBlockSet1(7,0),8));
            REQUIRE(moris::equal_to(tNodesInBlockSet1(8,0),9));

            // Verify that side sets contained the expected edges (no faces in 2D)
            REQUIRE(moris::equal_to(tEdgesInSideSet1(0,0),13));
            REQUIRE(moris::equal_to(tEdgesInSideSet1(1,0),25));
            REQUIRE(moris::equal_to(tEdgesInSideSet1(2,0),33));
            REQUIRE(moris::equal_to(tEdgesInSideSet1(3,0),45));
        }
    }

    SECTION("Compute centroids of given block and sidesets" )
    {
        if ( p_rank == 0 && p_size == 1 )
        {
            // In this section, two sets (one of elements and another one of sides) are accessed
            // and the centroid of each entity is computed. On top of that, when accessing sidesets,
            // also request the centroids of adjacent elements and compute the mean between them.

            bool printCentroidFlag = false;

            // String for MTK to generate mesh with sideset
            const std::string meshName = "generated:2x2x2|sideset:X";    // 512 elements, 729 nodes, 1944 edges, 1728 faces

            // Create MORIS mesh using MTK database
            moris::mesh meshWithSideSet( MeshType::MTK, meshName );

            // ---------------------------------------------------
            // Testing block set (universal, all elements in mesh)
            // ---------------------------------------------------

            uint tNumElems          = meshWithSideSet.get_num_elems();
            Mat< uint > tElemsOwned = meshWithSideSet.get_entities_owned_current_proc(EntityRank::ELEMENT);
            uint tElemTopoNodes     = meshWithSideSet.get_elem_topology_num_nodes( tElemsOwned(0) );

            // Initialize matrix where the centroids will be stored
            uint tNumDims = 3;
            Mat< real > tCentroid ( tNumDims, 1, 0.0 );

            // Get the centroids of elements by accessing their node connectivity and coordinates of each node
            for ( uint iElem = 0; iElem < tNumElems; ++iElem )
            {
                Mat< uint > tNodesInElem = meshWithSideSet.get_nodes_connected_to_element( tElemsOwned( iElem ) );
                Mat< real > tNodeCoords  = meshWithSideSet.get_selected_nodes_coords( tNodesInElem );

                // Centroid is defined as the average of the coordinates of all the nodes
                // Sum all coordinates components divided by number of vertices in polygon
                for ( uint iNode = 0; iNode < tElemTopoNodes; ++iNode )
                {
                    tCentroid( 0 ) += tNodeCoords( iNode, 0 ) / tElemTopoNodes;
                    tCentroid( 1 ) += tNodeCoords( iNode, 1 ) / tElemTopoNodes;
                    tCentroid( 2 ) += tNodeCoords( iNode, 2 ) / tElemTopoNodes;
                }

                if ( printCentroidFlag )
                {
                    std::cout << "Centroid of element "<< tElemsOwned( iElem ) << " is ("
                            << tCentroid( 0 ) << ", "
                            << tCentroid( 1 ) << ", "
                            << tCentroid( 2 ) << ")" << std::endl;
                }

                // Check centroid values for element with ID = 4
                if ( tElemsOwned( iElem ) == 4 )
                {
                    REQUIRE(moris::equal_to(tCentroid( 0 ),1.5));
                    REQUIRE(moris::equal_to(tCentroid( 1 ),1.5));
                    REQUIRE(moris::equal_to(tCentroid( 2 ),0.5));
                }

                tCentroid.fill( 0.0 );
            }

            // ---------------------------------------------
            // Testing side set (only one side set provided)
            // ---------------------------------------------

            Mat< uint > tFacesInSideSet = meshWithSideSet.get_faces_in_side_set( 1 );
            uint tNumFaces              = tFacesInSideSet.length();

            // Initialize matrix where the centroids will be stored
            Cell< Mat< real > > tAllCentroids( 2 );
            tAllCentroids ( 0 ) = tCentroid;
            tAllCentroids ( 1 ) = tCentroid;

            // Get the centroids of elements by accessing their node connectivity and coordinates of each node
            for ( uint iFace = 0; iFace < tNumFaces; ++iFace )
            {
                Mat< uint > tElemsSharingFace = meshWithSideSet.get_elements_connected_to_face( tFacesInSideSet( iFace ) );
                uint tNumElemsSharingFace     = tElemsSharingFace.length();

                for ( uint iElem = 0; iElem < tNumElemsSharingFace; ++iElem )
                {
                    Mat< uint > tNodesInElem      = meshWithSideSet.get_nodes_connected_to_element( tElemsSharingFace( iElem ) );
                    Mat< real > tNodeCoords       = meshWithSideSet.get_selected_nodes_coords( tNodesInElem );

                    // Centroid is defined as the average of the coordinates of all the nodes

                    // Sum all coordinates components divided by number of vertices in polygon
                    for ( uint iNode = 0; iNode < tElemTopoNodes; ++iNode )
                    {
                        tCentroid( 0 ) += tNodeCoords( iNode, 0 ) / tElemTopoNodes;
                        tCentroid( 1 ) += tNodeCoords( iNode, 1 ) / tElemTopoNodes;
                        tCentroid( 2 ) += tNodeCoords( iNode, 2 ) / tElemTopoNodes;
                    }

                    // Collect centroid of current element shared
                    tAllCentroids (iElem ) = tCentroid;
                    tCentroid.fill( 0.0 );
                }

                if ( printCentroidFlag )
                {
                    std::cout<< "Average of centroids of shared elements for face "<< tFacesInSideSet( iFace ) << " is ("
                            << 0.5 * ( tAllCentroids( 0 )( 0 ) + tAllCentroids( 1 )( 0 ) ) << ", "
                            << 0.5 * ( tAllCentroids( 0 )( 1 ) + tAllCentroids( 1 )( 1 ) ) << ", "
                            << 0.5 * ( tAllCentroids( 0 )( 2 ) + tAllCentroids( 1 )( 2 ) ) << ")" << std::endl;
                }

                // Check average of centroids for elements adjacent to face with ID = 42
                if ( tFacesInSideSet( iFace )== 42 )
                {
                    REQUIRE(moris::equal_to(0.5 * ( tAllCentroids( 0 )( 0 ) + tAllCentroids( 1 )( 0 ) ),0.75));
                    REQUIRE(moris::equal_to(0.5 * ( tAllCentroids( 0 )( 1 ) + tAllCentroids( 1 )( 1 ) ),0.75));
                    REQUIRE(moris::equal_to(0.5 * ( tAllCentroids( 0 )( 2 ) + tAllCentroids( 1 )( 2 ) ),0.25));
                }
            }
        }
    }

    SECTION(" Verify appropriate reading of higher-order meshes " )
    {
        if ( p_rank == 0 && p_size == 1 )
        {
            // In this section higher order elements are tested. Basic information for 2D and 3D matrices are requested and
            // contrasted against expected topology and ids.

            //////////////
            //  QUAD_9  //
            //////////////

            // String for MTK to generate mesh with sideset
            const std::string meshName = "/test/src/mesh/MeshFiles/Quad9_4Elem.g";

            // Create MORIS mesh using MTK database
            moris::mesh meshHigherOrderElems0( MeshType::MTK, meshName );

            // Access basic mesh information
            Mat< uint > tNodesOwned = meshHigherOrderElems0.get_entities_owned_current_proc( EntityRank::NODE );
            Mat< real > tNodeCoords = meshHigherOrderElems0.get_selected_nodes_coords( tNodesOwned );
            uint tNumElems          = meshHigherOrderElems0.get_num_elems();
            uint tNumNodes          = meshHigherOrderElems0.get_num_nodes();

            // Testing random information in entities
            REQUIRE(moris::equal_to(tNumElems,4));
            REQUIRE(moris::equal_to(tNumNodes,25));
            REQUIRE(moris::equal_to(tNodeCoords(2,0),0.025));
            REQUIRE(moris::equal_to(tNodeCoords(13,1),0.3));

            // -----------------------------------------------------------------------

            //////////////
            //  TET_10  //
            //////////////

            // String for MTK to generate mesh with sideset
            const std::string meshName1 = "/test/src/mesh/MeshFiles/Tet10_12Elem.g";

            // Create MORIS mesh using MTK database
            moris::mesh meshHigherOrderElems1( MeshType::MTK, meshName1 );

            // Access basic mesh information
            tNodesOwned = meshHigherOrderElems1.get_entities_owned_current_proc( EntityRank::NODE );
            tNodeCoords = meshHigherOrderElems1.get_selected_nodes_coords( tNodesOwned );
            tNumElems   = meshHigherOrderElems1.get_num_elems();
            tNumNodes   = meshHigherOrderElems1.get_num_nodes();

            // Testing random information in entities
            REQUIRE(moris::equal_to(tNumElems,12));
            REQUIRE(moris::equal_to(tNumNodes,35));
            REQUIRE(moris::equal_to(tNodeCoords(0,0),0.625));
            REQUIRE(moris::equal_to(tNodeCoords(10,1),0.6));
            REQUIRE(moris::equal_to(tNodeCoords(17,2),-0.675));

            // -----------------------------------------------------------------------

            //////////////
            //  HEX_20  //
            //////////////

            // String for MTK to generate mesh with sideset
            const std::string meshName2 = "/test/src/mesh/MeshFiles/Hex20_2Elem.g";

            // Create MORIS mesh using MTK database
            moris::mesh meshHigherOrderElems2( MeshType::MTK, meshName2 );

            // Access basic mesh information
            tNodesOwned = meshHigherOrderElems2.get_entities_owned_current_proc( EntityRank::NODE );
            tNodeCoords = meshHigherOrderElems2.get_selected_nodes_coords( tNodesOwned );
            tNumElems   = meshHigherOrderElems2.get_num_elems();
            tNumNodes   = meshHigherOrderElems2.get_num_nodes();

            // Testing random information in entities
            REQUIRE(moris::equal_to(tNumElems,2));
            REQUIRE(moris::equal_to(tNumNodes,32));
            REQUIRE(moris::equal_to(tNodeCoords(0,0),1.625));
            REQUIRE(moris::equal_to(tNodeCoords(10,1),0.8));
            REQUIRE(moris::equal_to(tNodeCoords(17,2),0.875));

            // -----------------------------------------------------------------------

            //////////////
            //  HEX_27  //
            //////////////

            // String for MTK to generate mesh with sideset
            const std::string meshName3 = "/test/src/mesh/MeshFiles/Hex27_2Elem.g";

            // Create MORIS mesh using MTK database
            moris::mesh meshHigherOrderElems3( MeshType::MTK, meshName3 );

            // Access basic mesh information
            tNodesOwned = meshHigherOrderElems3.get_entities_owned_current_proc( EntityRank::NODE );
            tNodeCoords = meshHigherOrderElems3.get_selected_nodes_coords( tNodesOwned );
            tNumElems   = meshHigherOrderElems3.get_num_elems();
            tNumNodes   = meshHigherOrderElems3.get_num_nodes();

            // Testing random information in entities
            REQUIRE(moris::equal_to(tNumElems,2));
            REQUIRE(moris::equal_to(tNumNodes,45));
            REQUIRE(moris::equal_to(tNodeCoords(0,0),1.625));
            REQUIRE(moris::equal_to(tNodeCoords(39,1),0.8));
            REQUIRE(moris::equal_to(tNodeCoords(8,2),-0.875));

            // -----------------------------------------------------------------------
        }
    }

    SECTION( "Writing Mesh to an Exodus File")
    {
        // In this test local to global maps are accessed. The goal is to check the member variables
        // that were created in mtk to emulate the data stored in stk. Using this alternative can be
        // translated into additional memory requirements, but it gives the user access to the entirety
        // of the communication maps without asking on the fly to stk.
        if( p_size == 2)
        {
            const std::string fileName = "generated:1x1x2";

            // Create MORIS mesh using MTK database
            moris::mesh tSimpleMesh( MeshType::MTK, fileName );

            // Access maps (all of the local to global maps for all entity types).
            Mat < uint > lcl2GblNodeMap = tSimpleMesh.get_nodal_local_map();
            Mat < uint > lcl2GblElemMap = tSimpleMesh.get_elemental_local_map();
            Mat < uint > lcl2GblEdgeMap = tSimpleMesh.get_edge_local_map();
            Mat < uint > lcl2GblFaceMap = tSimpleMesh.get_face_local_map();

            Mat < uint > lcl2GblNodeOwnerProc = tSimpleMesh.get_nodal_owner_proc_map();
            Mat < uint > lcl2GblElemOwnerProc = tSimpleMesh.get_elemental_owner_proc_map();
            Mat < uint > lcl2GblEdgeOwnerProc = tSimpleMesh.get_edge_owner_proc_map();
            Mat < uint > lcl2GblFaceOwnerProc = tSimpleMesh.get_face_owner_proc_map();

            // Testing parts of maps
            if ( p_rank == 0 )
            {
                // Nodal maps
                REQUIRE(moris::equal_to(lcl2GblNodeMap(3),4));
                REQUIRE(moris::equal_to(lcl2GblNodeMap(5),6));
                REQUIRE(moris::equal_to(lcl2GblNodeOwnerProc(7),0));
                REQUIRE(moris::equal_to(lcl2GblNodeOwnerProc(8),1));

                // Elemental maps
                REQUIRE(moris::equal_to(lcl2GblElemMap(0),1));
                REQUIRE(moris::equal_to(lcl2GblElemOwnerProc(0),0));
                REQUIRE(moris::equal_to(lcl2GblElemOwnerProc(1),1));

                // Edge maps
                REQUIRE(moris::equal_to(lcl2GblEdgeMap(12),29));
                REQUIRE(moris::equal_to(lcl2GblEdgeOwnerProc(10),0));
                REQUIRE(moris::equal_to(lcl2GblEdgeOwnerProc(15),1));

                // Face maps
                REQUIRE(moris::equal_to(lcl2GblFaceMap(7),14));
                REQUIRE(moris::equal_to(lcl2GblFaceOwnerProc(3),0));
                REQUIRE(moris::equal_to(lcl2GblFaceOwnerProc(9),1));
            }
            else if ( p_rank == 1 )
            {
                // Nodal map
                REQUIRE(moris::equal_to(lcl2GblNodeMap(7),8));
                REQUIRE(moris::equal_to(lcl2GblNodeMap(9),10));
                REQUIRE(moris::equal_to(lcl2GblNodeOwnerProc(2),0));
                REQUIRE(moris::equal_to(lcl2GblNodeOwnerProc(10),1));

                // Elemental map
                REQUIRE(moris::equal_to(lcl2GblElemMap(1),2));
                REQUIRE(moris::equal_to(lcl2GblElemOwnerProc(0),0));
                REQUIRE(moris::equal_to(lcl2GblElemOwnerProc(1),1));

                // Edge map
                REQUIRE(moris::equal_to(lcl2GblEdgeMap(18),35));
                REQUIRE(moris::equal_to(lcl2GblEdgeOwnerProc(10),0));
                REQUIRE(moris::equal_to(lcl2GblEdgeOwnerProc(15),1));

                // Face map
                REQUIRE(moris::equal_to(lcl2GblFaceMap(10),18));
                REQUIRE(moris::equal_to(lcl2GblFaceOwnerProc(3),0));
                REQUIRE(moris::equal_to(lcl2GblFaceOwnerProc(9),1));
            }
        }
    }

//    SECTION( "Reading Mesh and Check for Duplicate Coordinates and IDs")
//    {
//        const std::string fileName = "/test/src/mesh/MeshFiles/2cubes_coarse.g";
//
//        if(p_size == 1) // Check coordinates and IDs in serial
//        {
//            Mat< uint > id(54,1,UINT_MAX);
//            for (moris::uint  i = 0; i<54;i++) id(i) = i;   // dummy id-field for the nodes
//            id(28)=2; // Modify ID field, that one coordinate and id have a duplicate
//            id(27)=3; // Modify ID field, that one coordinate and id have a duplicate
//
//            moris::mesh Mesh( MeshType::MTK, fileName );
//            Mat< real > tCoord = Mesh.get_all_nodes_coords();
//
//            Mat< uint > duplicate_list_check1 = Mesh.duplicate_node_coord_check();
//            Mat< uint > duplicate_list_check2 = Mesh.duplicate_node_coord_and_id_check(tCoord,id);
//        }
//        if(p_size == 2) // Check coordinates and IDs in parallel
//        {
//            moris::mesh Mesh( MeshType::MTK, fileName );
//            Mat< real > tCoord = Mesh.get_all_nodes_coords();
//
//            Mat< uint > id(27,1,UINT_MAX);
//            for ( uint  i = 0; i<27;i++) id(i) = i;
//
//            if(p_size > 0)
//            {
//                moris::MPITools tWorld;
//                if(p_rank == 0)
//                {
//                    Cell< Mat< uint > > tGatheredId = tWorld.gather(id);
//                    Cell< Mat< real > > tGatheredCoord = tWorld.gather(tCoord);
//
//                    Mat< uint > duplicate_list = Mesh.duplicate_node_coord_and_id_check(tGatheredCoord,tGatheredId);
//                }
//                if(p_rank > 0)
//                {
//                    id(1)=3; // Modify ID field, that one coordinate and id have a duplicate
//                    Cell< Mat< uint > > tGatheredId = tWorld.gather(id);
//                    Cell< Mat< real > > tGatheredCoord = tWorld.gather(tCoord);
//                }
//            }
//        }
//        if(p_size == 3) // Check coordinates and IDs in parallel
//        {
//            const std::string fileNameP3 = "/test/src/mesh/MeshFiles/3cubes_coarse.g";
//            moris::mesh Mesh( MeshType::MTK, fileNameP3 );
//            Mat< real > tCoord = Mesh.get_all_nodes_coords();
//            Mat< uint >  id(27,1,UINT_MAX);
//
//            for ( uint  i = 0; i<27;i++) id(i) = i;
//
//            if(p_size > 0)
//            {
//                moris::MPITools tWorld;
//                if(p_rank == 0)
//                {
//                    Cell< Mat< uint > > tGatheredId = tWorld.gather(id);
//                    Cell< Mat< real > > tGatheredCoord = tWorld.gather(tCoord);
//
//                    Mat< uint > duplicate_list = Mesh.duplicate_node_coord_and_id_check(tGatheredCoord,tGatheredId);
//                }
//                if(p_rank > 0)
//                {
//                    id(5)=0; // Modify ID field, that one coordinate and id have a duplicate
//                    Cell< Mat< uint > > tGatheredId = tWorld.gather(id);
//                    Cell< Mat< real > > tGatheredCoord = tWorld.gather(tCoord);
//                }
//            }
//        }
//    }

//    SECTION( "Reading Mesh and Check for Duplicate Coordinates and IDs and Give only Duplicates with Problems")
//    {
//        if(p_size == 1)
//        {
//            Mat< uint >  id(54,1,UINT_MAX);
//            for ( uint  i = 0; i<54;i++) id(i) = i;   // dummy id-field for the nodes
//            id(28)=2; // Modify ID field, that one coordinate and id have a duplicate
//            id(27)=3; // Modify ID field, that one coordinate and id have a duplicate
//            const std::string fileName = "/test/src/mesh/MeshFiles/2cubes_coarse.g";
//
//            moris::mesh Mesh( MeshType::MTK, fileName );
//
//            Mat< real > tCoord = Mesh.get_all_nodes_coords();
//            Mat< uint > problem_list = Mesh.duplicate_node_coord_and_id_check_problems(tCoord,id);
//        }
//
//        if( p_size == 3 )
//        {
//            const std::string fileName = "/test/src/mesh/MeshFiles/3cubes_coarse.g";    // 8 elements, 27 nodes
//            moris::mesh Mesh( MeshType::MTK, fileName );
//            Mat< real > tCoord = Mesh.get_all_nodes_coords();
//            Mat< uint > id(27,1,UINT_MAX);
//
//            for ( uint  i = 0; i<27;i++) id(i) = i;
//
//            if(p_size >1)
//            {
//                moris::MPITools tWorld;
//                if(p_rank == 0)
//                {
//                    Cell<Mat< uint >> tGatheredId = tWorld.gather(id);
//                    Cell<Mat< real >> tGatheredCoord = tWorld.gather(tCoord);
//
//                    Mat< uint > problem_list = Mesh.duplicate_node_coord_and_id_check_problems(tGatheredCoord,tGatheredId);
//                }
//
//                if(p_rank > 0)
//                {
//                    id(5)=0; // Modify ID field, that one coordinate and id have a duplicate
//                    Cell< Mat< uint > > tGatheredId = tWorld.gather(id);
//                    Cell< Mat< real > > tGatheredCoord = tWorld.gather(tCoord);
//                }
//            }
//        }
//    }
}

