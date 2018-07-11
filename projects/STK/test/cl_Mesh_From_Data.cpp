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

TEST_CASE("MtkMeshFromData", "[moris],[mesh],[cl_Mesh],[Mesh]")
{
    // Parallel
    uint p_rank = 0;
    uint p_size = 1;
#ifdef MORIS_HAVE_PARALLEL
    p_rank = moris::par_rank();
    p_size = moris::par_size();
#endif

    SECTION( "Creating a 2D mesh from data in serial" )
    {
        if(p_rank == 0 && p_size == 1 ) // specify it is a serial test only
        {
            // Generate data for test
            uint aNumDim = 2;
            Mat< real > aCoords(6,2);
            aCoords(0,0) = 0.0, aCoords(0,1) = 0.0;
            aCoords(1,0) = 1.0, aCoords(1,1) = 0.0;
            aCoords(2,0) = 1.0, aCoords(2,1) = 1.0;
            aCoords(3,0) = 0.0, aCoords(3,1) = 1.0;
            aCoords(4,0) = 2.0, aCoords(4,1) = 0.0;
            aCoords(5,0) = 2.0, aCoords(5,1) = 1.0;
            Mat< uint > aElemConn( 2, 4 );

            SECTION( "using consecutive node and element ids" )
            {
                // 0D to 3D connectivity (node to element)
                aElemConn( 0, 0 ) = 1; aElemConn( 0, 1 ) = 2; aElemConn( 0, 2 ) = 3; aElemConn( 0, 3 ) = 4;
                aElemConn( 1, 0 ) = 2; aElemConn( 1, 1 ) = 5; aElemConn( 1, 2 ) = 6; aElemConn( 1, 3 ) = 3;

                // No need of an element map since elements in connectivity table are assumed to be contiguous

                // Create MORIS mesh using MTK database
                moris::MtkMeshData aMeshData;
                aMeshData.SpatialDim = &aNumDim;
                aMeshData.ElemConn   = &aElemConn;
                aMeshData.NodeCoords = &aCoords;
                moris::mesh tMesh2D_QUADs( MeshType::MTK, aMeshData );

                // =============================
                // Testing basic functionalities
                // =============================

                uint tNumElements = tMesh2D_QUADs.get_num_elems();
                REQUIRE(moris::equal_to(tNumElements,2));

                uint tNumNodes    = tMesh2D_QUADs.get_num_nodes();
                REQUIRE(moris::equal_to(tNumNodes,6));

                Mat< uint > tAvailableNodeIDs = tMesh2D_QUADs.generate_unique_node_ids(4);
                REQUIRE(moris::equal_to(tAvailableNodeIDs(0,0),7));
                REQUIRE(moris::equal_to(tAvailableNodeIDs(1,0),8));
                REQUIRE(moris::equal_to(tAvailableNodeIDs(2,0),9));
                REQUIRE(moris::equal_to(tAvailableNodeIDs(3,0),10));
            }

            SECTION( "using non-consecutive node and element ids" )
            {
                // 0D to 3D connectivity (node to element)
                aElemConn( 0, 0 ) = 10; aElemConn( 0, 1 ) = 2; aElemConn( 0, 2 ) = 23; aElemConn( 0, 3 ) = 4;
                aElemConn( 1, 0 ) = 2; aElemConn( 1, 1 ) = 35; aElemConn( 1, 2 ) = 6; aElemConn( 1, 3 ) = 23;

                // Local to global maps to tell the mesh the corresponding non-consecutive entity ids
                Mat< uint > aNodeLocaltoGlobal = {{10},{2},{23},{4},{35},{6}};
                Mat< uint > aElemLocaltoGlobal = {{31},{18}};

                // Create MORIS mesh using MTK database
                moris::MtkMeshData aMeshData;
                aMeshData.SpatialDim = &aNumDim;
                aMeshData.ElemConn   = &aElemConn;
                aMeshData.NodeCoords = &aCoords;
                aMeshData.LocaltoGlobalElemMap = &aElemLocaltoGlobal; // nodal map for non-consecutive ids
                aMeshData.LocaltoGlobalNodeMap = &aNodeLocaltoGlobal; // elemental map for non-consecutive ids

                moris::mesh tMesh2D_QUADs( MeshType::MTK, aMeshData );

                uint tNumElems = tMesh2D_QUADs.get_num_elems();
                uint tNumNodes = tMesh2D_QUADs.get_num_nodes();

                // =============================
                // Testing basic functionalities
                // =============================

                Mat< uint > tNodeIDs = tMesh2D_QUADs.get_entities_owned_current_proc(EntityRank::NODE);
                Mat< uint > tElemIDs = tMesh2D_QUADs.get_entities_owned_current_proc(EntityRank::ELEMENT);

                REQUIRE(moris::equal_to(tNodeIDs(0,0),2));
                REQUIRE(moris::equal_to(tNodeIDs(1,0),4));
                REQUIRE(moris::equal_to(tNodeIDs(2,0),6));
                REQUIRE(moris::equal_to(tNodeIDs(3,0),10));
                REQUIRE(moris::equal_to(tNodeIDs(4,0),23));
                REQUIRE(moris::equal_to(tNodeIDs(5,0),35));

                REQUIRE(moris::equal_to(tElemIDs(0,0),18));
                REQUIRE(moris::equal_to(tElemIDs(1,0),31));

                // Entity keys
                uint tTestId = tMesh2D_QUADs.get_entity_key_from_entity_id(EntityRank::NODE,23);
                uint tTestKey = tMesh2D_QUADs.get_entity_id_from_entity_key(EntityRank::NODE,tTestId);
                REQUIRE(moris::equal_to(tTestKey,23));
            }
        }
    }

    SECTION( "Creating a 3D 2 element mesh from data in serial using non-consecutive node and element id maps")
    {
        if(p_rank == 0 && p_size == 1 ) // specify it is a serial test only
        {
            uint aNumDim = 3;
            Mat< uint > aElemConn = {{1000,2,4,38,543,6,8,77},{543,6,8,77,93,10,12,111}};
            Mat< uint > aNodeLocaltoGlobalNC = {{1000},{2},{38},{4},{543},{6},{77},{8},{93},{10},{111},{12}};
            Mat< real > aCoords   = {{0.0, 0.0, 0.0},
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
            Mat< uint > aNodeLocaltoGlobalnc = {{1000},{2},{38},{4},{543},{6},{77},{8},{93},{10},{111},{12}};
            Mat< uint > aElemLocaltoGlobalNC = {{51},{32}};

            // Create MORIS mesh using MTK database
            moris::MtkMeshData aMeshData;
            aMeshData.SpatialDim = &aNumDim;
            aMeshData.ElemConn   = &aElemConn;
            aMeshData.NodeCoords = &aCoords;
            aMeshData.LocaltoGlobalElemMap = &aElemLocaltoGlobalNC;
            aMeshData.LocaltoGlobalNodeMap = &aNodeLocaltoGlobalNC;

            moris::mesh tMesh3D_HEXs( MeshType::MTK, aMeshData );

            // =============================
            // Testing basic functionalities
            // =============================
            uint tNumElems = tMesh3D_HEXs.get_num_elems();

            Mat< uint > tNodeIDs = tMesh3D_HEXs.get_entities_owned_current_proc(EntityRank::NODE);
            Mat< uint > tElemIDs = tMesh3D_HEXs.get_entities_owned_current_proc(EntityRank::ELEMENT);

            REQUIRE(moris::equal_to(tNumElems,2));
            REQUIRE(moris::equal_to(tNodeIDs(0),2));
            REQUIRE(moris::equal_to(tNodeIDs(1),4));
            REQUIRE(moris::equal_to(tNodeIDs(2),6));
            REQUIRE(moris::equal_to(tNodeIDs(3),8));
            REQUIRE(moris::equal_to(tNodeIDs(4),10));
            REQUIRE(moris::equal_to(tNodeIDs(5),12));
            REQUIRE(moris::equal_to(tNodeIDs(6),38));
            REQUIRE(moris::equal_to(tNodeIDs(7),77));
            REQUIRE(moris::equal_to(tNodeIDs(8),93));
            REQUIRE(moris::equal_to(tNodeIDs(9),111));
            REQUIRE(moris::equal_to(tNodeIDs(10),543));
            REQUIRE(moris::equal_to(tNodeIDs(11),1000));

            REQUIRE(moris::equal_to(tElemIDs(0),32));
            REQUIRE(moris::equal_to(tElemIDs(1),51));
        }
    }

    SECTION( "Creating a 3D 2 element mesh from data in serial ")
    {
        if(p_rank == 0 && p_size == 1 ) // specify it is a serial test only
        {
            // Generate data for test
            uint aNumDim = 3;
            Mat< uint > aElemConn = {{1,2,4,3,5,6,8,7},{5,6,8,7,9,10,12,11}};
            Mat< real > aCoords   = {{0.0, 0.0, 0.0},
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
            Mat< uint > aNodeLocaltoGlobal = {{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12}};
            Mat< uint > aElemLocaltoGlobal = {{1},{2}};

            SECTION("with multiple scalar/vector nodal/elemental fields")
            {
                // Create nodal and elemental fields of reals
                Mat< real > tNodalFieldData_1 = {{10.0},{10.0},{10.0},{10.0},
                        {20.0},{20.0},{20.0},{20.0},
                        {30.0},{35.0},{35.0},{30.0}};
                Mat< real > tNodalFieldData_2 = {{-1.2},{-2.3},{-3.4},{-4.5},
                        {-5.6},{-6.7},{-7.8},{-8.9},
                        {-1.2},{-2.3},{-3.4},{-4.5}};
                Mat< real > aElemFieldData_1 = {{111.11},{222.22}};
                Mat< real > aElemFieldData_2 = {{0.253, -0.253},{0.976,0.0}};
                Mat< real > aElemFieldData_3 = {{8.2, 1.1, 1.1, 8.2},{-3.45, 2.2, 2.2, -3.45}};

                Cell < Mat< real > > aFieldData      = { tNodalFieldData_1, tNodalFieldData_2, aElemFieldData_1, aElemFieldData_2, aElemFieldData_3 };
                Cell < std::string > aFieldName      = { "nodalField_1", "nodalField_2", "elementField_1", "elementField_2" , "elementField_3" };
                Cell < enum EntityRank > aFieldRanks = { EntityRank::NODE, EntityRank::NODE, EntityRank::ELEMENT, EntityRank::ELEMENT, EntityRank::ELEMENT };

                // Create MORIS mesh using MTK database
                moris::MtkFieldsInfo aFieldsInfo;
                aFieldsInfo.FieldsData = &aFieldData;
                aFieldsInfo.FieldsName = aFieldName;
                aFieldsInfo.FieldsRank = aFieldRanks;

                moris::MtkMeshData aMeshData;
                aMeshData.SpatialDim = &aNumDim;
                aMeshData.ElemConn   = &aElemConn;
                aMeshData.NodeCoords = &aCoords;
                aMeshData.FieldsInfo  = &aFieldsInfo;
                aMeshData.CreateAllEdgesAndFaces = true;
                aMeshData.LocaltoGlobalElemMap = &aElemLocaltoGlobal;
                aMeshData.LocaltoGlobalNodeMap = &aNodeLocaltoGlobal;

                moris::mesh tMesh( MeshType::MTK, aMeshData );

                // Output mesh
                std::string OutputFileName1 = "meshMultFields.exo";  //output name
                tMesh.create_output_mesh(OutputFileName1);

                // ========================
                // Testing outputting field
                // ========================
                Mat< real > tNodalField_1 = tMesh.get_field_values( EntityRank::NODE, "nodalField_1" );
                Mat< real > tNodalField_2 = tMesh.get_field_values( EntityRank::NODE, "nodalField_2" );
                Mat< real > tElemField_1  = tMesh.get_field_values( EntityRank::ELEMENT, "elementField_1" );
                Mat< real > tElemField_2  = tMesh.get_field_values( EntityRank::ELEMENT, "elementField_2" );
                Mat< real > tElemField_3  = tMesh.get_field_values( EntityRank::ELEMENT, "elementField_3" );

                // Check fields values
                REQUIRE(moris::equal_to(tNodalField_1(0,0),10.0));
                REQUIRE(moris::equal_to(tNodalField_1(1,0),10.0));
                REQUIRE(moris::equal_to(tNodalField_1(2,0),10.0));
                REQUIRE(moris::equal_to(tNodalField_1(3,0),10.0));
                REQUIRE(moris::equal_to(tNodalField_1(4,0),20.0));
                REQUIRE(moris::equal_to(tNodalField_1(5,0),20.0));
                REQUIRE(moris::equal_to(tNodalField_1(6,0),20.0));
                REQUIRE(moris::equal_to(tNodalField_1(7,0),20.0));
                REQUIRE(moris::equal_to(tNodalField_1(8,0),30.0));
                REQUIRE(moris::equal_to(tNodalField_1(9,0),35.0));
                REQUIRE(moris::equal_to(tNodalField_1(10,0),35.0));
                REQUIRE(moris::equal_to(tNodalField_1(11,0),30.0));

                REQUIRE(moris::equal_to(tNodalField_2(0,0),-1.2));
                REQUIRE(moris::equal_to(tNodalField_2(1,0),-2.3));
                REQUIRE(moris::equal_to(tNodalField_2(2,0),-3.4));
                REQUIRE(moris::equal_to(tNodalField_2(3,0),-4.5));
                REQUIRE(moris::equal_to(tNodalField_2(4,0),-5.6));
                REQUIRE(moris::equal_to(tNodalField_2(5,0),-6.7));
                REQUIRE(moris::equal_to(tNodalField_2(6,0),-7.8));
                REQUIRE(moris::equal_to(tNodalField_2(7,0),-8.9));
                REQUIRE(moris::equal_to(tNodalField_2(8,0),-1.2));
                REQUIRE(moris::equal_to(tNodalField_2(9,0),-2.3));
                REQUIRE(moris::equal_to(tNodalField_2(10,0),-3.4));
                REQUIRE(moris::equal_to(tNodalField_2(11,0),-4.5));

                REQUIRE(moris::equal_to(tElemField_1(0,0),111.11));
                REQUIRE(moris::equal_to(tElemField_1(1,0),222.22));

                REQUIRE(moris::equal_to(tElemField_2(0,0),0.253));
                REQUIRE(moris::equal_to(tElemField_2(0,1),-0.253));
                REQUIRE(moris::equal_to(tElemField_2(1,0),0.976));
                REQUIRE(moris::equal_to(tElemField_2(1,1),0.0));

                REQUIRE(moris::equal_to(tElemField_3(0,0),8.2));
                REQUIRE(moris::equal_to(tElemField_3(0,1),1.1));
                REQUIRE(moris::equal_to(tElemField_3(0,2),1.1));
                REQUIRE(moris::equal_to(tElemField_3(0,3),8.2));
                REQUIRE(moris::equal_to(tElemField_3(1,0),-3.45));
                REQUIRE(moris::equal_to(tElemField_3(1,1),2.2));
                REQUIRE(moris::equal_to(tElemField_3(1,2),2.2));
                REQUIRE(moris::equal_to(tElemField_3(1,3),-3.45));

                // ===================================================
                // Testing entities connected to element with ID = 7
                // ===================================================
                uint elementID = 2;

                Mat< uint > elemsConnectedToElement = tMesh.get_elements_connected_to_element(elementID);
                uint NumberOfElemsConnectedToElement = elemsConnectedToElement.length();

                REQUIRE(moris::equal_to(NumberOfElemsConnectedToElement,1));
                REQUIRE(moris::equal_to(elemsConnectedToElement(0),1));
            }

            SECTION( "with 2 block sets, 1 node set, and 1 side set" )
            {
                // Create 2  block sets (one over each element) and a node set that contains only 4 nodes
                // NOTE: A side set requires a two column matrix. The first column contains the ids of the elements to which
                // the faces are associated (for reference), and the second column represents the ordinal of the face,
                // which should go from 0 to n-1, being n the number of faces. numbering goes counterclockwise.

                // Declare block sets
                //////////////////////
                Mat<uint> tBlockSetsPartOwners = { {0},{1} };

                MtkBlockSetsInfo tBlockSetStruc;
                tBlockSetStruc.BSetInds = &tBlockSetsPartOwners;
                tBlockSetStruc.BSetNames   = { "blockset_1", "blockset_2" };

                // Declare side sets
                /////////////////////
                Mat<uint> tSideset_1  = { {1, 3}, {1, 4}, {1, 5}, {2, 1}, {2, 2} };
                Cell< Mat< uint > > tSideSetsInfo = { tSideset_1 };

                MtkSideSetsInfo tSideSetStruc;
                tSideSetStruc.ElemIdsAndSideOrds = &tSideSetsInfo;
                tSideSetStruc.SSetNames   = { "Sideset_1" };

                // Declare node sets
                /////////////////////
                Mat<uint> tNodeSet_1  = { {1}, {3}, {5}, {6} };
                Cell< Mat< uint > > tNodeSetsEntIds = { tNodeSet_1 };

                MtkNodeSetsInfo tNodeSetStruc;
                tNodeSetStruc.EntIds = &tNodeSetsEntIds;
                tNodeSetStruc.NSetNames = { "Nodeset_1" };

                // Create MORIS mesh using MTK database
                ///////////////////////////////////////
                MtkSetsInfo aMeshSets;
                aMeshSets.NodeSetsInfo = &tNodeSetStruc;
                aMeshSets.SideSetsInfo   = &tSideSetStruc;
                aMeshSets.BlockSetsInfo   = &tBlockSetStruc;

                moris::MtkMeshData aMeshData;
                aMeshData.SpatialDim = &aNumDim;
                aMeshData.ElemConn   = &aElemConn;
                aMeshData.NodeCoords = &aCoords;
                aMeshData.SetsInfo   = &aMeshSets;
                aMeshData.LocaltoGlobalElemMap   = &aElemLocaltoGlobal;
                aMeshData.LocaltoGlobalNodeMap   = &aNodeLocaltoGlobal;

                moris::mesh tMesh( MeshType::MTK, aMeshData );

                // ========================
                // Testing sets information
                // ========================

                Mat< uint > tBlock1 = tMesh.get_set_entity_ids( EntityRank::ELEMENT, "blockset_1" );
                Mat< uint > tBlock2 = tMesh.get_set_entity_ids( EntityRank::ELEMENT, "blockset_2" );
                Mat< uint > tNodeSet1 = tMesh.get_set_entity_ids( EntityRank::NODE, "Nodeset_1" );
                Mat< uint > tSideSet1 = tMesh.get_set_entity_ids( EntityRank::FACE, "Sideset_1" );

                Mat< uint > tNodesInBlockSet1 = tMesh.get_nodes_in_block_set( 1 );
                Mat< uint > tNodesInBlockSet2 = tMesh.get_nodes_in_block_set( 2 );
                Mat< uint > tNodesInNodeSet1 = tMesh.get_nodes_in_node_set( 1 );
                Mat< uint > tFacesInFaceSet1 = tMesh.get_faces_in_side_set( 1 );

                Mat< uint > tFacesInFaceSet1AndBlockSet1 = tMesh.get_intersected_entities_field_set(EntityRank::FACE, "Sideset_1", "blockset_1");
                Mat< uint > tFacesInFaceSet1AndBlockSet2 = tMesh.get_intersected_entities_field_set(EntityRank::FACE, "Sideset_1", "blockset_2");
                Mat< uint > tNodesInFaceSet1AndBlockSet2 = tMesh.get_intersected_entities_field_set(EntityRank::NODE, "Sideset_1", "blockset_2");

                uint tNumFaces = tMesh.get_num_faces();
                uint tNumEdges = tMesh.get_num_edges();

                REQUIRE(moris::equal_to(tBlock1(0,0),1));
                REQUIRE(moris::equal_to(tBlock2(0,0),2));

                REQUIRE(moris::equal_to(tNodeSet1(0,0),1));
                REQUIRE(moris::equal_to(tNodeSet1(1,0),3));
                REQUIRE(moris::equal_to(tNodeSet1(2,0),5));
                REQUIRE(moris::equal_to(tNodeSet1(3,0),6));

                REQUIRE(moris::equal_to(tNodesInNodeSet1(0,0),1));
                REQUIRE(moris::equal_to(tNodesInNodeSet1(1,0),3));
                REQUIRE(moris::equal_to(tNodesInNodeSet1(2,0),5));
                REQUIRE(moris::equal_to(tNodesInNodeSet1(3,0),6));

                REQUIRE(moris::equal_to(tNodesInBlockSet1(0,0),1));
                REQUIRE(moris::equal_to(tNodesInBlockSet1(1,0),2));
                REQUIRE(moris::equal_to(tNodesInBlockSet1(2,0),3));
                REQUIRE(moris::equal_to(tNodesInBlockSet1(3,0),4));
                REQUIRE(moris::equal_to(tNodesInBlockSet1(4,0),5));
                REQUIRE(moris::equal_to(tNodesInBlockSet1(5,0),6));
                REQUIRE(moris::equal_to(tNodesInBlockSet1(6,0),7));
                REQUIRE(moris::equal_to(tNodesInBlockSet1(7,0),8));

                REQUIRE(moris::equal_to(tNodesInBlockSet2(0,0),5));
                REQUIRE(moris::equal_to(tNodesInBlockSet2(1,0),6));
                REQUIRE(moris::equal_to(tNodesInBlockSet2(2,0),7));
                REQUIRE(moris::equal_to(tNodesInBlockSet2(3,0),8));
                REQUIRE(moris::equal_to(tNodesInBlockSet2(4,0),9));
                REQUIRE(moris::equal_to(tNodesInBlockSet2(5,0),10));
                REQUIRE(moris::equal_to(tNodesInBlockSet2(6,0),11));
                REQUIRE(moris::equal_to(tNodesInBlockSet2(7,0),12));

                REQUIRE(moris::equal_to(tSideSet1(0),14));
                REQUIRE(moris::equal_to(tSideSet1(1),15));
                REQUIRE(moris::equal_to(tFacesInFaceSet1(2),16));
                REQUIRE(moris::equal_to(tFacesInFaceSet1(3),22));
                REQUIRE(moris::equal_to(tFacesInFaceSet1(4),23));

                REQUIRE(moris::equal_to(tNumFaces,5));
                REQUIRE(moris::equal_to(tNumEdges,0));

                REQUIRE(moris::equal_to(tFacesInFaceSet1AndBlockSet1(0),14));
                REQUIRE(moris::equal_to(tFacesInFaceSet1AndBlockSet1(1),15));
                REQUIRE(moris::equal_to(tFacesInFaceSet1AndBlockSet1(2),16));

                REQUIRE(moris::equal_to(tFacesInFaceSet1AndBlockSet2(0),22));
                REQUIRE(moris::equal_to(tFacesInFaceSet1AndBlockSet2(1),23));

                REQUIRE(moris::equal_to(tNodesInFaceSet1AndBlockSet2(0),5));
                REQUIRE(moris::equal_to(tNodesInFaceSet1AndBlockSet2(2),7));
                REQUIRE(moris::equal_to(tNodesInFaceSet1AndBlockSet2(3),8));
                REQUIRE(moris::equal_to(tNodesInFaceSet1AndBlockSet2(5),11));
            }

            SECTION( "with multiple sets and fields allocated to sets instead of the entire mesh")
            {
                // Create 2  block sets (one over each element) and a node set that contains only 4 nodes
                // NOTE: A side set requires a two column matrix. The first column contains the ids of the elements to which
                // the faces are associated (for reference), and the second column represents the ordinal of the face,
                // which should go from 0 to n-1, being n the number of faces. numbering goes counterclockwise.

                // Declare side sets
                /////////////////////
                Mat<uint> tSideset_1  = { {1, 5}, {2, 0} };
                Cell< Mat< uint > > tSideSetsInfo = { tSideset_1 };

                MtkSideSetsInfo tSideSetStruc;
                tSideSetStruc.ElemIdsAndSideOrds = &tSideSetsInfo;
                tSideSetStruc.SSetNames   = { "Sideset_1" };

                // Declare node sets
                /////////////////////
                Mat<uint> tNodeSet_1  = { {4}, {5}, {7} };
                Cell< Mat< uint > > tNodeSetsEntIds = { tNodeSet_1 };

                MtkNodeSetsInfo tNodeSetStruc;
                tNodeSetStruc.EntIds = &tNodeSetsEntIds;
                tNodeSetStruc.NSetNames = { "Nodeset_1" };

                // Declare fields
                /////////////////

                // Create nodal and elemental fields of reals
                Mat< real > tElemFieldData_1      = {{123.45},{678.9}};
                Mat< real > tNodalFieldData_1     = {{-1.2},{-2.3},{-3.4}};
                Mat< real > tSideFieldNodalData_1 = {{555.555}};
                Mat< real > tSideFieldFaceData_1  = {{111.11},{222.22}};

                Cell < Mat< real > > aFieldData      = { tElemFieldData_1, tNodalFieldData_1, tSideFieldNodalData_1, tSideFieldFaceData_1 };
                // Field set owner. If not provided, it will use the universal part.
                Cell < std::string > aSetOwnerName   = { "", "NodeSet_1", "Sideset_1", "Sideset_1" };
                Cell < std::string > aFieldName      = { "ElementField_1", "NodesetField_1", "SidesetField_1", "SidesetField_2" };
                Cell < enum EntityRank > aFieldRanks = { EntityRank::ELEMENT, EntityRank::NODE, EntityRank::NODE, EntityRank::FACE };

                // Create MORIS mesh using MTK database
                ///////////////////////////////////////
                MtkSetsInfo aMeshSets;
                aMeshSets.NodeSetsInfo = &tNodeSetStruc;
                aMeshSets.SideSetsInfo   = &tSideSetStruc;

                moris::MtkFieldsInfo aFieldsInfo;
                aFieldsInfo.FieldsData = &aFieldData;
                aFieldsInfo.FieldsName = aFieldName;
                aFieldsInfo.FieldsRank = aFieldRanks;
                aFieldsInfo.SetsOwner  = &aSetOwnerName;

                moris::MtkMeshData aMeshData;
                aMeshData.SpatialDim = &aNumDim;
                aMeshData.ElemConn   = &aElemConn;
                aMeshData.NodeCoords = &aCoords;
                aMeshData.SetsInfo   = &aMeshSets;
                aMeshData.FieldsInfo = &aFieldsInfo;
                aMeshData.LocaltoGlobalElemMap = &aElemLocaltoGlobal;
                aMeshData.LocaltoGlobalNodeMap = &aNodeLocaltoGlobal;

                moris::mesh tMesh( MeshType::MTK, aMeshData );

                // ========================
                // Testing outputting field
                // ========================

                Mat< real > tField_1  = tMesh.get_field_values( EntityRank::ELEMENT, "ElementField_1" );
                Mat< real > tField_2  = tMesh.get_field_values( EntityRank::NODE, "NodesetField_1" );
                Mat< real > tField_3  = tMesh.get_field_values( EntityRank::NODE, "SidesetField_1" );
                Mat< real > tField_4  = tMesh.get_field_values( EntityRank::FACE, "SidesetField_2" );

                REQUIRE(moris::equal_to(tField_1(0,0),123.45));
                REQUIRE(moris::equal_to(tField_1(1,0),678.9));

                REQUIRE(moris::equal_to(tField_2(0,0),-1.2));
                REQUIRE(moris::equal_to(tField_2(1,0),-2.3));
                REQUIRE(moris::equal_to(tField_2(2,0),-3.4));

                REQUIRE(moris::equal_to(tField_3(0,0),555.555));
                REQUIRE(moris::equal_to(tField_3(1,0),555.555));
                REQUIRE(moris::equal_to(tField_3(2,0),555.555));
                REQUIRE(moris::equal_to(tField_3(3,0),555.555));
                REQUIRE(moris::equal_to(tField_3(4,0),555.555));

                REQUIRE(moris::equal_to(tField_4(0,0),111.11));
                REQUIRE(moris::equal_to(tField_4(1,0),222.22));
            }
        }
    }

    SECTION( "Creating a 3D 3 element mesh in parallel using element processor owner list")
    {
        if( p_size == 2 ) // specify it is a 2 processor test
        {
            uint aNumDim = 3;
            Mat<uint >  aElemProcs;
            Mat<uint >  aElemConn;
            Mat<real >  aNodeCoords;
            Mat<uint>    aElemLocaltoGlobal;
            Mat<uint>    aNodeLocaltoGlobal;

            if ( p_rank == 0 )
            {
                // Generate data for test
                Mat< uint>  tElemProcsDummy  = {{0}, {0}, {1}};
                Mat< uint > tElemConnDummy   = {{1000,2,4,3000,5,6,8,7000},
                        {5,6,8,7000,9000,10,12,11000},
                        {9000,10,12,11000,13,14,16,15000}};
                Mat< real > tNodeCoordsDummy = {{0.0, 0.0, 0.0},
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
                        {1.0, 1.0, 2.0},
                        {0.0, 0.0, 3.0},
                        {1.0, 0.0, 3.0},
                        {0.0, 1.0, 3.0},
                        {1.0, 1.0, 3.0}};
                Mat< uint > tElemLocaltoGlobalDummy = {{15},{7},{31}};
                Mat< uint > tNodeLocaltoGlobalDummy = {{1000},{2},{3000},{4},{5},
                        {6},{7000},{8},{9000},{10},
                        {11000},{12},{13},{14},{15000},{16}};
                aElemProcs         = tElemProcsDummy;
                aElemConn          = tElemConnDummy;
                aNodeCoords        = tNodeCoordsDummy;
                aElemLocaltoGlobal = tElemLocaltoGlobalDummy;
                aNodeLocaltoGlobal = tNodeLocaltoGlobalDummy;

            }
            else
            {

                // Generate data for test
                moris::Mat< moris::uint>  tElemProcsDummy  = {{0}, {1}};
                moris::Mat< moris::uint > tElemConnDummy = {{5,6,8,7000,9000,10,12,11000},
                        {9000,10,12,11000,13,14,16,15000}};
                moris::Mat< moris::real > tNodeCoordsDummy = {{0.0, 0.0, 1.0},
                        {1.0, 0.0, 1.0},
                        {0.0, 1.0, 1.0},
                        {1.0, 1.0, 1.0},
                        {0.0, 0.0, 2.0},
                        {1.0, 0.0, 2.0},
                        {0.0, 1.0, 2.0},
                        {1.0, 1.0, 2.0},
                        {0.0, 0.0, 3.0},
                        {1.0, 0.0, 3.0},
                        {0.0, 1.0, 3.0},
                        {1.0, 1.0, 3.0}};
                moris::Mat< moris::uint > tElemLocaltoGlobalDummy = {{7},{31}};
                moris::Mat< moris::uint > tNodeLocaltoGlobalDummy = {{5},{6},{7000},{8},{9000},{10},
                        {11000},{12},{13},{14},{15000},{16}};

                aElemProcs         = tElemProcsDummy;
                aElemConn          = tElemConnDummy;
                aNodeCoords        = tNodeCoordsDummy;
                aElemLocaltoGlobal = tElemLocaltoGlobalDummy;
                aNodeLocaltoGlobal = tNodeLocaltoGlobalDummy;
            }

            SECTION( "mesh without sets or fields")
            {
                moris::MtkMeshData aMeshData;
                aMeshData.SpatialDim   = &aNumDim;
                aMeshData.ElemConn     = &aElemConn;
                aMeshData.NodeCoords   = &aNodeCoords;
                aMeshData.EntProcOwner = &aElemProcs;
                aMeshData.CreateAllEdgesAndFaces = true;
                aMeshData.LocaltoGlobalElemMap   = &aElemLocaltoGlobal;
                aMeshData.LocaltoGlobalNodeMap   = &aNodeLocaltoGlobal;

                moris::mesh tParallelMesh( MeshType::MTK, aMeshData );

                // ========================================
                // Testing get_entities_owned_current_proc
                // ========================================

                Mat< uint > tElemsOwned  = tParallelMesh.get_entities_owned_current_proc(EntityRank::ELEMENT);
                Mat< uint > tNodesOwned  = tParallelMesh.get_entities_owned_current_proc(EntityRank::NODE);
                Mat< uint > tProcsSharedByEntityId  = tParallelMesh.get_procs_sharing_entity_by_id(9000,EntityRank::NODE);

                // Index = 5-> location of node with id 9000 for processor 1.
                Mat< uint > tProcsSharedByEntityInd = tParallelMesh.get_procs_sharing_entity_by_index(5,EntityRank::NODE);

                uint tEntityOwner = tParallelMesh.parallel_owner_rank_by_entity_index(5,EntityRank::NODE);

                // Access maps (all of the local to global maps for all entity types).
                Mat < uint > lcl2GblNodeMap = tParallelMesh.get_nodal_local_map();
                Mat < uint > lcl2GblElemMap = tParallelMesh.get_elemental_local_map();
                Mat < uint > lcl2GblEdgeMap = tParallelMesh.get_edge_local_map();
                Mat < uint > lcl2GblFaceMap = tParallelMesh.get_face_local_map();

                // Access maps (all of the local to global maps for all entity types).
                Mat < uint > lcl2GblNodeOwnerProc = tParallelMesh.get_nodal_owner_proc_map();
                Mat < uint > lcl2GblElemOwnerProc = tParallelMesh.get_elemental_owner_proc_map();
                Mat < uint > lcl2GblEdgeOwnerProc = tParallelMesh.get_edge_owner_proc_map();
                Mat < uint > lcl2GblFaceOwnerProc = tParallelMesh.get_face_owner_proc_map();

                Cell < Cell < uint > > tNodesSharedPerProc = tParallelMesh.get_nodes_shared_processors();
                Cell < Cell < uint > > tElemsSharedPerProc = tParallelMesh.get_elements_shared_processors();
                Cell < Cell < uint > > tEdgesSharedPerProc = tParallelMesh.get_edges_shared_processors();
                Cell < Cell < uint > > tFacesSharedPerProc = tParallelMesh.get_faces_shared_processors();

                //Check entity owner
                REQUIRE(moris::equal_to(tEntityOwner,0));

                if( p_rank == 0 )
                {
                    // Check elements locally owned
                    REQUIRE(moris::equal_to(tElemsOwned(0,0),7));
                    REQUIRE(moris::equal_to(tElemsOwned(1,0),15));

                    // Check nodes locally owned
                    REQUIRE(moris::equal_to(tNodesOwned(0,0),2));
                    REQUIRE(moris::equal_to(tNodesOwned(1,0),4));
                    REQUIRE(moris::equal_to(tNodesOwned(2,0),5));
                    REQUIRE(moris::equal_to(tNodesOwned(3,0),6));
                    REQUIRE(moris::equal_to(tNodesOwned(4,0),8));
                    REQUIRE(moris::equal_to(tNodesOwned(5,0),10));
                    REQUIRE(moris::equal_to(tNodesOwned(6,0),12));
                    REQUIRE(moris::equal_to(tNodesOwned(7,0),1000));
                    REQUIRE(moris::equal_to(tNodesOwned(8,0),3000));
                    REQUIRE(moris::equal_to(tNodesOwned(9,0),7000));
                    REQUIRE(moris::equal_to(tNodesOwned(10,0),9000));
                    REQUIRE(moris::equal_to(tNodesOwned(11,0),11000));

                    //Check processors sharing inquired node
                    REQUIRE(moris::equal_to(tProcsSharedByEntityId(0,0),1));
                    REQUIRE(moris::equal_to(tProcsSharedByEntityInd(0,0),UINT_MAX)); // Node for local index 5 not shared

                    // Nodal maps
                    REQUIRE(moris::equal_to(lcl2GblNodeMap(0),1000));
                    REQUIRE(moris::equal_to(lcl2GblNodeMap(1),2));
                    REQUIRE(moris::equal_to(lcl2GblNodeOwnerProc(7),0));
                    REQUIRE(moris::equal_to(lcl2GblNodeOwnerProc(14),1));

                    // Elemental maps
                    REQUIRE(moris::equal_to(lcl2GblElemMap(0),15));
                    REQUIRE(moris::equal_to(lcl2GblElemOwnerProc(0),0));
                    REQUIRE(moris::equal_to(lcl2GblElemOwnerProc(2),1));

                    // Edge maps
                    REQUIRE(moris::equal_to(lcl2GblEdgeMap(19),20));
                    REQUIRE(moris::equal_to(lcl2GblEdgeMap(24),45));
                    REQUIRE(moris::equal_to(lcl2GblEdgeOwnerProc(10),0));
                    REQUIRE(moris::equal_to(lcl2GblEdgeOwnerProc(25),1));

                    // Face maps
                    REQUIRE(moris::equal_to(lcl2GblFaceMap(11),19));
                    REQUIRE(moris::equal_to(lcl2GblFaceMap(15),24));
                    REQUIRE(moris::equal_to(lcl2GblFaceOwnerProc(3),0));
                    REQUIRE(moris::equal_to(lcl2GblFaceOwnerProc(14),1));

                    // Check shared structure (same as proc 1 since it is shared info)
                    REQUIRE(moris::equal_to(tNodesSharedPerProc.size(),1));
                    REQUIRE(moris::equal_to(tNodesSharedPerProc(0).size(),4));
                    REQUIRE(moris::equal_to(tNodesSharedPerProc(0)(0),10));
                    REQUIRE(moris::equal_to(tNodesSharedPerProc(0)(1),12));
                    REQUIRE(moris::equal_to(tNodesSharedPerProc(0)(2),9000));
                    REQUIRE(moris::equal_to(tNodesSharedPerProc(0)(3),11000));
                    REQUIRE(moris::equal_to(tElemsSharedPerProc(0).size(),0)); // no aura, no shared elements
                    REQUIRE(moris::equal_to(tEdgesSharedPerProc.size(),1));
                    REQUIRE(moris::equal_to(tEdgesSharedPerProc(0).size(),4));
                    REQUIRE(moris::equal_to(tEdgesSharedPerProc(0)(0),5));
                    REQUIRE(moris::equal_to(tEdgesSharedPerProc(0)(1),6));
                    REQUIRE(moris::equal_to(tEdgesSharedPerProc(0)(2),7));
                    REQUIRE(moris::equal_to(tEdgesSharedPerProc(0)(3),8)); //
                    REQUIRE(moris::equal_to(tFacesSharedPerProc.size(),1));
                    REQUIRE(moris::equal_to(tFacesSharedPerProc(0).size(),1));
                    REQUIRE(moris::equal_to(tFacesSharedPerProc(0)(0),6));

                    // ===================================================
                    // Testing entities connected to element with ID = 7
                    // ===================================================
                    uint elementID = 7;

                    // Initialize and fill cells to store IDs of faces, edges and nodes connected to current element (elementID = 1)

                    Mat< uint > elemsConnectedToElement = tParallelMesh.get_elements_connected_to_element(elementID);

                    // Get number of elements, faces and edges connected to current node
                    uint NumberOfElemsConnectedToElement = elemsConnectedToElement.length();

                    // Check the number of elements and its IDs connected to current element
                    REQUIRE( moris::equal_to(NumberOfElemsConnectedToElement, 2) );
                    REQUIRE( moris::equal_to(elemsConnectedToElement(0), 15) );
                    REQUIRE( moris::equal_to(elemsConnectedToElement(1), 31) );
                }
                else
                {
                    // Check elements locally owned
                    REQUIRE(moris::equal_to(tElemsOwned(0,0),31));

                    // Check nodes locally owned
                    REQUIRE(moris::equal_to(tNodesOwned(0,0),13));
                    REQUIRE(moris::equal_to(tNodesOwned(1,0),14));
                    REQUIRE(moris::equal_to(tNodesOwned(2,0),16));
                    REQUIRE(moris::equal_to(tNodesOwned(3,0),15000));

                    //Check processors sharing inquired node
                    REQUIRE(moris::equal_to(tProcsSharedByEntityId(0,0),0));
                    REQUIRE(moris::equal_to(tProcsSharedByEntityInd(0,0),0));

                    // Nodal map
                    REQUIRE(moris::equal_to(lcl2GblNodeMap(2),7000));
                    REQUIRE(moris::equal_to(lcl2GblNodeMap(5),10));
                    REQUIRE(moris::equal_to(lcl2GblNodeOwnerProc(7),0));
                    REQUIRE(moris::equal_to(lcl2GblNodeOwnerProc(9),1));

                    // Elemental map
                    REQUIRE(moris::equal_to(lcl2GblElemMap(0),7));
                    REQUIRE(moris::equal_to(lcl2GblElemOwnerProc(0),0));
                    REQUIRE(moris::equal_to(lcl2GblElemOwnerProc(1),1));

                    // Edge map
                    REQUIRE(moris::equal_to(lcl2GblEdgeMap(12),41));
                    REQUIRE(moris::equal_to(lcl2GblEdgeMap(17),46));
                    REQUIRE(moris::equal_to(lcl2GblEdgeOwnerProc(10),0));
                    REQUIRE(moris::equal_to(lcl2GblEdgeOwnerProc(15),1));

                    // Face map
                    REQUIRE(moris::equal_to(lcl2GblFaceMap(6),19));
                    REQUIRE(moris::equal_to(lcl2GblFaceMap(2),3));
                    REQUIRE(moris::equal_to(lcl2GblFaceOwnerProc(3),0));
                    REQUIRE(moris::equal_to(lcl2GblFaceOwnerProc(9),1));

                    // Check shared structure (same as proc 0 since it is shared info)
                    REQUIRE(moris::equal_to(tNodesSharedPerProc.size(),1));
                    REQUIRE(moris::equal_to(tNodesSharedPerProc(0).size(),4));
                    REQUIRE(moris::equal_to(tNodesSharedPerProc(0)(0),10));
                    REQUIRE(moris::equal_to(tNodesSharedPerProc(0)(1),12));
                    REQUIRE(moris::equal_to(tNodesSharedPerProc(0)(2),9000));
                    REQUIRE(moris::equal_to(tNodesSharedPerProc(0)(3),11000));
                    REQUIRE(moris::equal_to(tElemsSharedPerProc(0).size(),0)); // no aura, no shared elements
                    REQUIRE(moris::equal_to(tEdgesSharedPerProc.size(),1));
                    REQUIRE(moris::equal_to(tEdgesSharedPerProc(0).size(),4));
                    REQUIRE(moris::equal_to(tEdgesSharedPerProc(0)(0),5));
                    REQUIRE(moris::equal_to(tEdgesSharedPerProc(0)(1),6));
                    REQUIRE(moris::equal_to(tEdgesSharedPerProc(0)(2),7));
                    REQUIRE(moris::equal_to(tEdgesSharedPerProc(0)(3),8)); //
                    REQUIRE(moris::equal_to(tFacesSharedPerProc.size(),1));
                    REQUIRE(moris::equal_to(tFacesSharedPerProc(0).size(),1));
                    REQUIRE(moris::equal_to(tFacesSharedPerProc(0)(0),6));
                }
            }
        }
    }

    SECTION( "Creating a 3D 3 element mesh in parallel using nodal processor owner list in parallel")
    {
        if( p_size == 2 ) // specify it is a 2 processor test
        {
            uint aNumDim = 3;
            Mat< moris::uint >  aNodeProcs;
            Mat< moris::uint >  aElemConn;
            Mat< moris::real >  aNodeCoords;
            Mat<moris::uint>    aElemLocaltoGlobal;
            Mat<moris::uint>    aNodeLocaltoGlobal;

            Cell < Mat< real > > aFieldData;
            Cell < std::string > aFieldName;
            Cell < enum EntityRank > aFieldRanks ;

            if ( p_rank == 0 )
            {
                // Mesh information
                // Generate data for test
                Mat< moris::uint>  tNodeProcsDummy  = {{0}, {0}, {0}, {0},
                        {0}, {0}, {0}, {0},
                        {0}, {0}, {0}, {0}};
                Mat< moris::uint > tElemConnDummy   = {{1000,2,4,3000,5,6,8,7000},
                        {5,6,8,7000,9000,10,12,11000}};
                Mat< moris::real > tNodeCoordsDummy = {{0.0, 0.0, 0.0},
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
                Mat< moris::uint > tElemLocaltoGlobalDummy = {{15},{7}};
                Mat< moris::uint > tNodeLocaltoGlobalDummy = {{1000},{2},{3000},{4},{5},{6},{7000},{8},{9000},
                        {10},{11000},{12}};
                aNodeProcs         = tNodeProcsDummy;
                aElemConn          = tElemConnDummy;
                aNodeCoords        = tNodeCoordsDummy;
                aElemLocaltoGlobal = tElemLocaltoGlobalDummy;
                aNodeLocaltoGlobal = tNodeLocaltoGlobalDummy;

                // Field information
                // Create nodal and elemental fields of reals
                Mat< real > tNodalFieldData_1 = {{10.0},{10.0},{10.0},{10.0},
                        {50.0},{50.0},{35.0},{35.0},
                        {40.0},{45.0},{45.0},{40.0}};
                Mat< real > aElemFieldData_1 = {{111.11, 222.22, 333.33},{0.253, -0.253, 0.976}};

                aFieldData  = { tNodalFieldData_1, aElemFieldData_1};
                aFieldName  = { "nodalField_1", "elementField_1"};
                aFieldRanks = { EntityRank::NODE, EntityRank::ELEMENT };

            } else {
                // Mesh information
                // Generate data for test
                Mat< moris::uint>  tNodeProcsDummy  = {{0}, {0}, {0}, {0},
                        {1}, {1}, {1}, {1}};
                Mat< moris::uint > tElemConnDummy   = {{9000,10,12,11000,13,14,16,15000}};
                Mat< moris::real > tNodeCoordsDummy = {{0.0, 0.0, 2.0},
                        {1.0, 0.0, 2.0},
                        {0.0, 1.0, 2.0},
                        {1.0, 1.0, 2.0},
                        {0.0, 0.0, 3.0},
                        {1.0, 0.0, 3.0},
                        {0.0, 1.0, 3.0},
                        {1.0, 1.0, 3.0}};
                Mat< moris::uint > tElemLocaltoGlobalDummy = {{31}};
                Mat< moris::uint > tNodeLocaltoGlobalDummy = {{9000},{10},{11000},{12},{13},{14},{15000},{16}};

                aNodeProcs         = tNodeProcsDummy;
                aElemConn          = tElemConnDummy;
                aNodeCoords        = tNodeCoordsDummy;
                aElemLocaltoGlobal = tElemLocaltoGlobalDummy;
                aNodeLocaltoGlobal = tNodeLocaltoGlobalDummy;

                // Field information
                // Create nodal and elemental fields of reals
                Mat< real > tNodalFieldData_1 = {{40.0},{45.0},{45.0},{40.0},
                        {60.0},{65.0},{65.0},{60.0}};
                Mat< real > aElemFieldData_1 = {{-1.0, -3.0, -4.0}};

                aFieldData  = { tNodalFieldData_1, aElemFieldData_1};
                aFieldName  = { "nodalField_1", "elementField_1"};
                aFieldRanks = { EntityRank::NODE, EntityRank::ELEMENT };
            }

            SECTION( "with multiple fields")
            {
                // Create MORIS mesh using MTK database
                moris::MtkFieldsInfo aFieldsInfo;
                aFieldsInfo.FieldsData = &aFieldData;
                aFieldsInfo.FieldsName = aFieldName;
                aFieldsInfo.FieldsRank = aFieldRanks;

                moris::MtkMeshData aMeshData;
                aMeshData.SpatialDim    = &aNumDim;
                aMeshData.ElemConn      = &aElemConn;
                aMeshData.NodeCoords    = &aNodeCoords;
                aMeshData.EntProcOwner  = &aNodeProcs;
                aMeshData.FieldsInfo     = &aFieldsInfo;
                aMeshData.LocaltoGlobalElemMap = &aElemLocaltoGlobal;
                aMeshData.LocaltoGlobalNodeMap = &aNodeLocaltoGlobal;

                moris::mesh tParallelMesh( MeshType::MTK, aMeshData );

                // ========================================
                // Testing get_entities_owned_current_proc
                // ========================================

                Mat< uint > tElemsOwned  = tParallelMesh.get_entities_owned_current_proc(EntityRank::ELEMENT);
                Mat< uint > tNodesOwned  = tParallelMesh.get_entities_owned_current_proc(EntityRank::NODE);
                Mat< real > tNodalField_1 = tParallelMesh.get_field_values( EntityRank::NODE, "nodalField_1" );
                Mat< real > tElemField_1  = tParallelMesh.get_field_values( EntityRank::ELEMENT, "elementField_1" );

                if( p_rank == 0 )
                {
                    // Check some entities locally owned
                    REQUIRE(moris::equal_to(tNodesOwned(0,0),2));
                    REQUIRE(moris::equal_to(tNodesOwned(2,0),5));
                    REQUIRE(moris::equal_to(tNodesOwned(5,0),10));
                    REQUIRE(moris::equal_to(tElemsOwned(0,0),7));
                    REQUIRE(moris::equal_to(tElemsOwned(1,0),15));

                    // Check some field values
                    REQUIRE(moris::equal_to(tNodalField_1(0,0),10.0));
                    REQUIRE(moris::equal_to(tNodalField_1(2,0),50.0));
                    REQUIRE(moris::equal_to(tNodalField_1(4,0),35.0));
                    REQUIRE(moris::equal_to(tNodalField_1(10,0),40.0));
                    REQUIRE(moris::equal_to(tElemField_1(0,0),0.253));
                    REQUIRE(moris::equal_to(tElemField_1(1,1),222.22));
                }
                else
                {
                    // Check some entities locally owned
                    REQUIRE(moris::equal_to(tNodesOwned(0,0),13));
                    REQUIRE(moris::equal_to(tNodesOwned(3,0),15000));
                    REQUIRE(moris::equal_to(tElemsOwned(0,0),31));

                    // Check some field values
                    REQUIRE(moris::equal_to(tNodalField_1(0,0),45.0));
                    REQUIRE(moris::equal_to(tNodalField_1(2,0),60.0));
                    REQUIRE(moris::equal_to(tNodalField_1(4,0),60.0));
                    REQUIRE(moris::equal_to(tElemField_1(0,1),-3.0));
                    REQUIRE(moris::equal_to(tElemField_1(0,2),-4.0));
                }
            }
        }
    }

    SECTION( "Creating a 3D 3 element mesh in parallel")
    {
        if( p_size == 2 ) // specify it is a serial test only
        {
            uint aNumDim = 3;
            Mat< uint > aNodeProcs;
            Mat< uint > aElemConn;
            Mat< real > aNodeCoords;
            Mat< uint > aElemLocaltoGlobal;
            Mat< uint > aNodeLocaltoGlobal;

            if ( p_rank == 0 )
            {
                // Mesh information
                // Generate data for test
                Mat< moris::uint>  tNodeProcsDummy  = {{0}, {0}, {0}, {0},
                        {0}, {0}, {0}, {0},
                        {0}, {0}, {0}, {0}};
                Mat< moris::uint > tElemConnDummy   = {{1000,2,4,3000,5,6,8,7000},
                        {5,6,8,7000,9000,10,12,11000}};
                Mat< moris::real > tNodeCoordsDummy = {{0.0, 0.0, 0.0},
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
                Mat< moris::uint > tElemLocaltoGlobalDummy = {{15},{7}};
                Mat< moris::uint > tNodeLocaltoGlobalDummy = {{1000},{2},{3000},{4},{5},{6},{7000},{8},{9000},
                        {10},{11000},{12}};
                aNodeProcs         = tNodeProcsDummy;
                aElemConn          = tElemConnDummy;
                aNodeCoords        = tNodeCoordsDummy;
                aElemLocaltoGlobal = tElemLocaltoGlobalDummy;
                aNodeLocaltoGlobal = tNodeLocaltoGlobalDummy;
            }
            else
            {
                // Mesh information
                // Generate data for test
                Mat< moris::uint>  tNodeProcsDummy  = {{0}, {0}, {0}, {0},
                        {1}, {1}, {1}, {1}};
                Mat< moris::uint > tElemConnDummy   = {{9000,10,12,11000,13,14,16,15000}};
                Mat< moris::real > tNodeCoordsDummy = {{0.0, 0.0, 2.0},
                        {1.0, 0.0, 2.0},
                        {0.0, 1.0, 2.0},
                        {1.0, 1.0, 2.0},
                        {0.0, 0.0, 3.0},
                        {1.0, 0.0, 3.0},
                        {0.0, 1.0, 3.0},
                        {1.0, 1.0, 3.0}};
                Mat< moris::uint > tElemLocaltoGlobalDummy = {{31}};
                Mat< moris::uint > tNodeLocaltoGlobalDummy = {{9000},{10},{11000},{12},{13},{14},{15000},{16}};

                aNodeProcs         = tNodeProcsDummy;
                aElemConn          = tElemConnDummy;
                aNodeCoords        = tNodeCoordsDummy;
                aElemLocaltoGlobal = tElemLocaltoGlobalDummy;
                aNodeLocaltoGlobal = tNodeLocaltoGlobalDummy;
            }

            SECTION( "without fields or sets for parallel data check")
            {
                moris::MtkMeshData aMeshData;
                aMeshData.SpatialDim    = &aNumDim;
                aMeshData.ElemConn      = &aElemConn;
                aMeshData.NodeCoords    = &aNodeCoords;
                aMeshData.EntProcOwner  = &aNodeProcs;
                aMeshData.LocaltoGlobalElemMap = &aElemLocaltoGlobal;
                aMeshData.LocaltoGlobalNodeMap = &aNodeLocaltoGlobal;

                moris::mesh tParallelMesh( MeshType::MTK, aMeshData );

                // ========================================
                // Testing get_entities_owned_current_proc
                // ========================================

                Mat< uint > tNodeOwned  = tParallelMesh.get_entities_owned_current_proc(EntityRank::NODE);
                Mat< uint > tNodeShared = tParallelMesh.get_entities_glb_shared_current_proc(EntityRank::NODE);
                Mat< uint > tNodeAura   = tParallelMesh.get_entities_in_aura(EntityRank::NODE);
                Mat< uint > tNodeUniv   = tParallelMesh.get_entities_universal(EntityRank::NODE);
                Mat< uint > tNodeOwnedAndShared  = tParallelMesh.get_entities_owned_and_shared_by_current_proc(EntityRank::NODE);

                Mat< uint > tElemOwned = tParallelMesh.get_entities_owned_current_proc(EntityRank::ELEMENT);
                Mat< uint > tElemAura  = tParallelMesh.get_entities_in_aura(EntityRank::ELEMENT);
                Mat< uint > tElemUniv  = tParallelMesh.get_entities_universal(EntityRank::ELEMENT);
                Mat< uint > tElemOwnedAndShared  = tParallelMesh.get_entities_owned_and_shared_by_current_proc(EntityRank::ELEMENT);

                if( p_rank == 0 )
                {
                    // Check some of the requested data
                    // For this processor we don't have shared nodes
                    REQUIRE(moris::equal_to(tNodeOwned(3),6));
                    REQUIRE(moris::equal_to(tNodeOwnedAndShared(10),11000));
                    REQUIRE(moris::equal_to(tNodeAura(0),13));
                    REQUIRE(moris::equal_to(tNodeUniv(15),15000));

                    // Elements don't have shared entities
                    REQUIRE(moris::equal_to(tElemOwned(1),15));
                    REQUIRE(moris::equal_to(tElemAura(0),31));
                    REQUIRE(moris::equal_to(tElemUniv(2),31));
                }
                else
                {
                    // Check some of the requested data
                    // For this processor we don't have shared nodes
                    REQUIRE(moris::equal_to(tNodeOwned(1),14));
                    REQUIRE(moris::equal_to(tNodeShared(2),9000));
                    REQUIRE(moris::equal_to(tNodeOwnedAndShared(7),16));
                    REQUIRE(moris::equal_to(tNodeUniv(7),15000));

                    // Elements don't have shared entities
                    REQUIRE(moris::equal_to(tElemOwned(0),31));
                    REQUIRE(moris::equal_to(tElemUniv(0),31));
                }
            }

            SECTION( "with multiple fields on the entire domain")
            {
                // Test for showing how to provide information of multiple fields and output them into
                // a mesh file. In this example 1 nodal and 1 element fields of reals are provided for
                // a three element (HEX8) mesh generated from data.

                Cell < Mat< real > > aFieldData;
                Cell < std::string > aFieldName;
                Cell < enum EntityRank > aFieldRanks ;

                if ( p_rank == 0 )
                {
                    // Field information
                    // Create nodal and elemental fields of reals
                    Mat< real > tNodalFieldData_1 = {{10.0},{10.0},{10.0},{10.0},
                            {50.0},{50.0},{35.0},{35.0},
                            {40.0},{45.0},{45.0},{40.0}};
                    Mat< real > aElemFieldData_1 = {{111.11, 222.22, 333.33},{0.253, -0.253, 0.976}};

                    aFieldData  = { tNodalFieldData_1, aElemFieldData_1};
                    aFieldName  = { "nodalField_1", "elementField_1"};
                    aFieldRanks = { EntityRank::NODE, EntityRank::ELEMENT };
                }
                else
                {
                    // Field information
                    // Create nodal and elemental fields of reals
                    Mat< real > tNodalFieldData_1 = {{40.0},{45.0},{45.0},{40.0},
                            {60.0},{65.0},{65.0},{60.0}};
                    Mat< real > aElemFieldData_1 = {{-1.0, -3.0, -4.0}};

                    aFieldData  = { tNodalFieldData_1, aElemFieldData_1};
                    aFieldName  = { "nodalField_1", "elementField_1"};
                    aFieldRanks = { EntityRank::NODE, EntityRank::ELEMENT };
                }

                // Create MORIS mesh using MTK database
                moris::MtkFieldsInfo aFieldsInfo;
                aFieldsInfo.FieldsData = &aFieldData;
                aFieldsInfo.FieldsName = aFieldName;
                aFieldsInfo.FieldsRank = aFieldRanks;

                moris::MtkMeshData aMeshData;
                aMeshData.SpatialDim    = &aNumDim;
                aMeshData.ElemConn      = &aElemConn;
                aMeshData.NodeCoords    = &aNodeCoords;
                aMeshData.EntProcOwner  = &aNodeProcs;
                aMeshData.FieldsInfo     = &aFieldsInfo;
                aMeshData.LocaltoGlobalElemMap = &aElemLocaltoGlobal;
                aMeshData.LocaltoGlobalNodeMap = &aNodeLocaltoGlobal;

                moris::mesh tParallelMesh( MeshType::MTK, aMeshData );

                // ========================================
                // Testing get_entities_owned_current_proc
                // ========================================

                Mat< uint > tElemsOwned   = tParallelMesh.get_entities_owned_current_proc(EntityRank::ELEMENT);
                Mat< uint > tNodesOwned   = tParallelMesh.get_entities_owned_current_proc(EntityRank::NODE);
                Mat< real > tNodalField_1 = tParallelMesh.get_field_values( EntityRank::NODE, "nodalField_1" );
                Mat< real > tElemField_1  = tParallelMesh.get_field_values( EntityRank::ELEMENT, "elementField_1" );
                Mat< uint > tFieldEnts_1  = tParallelMesh.get_field_entities( EntityRank::NODE, "nodalField_1" );
                Mat< uint > tFieldEnts_2  = tParallelMesh.get_field_entities( EntityRank::ELEMENT, "elementField_1" );

                if( p_rank == 0 )
                {
                    // Check some entities owned
                    REQUIRE(moris::equal_to(tFieldEnts_1(0),2));
                    REQUIRE(moris::equal_to(tFieldEnts_1(2),5));
                    REQUIRE(moris::equal_to(tFieldEnts_1(4),8));
                    REQUIRE(moris::equal_to(tFieldEnts_1(8),3000));
                    REQUIRE(moris::equal_to(tFieldEnts_1(10),9000));
                    REQUIRE(moris::equal_to(tFieldEnts_2(0),7));
                    REQUIRE(moris::equal_to(tFieldEnts_2(1),15));

                    // Check some field values
                    REQUIRE(moris::equal_to(tNodalField_1(0,0),10.0));
                    REQUIRE(moris::equal_to(tNodalField_1(4,0),35.0));
                    REQUIRE(moris::equal_to(tNodalField_1(6,0),40.0));
                    REQUIRE(moris::equal_to(tNodalField_1(8,0),10.0));
                    REQUIRE(moris::equal_to(tElemField_1(1,0),111.11));
                    REQUIRE(moris::equal_to(tElemField_1(0,1),-0.253));

                }
                else
                {
                    // Check some entities owned
                    REQUIRE(moris::equal_to(tFieldEnts_1(0),10));
                    REQUIRE(moris::equal_to(tFieldEnts_1(2),13));
                    REQUIRE(moris::equal_to(tFieldEnts_1(7),15000));
                    REQUIRE(moris::equal_to(tFieldEnts_2(0),31));

                    // Check some field values
                    REQUIRE(moris::equal_to(tNodalField_1(0,0),45.0));
                    REQUIRE(moris::equal_to(tNodalField_1(2,0),60.0));
                    REQUIRE(moris::equal_to(tNodalField_1(4,0),60.0));
                    REQUIRE(moris::equal_to(tElemField_1(0,2),-4.0));
                }
            }

            SECTION( "with sets allocated in separate processors")
            {
                Mat< uint > tBlockSetsPartOwners;
                Cell< Mat< uint > > tNodeSetsEntIds;
                Cell< Mat< uint > > tSideSetsInfo;

                if ( p_rank == 0 )
                {
                    // Populate block sets
                    Mat<uint> tBlockSetsInds = { {1}, {0} };
                    tBlockSetsPartOwners = tBlockSetsInds;

                    // Populate side sets
                    Mat< uint > tSideset_1 = { {15, 5}, {7, 0} };
                    tSideSetsInfo = { tSideset_1 };

                    // Populate node sets
                    Mat< uint > tNodeSet_1 = { {1000}, {4}, {6} };
                    tNodeSetsEntIds = { tNodeSet_1 };
                }
                else
                {
                    // Populate block sets
                    Mat<uint> tBlockSetsInds = { {0} };
                    tBlockSetsPartOwners = tBlockSetsInds;

                    // Populate side sets
                    Mat< uint > tSideset_1 = { {31, 3} };
                    tSideSetsInfo = { tSideset_1 };

                    // Populate node sets
                    Mat< uint > tNodeSet_1 = { };
                    tNodeSetsEntIds = { tNodeSet_1 };
                }

                // Declare block sets
                /////////////////////
                MtkBlockSetsInfo tBlockSetStruc;
                tBlockSetStruc.BSetInds = &tBlockSetsPartOwners;
                tBlockSetStruc.BSetNames   = { "blockset_001", "blockset_002" };

                // Declare side sets
                /////////////////////
                MtkSideSetsInfo tSideSetStruc;
                tSideSetStruc.ElemIdsAndSideOrds = &tSideSetsInfo;
                tSideSetStruc.SSetNames   = { "Sideset_001" };

                // Declare node sets
                /////////////////////
                MtkNodeSetsInfo tNodeSetStruc;
                tNodeSetStruc.EntIds = &tNodeSetsEntIds;
                tNodeSetStruc.NSetNames = { "Nodeset_001" };

                // Create MORIS mesh using MTK database
                ///////////////////////////////////////
                MtkSetsInfo aMeshSets;
                aMeshSets.NodeSetsInfo = &tNodeSetStruc;
                aMeshSets.SideSetsInfo   = &tSideSetStruc;
                aMeshSets.BlockSetsInfo   = &tBlockSetStruc;

                moris::MtkMeshData aMeshData;
                aMeshData.SpatialDim   = &aNumDim;
                aMeshData.ElemConn     = &aElemConn;
                aMeshData.NodeCoords   = &aNodeCoords;
                aMeshData.EntProcOwner = &aNodeProcs;
                aMeshData.SetsInfo     = &aMeshSets;
                aMeshData.LocaltoGlobalElemMap = &aElemLocaltoGlobal;
                aMeshData.LocaltoGlobalNodeMap = &aNodeLocaltoGlobal;

                moris::mesh tParallelMesh( MeshType::MTK, aMeshData );

                // ========================================
                // Testing get_entities_owned_current_proc
                // ========================================

                Mat< uint > tElemsOwned = tParallelMesh.get_entities_owned_current_proc(EntityRank::ELEMENT);
                Mat< uint > tNodesOwned = tParallelMesh.get_entities_owned_current_proc(EntityRank::NODE);

                Mat< uint > tNodeSet = tParallelMesh.get_set_entity_ids( EntityRank::NODE, "NodeSet_001" );
                Mat< uint > tSideSet = tParallelMesh.get_set_entity_ids( EntityRank::FACE, "Sideset_001" );
                Mat< uint > tBlockSet1 = tParallelMesh.get_set_entity_ids( EntityRank::ELEMENT, "Blockset_001" );
                Mat< uint > tBlockSet2 = tParallelMesh.get_set_entity_ids( EntityRank::ELEMENT, "Blockset_002" );

                Mat< uint > tNodesInNodeSet = tParallelMesh.get_nodes_in_node_set( 1 );
                Mat< uint > tNodesInSideSet = tParallelMesh.get_nodes_in_side_set( 1 );

                if( p_rank == 0 )
                {
                    // Check some entities owned
                    REQUIRE(moris::equal_to(tNodeSet(0),4));
                    REQUIRE(moris::equal_to(tNodeSet(1),6));
                    REQUIRE(moris::equal_to(tNodeSet(2),1000));
                    REQUIRE(moris::equal_to(tSideSet(0),71));
                    REQUIRE(moris::equal_to(tSideSet(1),156));
                    REQUIRE(moris::equal_to(tBlockSet2(0),15));
                    REQUIRE(moris::equal_to(tBlockSet1(0),7));
                }
                else
                {
                    // zero nodes in node set for this processor
                    // Check entity owned
                    REQUIRE(moris::equal_to(tSideSet(0),314));
                    REQUIRE(moris::equal_to(tBlockSet1(0),31));
                }

                // Output mesh
                std::string OutputFileName1 = "meshWMultipleSetsPara.exo";  //output name
                tParallelMesh.create_output_mesh(OutputFileName1);
            }
        }
    }

    SECTION( "Creating a 1D mesh from data in parallel" )
    {
        if( p_size ==2 ) // specify it is a serial test only
        {
            // Generate data for test
            uint aNumDim = 1;
            Mat<uint >  aNodeProcs;
            Mat<uint >  aElemConn;
            Mat<real >  aNodeCoords;
            Mat<uint>    aElemLocaltoGlobal;
            Mat<uint>    aNodeLocaltoGlobal;

            if ( p_rank == 0 )
            {
                // Generate data for test
                Mat< uint>  tNodeProcsDummy  = {{0}, {0}, {0}, {0}};
                Mat< uint > tElemConnDummy   = {{1,2},{2,5},{5,6}};
                Mat< real > tNodeCoordsDummy = {{0.0,0.0},{1.0,0.0},{2.0,0.0},{3.0,0.0}};
                Mat< uint > tElemLocaltoGlobalDummy = {{1},{3},{10}};
                Mat< uint > tNodeLocaltoGlobalDummy = {{1},{2},{5},{6}};
                aNodeProcs         = tNodeProcsDummy;
                aElemConn          = tElemConnDummy;
                aNodeCoords        = tNodeCoordsDummy;
                aElemLocaltoGlobal = tElemLocaltoGlobalDummy;
                aNodeLocaltoGlobal = tNodeLocaltoGlobalDummy;

            } else {

                // Generate data for test
                moris::Mat< uint>  tNodeProcsDummy  = {{0}, {1}};
                Mat< uint > tElemConnDummy   = {{6,3}};
                Mat< real > tNodeCoordsDummy = {{3.0,0.0},{4.0,0.0}};
                Mat< uint > tElemLocaltoGlobalDummy = {{2}};
                Mat< uint > tNodeLocaltoGlobalDummy = {{6},{3}};

                aNodeProcs         = tNodeProcsDummy;
                aElemConn          = tElemConnDummy;
                aNodeCoords        = tNodeCoordsDummy;
                aElemLocaltoGlobal = tElemLocaltoGlobalDummy;
                aNodeLocaltoGlobal = tNodeLocaltoGlobalDummy;
            }

            moris::MtkMeshData aMeshData;
            aMeshData.SpatialDim   = &aNumDim;
            aMeshData.ElemConn     = &aElemConn;
            aMeshData.NodeCoords   = &aNodeCoords;
            aMeshData.EntProcOwner = &aNodeProcs;
            aMeshData.LocaltoGlobalElemMap = &aElemLocaltoGlobal;
            aMeshData.LocaltoGlobalNodeMap = &aNodeLocaltoGlobal;

            moris::mesh tMesh1DParallel( MeshType::MTK, aMeshData );

            // =============================
            // Testing basic functionalities
            // =============================

            Mat< uint > tElemsOwned  = tMesh1DParallel.get_entities_owned_current_proc(EntityRank::ELEMENT);
            Mat< uint > tNodesOwned  = tMesh1DParallel.get_entities_owned_current_proc(EntityRank::NODE);

            if ( p_rank == 0 )
            {
                REQUIRE(moris::equal_to(tNodesOwned(0),1));
                REQUIRE(moris::equal_to(tNodesOwned(1),2));
                REQUIRE(moris::equal_to(tNodesOwned(2),5));
                REQUIRE(moris::equal_to(tNodesOwned(3),6));

                REQUIRE(moris::equal_to(tElemsOwned(0),1));
                REQUIRE(moris::equal_to(tElemsOwned(1),3));
                REQUIRE(moris::equal_to(tElemsOwned(2),10));
            }
            else
            {
                REQUIRE(moris::equal_to(tNodesOwned(0),3));

                REQUIRE(moris::equal_to(tElemsOwned(0),2));
            }


            // Output mesh
            std::string OutputFileName1 = "mesh1DProblemParallel.exo";  //output name
            tMesh1DParallel.create_output_mesh(OutputFileName1);
        }
    }


    SECTION( "Creating a 2D mesh with 'hanging nodes' in serial")
    {
        if( p_size == 1 ) // serial test
        {
            // Verification of mesh behavior when giving hanging nodes in mesh

            //   *-----*-----*  nodes: 1,2,3
            //   |     |     |  elements: 1,2
            //   |     |     |
            //   *-----*-----*  nodes: 4,5,6
            //   |           |  element: 3
            //   |           |
            //   *-----------*  nodes 7,8

            moris::uint aNumDim = 2;
            moris::Mat< moris::uint >  aElemConn   = { {1,2,5,4}, {2,3,6,5}, {4,6,8,7} };
            moris::Cell< std::string > aPartNames  = { "block_1" };
            moris::Mat< moris::real >  aNodeCoords = {{ 0.0, 2.0, 0.0 },
                                                      { 1.0, 2.0, 0.0 },
                                                      { 2.0, 2.0, 0.0 },
                                                      { 0.0, 1.0, 0.0 },
                                                      { 1.0, 1.0, 0.0 },
                                                      { 2.0, 1.0, 0.0 },
                                                      { 0.0, 0.0, 0.0 },
                                                      { 2.0, 0.0, 0.0 }};

            moris::MtkMeshData aMeshData;
            aMeshData.SpatialDim             = &aNumDim;
            aMeshData.ElemConn               = &aElemConn;
            aMeshData.NodeCoords             = &aNodeCoords;
            aMeshData.CreateAllEdgesAndFaces = true;

            // Creating MTK mesh
            moris::mesh tHangingNodesMesh( MeshType::MTK, aMeshData );

            // Get node and element entities own by current processor.
            Mat< uint > nodesInElem1 = tHangingNodesMesh.get_nodes_connected_to_element( 1 );
            Mat< uint > nodesInElem2 = tHangingNodesMesh.get_nodes_connected_to_element( 2 );
            Mat< uint > nodesInElem3 = tHangingNodesMesh.get_nodes_connected_to_element( 3 );

            Mat< uint > edgesInElem1 = tHangingNodesMesh.get_edges_connected_to_element( 1 );
            Mat< uint > edgesInElem2 = tHangingNodesMesh.get_edges_connected_to_element( 2 );
            Mat< uint > edgesInElem3 = tHangingNodesMesh.get_edges_connected_to_element( 3 );

            Mat< uint > elemsConnToHangingNodes = tHangingNodesMesh.get_elements_connected_to_node( 5 );

            
            Mat< uint > elemsInEdge1 = tHangingNodesMesh.get_elements_connected_to_edge( 8 );
            Mat< uint > elemsInEdge2 = tHangingNodesMesh.get_elements_connected_to_edge( 9 );
            Mat< uint > elemsInEdge3 = tHangingNodesMesh.get_elements_connected_to_edge( 10 );
            Mat< uint > elemsInEdge4 = tHangingNodesMesh.get_elements_connected_to_edge( 11 );

            uint numEdgesInMesh = tHangingNodesMesh.get_num_edges( );

            // ========================================
            // Testing output related to special mesh
            // ========================================

            REQUIRE(moris::equal_to( nodesInElem1( 0 ), 1 ) );
            REQUIRE(moris::equal_to( nodesInElem1( 1 ), 2 ) );
            REQUIRE(moris::equal_to( nodesInElem1( 2 ), 5 ) );
            REQUIRE(moris::equal_to( nodesInElem1( 3 ), 4 ) );

            REQUIRE(moris::equal_to( nodesInElem2( 0 ), 2 ) );
            REQUIRE(moris::equal_to( nodesInElem2( 1 ), 3 ) );
            REQUIRE(moris::equal_to( nodesInElem2( 2 ), 6 ) );
            REQUIRE(moris::equal_to( nodesInElem2( 3 ), 5 ) );

            REQUIRE(moris::equal_to( nodesInElem3( 0 ), 4 ) );
            REQUIRE(moris::equal_to( nodesInElem3( 1 ), 6 ) );
            REQUIRE(moris::equal_to( nodesInElem3( 2 ), 8 ) );
            REQUIRE(moris::equal_to( nodesInElem3( 3 ), 7 ) );

            REQUIRE(moris::equal_to( edgesInElem1( 0 ), 1 ) );
            REQUIRE(moris::equal_to( edgesInElem1( 1 ), 2 ) );
            REQUIRE(moris::equal_to( edgesInElem1( 2 ), 3 ) );
            REQUIRE(moris::equal_to( edgesInElem1( 3 ), 4 ) );

            REQUIRE(moris::equal_to( edgesInElem2( 0 ), 5 ) );
            REQUIRE(moris::equal_to( edgesInElem2( 1 ), 6 ) );
            REQUIRE(moris::equal_to( edgesInElem2( 2 ), 7 ) );
            REQUIRE(moris::equal_to( edgesInElem2( 3 ), 2 ) );

            REQUIRE(moris::equal_to( edgesInElem3( 0 ), 8  ) );
            REQUIRE(moris::equal_to( edgesInElem3( 1 ), 9  ) );
            REQUIRE(moris::equal_to( edgesInElem3( 2 ), 10 ) );
            REQUIRE(moris::equal_to( edgesInElem3( 3 ), 11 ) );

            REQUIRE(moris::equal_to( elemsInEdge1( 0 ), 3 ) );
            REQUIRE(moris::equal_to( elemsInEdge2( 0 ), 3 ) );
            REQUIRE(moris::equal_to( elemsInEdge3( 0 ), 3 ) );
            REQUIRE(moris::equal_to( elemsInEdge4( 0 ), 3 ) );

            REQUIRE(moris::equal_to( elemsConnToHangingNodes( 0 ), 1 ) );
            REQUIRE(moris::equal_to( elemsConnToHangingNodes( 1 ), 2 ) );

            REQUIRE(moris::equal_to( numEdgesInMesh, 11 ) );
        }
    }
}

