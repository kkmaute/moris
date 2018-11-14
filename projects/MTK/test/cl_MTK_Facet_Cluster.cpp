/*
 * cl_MTK_Face_Cluster.cpp
 *
 *  Created on: Nov 5, 2018
 *      Author: doble
 */


#include <catch.hpp>
#include <iostream>

#include "cl_MTK_Facet_Cluster.hpp"
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

TEST_CASE( "MTK Face Cluster API", "[MTK_Face_Cluster]" )
{
//    // Spatial dimension of the mesh
//    uint aNumDim = 3;
//
//    // Number of element types (hex8)
//    uint aNumElemTypes = 1;
//
//    // Specify the element to node connectivity of all hex 8's in the mesh
//    Matrix< IdMat >  aElemConnHex8  = {{  7,  8, 11, 10, 1, 2,  5,  4},
//                                       {  8,  9, 12, 11, 2, 3,  6,  5},
//                                       { 13, 14, 17, 16, 7, 8, 11, 10},
//                                       { 14, 15, 18, 17, 8, 9, 12, 11},
//                                       { 16, 18, 22, 21, 4, 6, 20, 19}};
//
//
//    // Node coordinate matrix
//    Matrix< DDRMat >  aCoords   = {{0.0, 0.0, 1.0},
//                                   {0.5, 0.0, 1.0},
//                                   {1.0, 0.0, 1.0},
//                                   {0.0, 0.5, 1.0},
//                                   {0.5, 0.5, 1.0},
//                                   {1.0, 0.5, 1.0},
//                                   {0.0, 0.0, 0.5},
//                                   {0.5, 0.0, 0.5},
//                                   {1.0, 0.0, 0.5},
//                                   {0.0, 0.5, 0.5},
//                                   {0.5, 0.5, 0.5},
//                                   {1.0, 0.5, 0.5},
//                                   {0.0, 0.0, 0.0},
//                                   {0.5, 0.0, 0.0},
//                                   {1.0, 0.0, 0.0},
//                                   {0.0, 0.5, 0.0},
//                                   {0.5, 0.5, 0.0},
//                                   {1.0, 0.5, 0.0},
//                                   {0.0, 1.5, 1.0},
//                                   {1.0, 1.5, 1.0},
//                                   {0.0, 1.5, 0.0},
//                                   {1.0, 1.5, 0.0}};
//
//    // Specify the local to global map
//    Matrix< IdMat >  aNodeLocaltoGlobal     = {{1 ,2 ,3, 4, 5, 6 ,7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}};
//
//    // Specify the local to global element map for all hex8s in the mesh
//    Matrix< IdMat >  aElemLocaltoGlobalHex8 = {{1,2,3,4,5}};
//
//    // Initialize Sets information structure
//    moris::mtk::MtkSetsInfo tMtkMeshSets;
//
//    // Cells in block set ids
//    Matrix< IdMat >tCellIdsBS1({{1,2,3,4}});
//
//    // Place all elements with child face in a block together
//    moris::mtk::MtkBlockSetInfo tBlockSet1;
//    tBlockSet1.mCellIdsInSet = &tCellIdsBS1;
//    tBlockSet1.mBlockSetName = "elems_w_child_faces";
//    tBlockSet1.mBlockSetTopo = CellTopology::HEX8;
//
//    // Place all elements with parent face in a block together
//    Matrix< IdMat >tCellIdsBS2({{5}});
//    moris::mtk::MtkBlockSetInfo tBlockSet2;
//    tBlockSet2.mCellIdsInSet = &tCellIdsBS2;
//    tBlockSet2.mBlockSetName = "elems_w_parent_faces";
//    tBlockSet2.mBlockSetTopo = CellTopology::HEX8;
//
//    // Add block sets to mtk mesh sets
//    tMtkMeshSets.add_block_set(&tBlockSet1);
//    tMtkMeshSets.add_block_set(&tBlockSet2);
//
//    // Place data in mesh data container for input
//    moris::mtk::MtkMeshData aMeshData(aNumElemTypes);
//    aMeshData.SpatialDim              = &aNumDim;
//    aMeshData.ElemConn(0)             = &aElemConnHex8;
//    aMeshData.NodeCoords              = &aCoords;
//    aMeshData.SetsInfo                = &tMtkMeshSets;
//    aMeshData.LocaltoGlobalElemMap(0) = &aElemLocaltoGlobalHex8;
//    aMeshData.LocaltoGlobalNodeMap    = &aNodeLocaltoGlobal;
//    aMeshData.CreateAllEdgesAndFaces  = true;
//
//    // Create mesh from data with the factory
//    moris::mtk::Mesh* tMesh = create_mesh( MeshType::STK, aMeshData );
//
//
//    uint tNumFaces = tMesh->get_num_entities(EntityRank::FACE);
//
//    // Output to exodus file
//    std::string tPrefix = std::getenv("MORISOUTPUT");
//    std::string tMeshOutputFile = tPrefix + "/Face_Cluster_Unit_Test.e";
//    tMesh->create_output_mesh(tMeshOutputFile);
//    std::cout<<"Mesh output file: " << tMeshOutputFile<<std::endl;
//
//    // Construct a face cluster
//    // Note: parent face was found by looking at the faces created by STK
//    //       since STK does not allow for user declaration of faces
//    moris_index  tParentFaceIndex = 20;
//    uint         tNumChildFaces   = 4;
//    Facet_Cluster tFacetCluster(tParentFaceIndex,tNumChildFaces);
//
//    // Child face Ordinal 0 (node ids 4,5,11,10)
//    moris_index tChildOrdinal0          = 0;
//    moris_index tChildFaceIndex0        = 2;
//    Matrix< DDRMat > tParametricBounds0 = {{1.0,1.0}, {0.0,1.0}, {0.0,0.0}, {1.0,0.0}};
//
//    tFacetCluster.set_child_face(tChildOrdinal0, tChildFaceIndex0, tParametricBounds0);
//
//    // Child Face Ordinal 1 (node ids 5,6,12,11)
//    moris_index tChildOrdinal1          = 1;
//    moris_index tChildFaceIndex1        = 8;
//    Matrix< DDRMat > tParametricBounds1 = {{0.0,1.0}, {-1.0,1.0}, {-1.0,0.0}, {0.0,0.0}};
//    tFacetCluster.set_child_face(tChildOrdinal1, tChildFaceIndex1, tParametricBounds1);
//
//    // Child Face Ordinal 2
//    moris_index tChildOrdinal2          =  2;
//    moris_index tChildFaceIndex2        = 18;
//    Matrix< DDRMat > tParametricBounds2 = {{0.0,0.0},{-1.0,0.0},{-1.0,-1.0},{0.0,-1.0}};
//    tFacetCluster.set_child_face(tChildOrdinal2, tChildFaceIndex2, tParametricBounds2);
//
//    // Child Face Ordinal 3
//    moris_index tChildOrdinal3          =  3;
//    moris_index tChildFaceIndex3        = 13;
//    Matrix< DDRMat > tParametricBounds3 = {{1.0,0.0},{0.0,0.0},{0.0,-1.0},{1.0,-1.0}};
//    tFacetCluster.set_child_face(tChildOrdinal3, tChildFaceIndex3, tParametricBounds3);

    // Register the face cluster in mesh
//    tMesh->register_face_cluster(tFacetCluster);

//    delete tMesh;

    }
}
}

