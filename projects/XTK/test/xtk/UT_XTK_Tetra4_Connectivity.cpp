/*
 * cl_XTK_Tetra4_Connectivity.cpp
 *
 *  Created on: Jan 16, 2019
 *      Author: doble
 */


#include "catch.hpp"

#include "cl_Matrix.hpp"
#include "fn_all_true.hpp"
#include "fn_cross.hpp"
#include "fn_norm.hpp"
#include "op_equal_equal.hpp"

// For outputting
#include "cl_MTK_Mesh.hpp" // MTK/src
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_Mesh_Factory.hpp"


#include "cl_XTK_Tetra4_Connectivity.hpp"

using namespace moris;
using namespace xtk;

TEST_CASE("Test Tetra 4 connectivity","[TETRA_4_CONN]")
{
    // Tetra4 coordinates
    Matrix< DDRMat > tNodeCoordinates(4,3);
    tNodeCoordinates(0,0) =  0.0;     tNodeCoordinates(0,1) =  0.0;     tNodeCoordinates(0,2) =  0.0;
    tNodeCoordinates(1,0) =  1.0;     tNodeCoordinates(1,1) =  0.0;     tNodeCoordinates(1,2) =  0.0;
    tNodeCoordinates(2,0) =  0.0;     tNodeCoordinates(2,1) =  1.0;     tNodeCoordinates(2,2) =  0.0;
    tNodeCoordinates(3,0) =  0.0;     tNodeCoordinates(3,1) =  0.0;     tNodeCoordinates(3,2) =  1.0;

    moris::Matrix< moris::IndexMat > tNodeIndex({{0,1,2,3}});
    moris::Matrix< moris::IdMat >   tNodeIds({{1,2,3,4}});

    // verify face connectivity map
    Matrix< IndexMat > tGoldFaceConnMap = {{0, 1, 3}, {2, 1, 3}, {0, 2, 3}, {0, 2, 1}};
    CHECK(all_true(tGoldFaceConnMap == Tetra4_Connectivity::get_node_to_face_map()));

    Matrix<F31RMat> tGoldOutwardNormal0 = {{0.0},{-1.0}, {0.0}};
    Matrix<F31RMat> tGoldOutwardNormal1 = {{+5.773502691896258e-01},
                                           {+5.773502691896258e-01},
                                           {+5.773502691896258e-01}};
    Matrix<F31RMat> tGoldOutwardNormal2 = {{-1.0}, {0.0}, { 0.0}};
    Matrix<F31RMat> tGoldOutwardNormal3 = {{ 0.0}, {0.0}, {-1.0}};

    // Compute outward normal of face 0
    Matrix<IndexMat> tEdgeNodesForNormal = Tetra4_Connectivity::get_node_map_outward_normal(0);
    Matrix<F31RMat> tEdge0Vector = tNodeCoordinates.get_row(tEdgeNodesForNormal(1,0)) - tNodeCoordinates.get_row(tEdgeNodesForNormal(0,0));
    Matrix<F31RMat> tEdge1Vector = tNodeCoordinates.get_row(tEdgeNodesForNormal(1,1)) - tNodeCoordinates.get_row(tEdgeNodesForNormal(0,1));
    // TODO: add cross product to linalg
    Matrix<F31RMat> tOutwardNormal0 = cross(tEdge0Vector,tEdge1Vector);



    // Compute outward normal of face 1
    tEdgeNodesForNormal = Tetra4_Connectivity::get_node_map_outward_normal(1);
    tEdge0Vector = tNodeCoordinates.get_row(tEdgeNodesForNormal(1,0)) - tNodeCoordinates.get_row(tEdgeNodesForNormal(0,0));
    tEdge1Vector = tNodeCoordinates.get_row(tEdgeNodesForNormal(1,1)) - tNodeCoordinates.get_row(tEdgeNodesForNormal(0,1));
    Matrix<F31RMat> tOutwardNormal1 = cross(tEdge0Vector,tEdge1Vector) / std::pow(3,0.5);


    // Compute outward normal of face 2
    tEdgeNodesForNormal = Tetra4_Connectivity::get_node_map_outward_normal(2);
    tEdge0Vector = tNodeCoordinates.get_row(tEdgeNodesForNormal(1,0)) - tNodeCoordinates.get_row(tEdgeNodesForNormal(0,0));
    tEdge1Vector = tNodeCoordinates.get_row(tEdgeNodesForNormal(1,1)) - tNodeCoordinates.get_row(tEdgeNodesForNormal(0,1));
    Matrix<F31RMat> tOutwardNormal2 = cross(tEdge0Vector,tEdge1Vector);


    // Compute outward normal of face 3
    tEdgeNodesForNormal = Tetra4_Connectivity::get_node_map_outward_normal(3);
    tEdge0Vector = tNodeCoordinates.get_row(tEdgeNodesForNormal(1,0)) - tNodeCoordinates.get_row(tEdgeNodesForNormal(0,0));
    tEdge1Vector = tNodeCoordinates.get_row(tEdgeNodesForNormal(1,1)) - tNodeCoordinates.get_row(tEdgeNodesForNormal(0,1));
    Matrix<F31RMat> tOutwardNormal3 = cross(tEdge0Vector,tEdge1Vector);



    CHECK(all_true(tOutwardNormal0 == tGoldOutwardNormal0));
    CHECK(all_true(tOutwardNormal1 == tGoldOutwardNormal1));
    CHECK(all_true(tOutwardNormal2 == tGoldOutwardNormal2));
    CHECK(all_true(tOutwardNormal3 == tGoldOutwardNormal3));


    // Things needed for mesh output
    // specify spatial dimension
    uint tSpatialDim = 3;
    moris::Matrix<moris::IdMat> tElemLocalToGlobal({{1}});


    mtk::MtkMeshData aMeshData;
    aMeshData.SpatialDim              = &tSpatialDim;
    aMeshData.ElemConn(0)             = &tNodeIds;
    aMeshData.NodeCoords              = &tNodeCoordinates;
    aMeshData.LocaltoGlobalElemMap(0) = &tElemLocalToGlobal;
    aMeshData.LocaltoGlobalNodeMap    = &tNodeIds;
    aMeshData.CreateAllEdgesAndFaces  = false;

    moris::mtk::Mesh* tMesh = create_mesh( MeshType::STK, aMeshData );
    std::string tMeshOutputFile = "./Tetra4_Conn.e";
    tMesh->create_output_mesh(tMeshOutputFile);

    delete tMesh;







}


