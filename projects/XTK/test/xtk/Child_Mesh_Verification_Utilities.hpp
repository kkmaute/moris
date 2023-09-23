/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Child_Mesh_Verification_Utilities.hpp
 *
 */

#ifndef PROJECTS_XTK_TEST_XTK_CHILD_MESH_VERIFICATION_UTILITIES_HPP_
#define PROJECTS_XTK_TEST_XTK_CHILD_MESH_VERIFICATION_UTILITIES_HPP_

#include "cl_Matrix.hpp"
#include "fn_norm.hpp"
#include "op_minus.hpp"
#include "cl_XTK_Child_Mesh.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_XTK_Background_Mesh.hpp"
#include "fn_equal_to.hpp"
using namespace moris;
namespace xtk
{
bool
node_is_on_face(Matrix<DDRMat> const & aFaceCoords,
                Row_View_Real  const & aNodeCoord)
{
    // compute bounding box
    Matrix<DDRMat> tFaceXCoords = aFaceCoords.get_column(0);
    Matrix<DDRMat> tFaceYCoords = aFaceCoords.get_column(1);
    Matrix<DDRMat> tFaceZCoords = aFaceCoords.get_column(2);
   moris:: real tXHigh = tFaceXCoords.max();
    moris::real tXLow  = tFaceXCoords.min();
    moris::real tYHigh = tFaceYCoords.max();
    moris::real tYLow  = tFaceYCoords.min();
   moris:: real tZHigh = tFaceZCoords.max();
    moris::real tZLow  = tFaceZCoords.min();

    // check x coordinate of node is within bounding box
    bool tXInBox = false;
    if(aNodeCoord(0) <= tXHigh && aNodeCoord(0) >= tXLow )
    {
        tXInBox = true;
    }
    if(moris::equal_to(aNodeCoord(0),tXHigh) ||moris:: equal_to(aNodeCoord(0),tXLow) )
    {
        tXInBox = true;
    }

    // check y coordinate of node is within bounding box
    bool tYInBox = false;
    if(aNodeCoord(1) < tYHigh && aNodeCoord(1) > tYLow )
    {
        tYInBox = true;
    }
    if(moris::equal_to(aNodeCoord(1),tYHigh) || moris::equal_to(aNodeCoord(1),tYLow) )
    {
        tYInBox = true;
    }

    // check z coordinate of node is within bounding box
    bool tZInBox = false;
    if(aNodeCoord(2) < tZHigh && aNodeCoord(2) >tZLow )
    {
        tZInBox = true;
    }

    if(moris::equal_to(aNodeCoord(2),tZHigh) ||moris::equal_to(aNodeCoord(2),tZLow) )
    {
        tZInBox = true;
    }

    return (tXInBox && tYInBox && tZInBox);
}

void
verify_child_mesh_ancestry(Background_Mesh const & aBackgroundMeshData,
                           Cut_Mesh const &  tCutMesh)
{

    moris::Matrix< moris::DDRMat > tNodeCoordinates = aBackgroundMeshData.get_all_node_coordinates_loc_inds();

    moris::mtk::Mesh const & tMeshData = aBackgroundMeshData.get_mesh_data();

    // Iterate over child meshes
    for(size_t iCM = 0; iCM < tCutMesh.get_num_child_meshes(); iCM++)
    {
        // Get reference to child mesh
        Child_Mesh const & tChildMesh = tCutMesh.get_child_mesh(iCM);
        // verify edge ancestry
        Matrix<IndexMat> const & tEdgeToNode = tChildMesh.get_edge_to_node();
        moris::Matrix< moris::IndexMat > const & tEdgeParentIndices = tChildMesh.get_edge_parent_inds();
        moris::Matrix< moris::DDSTMat >  const & tEdgeParentRanks   = tChildMesh.get_edge_parent_ranks();

        for(moris::uint i = 0; i <tChildMesh.get_num_entities(mtk::EntityRank::EDGE); i++)
        {
            if(tEdgeParentRanks(i) == 2)
            {
                moris_index tFaceIndex = tEdgeParentIndices(i);
                Matrix<IndexMat> tFaceNodes = tMeshData.get_entity_connected_to_entity_loc_inds(tFaceIndex, mtk::EntityRank::FACE, mtk::EntityRank::NODE);
                Matrix<DDRMat>   tFaceNodeCoords(4,3);
                tFaceNodeCoords.get_row(0) = tMeshData.get_node_coordinate(tFaceNodes(0)).get_row(0);
                tFaceNodeCoords.get_row(1) = tMeshData.get_node_coordinate(tFaceNodes(1)).get_row(0);
                tFaceNodeCoords.get_row(2) = tMeshData.get_node_coordinate(tFaceNodes(2)).get_row(0);
                tFaceNodeCoords.get_row(3) = tMeshData.get_node_coordinate(tFaceNodes(3)).get_row(0);

                Row_View_Real const & tEdgeNodeCoord0 = tNodeCoordinates.get_row(tEdgeToNode(i,0));
                Row_View_Real const & tEdgeNodeCoord1 = tNodeCoordinates.get_row(tEdgeToNode(i,1));

                CHECK((node_is_on_face(tFaceNodeCoords,tEdgeNodeCoord0) && node_is_on_face(tFaceNodeCoords,tEdgeNodeCoord1)));

            }
        }
    }
}
}

#endif /* PROJECTS_XTK_TEST_XTK_CHILD_MESH_VERIFICATION_UTILITIES_HPP_ */

