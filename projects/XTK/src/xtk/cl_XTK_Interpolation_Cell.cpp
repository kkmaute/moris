/*
 * cl_XTK_Interpolation_Cell.cpp
 *
 *  Created on: Jul 23, 2019
 *      Author: doble
 */

#include "cl_XTK_Interpolation_Cell.hpp"
#include "cl_MTK_Cell_Info_Tet4.hpp"
#include "cl_MTK_Cell_Info_Hex8.hpp"
#include "cl_MTK_Cell_Info_Hex27.hpp"
#include "cl_MTK_Cell_Info_Hex64.hpp"
#include "fn_cross.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"
#include "op_div.hpp"

using namespace moris;

namespace xtk
{
moris::Cell<moris::mtk::Vertex const *>
Interpolation_Cell::get_vertices_on_side_ordinal(moris::moris_index aSideOrdinal) const
{
    moris::Cell< moris::mtk::Vertex* > tVertices = this->get_vertex_pointers();

    moris::Matrix<moris::IndexMat> tNodeOrdsOnSide = mCellInfo->get_node_to_facet_map(aSideOrdinal);

    moris::Cell<moris::mtk::Vertex const *> tVerticesOnSide(tNodeOrdsOnSide.numel());
    for(moris::uint i = 0; i < tNodeOrdsOnSide.numel(); i++)
    {
        tVerticesOnSide(i) = tVertices(tNodeOrdsOnSide(i));
    }
    return tVerticesOnSide;
}
moris::Matrix<moris::DDRMat>
Interpolation_Cell::compute_outward_side_normal(moris::moris_index aSideOrdinal) const
{
    // get the vertex coordinates
    moris::Matrix<moris::DDRMat> tVertexCoords = this->get_vertex_coords();

    // Get vector along these edges
    moris::Matrix<moris::DDRMat> tEdge0Vector(tVertexCoords.numel(),1);
    moris::Matrix<moris::DDRMat> tEdge1Vector(tVertexCoords.numel(),1);

    // Get the nodes which need to be used to compute normal
    moris::Matrix<moris::IndexMat> tEdgeNodesForNormal = mCellInfo->get_node_map_outward_normal(aSideOrdinal);

    // Get vector along these edges
    tEdge0Vector = moris::linalg_internal::trans(tVertexCoords.get_row(tEdgeNodesForNormal(1,0)) - tVertexCoords.get_row(tEdgeNodesForNormal(0,0)));
    tEdge1Vector = moris::linalg_internal::trans(tVertexCoords.get_row(tEdgeNodesForNormal(1,1)) - tVertexCoords.get_row(tEdgeNodesForNormal(0,1)));

    // Take the cross product to get the normal
    Matrix<DDRMat> tOutwardNormal = moris::cross(tEdge0Vector,tEdge1Vector);

    // Normalize
    Matrix<DDRMat> tUnitOutwardNormal = tOutwardNormal / moris::norm(tOutwardNormal);


    return tUnitOutwardNormal;
}

void Interpolation_Cell::set_id(moris_id aId)
{
    mCellId = aId;
}

}

