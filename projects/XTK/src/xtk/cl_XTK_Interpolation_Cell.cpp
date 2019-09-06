/*
 * cl_XTK_Interpolation_Cell.cpp
 *
 *  Created on: Jul 23, 2019
 *      Author: doble
 */

#include "cl_XTK_Interpolation_Cell.hpp"
#include "cl_MTK_Tetra4_Connectivity.hpp"
#include "cl_MTK_Hex8_Connectivity.hpp"
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
    mtk::Geometry_Type tGeomType = this->get_geometry_type();

    moris::Cell< mtk::Vertex* > tVertices = this->get_vertex_pointers();

    switch(tGeomType)
    {
        case(mtk::Geometry_Type::TET):
        {
            if(this->get_number_of_vertices() == 4)
            {
                MORIS_ASSERT(aSideOrdinal<4,"Side ordinal out of bounds for cell type tet");
                moris::Matrix<moris::IndexMat> tNodeOrdsOnSide = moris::Tetra4_Connectivity::get_node_to_face_map(aSideOrdinal);

                moris::Cell<moris::mtk::Vertex const *> tVerticesOnSide(3);
                tVerticesOnSide(0) = tVertices(tNodeOrdsOnSide(0));
                tVerticesOnSide(1) = tVertices(tNodeOrdsOnSide(1));
                tVerticesOnSide(2) = tVertices(tNodeOrdsOnSide(2));
                return tVerticesOnSide;
                break;
            }
            else
            {
                MORIS_ERROR(0,"Invalid geometry type, currently only supports hex8 and tet4");
                return moris::Cell<moris::mtk::Vertex const *>(0);
                break;
            }

        }
        case(mtk::Geometry_Type::HEX):
        {

            if(this->get_number_of_vertices() == 8)
            {
                MORIS_ASSERT(aSideOrdinal<6,"Side ordinal out of bounds for cell type hex");
                moris::Matrix<moris::IndexMat> tNodeOrdsOnSide = moris::Hex8::get_node_to_face_map(aSideOrdinal);

                moris::Cell<moris::mtk::Vertex const *> tVerticesOnSide(4);
                tVerticesOnSide(0) = tVertices(tNodeOrdsOnSide(0));
                tVerticesOnSide(1) = tVertices(tNodeOrdsOnSide(1));
                tVerticesOnSide(2) = tVertices(tNodeOrdsOnSide(2));
                tVerticesOnSide(3) = tVertices(tNodeOrdsOnSide(3));
                return tVerticesOnSide;
                break;
            }
            else
            {
                MORIS_ERROR(0,"Invalid geometry type, currently only supports hex8 and tet4");
                return moris::Cell<moris::mtk::Vertex const *>(0);
                break;
            }

        }

        default:
            MORIS_ERROR(0,"Invalid geometry type, currently only supports hex8 and tet4");
            return moris::Cell<moris::mtk::Vertex const *>(0);
            break;
    }
}
moris::Matrix<moris::DDRMat>
Interpolation_Cell::compute_outward_side_normal(moris::moris_index aSideOrdinal) const
  {

    enum mtk::Geometry_Type tGeomType = this->get_geometry_type();

    // Get vector along these edges
    moris::Matrix<moris::DDRMat> tEdge0Vector(3,1);
    moris::Matrix<moris::DDRMat> tEdge1Vector(3,1);


    switch(tGeomType)
    {
        case(mtk::Geometry_Type::HEX):
              {
            MORIS_ERROR(aSideOrdinal<6,"Side ordinal out of bounds.");

#ifdef DEBUG
            if(this->get_vertex_pointers().size() > 8)
            {
                MORIS_LOG_DEBUG("Warning: this normal computation only valid for flat facets. Ensure your higher order element has flat facets");
            }
#endif

            // get the vertex coordinates
            moris::Matrix<moris::DDRMat> tVertexCoords = this->get_vertex_coords();

            // Get the nodes which need to be used to compute normal
            moris::Matrix<moris::IndexMat> tEdgeNodesForNormal = Hex8::get_node_map_outward_normal(aSideOrdinal);

            // Get vector along these edges
            tEdge0Vector = moris::linalg_internal::trans(tVertexCoords.get_row(tEdgeNodesForNormal(1,0)) - tVertexCoords.get_row(tEdgeNodesForNormal(0,0)));
            tEdge1Vector = moris::linalg_internal::trans(tVertexCoords.get_row(tEdgeNodesForNormal(1,1)) - tVertexCoords.get_row(tEdgeNodesForNormal(0,1)));
            break;
              }
        case(mtk::Geometry_Type::TET):
              {
            MORIS_ERROR(aSideOrdinal<4,"Side ordinal out of bounds.");

#ifdef DEBUG
            if(this->get_vertex_pointers().size() > 4)
            {
                MORIS_LOG_DEBUG("Warning: this normal computation only valid for flat facets. Ensure your higher order element has flat facets");
            }
#endif

            // get the vertex coordinates
            moris::Matrix<moris::DDRMat> tVertexCoords = this->get_vertex_coords();

            // Get the nodes which need to be used to compute normal
            moris::Matrix<moris::IndexMat> tEdgeNodesForNormal = Tetra4_Connectivity::get_node_map_outward_normal(aSideOrdinal);

            // Get vector along these edges
            tEdge0Vector = moris::linalg_internal::trans(tVertexCoords.get_row(tEdgeNodesForNormal(1,0)) - tVertexCoords.get_row(tEdgeNodesForNormal(0,0)));
            tEdge1Vector = moris::linalg_internal::trans(tVertexCoords.get_row(tEdgeNodesForNormal(1,1)) - tVertexCoords.get_row(tEdgeNodesForNormal(0,1)));

            break;
              }

        default:
            MORIS_ERROR(0,"Only implemented for hex8 compute_outward_side_normal");
            break;
    }

    // Take the cross product to get the normal
    Matrix<DDRMat> tOutwardNormal = moris::cross(tEdge0Vector,tEdge1Vector);

    // Normalize
    Matrix<DDRMat> tUnitOutwardNormal = tOutwardNormal / moris::norm(tOutwardNormal);


    return tUnitOutwardNormal;

  }

}

