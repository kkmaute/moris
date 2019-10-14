/*
 * cl_MTK_Cell_Info_Quad4.cpp
 *
 *  Created on: Sep 26, 2019
 *      Author: doble
 */


#include "cl_MTK_Cell_Info_Quad4.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "fn_det.hpp"
#include "fn_dot.hpp"
#include "fn_sum.hpp"
#include "fn_trans.hpp"
#include "op_div.hpp"
#include "op_times.hpp"
namespace moris
{
namespace mtk
{
// ----------------------------------------------------------------------------------
enum Geometry_Type
Cell_Info_Quad4::get_cell_geometry() const
{
    return Geometry_Type::QUAD;
}
// ----------------------------------------------------------------------------------
enum Interpolation_Order
Cell_Info_Quad4::get_cell_interpolation_order() const
{
    return Interpolation_Order::LINEAR;
}
// ----------------------------------------------------------------------------------
uint
Cell_Info_Quad4::get_num_verts() const
{
    return 4;
}
// ----------------------------------------------------------------------------------
uint
Cell_Info_Quad4::get_num_facets() const
{
    return 4;
}
// ----------------------------------------------------------------------------------
uint
Cell_Info_Quad4::get_num_verts_per_facet() const
{
    return 2;
}
// ----------------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Quad4::get_node_to_face_map() const
{
    MORIS_ERROR(0,"Elements have no faces in 2D. Check the MTK mesh class to get nodes connected to an element.");
    return moris::Matrix<moris::IndexMat>(0,0);
}
// ----------------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Quad4::get_node_to_edge_map() const
{
    return {{0,1}, {1,2}, {2,3}, {3,0}};
}
// ----------------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Quad4::get_node_to_facet_map() const
{
    return this->get_node_to_edge_map();
}
// ----------------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Quad4::get_node_to_face_map(moris::uint aSideOrdinal) const
{
    MORIS_ERROR(0,"Elements have no faces in 2D. Check the MTK mesh class to get nodes connected to an element.");
    return moris::Matrix<moris::IndexMat>(0,0);
}
// ----------------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Quad4::get_node_to_edge_map(moris::uint aEdgeOrdinal) const
{
    switch (aEdgeOrdinal)
    {
        case(0):{ return {{0, 1}}; break; }
        case(1):{ return {{1, 2}}; break; }
        case(2):{ return {{2, 3}}; break; }
        case(3):{ return {{3, 0}}; break; }
        default:
        {
            MORIS_ASSERT(0,"Invalid edge ordinal specified");
            return moris::Matrix<moris::IndexMat>(0,0);
            break;
        }
    }
}
// ----------------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Quad4::get_node_to_facet_map(moris::uint aSideOrdinal) const
{
    return this->get_node_to_edge_map(aSideOrdinal);
}
// ----------------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Quad4::get_node_map_outward_normal(moris::uint aSideOrdinal) const
{
    switch (aSideOrdinal)
    {
        case(0):{ return {{1,0}}; break; }
        case(1):{ return {{2,1}}; break; }
        case(2):{ return {{3,2}}; break; }
        case(3):{ return {{0,3}}; break; }
        default:
        {
            MORIS_ERROR(0,"Invalid side ordinal specified");
            return moris::Matrix<moris::IndexMat>(0,0);
            break;
        }
    }
}
// ----------------------------------------------------------------------------------
moris::real
Cell_Info_Quad4::compute_cell_size( moris::mtk::Cell const * aCell ) const
{
    moris::Cell< Vertex* > tVertices = aCell->get_vertex_pointers();

    Matrix<DDRMat> tNode1Coords0 = tVertices(0)->get_coords();
    Matrix<DDRMat> tNodeCoords2 = tVertices(2)->get_coords();

    real tLx = std::abs(tNode1Coords0(0) - tNodeCoords2(0));
    real tLy = std::abs(tNode1Coords0(1) - tNodeCoords2(1));

    return tLx*tLy;
}
// ----------------------------------------------------------------------------------

}
}


