/*
 * cl_MTK_Cell_Info_Tet4.cpp
 *
 *  Created on: Sep 26, 2019
 *      Author: doble
 */


#include "cl_MTK_Cell_Info_Tet4.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"

namespace moris
{
namespace mtk
{
// ----------------------------------------------------------------------------------
enum Geometry_Type
Cell_Info_Tet4::get_cell_geometry() const
{
    return Geometry_Type::TET;
}
// ----------------------------------------------------------------------------------
enum Interpolation_Order
Cell_Info_Tet4::get_cell_interpolation_order() const
{
    return Interpolation_Order::LINEAR;
}
// ----------------------------------------------------------------------------------
uint
Cell_Info_Tet4::get_num_verts() const
{
    return 4;
}
// ----------------------------------------------------------------------------------
uint
Cell_Info_Tet4::get_num_facets() const
{
    return 4;
}
// ----------------------------------------------------------------------------------
uint
Cell_Info_Tet4::get_num_verts_per_facet() const
{
    return 3;
}
// ----------------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Tet4::get_node_to_face_map() const
{
    return {{ 0, 1, 3 },{ 1, 2, 3 }, { 0, 3, 2 }, { 0, 2, 1 }};
}
// ----------------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Tet4::get_node_to_facet_map() const
{
    return this->get_node_to_face_map();
}
// ----------------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Tet4::get_node_to_facet_map(moris::uint aSideOrdinal) const
{
    return this->get_node_to_face_map(aSideOrdinal);
}

// ----------------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Tet4::get_node_to_edge_map() const
{
    return {{0, 1},{1, 2},{0, 2},{0, 3},{1, 3},{2, 3}};
}
// ----------------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Tet4::get_node_to_edge_map(moris::uint aEdgeOrdinal) const
{
    switch (aEdgeOrdinal)
    {
        case(0):{ return {{0,1}}; break; }
        case(1):{ return {{1,2}}; break; }
        case(2):{ return {{0,2}}; break; }
        case(3):{ return {{0,3}}; break; }
        case(4):{ return {{1,3}}; break; }
        case(5):{ return {{2,3}}; break; }
        default:{ MORIS_ASSERT(0,"Invalid edge ordinal specified"); return moris::Matrix<moris::IndexMat>(0,0); break; }
    }
}
// ----------------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Tet4::get_node_to_face_map(moris::uint aSideOrdinal) const
{
    switch (aSideOrdinal)
    {
        case(0):{ return {{0, 1, 3}}; break; }
        case(1):{ return {{1, 2, 3}}; break; }
        case(2):{ return {{0, 3, 2}}; break; }
        case(3):{ return {{0, 2, 1}}; break; }
        default:{ MORIS_ASSERT(0,"Invalid side ordinal specified"); return moris::Matrix<moris::IndexMat>(0,0); break;}
    }
}
// ----------------------------------------------------------------------------------
inline
moris::Matrix<moris::IndexMat>
Cell_Info_Tet4::get_node_map_outward_normal(moris::uint aSideOrdinal) const
{
    switch (aSideOrdinal)
    {
        case(0):{ return {{0,1},{1,3}}; break; }
        case(1):{ return {{1,2},{2,3}}; break; }
        case(2):{ return {{0,3},{3,2}}; break; }
        case(3):{ return {{0,2},{2,1}}; break; }
        default:{ MORIS_ERROR(0,"Invalid side ordinal specified"); return moris::Matrix<moris::IndexMat>(0,0); break;}
    }
}
// ----------------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Tet4::get_edge_to_face_map() const
{
    return  {{0,3},{1,3},{2,3},{0,2},{0,1},{1,2}};
}
// ----------------------------------------------------------------------------------
moris::real
Cell_Info_Tet4::compute_cell_size( moris::mtk::Cell const * aCell ) const
{
    moris::Cell< Vertex* > tVertices = aCell->get_vertex_pointers();

    Matrix<DDRMat> tNode0Coords = tVertices(0)->get_coords();
    Matrix<DDRMat> tNode1Coords = tVertices(1)->get_coords();
    Matrix<DDRMat> tNode2Coords = tVertices(2)->get_coords();
    Matrix<DDRMat> tNode3Coords = tVertices(3)->get_coords();

    moris::Matrix< moris::DDRMat > J(4,4,1.0);

    J( 1, 0 ) = tNode0Coords(0); J( 1, 1 ) = tNode1Coords(0); J( 1, 2 ) = tNode2Coords(0); J( 1, 3 ) = tNode3Coords(0);
    J( 2, 0 ) = tNode0Coords(1); J( 2, 1 ) = tNode1Coords(1); J( 2, 2 ) = tNode2Coords(1); J( 2, 3 ) = tNode3Coords(1);
    J( 3, 0 ) = tNode0Coords(2); J( 3, 1 ) = tNode1Coords(2); J( 3, 2 ) = tNode2Coords(2); J( 3, 3 ) = tNode3Coords(2);

    moris::real tVolume = (moris::det(J))/6; // Volume = 1/6*det(J)
    return tVolume;
}
// ----------------------------------------------------------------------------------
}
}


