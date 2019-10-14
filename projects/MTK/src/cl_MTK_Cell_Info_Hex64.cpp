/*
 * cl_MTK_Cell_Info_Hex64.cpp
 *
 *  Created on: Sep 26, 2019
 *      Author: doble
 */


#include "cl_MTK_Cell_Info_Hex64.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"

namespace moris
{
namespace mtk
{
//-----------------------------------------------------------------------------
enum Geometry_Type
Cell_Info_Hex64::get_cell_geometry() const
{
    return Geometry_Type::HEX;
}
//-----------------------------------------------------------------------------
enum Interpolation_Order
Cell_Info_Hex64::get_cell_interpolation_order() const
{
    return Interpolation_Order::CUBIC;
}
//-----------------------------------------------------------------------------
uint
Cell_Info_Hex64::get_num_verts() const
{
    return 64;
}
//-----------------------------------------------------------------------------
uint
Cell_Info_Hex64::get_num_facets() const
{
    return 6;
}
//-----------------------------------------------------------------------------
uint
Cell_Info_Hex64::get_num_verts_per_facet() const
{
    return 16;
}
//-----------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Hex64::get_node_to_face_map() const
{
    return {{ 0, 1, 5, 4,  8,  9, 16, 17, 25, 24, 13, 12, 36, 37, 38, 39 },
        { 1, 2, 6, 5, 14, 15, 20, 21, 29, 28, 17, 16, 44, 45, 46, 47 },
        { 2, 3, 7, 6, 18, 19, 22, 23, 31, 30, 21, 20, 48, 49, 50, 51 },
        { 0, 4, 7, 3, 12, 13, 26, 27, 23, 22, 11, 10, 40, 41, 42, 43 },
        { 0, 3, 2, 1, 10, 11, 19, 18, 15, 14,  9,  8, 32, 33, 34, 35 },
        { 4, 5, 6, 7, 24, 25, 28, 29, 30, 31, 27, 26, 52, 53, 54, 55 }};
}
//-----------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Hex64::get_node_to_face_map(moris::uint aSideOrdinal) const
{
    switch (aSideOrdinal)
    {
        case(0):{ return {{ 0, 1, 5, 4,  8,  9, 16, 17, 25, 24, 13, 12, 36, 37, 38, 39 }}; break; }
        case(1):{ return {{ 1, 2, 6, 5, 14, 15, 20, 21, 29, 28, 17, 16, 44, 45, 46, 47 }}; break; }
        case(2):{ return {{ 2, 3, 7, 6, 18, 19, 22, 23, 31, 30, 21, 20, 48, 49, 50, 51 }}; break; }
        case(3):{ return {{ 0, 4, 7, 3, 12, 13, 26, 27, 23, 22, 11, 10, 40, 41, 42, 43 }}; break; }
        case(4):{ return {{ 0, 3, 2, 1, 10, 11, 19, 18, 15, 14,  9,  8, 32, 33, 34, 35 }}; break; }
        case(5):{ return {{ 4, 5, 6, 7, 24, 25, 28, 29, 30, 31, 27, 26, 52, 53, 54, 55 }}; break; }
        default:
            MORIS_ERROR(0,"Invalid side ordinal specified");
            return moris::Matrix<moris::IndexMat>(0,0);
            break;
    }
}
//-----------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Hex64::get_node_to_edge_map() const
{
    return {{0,1}, {1,2}, {2,3}, {3,0}, {4,5}, {5,6}, {6,7}, {7,4}, {0,4}, {1,5}, {2,6}, {3,7}};
}
//-----------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Hex64::get_node_to_edge_map(moris::uint aEdgeOrdinal) const
 {
     switch (aEdgeOrdinal)
     {
         case( 0):{ return {{0, 1}}; break; }
         case( 1):{ return {{1, 2}}; break; }
         case( 2):{ return {{2, 3}}; break; }
         case( 3):{ return {{3, 0}}; break; }
         case( 4):{ return {{4, 5}}; break; }
         case( 5):{ return {{5, 6}}; break; }
         case( 6):{ return {{6, 7}}; break; }
         case( 7):{ return {{7, 4}}; break; }
         case( 8):{ return {{0, 4}}; break; }
         case( 9):{ return {{1, 5}}; break; }
         case(10):{ return {{2, 6}}; break; }
         case(11):{ return {{3, 7}}; break; }
         default:
             MORIS_ASSERT(0,"Invalid edge ordinal specified");
             return moris::Matrix<moris::IndexMat>(0,0);
             break;
     }
 }
//-----------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Hex64::get_node_to_facet_map() const
{
    return this->get_node_to_face_map();
}
//-----------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Hex64::get_node_to_facet_map(moris::uint aSideOrdinal) const
{
    return this->get_node_to_face_map(aSideOrdinal);
}
//-----------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Hex64::get_node_map_outward_normal(moris::uint aSideOrdinal)  const
{
    switch (aSideOrdinal)
    {
        case(0):{ return {{0,1},{1,5}}; break; }
        case(1):{ return {{1,2},{2,6}}; break; }
        case(2):{ return {{2,3},{3,7}}; break; }
        case(3):{ return {{7,3},{3,0}}; break; }
        case(4):{ return {{0,3},{3,2}}; break; }
        case(5):{ return {{4,5},{5,6}}; break; }
        default:
            MORIS_ERROR(0,"Invalid side ordinal specified");
            return moris::Matrix<moris::IndexMat>(0,0);
            break;
    }
}
//-----------------------------------------------------------------------------
moris::real
Cell_Info_Hex64::compute_cell_size( moris::mtk::Cell const * aCell ) const
{
    moris::Cell< Vertex* > tVertices = aCell->get_vertex_pointers();

    Matrix<DDRMat> tNode0Coords = tVertices(0)->get_coords();
    Matrix<DDRMat> tNode6Coords = tVertices(6)->get_coords();

    real tLx = std::abs(tNode0Coords(0) - tNode6Coords(0));
    real tLy = std::abs(tNode0Coords(1) - tNode6Coords(1));
    real tLz = std::abs(tNode0Coords(2) - tNode6Coords(2));

    return tLx*tLy*tLz;
}
//-----------------------------------------------------------------------------


}
}

