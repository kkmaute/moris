/*
 * cl_MTK_Cell_Info_Hex8.cpp
 *
 *  Created on: Sep 26, 2019
 *      Author: doble
 */

#include "cl_MTK_Cell_Info_Hex8.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "fn_cross.hpp"

namespace moris
{
namespace mtk
{
// ----------------------------------------------------------------------------------
enum Geometry_Type
Cell_Info_Hex8::get_cell_geometry() const
{
    return Geometry_Type::HEX;
}
// ----------------------------------------------------------------------------------
enum Interpolation_Order
Cell_Info_Hex8::get_cell_interpolation_order() const
{
    return Interpolation_Order::LINEAR;
}
// ----------------------------------------------------------------------------------
uint
Cell_Info_Hex8::get_num_verts() const
{
    return 8;
}
// ----------------------------------------------------------------------------------
uint
Cell_Info_Hex8::get_num_facets() const
{
    return 6;
}
// ----------------------------------------------------------------------------------
uint
Cell_Info_Hex8::get_num_verts_per_facet() const
{
    return 4;
}

// ----------------------------------------------------------------------------------

uint
Cell_Info_Hex8::get_loc_coord_dim() const
{
    return 3;
}

// ----------------------------------------------------------------------------------

moris::Matrix<moris::IndexMat>
Cell_Info_Hex8::get_node_to_face_map() const
{
    return {{1, 5, 4, 0}, {1,2,6,5}, {3, 7, 6, 2}, {0,4,7,3}, {0,3,2,1}, {4,5,6,7}};
}
// ----------------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Hex8::get_node_to_edge_map() const
{
    return {{0,1}, {1,2}, {2,3}, {3,0}, {4,5}, {5,6}, {6,7}, {7,4}, {0,4}, {1,5}, {2,6}, {3,7}};
}
// ----------------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Hex8::get_node_to_face_map(moris::uint aSideOrdinal) const
{
    switch (aSideOrdinal)
    {
        case(0):{ return {{1, 5, 4, 0}}; break; }
        case(1):{ return {{1, 2, 6, 5}}; break; }
        case(2):{ return {{3, 7, 6, 2}}; break; }
        case(3):{ return {{0, 4, 7, 3}}; break; }
        case(4):{ return {{0, 3, 2, 1}}; break; }
        case(5):{ return {{4, 5, 6, 7}}; break; }
        default:
            MORIS_ERROR(0,"Invalid side ordinal specified");
            return moris::Matrix<moris::IndexMat>(0,0);
            break;
    }
}
// ----------------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Hex8::get_node_to_edge_map(moris::uint aEdgeOrdinal) const
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
// ----------------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Hex8::get_node_to_facet_map() const
{
    return this->get_node_to_face_map();
}
// ----------------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Hex8::get_node_to_facet_map(moris::uint aSideOrdinal) const
{
    return this->get_node_to_face_map(aSideOrdinal);
}

// ----------------------------------------------------------------------------------

moris::Matrix<moris::IndexMat>
Cell_Info_Hex8::get_geometric_node_to_facet_map() const
{
    return this->get_node_to_face_map();
}

// ----------------------------------------------------------------------------------

Matrix<DDRMat>
Cell_Info_Hex8::get_vertex_loc_coord(moris_index const & aVertexOrdinal) const
{
    switch (aVertexOrdinal)
    {
        case  0: { return {{ -1.0, -1.0, -1.0 }}; break; }
        case  1: { return {{ +1.0, -1.0, -1.0 }}; break; }
        case  2: { return {{ +1.0, +1.0, -1.0 }}; break; }
        case  3: { return {{ -1.0, +1.0, -1.0 }}; break; }
        case  4: { return {{ -1.0, -1.0, +1.0 }}; break; }
        case  5: { return {{ +1.0, -1.0, +1.0 }}; break; }
        case  6: { return {{ +1.0, +1.0, +1.0 }}; break; }
        case  7: { return {{ -1.0, +1.0, +1.0 }}; break; }
        default:
        {
            MORIS_ERROR(0,"Invalid vertex ordinal specified");
            return moris::Matrix<moris::DDRMat>(0,0);
            break;
        }
    }
}

// ----------------------------------------------------------------------------------

moris::Matrix<moris::IndexMat>
Cell_Info_Hex8::get_geometric_node_to_facet_map(moris::uint aSideOrdinal) const
{
    return this->get_node_to_face_map(aSideOrdinal);
}
// ----------------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Cell_Info_Hex8::get_node_map_outward_normal(moris::uint aSideOrdinal) const
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
// ----------------------------------------------------------------------------------
moris::real
Cell_Info_Hex8::compute_cell_size( moris::mtk::Cell const * aCell ) const
{
    moris::Cell< Vertex* > tVertices = aCell->get_vertex_pointers();

    Matrix<DDRMat> tNode0Coords = tVertices(0)->get_coords();
    Matrix<DDRMat> tNode6Coords = tVertices(6)->get_coords();

    real tLx = std::abs(tNode0Coords(0) - tNode6Coords(0));
    real tLy = std::abs(tNode0Coords(1) - tNode6Coords(1));
    real tLz = std::abs(tNode0Coords(2) - tNode6Coords(2));

    return tLx*tLy*tLz;
}

// ----------------------------------------------------------------------------------
moris::real
Cell_Info_Hex8::compute_cell_side_size( moris::mtk::Cell const * aCell ,
                                        moris_index const & aSideOrd) const
{
    moris::Cell< mtk::Vertex const* > tVertices = aCell->get_vertices_on_side_ordinal(aSideOrd);

    Matrix<DDRMat> tNodeCoords0 = tVertices(0)->get_coords();
    Matrix<DDRMat> tNodeCoords1 = tVertices(1)->get_coords();
    Matrix<DDRMat> tNodeCoords2 = tVertices(3)->get_coords();

    return norm( cross( tNodeCoords1 - tNodeCoords0, tNodeCoords2 - tNodeCoords0 ) );
}
// ----------------------------------------------------------------------------------
moris::uint
Cell_Info_Hex8::get_adjacent_side_ordinal(moris::uint aSideOrdinal) const
{
    switch (aSideOrdinal)
    {
        case(0):{ return 2; break; }
        case(1):{ return 3; break; }
        case(2):{ return 0; break; }
        case(3):{ return 1; break; }
        case(4):{ return 5; break; }
        case(5):{ return 4; break; }
        default:
        {
            MORIS_ERROR(0,"Invalid side ordinal specified");
            return MORIS_UINT_MAX;
            break;
        }
    }
}
// ----------------------------------------------------------------------------------
void
Cell_Info_Hex8::eval_N( const Matrix< DDRMat > & aXi,
                              Matrix< DDRMat > & aNXi ) const
{
    // make sure that input is correct
    MORIS_ASSERT( aXi.length() >= 3, "HEX8 - eval_N: aXi not allocated or hat wrong size." );

    // unpack xi and eta from input vector
    moris::real    xi = aXi( 0 );
    moris::real   eta = aXi( 1 );
    moris::real  zeta = aXi( 2 );

    // populate output matrix
    aNXi.set_size(1,8);
    aNXi( 0 ) =  - ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * 0.125;
    aNXi( 1 ) =    ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * 0.125;
    aNXi( 2 ) =  - ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * 0.125;
    aNXi( 3 ) =    ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * 0.125;
    aNXi( 4 ) =    ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * 0.125;
    aNXi( 5 ) =  - ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * 0.125;
    aNXi( 6 ) =    ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * 0.125;
    aNXi( 7 ) =  - ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * 0.125;
}

// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------

}
}
