/*
 * cl_MTK_HEX27_CELL_INFO.HPP
 *
 *  Created on: May 14, 2019
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_HEX27_CELL_INFO_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_HEX27_CELL_INFO_HPP_


#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_MTK_Cell_Info.hpp"

namespace moris
{
namespace mtk
{
class Cell_Info_Hex27 : public mtk::Cell_Info
{
public:
    enum Geometry_Type
    get_cell_geometry() const
    {
        return Geometry_Type::HEX;
    }

    enum Interpolation_Order
    get_cell_interpolation_order() const
    {
        return Interpolation_Order::QUADRATIC;
    }

    uint
    get_num_verts() const
    {
        return 27;
    }

    uint
    get_num_facets() const
    {
        return 6;
    }

    uint
    get_num_verts_per_facet() const
    {
        return 9;
    }
    //-----------------------------------------------------------------------------
    /*!
     * Node ordinal to face map
     */
    inline
    moris::Matrix<moris::IndexMat>
    get_node_to_face_map() const
    {
        return {{0, 1, 5, 4,  8, 13, 16, 12, 25},
                {1, 2, 6, 5,  9, 14, 17, 13, 24},
                {2, 3, 7, 6, 10, 15, 18, 14, 26},
                {0, 4, 7, 3, 12, 19, 15, 11, 23},
                {0, 3, 2, 1, 11, 10,  9,  8, 21},
                {4, 5, 6, 7, 16, 17, 18, 19, 22}};
    }

    moris::Matrix<moris::IndexMat>
    get_node_to_edge_map() const
    {
        return {{0,1}, {1,2}, {2,3}, {3,0}, {4,5}, {5,6}, {6,7}, {7,4}, {0,4}, {1,5}, {2,6}, {3,7}};
    }

    //-----------------------------------------------------------------------------
    inline
    moris::Matrix<moris::IndexMat>
    get_node_to_face_map(moris::uint aSideOrdinal) const
    {
        switch (aSideOrdinal)
        {
            case(0):{ return {{0, 1, 5, 4,  8, 13, 16, 12, 25}}; break; }
            case(1):{ return {{1, 2, 6, 5,  9, 14, 17, 13, 24}}; break; }
            case(2):{ return {{2, 3, 7, 6, 10, 15, 18, 14, 26}}; break; }
            case(3):{ return {{0, 4, 7, 3, 12, 19, 15, 11, 23}}; break; }
            case(4):{ return {{0, 3, 2, 1, 11, 10,  9,  8, 21}}; break; }
            case(5):{ return {{4, 5, 6, 7, 16, 17, 18, 19, 22}}; break; }
            default:
                MORIS_ERROR(0,"Invalid side ordinal specified");
                return moris::Matrix<moris::IndexMat>(0,0);
                break;
        }
    }
    //-----------------------------------------------------------------------------
    moris::Matrix<moris::IndexMat>
    get_node_to_edge_map(moris::uint aEdgeOrdinal) const
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
    get_node_to_facet_map() const
    {
        return this->get_node_to_face_map();
    }

    //-----------------------------------------------------------------------------
    moris::Matrix<moris::IndexMat>
    get_node_to_facet_map(moris::uint aSideOrdinal) const
    {
        return this->get_node_to_face_map(aSideOrdinal);
    }
    //-----------------------------------------------------------------------------
    inline
    moris::Matrix<moris::IndexMat>
    get_node_map_outward_normal(moris::uint aSideOrdinal) const
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

};
}
}
#endif /* PROJECTS_MTK_SRC_CL_MTK_Hex27_CELL_INFO_HPP_ */
