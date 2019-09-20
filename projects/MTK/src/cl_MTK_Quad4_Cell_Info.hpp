/*
 * cl_MTK_Quad4_CELL_INFO.hpp
 *
 *  Created on: Aug 5, 2019
 *      Author: ryan
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_QUAD4_CELL_INFO_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_QUAD4_CELL_INFO_HPP_

#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Cell_Info.hpp"



namespace moris
{
namespace mtk
{
class Cell_Info_Quad4 : public mtk::Cell_Info
{
public:

    enum Geometry_Type
    get_cell_geometry() const
    {
        return Geometry_Type::QUAD;
    }
    enum Interpolation_Order
    get_cell_interpolation_order() const
    {
        return Interpolation_Order::LINEAR;
    }
    uint
    get_num_verts() const
    {
        return 4;
    }
    uint
    get_num_facets() const
    {
        return 4;
    }
    uint
    get_num_verts_per_facet() const
    {
        return 2;
    }
    /*!
     * Node ordinal to face map
     */

    moris::Matrix<moris::IndexMat>
    get_node_to_face_map() const
    {
        MORIS_ERROR(0,"Elements have no faces in 2D. Check the MTK mesh class to get nodes connected to an element.");
        return moris::Matrix<moris::IndexMat>(0,0);
    }

    moris::Matrix<moris::IndexMat>
    get_node_to_edge_map() const
    {
        return {{0,1}, {1,2}, {2,3}, {3,0}};
    }

    moris::Matrix<moris::IndexMat>
    get_node_to_facet_map() const
    {
        return this->get_node_to_edge_map();
    }

    moris::Matrix<moris::IndexMat>
    get_node_to_face_map(moris::uint aSideOrdinal) const
    {
        MORIS_ERROR(0,"Elements have no faces in 2D. Check the MTK mesh class to get nodes connected to an element.");
        return moris::Matrix<moris::IndexMat>(0,0);
    }

    moris::Matrix<moris::IndexMat>
    get_node_to_edge_map(moris::uint aEdgeOrdinal) const
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

    moris::Matrix<moris::IndexMat>
    get_node_to_facet_map(moris::uint aSideOrdinal) const
    {
        return this->get_node_to_edge_map(aSideOrdinal);
    }


    moris::Matrix<moris::IndexMat>
    get_node_map_outward_normal(moris::uint aSideOrdinal) const
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
};
}
}

#endif /* PROJECTS_MTK_SRC_CL_MTK_QUAD4_CELL_INFO_HPP_ */


