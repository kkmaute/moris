/*
 * cl_MTK_Cell_Info_Tri3.cpp
 *
 *  Created on: Sep 26, 2019
 *      Author: doble
 */

#include "cl_MTK_Cell_Info_Tri3.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "fn_det.hpp"
#include "fn_dot.hpp"
#include "fn_sum.hpp"
#include "fn_trans.hpp"
#include "op_div.hpp"
#include "op_times.hpp"
#include "fn_norm.hpp"

namespace moris
{
    namespace mtk
    {
        // ----------------------------------------------------------------------------------

        enum Geometry_Type
        Cell_Info_Tri3::get_cell_geometry() const
        {
            return Geometry_Type::TRI;
        }

        // ----------------------------------------------------------------------------------

        enum Interpolation_Order
        Cell_Info_Tri3::get_cell_interpolation_order() const
        {
            return Interpolation_Order::LINEAR;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Tri3::get_num_verts() const
        {
            return 3;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Tri3::get_num_facets() const
        {
            return 3;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Tri3::get_num_verts_per_facet() const
        {
            return 2;
        }

        // ----------------------------------------------------------------------------------

        uint
        Cell_Info_Tri3::get_loc_coord_dim() const
        {
            return 3;
        }
        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Tri3::get_node_to_face_map() const
        {
            MORIS_ERROR(0,"Elements have no faces in 2D. Check the MTK mesh class to get nodes connected to an element.");
            return moris::Matrix<moris::IndexMat>(0,0);
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Tri3::get_node_to_edge_map() const
        {
            return {{0,1}, {1,2}, {2,0}};
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Tri3::get_node_to_facet_map() const
        {
            return this->get_node_to_edge_map();
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Tri3::get_node_to_face_map(moris::uint aSideOrdinal) const
        {
            MORIS_ERROR(0,"Elements have no faces in 2D. Check the MTK mesh class to get nodes connected to an element.");
            return moris::Matrix<moris::IndexMat>(0,0);
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Tri3::get_node_to_edge_map(moris::uint aEdgeOrdinal) const
        {
            switch (aEdgeOrdinal)
            {
                case(0):{ return {{0, 1}}; break; }
                case(1):{ return {{1, 2}}; break; }
                case(2):{ return {{2, 0}}; break; }
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
        Cell_Info_Tri3::get_node_to_facet_map(moris::uint aSideOrdinal) const
        {
            return this->get_node_to_edge_map(aSideOrdinal);
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info_Tri3::get_node_map_outward_normal(moris::uint aSideOrdinal) const
        {
            switch (aSideOrdinal)
            {
                case(0):{ return {{1,0}}; break; }
                case(1):{ return {{2,1}}; break; }
                case(2):{ return {{0,2}}; break; }
                default:
                {
                    MORIS_ERROR(0,"Invalid side ordinal specified");
                    return moris::Matrix<moris::IndexMat>(0,0);
                    break;
                }
            }
        }

        // ----------------------------------------------------------------------------------

        moris::real Cell_Info_Tri3::compute_cell_size( moris::mtk::Cell const * aCell ) const
        {
            // cell coordinates
            moris::Cell< Vertex* > tVertices = aCell->get_vertex_pointers();

            const Matrix<DDRMat> tNodeCoords0 = tVertices(0)->get_coords();

            MORIS_ASSERT(tNodeCoords0.numel() == 2,"Cell_Info_Tri3::compute_cell_size only works in 2D.\n");

            const Matrix<DDRMat> tNodeCoords10 = tVertices(1)->get_coords() - tNodeCoords0;
            const Matrix<DDRMat> tNodeCoords20 = tVertices(2)->get_coords() - tNodeCoords0;

            real tArea = 0.5 * std::abs( tNodeCoords10(0)*tNodeCoords20(1) - tNodeCoords20(0)*tNodeCoords10(1) );

            return tArea;
        }

        // ----------------------------------------------------------------------------------

        moris::real
        Cell_Info_Tri3::compute_cell_side_size(
                moris::mtk::Cell const * aCell ,
                moris_index      const & aSideOrd) const
        {
            moris::Cell< mtk::Vertex const* > tVertices = aCell->get_vertices_on_side_ordinal(aSideOrd);

            Matrix<DDRMat> tLVec = tVertices(1)->get_coords() - tVertices(0)->get_coords();

            return moris::norm(tLVec);
        }
    }
}
