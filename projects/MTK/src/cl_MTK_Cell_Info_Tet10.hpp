/*
 * cl_MTK_Cell_Info_Tet10.hpp
 *
 *  Created on: Oct 19, 2021
 *      Author: wunsch
 */

#ifndef PROJECTS_MTK_SRC_MTK_CL_MTK_CELL_INFO_TET10_HPP_
#define PROJECTS_MTK_SRC_MTK_CL_MTK_CELL_INFO_TET10_HPP_

#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "fn_det.hpp"

namespace moris
{
    namespace mtk
    {

        class Cell;

        class Cell_Info_Tet10 : public mtk::Cell_Info
        {
        public:
            enum Geometry_Type
            get_cell_geometry() const;

            // ----------------------------------------------------------------------------------

            enum Interpolation_Order
            get_cell_interpolation_order() const;

            //-----------------------------------------------------------------------------

            enum Integration_Order
            get_cell_integration_order() const;

            //-----------------------------------------------------------------------------

            enum CellShape
            compute_cell_shape(moris::mtk::Cell const *aCell) const;

            // ----------------------------------------------------------------------------------

            uint
            get_num_verts() const;

            // ----------------------------------------------------------------------------------

            uint
            get_num_facets() const;

            // ----------------------------------------------------------------------------------

            uint
            get_num_edges() const;

            // ----------------------------------------------------------------------------------

            uint
            get_num_verts_per_facet() const;

            // ----------------------------------------------------------------------------------

            uint
            get_loc_coord_dim() const;

            // ----------------------------------------------------------------------------------

            moris::Matrix<moris::IndexMat>
            get_node_to_face_map() const;

            // ----------------------------------------------------------------------------------

            moris::Matrix<moris::IndexMat>
            get_node_to_facet_map() const;

            // ----------------------------------------------------------------------------------

            moris::Matrix<moris::IndexMat>
            get_node_to_facet_map(moris::uint aSideOrdinal) const;

            // ----------------------------------------------------------------------------------

            moris::Matrix<moris::IndexMat>
            get_node_to_edge_map() const;

            // ----------------------------------------------------------------------------------

            moris::Matrix<moris::IndexMat>
            get_node_to_edge_map(moris::uint aEdgeOrdinal) const;

            // ----------------------------------------------------------------------------------

            moris::Matrix<moris::IndexMat>
            get_node_to_face_map(moris::uint aSideOrdinal) const;

            // ----------------------------------------------------------------------------------

            moris::Matrix<moris::IndexMat>
            get_node_map_outward_normal(moris::uint aSideOrdinal) const;

            // ----------------------------------------------------------------------------------

            moris::Matrix< moris::DDRMat >
            get_vertex_loc_coord(moris_index const & aVertexOrdinal) const;

            // ----------------------------------------------------------------------------------

            moris::Matrix<moris::IndexMat>
            get_edge_to_face_map() const;

            // ----------------------------------------------------------------------------------

            moris_index
            get_shared_vertex_ordinal_between_edges(
                moris_index aEdgeOrdinal1,
                moris_index aEdgeOrdinal2 ) const;

            moris::real
            compute_cell_size_special(moris::mtk::Cell const *aCell) const;

            // ----------------------------------------------------------------------------------

            /**
             * Computes the cell size if this isn't a rectangular cell
             * @param[in] aCell          MTK cell to compute size of.
             * 
             * @return return the cell size.
            */
            moris::real
            compute_cell_size_straight(moris::mtk::Cell const *aCell) const;

            // ----------------------------------------------------------------------------------

            /**
             * Computes the cell size derivative wrt to a single dof
             * @param[in] aCell           MTK cell to compute size of.
             * @param[in] aLocalVertexID  Local ID of vertex to use (0, 1, 2, or 3).
             * @param[in] aDirection      Direction to take derivative (0,1, or 2).
             * 
             * @return return the cell size.
            */
            moris::real
            compute_cell_size_deriv(moris::mtk::Cell const *aCell, uint aLocalVertexID, uint aDirection) const;

            // ----------------------------------------------------------------------------------

            moris::real
            compute_cell_side_size(moris::mtk::Cell const *aCell,
                                   moris_index const &aSideOrd) const;

            // ----------------------------------------------------------------------------------

            moris::real
            compute_cell_side_size_deriv(
                    moris::mtk::Cell const *aCell,
                    moris_index const      &aSideOrd,
                    uint                    aLocalVertexID,
                    uint                    aDirection ) const;

            // ----------------------------------------------------------------------------------

        };
    } // namespace mtk
} // namespace moris

#endif /* PROJECTS_MTK_SRC_MTK_CL_MTK_CELL_INFO_TET4_HPP_ */
