/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Cell_Info_Hex8.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_HEX8_CELL_INFO_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_HEX8_CELL_INFO_HPP_

#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Cell.hpp"

namespace moris
{
    namespace mtk
    {
        class Cell_Info_Hex8 : public Cell_Info
        {
          public:
            // ----------------------------------------------------------------------------------

            enum Geometry_Type
            get_cell_geometry() const override;

            // ----------------------------------------------------------------------------------

            enum CellTopology
            get_cell_topology() const override;

            // ----------------------------------------------------------------------------------

            enum Interpolation_Order
            get_cell_interpolation_order() const override;

            //-----------------------------------------------------------------------------

            enum Integration_Order
            get_cell_integration_order() const override;

            //-----------------------------------------------------------------------------

            enum CellShape
            compute_cell_shape( moris::mtk::Cell const *aCell ) const override;

            // ----------------------------------------------------------------------------------

            uint
            get_num_verts() const override;

            // ----------------------------------------------------------------------------------

            uint
            get_num_facets() const override;

            // ----------------------------------------------------------------------------------

            uint
            get_num_edges() const override;

            // ----------------------------------------------------------------------------------

            uint
            get_num_verts_per_facet() const override;

            // ----------------------------------------------------------------------------------

            uint
            get_loc_coord_dim() const override;

            //----------------------------------------------------------------------------

            moris::Matrix< moris::IndexMat >
            get_node_to_face_map() const override;

            // ----------------------------------------------------------------------------------

            moris::Matrix< moris::IndexMat >
            get_node_to_edge_map() const override;

            // ----------------------------------------------------------------------------------

            moris::Matrix< moris::IndexMat >
            get_node_to_face_map( moris::uint aSideOrdinal ) const override;

            // ----------------------------------------------------------------------------------

            moris::Matrix< moris::IndexMat >
            get_node_to_edge_map( moris::uint aEdgeOrdinal ) const override;

            // ----------------------------------------------------------------------------------

            moris::Matrix< moris::IndexMat >
            get_node_to_facet_map() const override;

            // ----------------------------------------------------------------------------------

            moris::Matrix< moris::IndexMat >
            get_node_to_facet_map( moris::uint aSideOrdinal ) const override;

            // ----------------------------------------------------------------------------------

            moris::Matrix< moris::IndexMat >
            get_geometric_node_to_facet_map() const override;

            //-----------------------------------------------------------------------------

            moris::Matrix< moris::IndexMat >
            get_geometric_node_to_facet_map( moris::uint aSideOrdinal ) const override;

            // ----------------------------------------------------------------------------------

            moris::uint
            get_adjacent_side_ordinal( moris::uint aSideOrdinal ) const override;

            // ----------------------------------------------------------------------------------

            Vector< moris_index >
            get_vertex_path_to_entity_rank_and_ordinal(
                    moris_index aVertexOrdinal,
                    moris_index aOtherEntityOrdinal,
                    moris_index aOtherEntityRank ) const override;

            Vector< moris_index >
            get_edge_path_to_entity_rank_and_ordinal(
                    moris_index aEdgeOrdinal,
                    moris_index aOtherEntityOrdinal,
                    moris_index aOtherEntityRank ) const override;

            bool
            is_entity_connected_to_facet(
                    moris_index aFacetOrdinal,
                    moris_index aOtherEntityOrdinal,
                    moris_index aOtherEntityRank ) const override;

            // ----------------------------------------------------------------------------------

            Matrix< DDRMat >
            get_vertex_loc_coord( moris_index const &aVertexOrdinal ) const override;

            // ----------------------------------------------------------------------------------

            moris::Matrix< moris::IndexMat >
            get_node_map_outward_normal( moris::uint aSideOrdinal ) const override;

            // ----------------------------------------------------------------------------------

            moris::real
            compute_cell_size_special( moris::mtk::Cell const *aCell ) const override;

            // ----------------------------------------------------------------------------------

            /**
             * Computes the cell size if this isn't a rectangular cell
             * @param[in] aCell          MTK cell to compute size of.
             *
             * @return return the cell size.
             */
            moris::real
            compute_cell_size_straight( moris::mtk::Cell const *aCell ) const override;

            // ----------------------------------------------------------------------------------
            moris::real
            compute_cell_side_size(
                    moris::mtk::Cell const *aCell,
                    moris_index const      &aSideOrd ) const override;
            // ----------------------------------------------------------------------------------

            void
            eval_N(
                    const Matrix< DDRMat > &aXi,
                    Matrix< DDRMat >       &aNXi ) const override;
        };
    }    // namespace mtk
}    // namespace moris
#endif /* PROJECTS_MTK_SRC_CL_MTK_HEX8_CELL_INFO_HPP_ */

