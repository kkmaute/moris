/*
 * cl_MTK_CELL_INFO.hpp
 *
 *  Created on: Jul 25, 2019
 *      Author: ryan
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_CELL_INFO_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_CELL_INFO_HPP_

#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "assert.hpp"
#include "cl_MTK_Enums.hpp"
namespace moris
{
    namespace mtk
    {

        class Cell;

        class Cell_Info
        {
            public:
                Cell_Info(){}

                virtual
                ~Cell_Info(){}
                // ---------------------------------------------------------------------------------
                virtual
                enum Geometry_Type
                get_cell_geometry() const = 0;
                // ---------------------------------------------------------------------------------
                virtual
                enum Interpolation_Order
                get_cell_interpolation_order() const = 0;
                // ---------------------------------------------------------------------------------
                virtual
                uint
                get_num_verts() const  = 0;
                // ---------------------------------------------------------------------------------
                virtual
                uint
                get_num_verts_per_facet() const = 0;
                // ---------------------------------------------------------------------------------
                virtual
                uint
                get_num_facets() const = 0;
                // ---------------------------------------------------------------------------------
                virtual
                moris::Matrix<moris::IndexMat>
                get_node_to_face_map() const = 0;
                // ---------------------------------------------------------------------------------
                virtual
                moris::Matrix<moris::IndexMat>
                get_node_to_face_map(moris::uint aSideOrdinal) const  = 0;
                // ---------------------------------------------------------------------------------
                virtual
                moris::Matrix<moris::IndexMat>
                get_node_to_facet_map() const = 0;
                // ---------------------------------------------------------------------------------
                virtual
                moris::Matrix<moris::IndexMat>
                get_geometric_node_to_facet_map() const
                {
                    MORIS_ERROR(0,"ERROR");
                    return moris::Matrix<moris::IndexMat>();
                }

                virtual
                moris::Matrix<moris::IndexMat>
                get_geometric_node_to_facet_map(moris::uint aSideOrdinal) const
                {
                    MORIS_ERROR(0,"ERROR");
                    return moris::Matrix<moris::IndexMat>();
                }

                // ---------------------------------------------------------------------------------
                virtual
                moris::Matrix<moris::IndexMat>
                get_node_to_facet_map(moris::uint aSideOrdinal) const = 0;
                // ---------------------------------------------------------------------------------
                virtual
                moris::Matrix<moris::IndexMat>
                get_node_to_edge_map() const  = 0;
                // ---------------------------------------------------------------------------------
                virtual
                moris::Matrix<moris::IndexMat>
                get_node_to_edge_map(moris::uint aEdgeOrdinal) const  = 0;
                // ---------------------------------------------------------------------------------
                virtual
                moris::Matrix<moris::IndexMat>
                get_node_map_outward_normal(moris::uint aSideOrdinal) const  = 0;

                // ---------------------------------------------------------------------------------
                /*
                 * Returns the adjacent facet of the given
                 */
                virtual
                moris::uint
                get_adjacent_side_ordinal(moris::uint aSideOrdinal) const
                {
                    MORIS_ERROR(0,"get_adjacent_side_ordinal only makes sense for hex/quad type cells");
                    return MORIS_UINT_MAX;
                }

                // ---------------------------------------------------------------------------------
                virtual
                Matrix<DDRMat>
                get_loc_coord_on_side_ordinal(moris::uint aSideOrdinal) const
                {
                    MORIS_ERROR(0,"get_loc_coord_on_side_ordinal not implemented for given cell info type");
                    return Matrix<DDRMat>(0,0);
                }
                // ---------------------------------------------------------------------------------
                virtual
                moris::Matrix<moris::IndexMat>
                get_edge_to_face_map() const
                {
                    MORIS_ERROR(0,"get edge to face map not implemented for this CELL_INFO, this CELL_INFO currently used in XTK only for tets/tris");
                    return moris::Matrix<moris::IndexMat>(0,0);
                }
                // ---------------------------------------------------------------------------------
                /*!
                 * Compute the volume of 3D cell or the surface area of 2d cell
                 */
                virtual
                moris::real
                compute_cell_size( moris::mtk::Cell const * aCell ) const = 0;
                // ---------------------------------------------------------------------------------
                /*!
                 * Compute the side surface area of 3D cell or the side length of 2d cell
                 */
                virtual
                moris::real
                compute_cell_side_size( moris::mtk::Cell const * aCell ,
                        moris_index const & aSideOrd) const = 0;
                // ---------------------------------------------------------------------------------

                virtual
                void
                eval_N( const Matrix< DDRMat > & aXi,
                        Matrix< DDRMat > & aNXi ) const
                {
                    MORIS_ERROR(0,"eval_N not implemented for this type of cell info");
                }

        };

    }

}



#endif /* PROJECTS_MTK_SRC_CL_MTK_CELL_INFO_HPP_ */
