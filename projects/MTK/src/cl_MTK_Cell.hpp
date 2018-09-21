/*
 * cl_MTK_Cell.hpp
 *
 *  Created on: Jul 23, 2018
 *      Author: messe
 */

#ifndef SRC_MESH_CL_MTK_CELL_HPP_
#define SRC_MESH_CL_MTK_CELL_HPP_

#include "typedefs.hpp" //MRS/COR/src
#include "cl_Cell.hpp" //MRS/CON/src
#include "cl_Matrix.hpp" //LNA/src
#include "linalg_typedefs.hpp"
#include "cl_MTK_Vertex.hpp" //MTK/src
#include "cl_MTK_Enums.hpp" //MTK/src

//------------------------------------------------------------------------------
namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------
        /**
         * \brief the mtk::Cell class provides the cell information that is
         * provided by the mesh.
         */

        class Cell
        {
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Cell(){};

//------------------------------------------------------------------------------

            /**
             * Destructor. Must be virtual.
             */
            virtual
            ~Cell(){};

//------------------------------------------------------------------------------

            /**
             * returns the domain wide id of the cell
             *
             * @return luint ID
             */
            virtual moris_id
            get_id() const = 0;

            /**
             * returns the domain wide id of the cell
             *
             * @return luint ID
             */
            virtual moris_index
            get_index() const = 0;

//------------------------------------------------------------------------------

            /**
             * tells how many vertices are connected to this cell
             */
            virtual uint
            get_number_of_vertices() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns the proc id of the owner of this cell
             * ( this information is needed for STK )
             */
            virtual moris_id
            get_owner() const = 0;

//------------------------------------------------------------------------------

            /**
             * fills a moris::cell with pointers to connected vertices
             */
            virtual moris::Cell< Vertex* >
            get_vertex_pointers() const = 0;


//------------------------------------------------------------------------------

            /**
             * returns a Mat with IDs of connected vertices
             */
            virtual Matrix< IdMat >
            get_vertex_ids() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns a Mat with indices of connected vertices
             */
            virtual Matrix< IndexMat >
            get_vertex_inds() const = 0;


//------------------------------------------------------------------------------

            /**
             * returns a Mat of dimension
             * < number of vertices * number of dimensions >
             */
            virtual Matrix< DDRMat >
            get_vertex_coords() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns an enum that defines the geometry type of the element
             */
            virtual Geometry_Type
            get_geometry_type() const = 0;

//------------------------------------------------------------------------------

            /**
             * T-Matrix is calculated for vertices if this flag is set
             */
            virtual void
            set_t_matrix_flag() = 0;

//------------------------------------------------------------------------------

            /**
             * T-Matrix is not calculated for vertices if this flag is not set
             */
            virtual void
            unset_t_matrix_flag() = 0;

//------------------------------------------------------------------------------

            /**
             * query if the T-Matrix is to be calculated
             */
            virtual bool
            get_t_matrix_flag() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns the order of the element
             */
            virtual Interpolation_Order
            get_interpolation_order() const = 0;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
//------------------------------------------------------------------------------

#endif /* SRC_MESH_CL_MTK_CELL_HPP_ */
