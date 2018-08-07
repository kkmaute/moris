/*
 * cl_MTK_Vertex.hpp
 *
 *  Created on: Jul 23, 2018
 *      Author: messe
 */

#ifndef SRC_MESH_CL_MTK_VERTEX_HPP_
#define SRC_MESH_CL_MTK_VERTEX_HPP_



#include "typedefs.hpp" //MRS/COR/src
#include "cl_Cell.hpp" //MRS/CON/src
#include "cl_Mat.hpp" //LNA/src

//------------------------------------------------------------------------------
namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------
        class Vertex
        {
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Vertex(){};

//------------------------------------------------------------------------------

            /**
             * Destructor, virtual
             */
            virtual
            ~Vertex(){};

//------------------------------------------------------------------------------

            /**
             * returns a moris::Mat with node coordinates
             */
            virtual Mat< real >
            get_coords() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns the domain wide id of this vertex
             */
            virtual luint
            get_id() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns the B-Spline IDs of this vertex
             */
            virtual Mat< luint >
            get_bspline_ids() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns the proc owners of the IDs of this vertex
             */
            virtual Mat< uint >
            get_bspline_owners() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns the B-Spline IDs of this vertex
             */
            virtual moris::Cell< Vertex* >
            get_bspline_pointers() = 0;

//------------------------------------------------------------------------------

            /**
             * returns the T-Matrix of this vertex
             */
            virtual Mat< real >
            get_t_matrix() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns the id of the proc that owns this vertex
             */
            virtual uint
            get_owner() const = 0;

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
    //------------------------------------------------------------------------------
#endif /* SRC_MESH_CL_MTK_VERTEX_HPP_ */
