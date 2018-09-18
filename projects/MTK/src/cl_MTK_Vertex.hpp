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
#include "cl_MTK_Vertex_Interpolation.hpp" //LNA/src

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
            virtual moris_id
            get_id() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns the domain wide id of this vertex
             */
            virtual moris_index
            get_index() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns the id of the proc that owns this vertex
             */
            virtual uint
            get_owner() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns a pointer to the interpolation object
             */
            virtual Vertex_Interpolation *
            get_interpolation() = 0;

//------------------------------------------------------------------------------

            /**
             * returns a pointer to the interpolation object ( const version )
             */
            virtual const Vertex_Interpolation *
            get_interpolation() const = 0;

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
    //------------------------------------------------------------------------------
#endif /* SRC_MESH_CL_MTK_VERTEX_HPP_ */
