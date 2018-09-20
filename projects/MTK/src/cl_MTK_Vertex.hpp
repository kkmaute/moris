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
#include "cl_Matrix.hpp" //LNA/src
#include "linalg_typedefs.hpp"
#include "fn_assert.hpp"
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
             * returns a moris::Matrix with node coordinates
             */
            virtual Matrix< DDRMat >
            get_coords() const
            {
                MORIS_ERROR(0,"Function not implemented in base vertex");
                return Matrix < DDRMat >(0,0);
            }

//------------------------------------------------------------------------------

            /**
             * returns the domain wide id of this vertex
             */
            virtual moris_id
            get_id() const
            {
                MORIS_ERROR(0,"Function not implemented in base vertex");
                return 0;
            }

//------------------------------------------------------------------------------

            /**
             * returns the domain wide id of this vertex
             */
            virtual moris_index
            get_index() const
            {
                MORIS_ERROR(0,"Function not implemented in base vertex");
                return 0;
            }


            virtual moris_index
            get_owner() const = 0;

            virtual Vertex_Interpolation *
            get_interpolation() = 0;

            virtual const Vertex_Interpolation *
            get_interpolation() const = 0;


//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
    //------------------------------------------------------------------------------
#endif /* SRC_MESH_CL_MTK_VERTEX_HPP_ */
