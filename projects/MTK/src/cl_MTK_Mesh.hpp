/*
 * cl_MTK_Mesh.hpp
 *
 *  Created on: Jul 24, 2018
 *      Author: messe
 */

#ifndef SRC_MESH_CL_MTK_MESH_HPP_
#define SRC_MESH_CL_MTK_MESH_HPP_

#include "typedefs.hpp" //MRS/COR/src
#include "cl_MTK_Block.hpp" //MTK/src

namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------
        class Mesh
        {
//------------------------------------------------------------------------------
        public :
//------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Mesh(){};

//------------------------------------------------------------------------------

            /**
             * virtual destructor
             */
            virtual
            ~Mesh(){};

//------------------------------------------------------------------------------

            /**
             * returns the number of blocks on this mesh
             */
            virtual uint
            get_number_of_blocks() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns a pointer to a block
             */
            virtual Block *
            get_block_by_index( const uint& aIndex ) = 0;

//------------------------------------------------------------------------------

            /**
             * popules the member variables of the relevant nodes
             * with their T-Matrices
             */
            virtual void
            finalize() = 0;

//------------------------------------------------------------------------------

        };
//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */



#endif /* SRC_MESH_CL_MTK_MESH_HPP_ */
