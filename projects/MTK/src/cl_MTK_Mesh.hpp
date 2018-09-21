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
        class Mesh : public Block
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
            get_block_by_index( const moris_index & aIndex ) = 0;

//------------------------------------------------------------------------------

            /**
             * returns a pointer to a block ( const version )
             */
            virtual const Block *
            get_block_by_index( const moris_index & aIndex ) const = 0;

//------------------------------------------------------------------------------

            //fixme: this function needs to go
            /**
             * populates the member variables of the relevant nodes
             * with their T-Matrices
             */
            virtual void
            finalize() = 0;

//------------------------------------------------------------------------------

            /**
             * provides a moris::Mat<uint> containing the IDs this mesh has
             * to communicate with
             */
            virtual Matrix< IdMat >
            get_communication_table() const = 0;

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */



#endif /* SRC_MESH_CL_MTK_MESH_HPP_ */
