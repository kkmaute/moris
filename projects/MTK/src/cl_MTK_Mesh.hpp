/*
 * cl_MTK_Mesh.hpp
 *
 *  Created on: Jul 24, 2018
 *      Author: messe
 */

#ifndef SRC_MESH_CL_MTK_MESH_HPP_
#define SRC_MESH_CL_MTK_MESH_HPP_

#include <MTK/src/cl_MTK_Blockset.hpp> //MTK/src
#include "typedefs.hpp" //MRS/COR/src

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
            get_number_of_blocksets() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns a pointer to a block
             */
            virtual Blockset *
            get_blockset_by_index( const uint& aIndex ) = 0;

//------------------------------------------------------------------------------

            /**
             * popules the member variables of the relevant nodes
             * with their T-Matrices
             */
            virtual void
            finalize() = 0;

//------------------------------------------------------------------------------

            /**
             * provides a moris::Mat<uint> containing the IDs this mesh has
             * to communicate with
             */
            virtual Mat< uint >
            get_communication_table() const = 0;

        };
//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */



#endif /* SRC_MESH_CL_MTK_MESH_HPP_ */
