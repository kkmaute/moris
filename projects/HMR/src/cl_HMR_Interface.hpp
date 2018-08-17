/*
 * cl_HMR_Interface.hpp
 *
 *  Created on: Jul 24, 2018
 *      Author: messe
 */

#ifndef SRC_HMR_CL_HMR_INTERFACE_HPP_
#define SRC_HMR_CL_HMR_INTERFACE_HPP_

#include "cl_Cell.hpp" //CON/src
#include "cl_MTK_Mesh.hpp" //MTK/src
#include "cl_HMR_Block.hpp" //HMR/src
namespace moris
{
    namespace hmr
    {
//-------------------------------------------------------------------------------

        // forward declaration of HMR MEsh object
        class HMR;

//-------------------------------------------------------------------------------
        /**
         * \brief mesh interface class
         */
        class Interface : public mtk::Mesh
        {

            //! ref to hmr object
            HMR & mHMR;

            //! cell of blocks
            Cell< Block* > mBlocks;

            //! tells how many blocks exist on this mesh
            uint mNumberOfBlocks;
//-------------------------------------------------------------------------------
        public:
//-------------------------------------------------------------------------------

            /**
             * mesh constructor, to be called from HMR
             */
            Interface( HMR & aHMR );

//-------------------------------------------------------------------------------

            // destructor
            ~Interface();

//------------------------------------------------------------------------------

            /**
             * returns the number of blocks on this mesh
             */
            uint
            get_number_of_blocks() const;

//-------------------------------------------------------------------------------

            /**
             * returns a pointer to a block
             */
            Block *
            get_block_by_index( const uint& aIndex );

//-------------------------------------------------------------------------------


            /**
             * popules the member variables of the relevant nodes
             * with their T-Matrices
             */
            void
            finalize();

//-------------------------------------------------------------------------------

            /**
             * provides a moris::Mat<uint> containing the IDs this mesh has
             * to communicate with
             */
            Mat< uint >
            get_communication_table() const ;

//-------------------------------------------------------------------------------
        };

    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_INTERFACE_HPP_ */
