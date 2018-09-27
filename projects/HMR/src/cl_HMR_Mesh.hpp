/*
 * cl_HMR_Mesh.hpp
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

        // forward declaration of HMR Field object
        class Field;

        // forward declaration of HMR Mesh object
        class HMR;

//-------------------------------------------------------------------------------
        /**
         * \brief mesh interface class
         */
        class Mesh : public mtk::Mesh
        {

            //! ref to hmr object
            HMR & mHMR;

            //! cell of blocks
            Cell< hmr::Block* > mBlocks;

            //! tells how many blocks exist on this mesh
            uint mNumberOfBlocks;

            //! describing label
            std::string mLabel;

//-------------------------------------------------------------------------------
        public:
//-------------------------------------------------------------------------------

            /**
             * mesh constructor, to be called from HMR
             */
            Mesh( HMR & aHMR, const uint & aActivationPattern );

//-------------------------------------------------------------------------------

            // destructor
            ~Mesh();

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
            mtk::Block *
            get_block_by_index( const moris_index & aIndex );

//-------------------------------------------------------------------------------

            /**
             * returns a pointer to a block ( const version )
             */
            const mtk::Block *
            get_block_by_index( const moris_index & aIndex ) const;

//-------------------------------------------------------------------------------

            /**
             * returns a pointer to a block ( const version )
             */
            hmr::Block *
            get_hmr_block_by_index( const moris_index & aIndex )
            {
                return mBlocks( aIndex );
            }

//-------------------------------------------------------------------------------

            /**
             * provides a moris::Matrix< DDUMat > containing the IDs this mesh has
             * to communicate with
             */
            Matrix< IdMat >
            get_communication_table() const ;

//-------------------------------------------------------------------------------
// Functions interited by Block
//-------------------------------------------------------------------------------

            /**
             * return a label that describes the block
             */
            std::string
            get_label() const
            {
                return mLabel;
            }

//-------------------------------------------------------------------------------

            /**
             * sets the name of a mesh
             */
            void
            set_label( const std::string & aLabel )
            {
                mLabel = aLabel;
            }

//-------------------------------------------------------------------------------

            mtk::Field *
            create_field( const std::string & aLabel );
        };

    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_INTERFACE_HPP_ */
