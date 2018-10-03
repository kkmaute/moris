/*
 * cl_MTK_Element_Block.hpp
 *
 *  Created on: Jul 24, 2018
 *      Author: messe
 */

#ifndef SRC_MESH_CL_MTK_BLOCK_HPP_
#define SRC_MESH_CL_MTK_BLOCK_HPP_

#include <string>

#include "typedefs.hpp" //MRS/COR/src
#include "cl_Map.hpp"
#include "cl_MTK_Vertex.hpp" //MTK/src
#include "cl_MTK_Cell.hpp" //MTK/src
#include "cl_MTK_Field.hpp" //MTK/src

namespace moris
{
    namespace mtk
    {

//------------------------------------------------------------------------------
        class Block
        {
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Block(){};

//------------------------------------------------------------------------------

            /**
             * virtual destructor
             */
            virtual
            ~Block(){};

//------------------------------------------------------------------------------

            /**
             * return a label that describes the block
             */
            virtual std::string
            get_label() const;

//------------------------------------------------------------------------------

            /**
             * sets the name of a block
             */
            virtual void
            set_label( const std::string & aLabel );

//------------------------------------------------------------------------------

            /**
             * returns the Id of this block
             */
            virtual moris_id
            get_id() const;

//------------------------------------------------------------------------------

            /**
             * returns the number of vertices owned and shared on this proc
             */
            virtual uint
            get_number_of_vertices() const;

//------------------------------------------------------------------------------

            /**
             * returns the number of element owned by this proc
             */
            virtual uint
            get_number_of_cells() const;

//------------------------------------------------------------------------------

            /**
             * returns a pointer to a vertex
             */
            virtual Vertex *
            get_vertex_by_index( const moris_index & aIndex );

//------------------------------------------------------------------------------

            /**
             * returns a pointer to a vertex ( const version )
             */
            virtual const Vertex *
            get_vertex_by_index( const moris_index & aIndex ) const;

//------------------------------------------------------------------------------

            /**
             * returns a pointer to a cell
             */
            virtual Cell *
            get_cell_by_index( const moris_index & aIndex );

//------------------------------------------------------------------------------

            /**
             * returns a pointer to a cell ( const version )
             */
            virtual const Cell *
            get_cell_by_index( const moris_index & aIndex ) const;

//------------------------------------------------------------------------------

            virtual sint
            get_number_of_adofs_used_by_proc() const;

//------------------------------------------------------------------------------

            virtual void
            get_adof_map( map< moris_id, moris_index > & aAdofMap ) const;

//------------------------------------------------------------------------------

            /**
             * return the interpolation order of the cells on the block
             */
            virtual uint
            get_interpolation_order() const;

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
//------------------------------------------------------------------------------
#endif /* SRC_MESH_CL_MTK_BLOCK_HPP_ */
