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

namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------
        class Blockset
        {
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Blockset(){};

//------------------------------------------------------------------------------

            /**
             * virtual destructor
             */
            virtual
            ~Blockset(){};

//------------------------------------------------------------------------------

            /**
             * return a label that describes the block
             */
            virtual std::string
            get_label() const = 0;

//------------------------------------------------------------------------------

            /**
             * sets the name of a block
             */
            virtual void
            set_label( const std::string & aLabel ) = 0;

//------------------------------------------------------------------------------

            /**
             * returns the Id of this block
             */
            virtual luint
            get_id() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns the number of vertices owned and shared on this proc
             */
            virtual luint
            get_number_of_vertices() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns the number of element owned by this proc
             */
            virtual luint
            get_number_of_cells() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns a pointer to a vertex
             */
            virtual Vertex *
            get_vertex_by_index( const luint & aIndex ) = 0;

//------------------------------------------------------------------------------

            /**
             * returns a pointer to a cell
             */
            virtual Cell *
            get_cell_by_index( const luint & aIndex ) = 0;

//------------------------------------------------------------------------------

            virtual sint
            get_number_of_adofs_used_by_proc() const = 0;

//------------------------------------------------------------------------------

            virtual void
            get_adof_map( map< moris_id, moris_index > & aAdofMap ) const = 0;

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
//------------------------------------------------------------------------------
#endif /* SRC_MESH_CL_MTK_BLOCK_HPP_ */
