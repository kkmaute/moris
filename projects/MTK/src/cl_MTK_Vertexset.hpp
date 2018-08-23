/*
 * cl_MTK_Vertexset.hpp
 *
 *  Created on: Aug 22, 2018
 *      Author: messe
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_VERTEXSET_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_VERTEXSET_HPP_

#include <string>

#include "typedefs.hpp"       //MRS/COR/src
#include "cl_MTK_Vertex.hpp"  //MTK/src
//------------------------------------------------------------------------------
namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------
        class Vertexset
        {
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Vertexset(){};

//------------------------------------------------------------------------------

            /**
             * Destructor, virtual
             */
            virtual
            ~Vertexset(){};

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
             * return the ID of this set
             */
            virtual uint
            get_id() = 0;

//------------------------------------------------------------------------------

            /**
             * returns the number of vertices in this set
             */
            virtual luint
            get_number_of_vertices() const = 0;

//------------------------------------------------------------------------------

            virtual Vertex *
            get_vertex( const luint & aVertexIndex ) const = 0;

//------------------------------------------------------------------------------
        } /* namespace mtk */
} /* namespace moris */

#endif /* PROJECTS_MTK_SRC_CL_MTK_VERTEXSET_HPP_ */
