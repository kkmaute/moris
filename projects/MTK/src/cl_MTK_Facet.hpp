/*
 * cl_MTK_Face.hpp
 *
 *  Created on: Sep 23, 2018
 *      Author: messe
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_FACET_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_FACET_HPP_

#include "cl_MTK_Cell.hpp"

namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------

        class Facet : public mtk::Cell
        {
//------------------------------------------------------------------------------
        public :
//------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Facet(){};

//------------------------------------------------------------------------------

            /**
             * trivial destructor
             */
            virtual ~Facet(){};

//------------------------------------------------------------------------------

            virtual
            mtk::Cell *
            get_master() = 0;

//------------------------------------------------------------------------------

            virtual
            const mtk::Cell *
            get_master() const = 0;

//------------------------------------------------------------------------------

            virtual
            mtk::Cell *
            get_slave() = 0;

//------------------------------------------------------------------------------

            virtual
            const mtk::Cell *
            get_slave() const = 0;


//-----------------------------------------------------------------------------

            // fixme add active flag
        };

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */



#endif /* PROJECTS_MTK_SRC_CL_MTK_FACET_HPP_ */
