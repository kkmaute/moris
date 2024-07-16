/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Facet.hpp
 *
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

          public:
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

            virtual mtk::Cell *
            get_leader() = 0;

            //------------------------------------------------------------------------------

            virtual const mtk::Cell *
            get_leader() const = 0;

            //------------------------------------------------------------------------------

            virtual mtk::Cell *
            get_follower() = 0;

            //------------------------------------------------------------------------------

            virtual const mtk::Cell *
            get_follower() const = 0;

            //-----------------------------------------------------------------------------

            // fixme add active flag
        };

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */

#endif /* PROJECTS_MTK_SRC_CL_MTK_FACET_HPP_ */
