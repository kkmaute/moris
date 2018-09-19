/*
 * cl_HMR_Field.hpp
 *
 *  Created on: Sep 13, 2018
 *      Author: messe
 */

#ifndef PROJECTS_HMR_SRC_CL_HMR_FIELD_HPP_
#define PROJECTS_HMR_SRC_CL_HMR_FIELD_HPP_

#include "typedefs.hpp"
#include "cl_MTK_Field.hpp"
#include "cl_HMR_Block.hpp"

namespace moris
{
    namespace hmr
    {
//------------------------------------------------------------------------------

        class Field : public mtk::Field
        {

//------------------------------------------------------------------------------
        public :
//------------------------------------------------------------------------------

            Field(  const std::string & aLabel,
                          hmr::Block *  aBlock );

//------------------------------------------------------------------------------

            ~Field();

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */



#endif /* PROJECTS_HMR_SRC_CL_HMR_FIELD_HPP_ */
