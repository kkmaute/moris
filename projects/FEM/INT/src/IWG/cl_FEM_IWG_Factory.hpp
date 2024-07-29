/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Factory.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_FACTORY_HPP_
#define SRC_FEM_CL_FEM_IWG_FACTORY_HPP_

#include "assert.hpp"

#include "cl_Matrix.hpp"
#include "cl_FEM_IWG.hpp"    //FEM/INT/src

namespace moris
{
    //------------------------------------------------------------------------------
    namespace fem
    {
        //------------------------------------------------------------------------------

        /**
         * \brief IWG factory
         */
        class IWG_Factory
        {

            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------
            /**
             * constructor
             */
            IWG_Factory() = default;

            //------------------------------------------------------------------------------
            /**
             * create IWGs
             */
            std::shared_ptr< IWG > create_IWG( IWG_Type aIWGType );
        };

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_FACTORY_HPP_ */
