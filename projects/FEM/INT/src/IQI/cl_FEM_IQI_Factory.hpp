/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Factory.hpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_FACTORY_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_FACTORY_HPP_

#include "assert.hpp"
#include "cl_Matrix.hpp"
#include "cl_FEM_IQI.hpp" //FEM/INT/src

namespace moris
{
    //------------------------------------------------------------------------------
    namespace fem
    {
        //------------------------------------------------------------------------------

        /**
         * IQI factory
         */
        class IQI_Factory
        {

                //------------------------------------------------------------------------------
            public:
                //------------------------------------------------------------------------------
                /**
                 * constructor
                 */
                IQI_Factory(){};

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IQI_Factory(){};

                //------------------------------------------------------------------------------
                /**
                 * create IQI
                 */
                std::shared_ptr< IQI > create_IQI( IQI_Type aIQIType );

        };

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_FACTORY_HPP_ */

