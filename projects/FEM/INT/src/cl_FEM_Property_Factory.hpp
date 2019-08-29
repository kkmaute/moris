/*
 * cl_FEM_Property_Factory.hpp
 *
 *  Created on: Jul 12, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_PROPERTY_FACTORY_HPP_
#define SRC_FEM_CL_FEM_PROPERTY_FACTORY_HPP_

#include "assert.hpp"
#include "cl_Matrix.hpp"
#include "cl_FEM_Property.hpp" //FEM/INT/src

namespace moris
{
//------------------------------------------------------------------------------
    namespace fem
    {
//------------------------------------------------------------------------------

    /**
     * \brief IWG factory
     */
    class Property_Factory
    {

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------
        /**
         * constructor
         */
        Property_Factory(){};

//------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~Property_Factory(){};

//------------------------------------------------------------------------------
        /**
         * create property
         */
        Property * create_property( Property_Type aPropertyType );

    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_Property_FACTORY_HPP_ */
