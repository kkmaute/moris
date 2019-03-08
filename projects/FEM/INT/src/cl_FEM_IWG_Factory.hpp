/*
 * cl_FEM_IWG_Factory.hpp
 *
 *  Created on: Feb 20, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_IWG_FACTORY_HPP_
#define SRC_FEM_CL_FEM_IWG_FACTORY_HPP_

#include "assert.hpp"

#include "cl_Matrix.hpp"
#include "cl_FEM_IWG.hpp" //FEM/INT/src

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
        IWG_Factory(){};

//------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IWG_Factory(){};

//------------------------------------------------------------------------------
        /**
         * create IWGs
         */
        IWG * create_IWGs( IWG_Type aIWGType );

    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_FACTORY_HPP_ */
