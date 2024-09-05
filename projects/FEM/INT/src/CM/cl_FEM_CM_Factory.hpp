/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Factory.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_CM_FACTORY_HPP_
#define SRC_FEM_CL_FEM_CM_FACTORY_HPP_

#include "assert.hpp"

#include "cl_Matrix.hpp"
#include "cl_FEM_Constitutive_Model.hpp" //FEM/INT/src

//------------------------------------------------------------------------------
namespace moris::fem
{
    //------------------------------------------------------------------------------

    /**
     * \brief CM factory
     */
    class CM_Factory
    {

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------
        /**
         * trivial constructor
         */
        CM_Factory(){};

//------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~CM_Factory(){};

//------------------------------------------------------------------------------
        /**
         * create a constitutive model
         */
        std::shared_ptr< Constitutive_Model > create_CM( Constitutive_Type aConstitutiveType );

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_CM_FACTORY_HPP_ */

