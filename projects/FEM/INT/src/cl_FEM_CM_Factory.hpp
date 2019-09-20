/*
 * cl_FEM_CM_Factory.hpp
 *
 *  Created on: Sep 17, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_CM_FACTORY_HPP_
#define SRC_FEM_CL_FEM_CM_FACTORY_HPP_

#include "assert.hpp"

#include "cl_Matrix.hpp"
#include "cl_FEM_Constitutive_Model.hpp" //FEM/INT/src

namespace moris
{
//------------------------------------------------------------------------------
    namespace fem
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
        Constitutive_Model * create_CM( Constitutive_Type aConstitutiveType );

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_CM_FACTORY_HPP_ */
