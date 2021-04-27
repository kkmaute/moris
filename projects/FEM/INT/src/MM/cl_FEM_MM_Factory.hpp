/*
 * cl_FEM_MM_Factory.hpp
 *
 *  Created on: Feb 1, 2021
 *      Author: wunsch
 */

#ifndef SRC_FEM_CL_FEM_MM_FACTORY_HPP_
#define SRC_FEM_CL_FEM_MM_FACTORY_HPP_

#include "assert.hpp"

#include "cl_Matrix.hpp"
#include "cl_FEM_Material_Model.hpp" //FEM/INT/src

namespace moris
{
//------------------------------------------------------------------------------
    namespace fem
    {
//------------------------------------------------------------------------------

    /**
     * \brief CM factory
     */
    class MM_Factory
    {

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------
        /**
         * trivial constructor
         */
        MM_Factory(){};

//------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~MM_Factory(){};

//------------------------------------------------------------------------------
        /**
         * create a material model
         */
        std::shared_ptr< Material_Model > create_MM( Material_Type aMaterialType );

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_MM_FACTORY_HPP_ */
