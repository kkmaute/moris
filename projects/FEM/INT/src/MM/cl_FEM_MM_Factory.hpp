/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_MM_Factory.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_MM_FACTORY_HPP_
#define SRC_FEM_CL_FEM_MM_FACTORY_HPP_

#include "assert.hpp"

#include "cl_Matrix.hpp"
#include "cl_FEM_Material_Model.hpp" //FEM/INT/src

//------------------------------------------------------------------------------
namespace moris::fem
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
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_MM_FACTORY_HPP_ */

