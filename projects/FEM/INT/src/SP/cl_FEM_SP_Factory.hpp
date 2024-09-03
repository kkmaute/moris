/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Factory.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_SP_FACTORY_HPP_
#define SRC_FEM_CL_FEM_SP_FACTORY_HPP_

#include "assert.hpp"
#include "cl_Matrix.hpp"
#include "cl_FEM_Stabilization_Parameter.hpp" //FEM/INT/src

    //------------------------------------------------------------------------------
namespace moris::fem
{
    //------------------------------------------------------------------------------

    /**
     * \brief PS factory
     */
    class SP_Factory
    {

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------
        /**
         * trivial constructor
         */
        SP_Factory(){};

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~SP_Factory(){};

        //------------------------------------------------------------------------------
        /**
         * create a stabilization parameter
         */
        std::shared_ptr< Stabilization_Parameter > create_SP( fem::Stabilization_Type aStabilizationType );

        //------------------------------------------------------------------------------
    };

    //------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_SP_FACTORY_HPP_ */

