/*
 * cl_FEM_SP_Factory.hpp
 *
 *  Created on: Nov 14, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_SP_FACTORY_HPP_
#define SRC_FEM_CL_FEM_SP_FACTORY_HPP_

#include "assert.hpp"
#include "cl_Matrix.hpp"
#include "cl_FEM_Stabilization_Parameter.hpp" //FEM/INT/src

namespace moris
{
    //------------------------------------------------------------------------------
    namespace fem
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
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SP_FACTORY_HPP_ */
