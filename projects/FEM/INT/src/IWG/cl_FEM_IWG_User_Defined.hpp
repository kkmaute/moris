/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_USER_DEFINED.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_USER_DEFINED_HPP_
#define SRC_FEM_CL_FEM_IWG_USER_DEFINED_HPP_

#include "moris_typedefs.hpp"    //MRS/COR/src

#include "cl_Vector.hpp"    //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_IWG.hpp"                   //FEM/INT/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IWG_User_Defined : public IWG
        {
            enum class IWG_Property_Type
            {
                IWG_FORMULATION,    // IWG formulation
                THICKNESS,          // thickness
                MAX_ENUM
            };

            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------
            // constructor
            IWG_User_Defined();

            //------------------------------------------------------------------------------
            // trivial destructor
            ~IWG_User_Defined(){};

            //------------------------------------------------------------------------------
            /**
             * compute jacobian and residual
             */
            void compute_jacobian_and_residual( real aWStar );

            //------------------------------------------------------------------------------
            /**
             * compute jacobian
             */
            void compute_jacobian( real aWStar );

            //------------------------------------------------------------------------------
            /**
             * compute residual
             */
            void compute_residual( real aWStar );

            //------------------------------------------------------------------------------
            /**
             * compute the derivative of the residual wrt design variables
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void compute_dRdp( real aWStar );
            //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_USER_DEFINED_HPP_ */