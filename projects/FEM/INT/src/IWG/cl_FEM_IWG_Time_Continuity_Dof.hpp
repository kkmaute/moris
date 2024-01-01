/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Time_Continuity_Dof.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_TIME_CONTINUITY_DOF_HPP_
#define SRC_FEM_CL_FEM_IWG_TIME_CONTINUITY_DOF_HPP_

#include <map>
#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"     //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_IWG.hpp"    //FEM/INT/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IWG_Time_Continuity_Dof : public IWG
        {

            //------------------------------------------------------------------------------

          public:
            enum class IWG_Property_Type
            {
                WEIGHT_CURRENT,
                WEIGHT_PREVIOUS,
                WEIGHT_RESIDUAL,
                INITIAL_CONDITION,
                THICKNESS,
                MAX_ENUM
            };

            //------------------------------------------------------------------------------
            /*
             *  constructor
             */
            IWG_Time_Continuity_Dof();

            //------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IWG_Time_Continuity_Dof(){};

            //------------------------------------------------------------------------------
            /**
             * compute the residual
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void compute_residual( real aWStar );

            //------------------------------------------------------------------------------
            /**
             * compute the jacobian
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void compute_jacobian( real aWStar );

            //------------------------------------------------------------------------------
            /**
             * compute the residual and the jacobian
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void compute_jacobian_and_residual( real aWStar );

            //------------------------------------------------------------------------------
            /**
             * compute dR(n)/du(n-1)
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void compute_jacobian_previous( real aWStar );

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

#endif /* SRC_FEM_CL_FEM_IWG_TIME_CONTINUITY_DOF_HPP_ */
