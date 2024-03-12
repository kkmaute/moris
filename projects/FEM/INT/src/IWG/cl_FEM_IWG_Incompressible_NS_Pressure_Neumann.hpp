/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Incompressible_NS_Pressure_Neumann.hpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IWG_INCOMPRESSIBLE_NS_PRESSURE_NEUMANN_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IWG_INCOMPRESSIBLE_NS_PRESSURE_NEUMANN_HPP_

#include <map>
#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_IWG.hpp"    //FEM/INT/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IWG_Incompressible_NS_Pressure_Neumann : public IWG
        {

            //------------------------------------------------------------------------------

          public:
            // local property enums
            enum class IWG_Property_Type
            {
                PRESSURE,
                TOTAL_PRESSURE,
                DENSITY,
                BACKFLOW_PREVENTION,
                SELECT,
                MAX_ENUM
            };

            //------------------------------------------------------------------------------
            /*
             *  constructor
             */
            IWG_Incompressible_NS_Pressure_Neumann();

            //------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IWG_Incompressible_NS_Pressure_Neumann(){};

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
             * compute the derivative of the residual wrt design variables
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void compute_dRdp( real aWStar );

            //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IWG_INCOMPRESSIBLE_NS_PRESSURE_BULK_HPP_ */
