/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_SP_Crosswind_SA.hpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_SP_CROSSWIND_SA_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_SP_CROSSWIND_SA_HPP_

#include <map>

#include "typedefs.hpp"    //MRS/COR/src
#include "cl_Cell.hpp"     //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_IQI.hpp"    //FEM/INT/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IQI_SP_Crosswind_SA : public IQI
        {
            //------------------------------------------------------------------------------

            // local constitutive enums
            enum class IQI_Constitutive_Type
            {
                SPALART_ALLMARAS_TURBULENCE,
                MAX_ENUM
            };

            // local stabilization enums
            enum class IQI_Stabilization_Type
            {
                CROSSWIND,
                MAX_ENUM
            };

          public:
            //------------------------------------------------------------------------------
            /*
             * constructor
             */
            IQI_SP_Crosswind_SA();

            //------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IQI_SP_Crosswind_SA(){};

            //------------------------------------------------------------------------------

          private:
            //------------------------------------------------------------------------------
            /**
             * compute the quantity of interest
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void compute_QI( real aWStar );

            //------------------------------------------------------------------------------
            /**
             * Evaluate the quantity of interest and fill aQI with value
             * @param[ in ] aQI IQI value at evaluation point
             */
            void compute_QI( Matrix< DDRMat >& aQI );

            //------------------------------------------------------------------------------
            /**
             * compute the derivative of the quantity of interest wrt dof types
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void compute_dQIdu( real aWStar );

            //------------------------------------------------------------------------------
            /**
             * compute the derivative of the quantity of interest wrt dof types
             * @param[ in ] aDofType group of dof types wrt which derivatives are evaluated
             * @param[ in ] adQIdu   derivative of quantity of interest matrix to fill
             */
            void compute_dQIdu(
                    moris::Cell< MSI::Dof_Type >& aDofType,
                    Matrix< DDRMat >&             adQIdu );

            //------------------------------------------------------------------------------
            /**
             * compute the residual strong form
             * @param[ in ] aR real to fill with R
             */
            void compute_residual_strong_form( Matrix< DDRMat > & aR );

            //------------------------------------------------------------------------------
        };
    } /* end namespace fem */
} /* end namespace moris */

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_SP_CROSSWIND_SA_HPP_ */

