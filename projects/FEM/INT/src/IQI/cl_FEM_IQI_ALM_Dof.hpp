/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Max_Dof.hpp
 *
 */

#pragma once

#include <map>

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Cell.hpp"     //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_IQI.hpp"    //FEM/INT/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IQI_ALM_Dof : public IQI
        {
            /**
             * @brief Implementation of Augmented Lagrange Multiplier Method defined by integral over space and time of
             *
             * lambda * gplus + 0.5 * cpen * gplus^2
             *
             * where gplus is the constraint (g = dof-value / mRefValue - mShift <= 0) and defined by
             *
             * gplus = max(g, - lambda/cpen)
             *
             * and lambda is the Lagrange multiplier field and cpen the penalty factor
             *
             */

            //------------------------------------------------------------------------------

            enum class IQI_Property_Type
            {
                LAGRANGE_MULTIPLIER,
                PENALTY_FACTOR,
                MAX_ENUM
            };

          public:
            //------------------------------------------------------------------------------
            /*
             * constructor
             */
            IQI_ALM_Dof();

            //------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IQI_ALM_Dof(){};

            //------------------------------------------------------------------------------

          private:
            //------------------------------------------------------------------------------

            //! initialization flag
            bool mIsInitialized = false;

            //! parameters: reference value and shift
            real mRefValue;
            real mShift;

            //------------------------------------------------------------------------------
            /**
             * initialize parameters
             */
            void initialize();

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
        };
    } /* end namespace fem */
} /* end namespace moris */
