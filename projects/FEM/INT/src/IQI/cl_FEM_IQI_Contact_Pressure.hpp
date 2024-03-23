/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Contact_Pressure.hpp
 *
 */

#pragma once
#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_IQI.hpp"    //FEM/INT/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IQI_Contact_Pressure : public IQI
        {
            bool mInReferenceConfiguration = false;

            enum class IQI_Constitutive_Type
            {
                TRACTION_CM,
                MAX_ENUM
            };

            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------
            /*
             * constructor
             */
            IQI_Contact_Pressure( bool aInReferenceConfiguration );

            //------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IQI_Contact_Pressure(){};

            //------------------------------------------------------------------------------

          private:
            //------------------------------------------------------------------------------
            /**
             * compute the quantity of interest
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void compute_QI( Matrix< DDRMat >& aQI );

            //------------------------------------------------------------------------------
            /**
             * Evaluate the quantity of interest and fill aQI with value
             * @param[ in ] aQI IQI value at evaluation point
             */
            void compute_QI( real aWStar );

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
                    Vector< MSI::Dof_Type >& aDofType,
                    Matrix< DDRMat >&        adQIdu );

            //------------------------------------------------------------------------------

        };    // class IQI_Contact_Pressure

        //------------------------------------------------------------------------------

    } /* end namespace fem */
} /* end namespace moris */