/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Max_Dof.hpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_MAX_DOF_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_MAX_DOF_HPP_

#include <map>

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"     //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_IQI.hpp"    //FEM/INT/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IQI_Max_Dof : public IQI
        {

            /**
             * @brief IQI to compute integral over space and time of
             *
             * ( ( dof-value / mRefValue - mShift)^pm )^mExponent
             *
             * where the integrand can be considered if strictly positive or negative or both
             */

            //------------------------------------------------------------------------------

            enum class IQI_Property_Type
            {
                WEIGHT,
                MAX_ENUM
            };

          public:
            //------------------------------------------------------------------------------
            /*
             * constructor
             */
            IQI_Max_Dof()
            {
                // set size for the property pointer cell
                mLeaderProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

                // populate the property map
                mPropertyMap[ "Weight" ] = static_cast< uint >( IQI_Property_Type::WEIGHT );

                // set FEM IQI type
                mFEMIQIType = fem::IQI_Type::MAX_DOF;
            }

            //------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IQI_Max_Dof(){};

            //------------------------------------------------------------------------------

          private:
            //------------------------------------------------------------------------------

            //! initialization flag
            bool mIsInitialized = false;

            //! parameters: reference value, exponent, shift, sign
            real mRefValue;
            real mExponent;
            real mShift;
            sint mSign;    /// flag to consider the integrand for
                           /// only positive values (1) or
                           /// only negative values (-1) or
                           /// all values (0) (default)

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
                    Vector< MSI::Dof_Type >& aDofType,
                    Matrix< DDRMat >&             adQIdu );

            //------------------------------------------------------------------------------
        };
    } /* end namespace fem */
} /* end namespace moris */

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_MAX_DOF_HPP_ */
