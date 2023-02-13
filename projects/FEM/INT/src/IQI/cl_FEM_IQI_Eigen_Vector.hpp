/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Eigen_Vector.hpp
 *
 */

#pragma once

#include <map>

#include "typedefs.hpp"    //MRS/COR/src
#include "cl_Cell.hpp"     //MRS/CON/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_IQI.hpp"    //FEM/INT/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IQI_Eigen_Vector : public IQI
        {
          private:
            //! index of eigen vector
            uint mEigenVectorIndex = 0;

          public:
            //------------------------------------------------------------------------------
            /*
             * constructor
             */
            IQI_Eigen_Vector();

            //------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IQI_Eigen_Vector(){};

            //------------------------------------------------------------------------------
            /**
             * child implementation of set_parameter function
             */
            void
            set_parameters( const moris::Cell< Matrix< DDRMat > >& aParameters );

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
             * evaluate the quantity of interest
             * @param[ in ] Matrix
             *
             * @return void
             */
            void evaluate_QI( Matrix< DDRMat >& aMat );

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
