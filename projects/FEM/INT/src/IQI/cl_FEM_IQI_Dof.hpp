/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Dof.hpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_DOF_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_DOF_HPP_

#include <map>

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_IQI.hpp"    //FEM/INT/src

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class IQI_Dof : public IQI
    {
      private:
        //! initialization flag
        bool mIsInitialized = false;

        //! spatial and time derivative information
        uint mSpatialDerivativeDirection = 0;
        uint mSpatialDerivativeOrder     = 0;
        uint mTimeDerivativeOrder        = 0;

      public:
        //------------------------------------------------------------------------------
        /*
         * constructor
         */
        IQI_Dof();

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IQI_Dof() override{};

        //------------------------------------------------------------------------------

      private:
        //------------------------------------------------------------------------------
        /**
         * initialize parameters
         */
        void initialize();

        //------------------------------------------------------------------------------
        /**
         * compute the quantity of interest
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_QI( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * Evaluate the quantity of interest and fill aQI with value
         * @param[ in ] aQI IQI value at evaluation point
         */
        void compute_QI( Matrix< DDRMat >& aQI ) override;

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
        void compute_dQIdu( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the quantity of interest wrt dof types
         * @param[ in ] aDofType group of dof types wrt which derivatives are evaluated
         * @param[ in ] adQIdu   derivative of quantity of interest matrix to fill
         */
        void compute_dQIdu(
                Vector< MSI::Dof_Type >& aDofType,
                Matrix< DDRMat >&        adQIdu ) override;

        //------------------------------------------------------------------------------
    };
}    // namespace moris::fem

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_DOF_HPP_ */

