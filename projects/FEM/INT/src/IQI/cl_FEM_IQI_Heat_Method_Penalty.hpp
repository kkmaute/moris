/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Heat_Method_Penalty.hpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_HEAT_METHOD_PENALTY__HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_HEAT_METHOD_PENALTY__HPP_

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

        class IQI_Heat_Method_Penalty : public IQI
        {
            //------------------------------------------------------------------------------

            /* Heat Method Penalty
             *
             * See https://doi.org/10.1007/s00158-019-02480-8
             *
             * out options controlled by vectorial index:
             *
             * 0: heat penalty (default)
             * 1: L2 contribution of heat penalty
             * 2: H1 contribution of heat penalty
             * 3: projected level set field
             * 4: norm of spatial gradients of projected level set field
             * 5: distance to interface measure
             * 6: L2 weight
             * 7: H1 weight
             * 8: intended sign of design level set function (optional)
             */

          public:
            //------------------------------------------------------------------------------
            /*
             * constructor
             */
            IQI_Heat_Method_Penalty();

            //------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IQI_Heat_Method_Penalty(){};

            //------------------------------------------------------------------------------

          private:
            enum class IQI_Property_Type
            {
                L2_REFERENCE_VALUE,
                H1S_REFERENCE_VALUE,
                SELECT,

                MAX_ENUM
            };

            //! initialization flag
            bool mIsInitialized = false;

            //! weight of L2 contribution
            real mL2Weight;

            //! weight of H1 semi-norm contribution
            real mH1SWeight;

            // Parameters from Markus,Coco paper
            real mPhiBound;
            real mPhiGamma;
            real mPhiGradient;
            real mWeightPhi1;
            real mWeightPhi2;
            real mWeightDelPhi1;
            real mWeightDelPhi2;

            // Sign of level set function
            real mLevelSetSign = 1.0;

            //! flag whether to skip computing dQIdu; skipping useful e.g. for level set regularization
            bool mSkipComputeDQIDU = false;

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

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_HEAT_METHOD_PENALTY__HPP_ */
