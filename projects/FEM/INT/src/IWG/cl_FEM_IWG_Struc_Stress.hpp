/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Struc_Stress.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_STRUC_STRESS_HPP_
#define SRC_FEM_CL_FEM_IWG_STRUC_STRESS_HPP_

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

        class IWG_Struc_Stress : public IWG
        {

            //------------------------------------------------------------------------------

          public:
            // stress type to evaluate
            enum Stress_Type mStressType = Stress_Type::END_STRESS_TYPE;

            enum class IWG_Property_Type
            {
                THICKNESS,
                MAX_ENUM
            };

            enum class IWG_Constitutive_Type
            {
                ELAST_LIN_ISO,
                MAX_ENUM
            };

            enum class IWG_Stabilization_Type
            {
                // GGLS_DIFFUSION,
                MAX_ENUM
            };

            //------------------------------------------------------------------------------
            /*
             *  constructor
             */
            IWG_Struc_Stress( enum Stress_Type aStressType );

            //------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IWG_Struc_Stress(){};

            //------------------------------------------------------------------------------
            /**
             * compute the residual
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void compute_residual( real tWStar );

            //------------------------------------------------------------------------------
            /**
             * compute the jacobian
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void compute_jacobian( real tWStar );

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

          private:
            //------------------------------------------------------------------------------

            void eval_stress_criterion(
                    moris::real&             aStressVal,
                    Matrix< DDRMat >&        aDStressVal,
                    Matrix< DDRMat > const & aStressTensor,
                    Matrix< DDRMat > const & aDStressTensor );

            //------------------------------------------------------------------------------

            void eval_Von_Mises_stress(
                    moris::real&             aStressVal,
                    Matrix< DDRMat >&        aDStressVal,
                    Matrix< DDRMat > const & aStressTensor,
                    Matrix< DDRMat > const & aDStressTensor );

            //------------------------------------------------------------------------------

            //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_STRUC_STRESS_HPP_ */
