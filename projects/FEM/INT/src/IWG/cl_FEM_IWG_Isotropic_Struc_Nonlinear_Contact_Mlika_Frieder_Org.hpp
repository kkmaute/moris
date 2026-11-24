/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Gap_Nitsche_Unbiased.hpp
 *
 */
#pragma once

#include <map>

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_IWG.hpp"                   //FEM/INT/src


namespace moris::fem
{
    //------------------------------------------------------------------------------

    class IWG_Isotropic_Struc_Nonlinear_Contact_Mlika_Frieder_Org : public IWG
    {
      private:
        // sign for symmetric/unsymmetric Nitsche
        sint mBeta = 1;

        // the literature used for this implementation uses a different parameter to switch between symmetric and unsymmetric Nitsche.
        // to avoid confusion with other IWGs, mBeta will be used to switch between symmetric and unsymmetric Nitsche in the factory, and
        // mTheta will be used to switch between the two different Nitsche formulations in the implementation.
        //      symmetric: mBeta = -1, mTheta =  1
        //    unsymmetric: mBeta =  0, mTheta =  0
        // skew-symmetric: mBeta =  1, mTheta = -1
        sint mTheta = 1;

      public:
        enum class IWG_Property_Type
        {
            MATERIAL,
            THICKNESS,
            GAP,
            MAX_ENUM
        };

        enum class IWG_Constitutive_Type
        {
            ELAST_LIN_ISO,
            MAX_ENUM
        };

        enum class IWG_Stabilization_Type
        {
            NITSCHE_INTERFACE,
            MAX_ENUM
        };

        //------------------------------------------------------------------------------
        /*
         * constructor
         */
        IWG_Isotropic_Struc_Nonlinear_Contact_Mlika_Frieder_Org( sint aBeta );

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IWG_Isotropic_Struc_Nonlinear_Contact_Mlika_Frieder_Org() override {};

        //------------------------------------------------------------------------------
        /**
         * compute the residual
         * @param[ in ] aResidual cell of residual vectors to fill
         */
        void compute_residual( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the jacobian
         * @param[ in ] aJacobians cell of cell of jacobian matrices to fill
         */
        void compute_jacobian( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the residual and the jacobian
         * @param[ in ] aJacobians cell of cell of jacobian matrices to fill
         * @param[ in ] aResidual  cell of residual vectors to fill
         */
        void compute_jacobian_and_residual( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the residual wrt design variables
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_dRdp( real aWStar ) override;

        //------------------------------------------------------------------------------
    };
    //------------------------------------------------------------------------------
}    // namespace moris::fem
