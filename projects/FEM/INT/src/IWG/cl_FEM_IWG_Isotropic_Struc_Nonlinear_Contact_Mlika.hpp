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

    class IWG_Isotropic_Struc_Nonlinear_Contact_Mlika : public IWG
    {
      private:
        // mTheta is used to switch between the following formulations of the Nitsche term:
        //      symmetric :mTheta =  1
        //        neutral: mTheta =  0
        //    unsymmetric: mTheta = -1
        sint mTheta = 1;

        // default for traction
        CM_Function_Type mCMFunctionType = CM_Function_Type::PK1;

      public:
        enum class IWG_Property_Type
        {
            MATERIAL,
            THICKNESS,
            SELECT,
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
        IWG_Isotropic_Struc_Nonlinear_Contact_Mlika(
                sint             aBeta,
                CM_Function_Type aCMFunctionType = CM_Function_Type::PK1 );

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IWG_Isotropic_Struc_Nonlinear_Contact_Mlika() override {};

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
        /**
         * set parameters
         */
        void set_parameters( const Vector< Matrix< DDRMat > >& aParameters ) override;

        //------------------------------------------------------------------------------
    };
    //------------------------------------------------------------------------------
}    // namespace moris::fem
