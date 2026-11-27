/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Nonlinear_Contact_Seitz.hpp
 *
 */
#pragma once

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_IWG.hpp"    //FEM/INT/src

namespace moris::fem
{
    //------------------------------------------------------------------------------

    /**
     * Seitz friction uses current constitutive traction in the trial projection
     * (Seitz, Wall and Popp, 2019, Eq. (33), DOI: 10.1007/s00466-018-1638-x).
     * Previous-slab displacements enter only the slip increment; no independent
     * traction history or previous-stress reconstruction is required.
     */
    class IWG_Isotropic_Struc_Nonlinear_Contact_Seitz : public IWG
    {
      private:
        // Retain the constructor-selected variant when penalty-only mode
        // temporarily disables the Nitsche terms during parameter replacement.
        const sint mNitscheTheta;

        // mTheta is used to switch between the following formulations of the Nitsche term:
        //      symmetric :mTheta =  1
        //        neutral: mTheta =  0
        //    unsymmetric: mTheta = -1
        sint mTheta = 1;

        // default for traction
        CM_Function_Type mCMFunctionType = CM_Function_Type::PK1;

        // scaling factor for the traction; set to zero for pure penalty formulation
        real mTractionScaling = 1.0;

        real mFrictionCoefficient = -1.0;    // Nonpositive values select frictionless contact.
        real mGammaTau            = 1.0;     // Nonpositive values use the Nitsche parameter.

        // flag to trigger debug output
        sint mDebugFlag = 0;

        // debug functions
        void debug_function();

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
        IWG_Isotropic_Struc_Nonlinear_Contact_Seitz(
                sint             aBeta,
                CM_Function_Type aCMFunctionType = CM_Function_Type::PK1 );

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IWG_Isotropic_Struc_Nonlinear_Contact_Seitz() override {};

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
         * compute dR(n)/du(n-1)
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_jacobian_previous( real aWStar ) override;

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
    };
    //------------------------------------------------------------------------------
}    // namespace moris::fem
