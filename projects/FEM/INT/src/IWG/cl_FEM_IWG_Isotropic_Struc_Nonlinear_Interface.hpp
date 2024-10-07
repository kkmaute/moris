/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Nonlinear_Interface.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_NONLINEAR_INTERFACE_HPP_
#define SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_NONLINEAR_INTERFACE_HPP_

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

    class IWG_Isotropic_Struc_Nonlinear_Interface : public IWG
    {

      public:
        // sint for symmetric/unsymmetric Nitsche
        sint mBeta = 1;

        // stress and strain type to evaluate the IWG
        enum CM_Function_Type mStressType;
        enum CM_Function_Type mStrainType;

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
            NITSCHE_INTERFACE,
            MAX_ENUM
        };

        //------------------------------------------------------------------------------
        /*
         * constructor
         * @param[ in ] aBeta +1 or -1 for symmetric/unsymmetric symmetric Nitsche
         */
        IWG_Isotropic_Struc_Nonlinear_Interface(
                enum CM_Function_Type aStressType,
                enum CM_Function_Type aStrainType,
                sint                  aBeta );

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IWG_Isotropic_Struc_Nonlinear_Interface() override{};

        //------------------------------------------------------------------------------
        /**
         * compute the residual
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_residual( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the jacobian
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_jacobian( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the residual and the jacobian
         * @param[ in ] aWStar weight associated to the evaluation point
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

#endif /* SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_NONLINEAR_INTERFACE_HPP_ */

