/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Spalart_Allmaras_Turbulence_Dirichlet.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_SPALART_ALLMARAS_TURBULENCE_DIRICHLET_HPP_
#define SRC_FEM_CL_FEM_IWG_SPALART_ALLMARAS_TURBULENCE_DIRICHLET_HPP_
//MRS/COR/src
#include <map>
#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"
//LINALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
//FEM/INT/src
#include "cl_FEM_IWG.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class IWG_Spalart_Allmaras_Turbulence_Dirichlet : public IWG
    {

        //------------------------------------------------------------------------------

      public:
        // sint for symmetric/unsymmetric Nitsche formulation
        sint mBeta = 1.0;

        // local property enums
        enum class IWG_Property_Type
        {
            DIRICHLET,
            SELECT,
            UPWIND,
            MAX_ENUM
        };

        // local constitutive enums
        enum class IWG_Constitutive_Type
        {
            SPALART_ALLMARAS_TURBULENCE,
            MAX_ENUM
        };

        // local stabilization enums
        enum class IWG_Stabilization_Type
        {
            NITSCHE,
            MAX_ENUM
        };

      private:
        // Spalart Allmaras model constants
        real mCb2   = 0.6220;
        real mSigma = 2.0 / 3.0;

      public:
        //------------------------------------------------------------------------------
        /*
         *  constructor
         *  aBeta signed int for symmetric/unsymmetric Nitsche
         */
        IWG_Spalart_Allmaras_Turbulence_Dirichlet( sint aBeta );

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IWG_Spalart_Allmaras_Turbulence_Dirichlet() override{};

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
    };
    //------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_IWG_SPALART_ALLMARAS_TURBULENCE_DIRICHLET_HPP_ */

