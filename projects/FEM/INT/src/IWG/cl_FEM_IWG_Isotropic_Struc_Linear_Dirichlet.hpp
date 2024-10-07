/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Linear_Dirichlet.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_DIRICHLET_HPP_
#define SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_DIRICHLET_HPP_

#include <map>

#include "moris_typedefs.hpp"                     //MRS/COR/src
#include "cl_Vector.hpp"                          //MRS/CNT/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_IWG.hpp"                   //FEM/INT/src

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class IWG_Isotropic_Struc_Linear_Dirichlet : public IWG
    {

        //------------------------------------------------------------------------------

      public:
        // sign for symmetric/unsymmetric Nitsche
        sint mBeta = 1;

        enum class IWG_Property_Type
        {
            DIRICHLET,
            SELECT,
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
            DIRICHLET_NITSCHE,
            MAX_ENUM
        };

        //------------------------------------------------------------------------------
        /*
         * constructor
         */
        IWG_Isotropic_Struc_Linear_Dirichlet( sint aBeta );

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IWG_Isotropic_Struc_Linear_Dirichlet() override{};

        //------------------------------------------------------------------------------
        /**
         * computes the residual
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_residual( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * computes the jacobian
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_jacobian( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * computes the residual and the jacobian
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

#endif /* SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_DIRICHLET_HPP_ */

