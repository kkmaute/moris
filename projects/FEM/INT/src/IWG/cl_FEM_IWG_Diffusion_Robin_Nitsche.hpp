/*
 * Copyright (c) 2022 University of Colorado
 *Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Diffusion_Robin_Nitsche.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_DIFFUSION_ROBIN_NITSCHE_HPP_
#define SRC_FEM_CL_FEM_IWG_DIFFUSION_ROBIN_NITSCHE_HPP_

#include <map>
// MRS/COR/src
#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"
// LINALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
// FEM/INT/src
#include "cl_FEM_IWG.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class IWG_Diffusion_Robin_Nitsche : public IWG
    {

        //------------------------------------------------------------------------------

      public:
        // sign for symmetric/unsymmetric Nitsche
        sint mBeta;

        enum class IWG_Property_Type
        {
            DIRICHLET,
            NEUMANN_PENALTY,
            TRACTION,
            MATERIAL_COEFFICIENT,
            MAX_ENUM
        };

        // local constitutive enums
        enum class IWG_Constitutive_Type
        {
            DIFFUSION,
            MAX_ENUM
        };

        // local stabilization enums
        enum class IWG_Stabilization_Type
        {
            ROBIN_NITSCHE,
            MAX_ENUM
        };

        //------------------------------------------------------------------------------
        /*
         *  constructor
         */
        IWG_Diffusion_Robin_Nitsche( sint aBeta );

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IWG_Diffusion_Robin_Nitsche() override{};

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

#endif /* SRC_FEM_CL_FEM_IWG_INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_NITSCHE_HPP_ */
