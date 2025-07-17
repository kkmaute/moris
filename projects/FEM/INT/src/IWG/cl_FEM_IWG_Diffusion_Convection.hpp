/*
 * Copyright (c) 2022 University of Colorado
 *Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Diffusion_Convection.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_DIFFUSION_CONVECTION_HPP_
#define SRC_FEM_CL_FEM_IWG_DIFFUSION_CONVECTION_HPP_

#include "moris_typedefs.hpp"    //MRS/COR/src

#include "cl_FEM_IWG.hpp"                   //FEM/INT/src

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class IWG_Diffusion_Convection : public IWG
    {
        //------------------------------------------------------------------------------

      public:
        enum class IWG_Property_Type
        {
            HEAT_TRANSFER_COEFFICIENT,
            AMBIENT_TEMP,
            THICKNESS,
            MAX_ENUM
        };

        //------------------------------------------------------------------------------
        /*
         * constructor
         */
        IWG_Diffusion_Convection();

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IWG_Diffusion_Convection() override{};

        //------------------------------------------------------------------------------
        /**
         * compute the residual
         * @param[ in ] aWStar weight associated with evaluation point
         */
        void compute_residual( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the jacobian
         * @param[ in ] aWStar weight associated with evaluation point
         */
        void compute_jacobian( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the residual and the jacobian
         * @param[ in ] aWStar weight associated with evaluation point
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

#endif /* SRC_FEM_CL_FEM_IWG_Diffusion_Neumann_HPP_ */
