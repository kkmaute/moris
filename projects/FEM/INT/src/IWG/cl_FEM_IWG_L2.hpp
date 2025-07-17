/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_L2.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_L2_HPP_
#define SRC_FEM_CL_FEM_IWG_L2_HPP_

#include "moris_typedefs.hpp"    //MRS/COR/src

#include "cl_FEM_IWG.hpp"                   //FEM/INT/src

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class IWG_L2 : public IWG
    {
        enum class IWG_Property_Type
        {
            L2COEFFICIENT,    // coefficient for L2 norm term
            H1COEFFICIENT,    // coefficient for H1 semi-norm term
            DIFFUSION,        // diffusion coefficient
            SOURCE,           // source term
            THICKNESS,        // thickness
            MAX_ENUM
        };

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------
        // constructor
        IWG_L2();

        //------------------------------------------------------------------------------
        // trivial destructor
        ~IWG_L2() override{};

        //------------------------------------------------------------------------------
        /**
         * compute jacobian and residual
         */
        void compute_jacobian_and_residual( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * compute jacobian
         */
        void compute_jacobian( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * compute residual
         */
        void compute_residual( real aWStar ) override;

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

#endif /* SRC_FEM_CL_FEM_IWG_L2_HPP_ */
