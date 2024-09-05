/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Diffusion_Bulk.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_DIFFUSION_BULK_HPP_
#define SRC_FEM_CL_FEM_IWG_DIFFUSION_BULK_HPP_

#include <map>
#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_IWG.hpp"    //FEM/INT/src

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class IWG_Diffusion_Bulk : public IWG
    {

        //------------------------------------------------------------------------------

      public:
        enum class IWG_Property_Type
        {
            BODY_LOAD,
            THICKNESS,
            H2_PENALTY,
            H3_PENALTY,
            Phase_Field,
            SELECT,
            MAX_ENUM
        };

        enum class IWG_Constitutive_Type
        {
            DIFFUSION,
            MAX_ENUM
        };

        enum class IWG_Stabilization_Type
        {
            GGLS_DIFFUSION,
            MAX_ENUM
        };

        //------------------------------------------------------------------------------
        /*
         *  constructor
         */
        IWG_Diffusion_Bulk();

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IWG_Diffusion_Bulk() override{};

        //------------------------------------------------------------------------------
        /**
         * compute the residual
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_residual( real tWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the jacobian
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_jacobian( real tWStar ) override;

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

#endif /* SRC_FEM_CL_FEM_IWG_DIFFUSION_BULK_HPP_ */
