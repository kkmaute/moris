/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Incompressible_NS_Convective_Velocity_Ghost.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_INCOMPRESSIBLE_NS_CONVECTIVE_VELOCITY_GHOST_HPP_
#define SRC_FEM_CL_FEM_IWG_INCOMPRESSIBLE_NS_CONVECTIVE_VELOCITY_GHOST_HPP_

#include <map>
//MRS/COR/src
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

    class IWG_Incompressible_NS_Convective_Velocity_Ghost : public IWG
    {

        //------------------------------------------------------------------------------

      public:
        // local stabilization enums
        enum class IWG_Stabilization_Type
        {
            CONVECTIVE_GHOST,
            MAX_ENUM
        };

        //------------------------------------------------------------------------------
        /*
         *  constructor
         */
        IWG_Incompressible_NS_Convective_Velocity_Ghost();

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IWG_Incompressible_NS_Convective_Velocity_Ghost() override{};

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

      private:
        //------------------------------------------------------------------------------
        /**
         * compute the block matrix for dnNdxn
         * @param[ in ] adNdx     matrix to fill with derivatives
         * @param[ in ] aIsLeader enum for leader or follower
         */
        void compute_dnNdxn(
                Matrix< DDRMat >    &adNdx,
                mtk::Leader_Follower aIsLeader );

        //------------------------------------------------------------------------------
    };
    //------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_IWG_INCOMPRESSIBLE_NS_CONVECTIVE_VELOCITY_GHOST_HPP_ */

