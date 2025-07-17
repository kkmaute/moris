/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Ghost_Normal_Field.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_Ghost_Normal_Field_HPP_
#define SRC_FEM_CL_FEM_IWG_Ghost_Normal_Field_HPP_

//MRS/COR/src
#include "moris_typedefs.hpp"
//LINALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
//FEM/INT/src
#include "cl_FEM_IWG.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class IWG_Ghost_Normal_Field : public IWG
    {
        //------------------------------------------------------------------------------

      private:
        enum class IWG_Stabilization_Type
        {
            GHOST_SP,
            MAX_ENUM
        };

      public:
        //------------------------------------------------------------------------------
        /*
         *  constructor
         */
        IWG_Ghost_Normal_Field();

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IWG_Ghost_Normal_Field() override{};

        //------------------------------------------------------------------------------
        /**
         * compute the residual
         * @param[ in ] aResidual cell of residual vectors to fill
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
         * get flattened normal matrix
         * @param[ in ] aNormal flat normal matrix to fill
         * @param[ in ] aOrder  interpolation order
         */
        void get_flat_normal_matrix(
                Matrix< DDRMat > &aFlatNormal,
                uint              aOrder );

        //------------------------------------------------------------------------------
        /**
         * compute the block matrix for directional derivative
         * (dnNdxn . normal)
         * @param[ in ] aFlatdnNdxn matrix to fill with derivatives
         * @param[ in ] aOrder      interpolation order for residual dof type
         * @param[ in ] aIsLeader   enum for leader or follower
         */
        void compute_flat_dnNdxn(
                Matrix< DDRMat >    &aFlatdnNdxn,
                uint                 aOrder,
                mtk::Leader_Follower aIsLeader );

        //------------------------------------------------------------------------------
    };
    //------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_IWG_Ghost_Normal_Field_HPP_ */

