/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Compressible_NS_Advective_Momentum_Flux_Boundary.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_ADVECTIVE_MOMENTUM_FLUX_BOUNDARY_HPP_
#define SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_ADVECTIVE_MOMENTUM_FLUX_BOUNDARY_HPP_

#include <map>
#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_IWG.hpp"    //FEM/INT/src

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class IWG_Compressible_NS_Advective_Momentum_Flux_Boundary : public IWG
    {

        //------------------------------------------------------------------------------

      private:
        // default dof types
        MSI::Dof_Type mDofDensity  = MSI::Dof_Type::RHO;
        MSI::Dof_Type mDofVelocity = MSI::Dof_Type::VX;

      public:
        // local constitutive enums
        enum class IWG_Constitutive_Type
        {
            FLUID,
            MAX_ENUM
        };

        //------------------------------------------------------------------------------
        /*
         *  constructor
         */
        IWG_Compressible_NS_Advective_Momentum_Flux_Boundary();

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IWG_Compressible_NS_Advective_Momentum_Flux_Boundary() override{};

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

      private:
        //------------------------------------------------------------------------------
        /**
         * compute the residual strong form
         * @param[ in ] aRM a matrix to fill with RM
         * @param[ in ] aRC a matrix to fill with RC
         */
        void compute_residual_strong_form(
                Matrix< DDRMat > &aRM,
                real             &aRC );

        //------------------------------------------------------------------------------
        /**
         * compute the residual strong form
         * @param[ in ] aDofTypes a list of dof type wrt which
         *                        the derivative is requested
         * @param[ in ] aJM       a matrix to fill with dRMdDof
         * @param[ in ] aJC       a matrix to fill with dRCdDof
         */
        void compute_jacobian_strong_form(
                const Vector< MSI::Dof_Type > &aDofTypes,
                Matrix< DDRMat >              &aJM,
                Matrix< DDRMat >              &aJC );

        //------------------------------------------------------------------------------
        /**
         * compute the term ui uj
         * @param[ in ] auiuj a matrix to fill with ui uj
         */
        void compute_uiuj( Matrix< DDRMat > &auiuj );

        //------------------------------------------------------------------------------
        /**
         * computes the term d(ui uj)/duhat
         * @param[ in ] aduiujdDOF a matrix to fill with result
         */
        void compute_duiujdDOF( Matrix< DDRMat > &aduiujdDOF );

        //------------------------------------------------------------------------------
        /**
         * returns a matrix containing the entries of the surface normal for multiplication with flattened tensors
         * @param[ in ] aduiujdDOF a matrix to fill with result
         */
        void compute_normal_matrix( Matrix< DDRMat > &aNormalMatrix );
    };
    //------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_ADVECTIVE_MOMENTUM_FLUX_BOUNDARY_HPP_ */
