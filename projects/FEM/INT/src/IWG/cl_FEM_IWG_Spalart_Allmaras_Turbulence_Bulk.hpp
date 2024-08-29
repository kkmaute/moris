/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Spalart_Allmaras_Turbulence_Bulk.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_SPALART_ALLMARAS_TURBULENCE_BULK_HPP_
#define SRC_FEM_CL_FEM_IWG_SPALART_ALLMARAS_TURBULENCE_BULK_HPP_
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

    class IWG_Spalart_Allmaras_Turbulence_Bulk : public IWG
    {

        //------------------------------------------------------------------------------

      public:
        // local constitutive enums
        enum class IWG_Constitutive_Type
        {
            SPALART_ALLMARAS_TURBULENCE,
            MAX_ENUM
        };

        // local stabilization enums
        enum class IWG_Stabilization_Type
        {
            SUPG,
            DIFFUSION_CROSSWIND,
            DIFFUSION_ISOTROPIC,
            MAX_ENUM
        };

      public:
        //------------------------------------------------------------------------------
        /*
         *  constructor
         */
        IWG_Spalart_Allmaras_Turbulence_Bulk();

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IWG_Spalart_Allmaras_Turbulence_Bulk() override{};

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
         * @param[ in ] aR real to fill with R
         */
        void compute_residual_strong_form( Matrix< DDRMat > &aR );

        //------------------------------------------------------------------------------
        /**
         * compute the jacobian strong form
         * @param[ in ] aDofTypes a list of dof type wrt which
         *                        the derivative is requested
         * @param[ in ] aJ        a matrix to fill with dRdDof
         */
        void compute_jacobian_strong_form(
                const Vector< MSI::Dof_Type > &aDofTypes,
                Matrix< DDRMat >              &aJ );

        //------------------------------------------------------------------------------
    };
    //------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_IWG_SPALART_ALLMARAS_TURBULENCE_BULK_HPP_ */

