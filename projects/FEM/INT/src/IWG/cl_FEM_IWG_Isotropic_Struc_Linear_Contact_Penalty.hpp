/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Penalty.hpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_CONTACT_PENALTY_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_CONTACT_PENALTY_HPP_

#include <map>

#include "moris_typedefs.hpp"                     //MRS/COR/src
#include "cl_Vector.hpp"                          //MRS/CNT/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_IWG.hpp"                   //FEM/INT/src

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class IWG_Isotropic_Struc_Linear_Contact_Penalty : public IWG
    {

      public:
        enum class IWG_Property_Type
        {
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
            PENALTY_CONTACT,
            STAB_PENALTY_CONTACT,
            MAX_ENUM
        };

        //------------------------------------------------------------------------------
        /*
         * constructor
         */
        IWG_Isotropic_Struc_Linear_Contact_Penalty();

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IWG_Isotropic_Struc_Linear_Contact_Penalty() override{};

        //------------------------------------------------------------------------------
        /**
         * compute the residual
         * @param[ in ] aResidual cell of residual vectors to fill
         */
        void compute_residual( real tWStar ) override;
        //------------------------------------------------------------------------------
        /**
         * compute the jacobian
         * @param[ in ] aJacobians cell of cell of jacobian matrices to fill
         */
        void compute_jacobian( real tWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the residual and the jacobian
         * @param[ in ] aJacobians cell of cell of jacobian matrices to fill
         * @param[ in ] aResidual  cell of residual vectors to fill
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

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_CONTACT_PENALTY_HPP_ */

