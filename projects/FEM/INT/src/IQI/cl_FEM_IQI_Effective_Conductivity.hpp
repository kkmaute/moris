/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Effective_Conductivity.hpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_EFFECTIVE_CONDUCTIVITY_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_EFFECTIVE_CONDUCTIVITY_HPP_

#include <map>

#include "moris_typedefs.hpp"                     //MRS/COR/src
#include "cl_Vector.hpp"                          //MRS/CNT/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_IQI.hpp"                   //FEM/INT/src
#include "cl_FEM_CM_Diffusion_Linear_Isotropic_Turbulence.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class IQI_Effective_Conductivity : public IQI
    {
      private:
        enum class IQI_Constitutive_Type
        {
            DIFFUSION_TURBULENCE,
            MAX_ENUM
        };

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------
        /*
         * constructor
         */
        IQI_Effective_Conductivity();

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IQI_Effective_Conductivity() override{};

      private:
        //------------------------------------------------------------------------------
        /**
         * compute the quantity of interest
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_QI( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * Evaluate the quantity of interest and fill aQI with value
         * @param[ in ] aQI IQI value at evaluation point
         */
        void compute_QI( Matrix< DDRMat > &aQI ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the quantity of interest wrt dof types
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_dQIdu( real aWStar ) override
        {
            MORIS_ERROR( false, "IQI_Effective_Conductivity::compute_dQIdu - not implemented." );
        }

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the quantity of interest wrt dof types
         * @param[ in ] aDofType group of dof types wrt which derivatives are evaluated
         * @param[ in ] adQIdu   derivative of quantity of interest matrix to fill
         */
        void compute_dQIdu(
                Vector< MSI::Dof_Type > &aDofType,
                Matrix< DDRMat >        &adQIdu ) override
        {
            MORIS_ERROR( false, "IQI_Effective_Conductivity::compute_dQIdu() - not implemented for a drag/lift coefficient IQI." );
        }
    };
}    // namespace moris::fem

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_EFFECTIVE_CONDUCTIVITY_HPP_ */

