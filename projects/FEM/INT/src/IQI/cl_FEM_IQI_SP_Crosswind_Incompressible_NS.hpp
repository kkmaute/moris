/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_SP_Crosswind_Incompressible_NS.hpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_SP_CROSSWIND_INCOMPRESSIBLE_NS_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_SP_CROSSWIND_INCOMPRESSIBLE_NS_HPP_

#include <map>

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_IQI.hpp"    //FEM/INT/src

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class IQI_SP_Crosswind_Incompressible_NS : public IQI
    {
        //------------------------------------------------------------------------------

        // local property enums
        enum class IQI_Property_Type
        {
            GRAVITY,
            THERMAL_EXPANSION,
            REF_TEMP,
            INV_PERMEABILITY,
            MASS_SOURCE,
            BODY_LOAD,
            MAX_ENUM
        };

        // local constitutive enums
        enum class IQI_Constitutive_Type
        {
            INCOMPRESSIBLE_FLUID,
            MAX_ENUM
        };

        // local stabilization enums
        enum class IQI_Stabilization_Type
        {
            CROSSWIND,
            MAX_ENUM
        };

      public:
        //------------------------------------------------------------------------------
        /*
         * constructor
         */
        IQI_SP_Crosswind_Incompressible_NS();

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IQI_SP_Crosswind_Incompressible_NS() override{};

        //------------------------------------------------------------------------------

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
        void compute_QI( Matrix< DDRMat >& aQI ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the quantity of interest wrt dof types
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_dQIdu( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the quantity of interest wrt dof types
         * @param[ in ] aDofType group of dof types wrt which derivatives are evaluated
         * @param[ in ] adQIdu   derivative of quantity of interest matrix to fill
         */
        void compute_dQIdu(
                Vector< MSI::Dof_Type >& aDofType,
                Matrix< DDRMat >&        adQIdu ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the residual strong form for momentum
         * @param[ in ] aRM a matrix to fill with RM
         */
        void compute_residual_strong_form_momentum(
                Matrix< DDRMat >& aRM );

        //------------------------------------------------------------------------------
    };
}    // namespace moris::fem

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_SP_CROSSWIND_INCOMPRESSIBLE_NS_HPP_ */

