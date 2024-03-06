/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Turbulent_Dynamic_Viscosity.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Turbulent_Dynamic_Viscosity.hpp"
#include "fn_norm.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Turbulent_Dynamic_Viscosity::IQI_Turbulent_Dynamic_Viscosity()
        {
            // set fem IQI type
            mFEMIQIType = fem::IQI_Type::TURBULENT_DYNAMIC_VISCOSITY;
            init_constitutive_model( "Fluid_Turbulence", IQI_Constitutive_Type::FLUID_TURBULENCE );
        }

        //------------------------------------------------------------------------------

        void IQI_Turbulent_Dynamic_Viscosity::compute_QI( Matrix< DDRMat > &aQI )
        {
            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > &tCMFluidTurbulence = get_leader_constitutive_model(IQI_Constitutive_Type::FLUID_TURBULENCE);

            // compute turbulent dynamic viscosity
            aQI = tCMFluidTurbulence->turbulent_dynamic_viscosity();
        }

        //------------------------------------------------------------------------------

        void IQI_Turbulent_Dynamic_Viscosity::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( get_name() );

            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > &tCMFluidTurbulence = get_leader_constitutive_model(IQI_Constitutive_Type::FLUID_TURBULENCE);

            // compute turbulent dynamic viscosity
            mSet->get_QI()( tQIIndex ) += aWStar * ( tCMFluidTurbulence->turbulent_dynamic_viscosity() );
        }
        //------------------------------------------------------------------------------
    }    // namespace fem
}    // namespace moris
