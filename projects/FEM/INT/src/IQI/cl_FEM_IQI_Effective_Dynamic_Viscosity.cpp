/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Effective_Dynamic_Viscosity.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Effective_Dynamic_Viscosity.hpp"
#include "fn_norm.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_Effective_Dynamic_Viscosity::IQI_Effective_Dynamic_Viscosity()
    {
        // set fem IQI type
        mFEMIQIType = fem::IQI_Type::EFFECTIVE_DYNAMIC_VISCOSITY;

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "Fluid_Turbulence" ] =
                static_cast< uint >( IQI_Constitutive_Type::FLUID_TURBULENCE );
    }

    //------------------------------------------------------------------------------

    void IQI_Effective_Dynamic_Viscosity::compute_QI( Matrix< DDRMat > &aQI )
    {
        // get the diffusion CM
        const std::shared_ptr< Constitutive_Model > &tCMFluidTurbulence =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::FLUID_TURBULENCE ) );

        // compute effective dynamic viscosity
        aQI = tCMFluidTurbulence->effective_dynamic_viscosity();
    }

    //------------------------------------------------------------------------------

    void IQI_Effective_Dynamic_Viscosity::compute_QI( real aWStar )
    {
        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // get the diffusion CM
        const std::shared_ptr< Constitutive_Model > &tCMFluidTurbulence =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::FLUID_TURBULENCE ) );

        // compute effective dynamic viscosity
        mSet->get_QI()( tQIIndex ) += aWStar * ( tCMFluidTurbulence->effective_dynamic_viscosity() );
    }
    //------------------------------------------------------------------------------
}    // namespace moris::fem
