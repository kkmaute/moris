/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Effective_Conductivity.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Effective_Conductivity.hpp"
#include "fn_norm.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Effective_Conductivity::IQI_Effective_Conductivity()
        {
            // set fem IQI type
            mFEMIQIType = fem::IQI_Type::EFFECTIVE_CONDUCTIVITY;
            init_constitutive_model("Diffusion_Turbulence", IQI_Constitutive_Type::DIFFUSION_TURBULENCE);
        }

        //------------------------------------------------------------------------------

        void IQI_Effective_Conductivity::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > & tCMDiffusionTurbulence = get_leader_constitutive_model(IQI_Constitutive_Type::DIFFUSION_TURBULENCE);

            // compute turbulent dynamic viscosity
            aQI = tCMDiffusionTurbulence->effective_conductivity();
        }

        //------------------------------------------------------------------------------

        void IQI_Effective_Conductivity::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( get_name() );

            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > & tCMDiffusionTurbulence = get_leader_constitutive_model(IQI_Constitutive_Type::DIFFUSION_TURBULENCE);

            // compute turbulent dynamic viscosity
            mSet->get_QI()( tQIIndex ) += aWStar * ( tCMDiffusionTurbulence->effective_conductivity() );
        }
        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */

