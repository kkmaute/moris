/*
 * cl_FEM_IQI_Turbulent_Kinematic_Viscosity.cpp
 *
 *  Created on: Jul 20, 2020
 *      Author: noel
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

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Fluid_Turbulence" ] =
                    static_cast< uint >( IQI_Constitutive_Type::FLUID_TURBULENCE );
        }
        
        //------------------------------------------------------------------------------

        void IQI_Turbulent_Dynamic_Viscosity::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > & tCMFluidTurbulence =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::FLUID_TURBULENCE ) );

            // compute turbulent dynamic viscosity
            aQI = tCMFluidTurbulence->turbulent_dynamic_viscosity();
        }
        
        //------------------------------------------------------------------------------

        void IQI_Turbulent_Dynamic_Viscosity::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > & tCMFluidTurbulence =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::FLUID_TURBULENCE ) );

            // compute turbulent dynamic viscosity
            mSet->get_QI()( tQIIndex ) += aWStar * ( tCMFluidTurbulence->turbulent_dynamic_viscosity() );
        }
        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */



