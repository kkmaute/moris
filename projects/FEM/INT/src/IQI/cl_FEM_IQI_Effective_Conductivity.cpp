
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

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Diffusion_Turbulence" ] =
                    static_cast< uint >( IQI_Constitutive_Type::DIFFUSION_TURBULENCE );
        }
        
        //------------------------------------------------------------------------------

        void IQI_Effective_Conductivity::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > & tCMDiffusionTurbulence =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::DIFFUSION_TURBULENCE ) );

            // compute turbulent dynamic viscosity
            aQI = tCMDiffusionTurbulence->effective_conductivity();
        }
        

        //------------------------------------------------------------------------------

        void IQI_Effective_Conductivity::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > & tCMDiffusionTurbulence =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::DIFFUSION_TURBULENCE ) );

            // compute turbulent dynamic viscosity
            mSet->get_QI()( tQIIndex ) += aWStar * ( tCMDiffusionTurbulence->effective_conductivity() );
        }
        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */



