
#include "cl_FEM_IQI_Strong_Residual_SA.hpp"

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Strong_Residual_SA::IQI_Strong_Residual_SA()
        {
            // set fem IQI type
            mFEMIQIType = fem::IQI_Type::STRONG_RESIDUAL_SA;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "SpalartAllmarasTurbulence" ] =
                    static_cast< uint >( IQI_Constitutive_Type::TURBULENCE );
        }
        
        //------------------------------------------------------------------------------

        void IQI_Strong_Residual_SA::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get the viscosity FI
            // FIXME protect dof type
            Field_Interpolator* tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VISCOSITY );

            // get the velocity FI
            // FIXME protect dof type
            Field_Interpolator* tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the SA turbulence CM
            const std::shared_ptr< Constitutive_Model >& tCMSATurbulence =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::TURBULENCE ) );

            // compute modified velocity
            Matrix< DDRMat > tModVelocity =
                    tFIVelocity->val() - mCb2 * tFIViscosity->gradx( 1 ) / mSigma;

            // compute strong form of residual
            Matrix< DDRMat> tR = tFIViscosity->gradt( 1 ) +
                    trans( tModVelocity ) * tFIViscosity->gradx( 1 ) -
                    tCMSATurbulence->production_term() +
                    tCMSATurbulence->wall_destruction_term() -
                    tCMSATurbulence->divflux();

            // evaluate the QI
            aQI = tR(0,0);
        }

        //------------------------------------------------------------------------------

        void IQI_Strong_Residual_SA::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // evaluate strong form
            Matrix< DDRMat > tQI(1,1);
            this->compute_QI( tQI );

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * tQI;
        }

        //------------------------------------------------------------------------------

        void IQI_Strong_Residual_SA::compute_dQIdu( real aWStar )
        {
            MORIS_ERROR(false,
                    "IQI_Strong_Residual_SA::compute_dQIdu - not implemented\n.");
        }

        //------------------------------------------------------------------------------

        void IQI_Strong_Residual_SA::compute_dQIdu(
                moris::Cell< MSI::Dof_Type > & aDofType,
                Matrix< DDRMat >             & adQIdu )
        {
            MORIS_ERROR(false,
                    "IQI_Strong_Residual_SA::compute_dQIdu - not implemented\n.");
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */
