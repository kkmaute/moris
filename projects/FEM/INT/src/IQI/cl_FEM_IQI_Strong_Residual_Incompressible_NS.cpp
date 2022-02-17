
#include "cl_FEM_IQI_Strong_Residual_Incompressible_NS.hpp"

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Strong_Residual_Incompressible_NS::IQI_Strong_Residual_Incompressible_NS()
        {
            // set fem IQI type
            mFEMIQIType = fem::IQI_Type::STRONG_RESIDUAL_INCOMPRESSIBLE_NS;

            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Gravity" ]          = static_cast< uint >( IQI_Property_Type::GRAVITY );
            mPropertyMap[ "ThermalExpansion" ] = static_cast< uint >( IQI_Property_Type::THERMAL_EXPANSION );
            mPropertyMap[ "ReferenceTemp" ]    = static_cast< uint >( IQI_Property_Type::REF_TEMP );
            mPropertyMap[ "InvPermeability" ]  = static_cast< uint >( IQI_Property_Type::INV_PERMEABILITY );
            mPropertyMap[ "MassSource" ]       = static_cast< uint >( IQI_Property_Type::MASS_SOURCE );
            mPropertyMap[ "Load" ]             = static_cast< uint >( IQI_Property_Type::BODY_LOAD );
            mPropertyMap[ "Homotopy" ]         = static_cast< uint >( IQI_Property_Type::HOMOTOPY );

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "IncompressibleFluid" ] =
                    static_cast< uint >( IQI_Constitutive_Type::INCOMPRESSIBLE_FLUID );
        }
        
        //------------------------------------------------------------------------------

        void IQI_Strong_Residual_Incompressible_NS::compute_QI( Matrix< DDRMat > & aQI )
        {
            // check if index was set
            if( mQuantityDofType.size() > 1 )
            {
                MORIS_ERROR( mIQITypeIndex != -1, "IQI_Strong_Residual_Incompressible_NS::compute_QI - mIQITypeIndex not set." );
            }
            else
            {
                mIQITypeIndex = 0;
            }

            // get the velocity and pressure FIs
            Field_Interpolator * tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the properties
            const std::shared_ptr< Property > & tGravityProp =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::GRAVITY ) );

            const std::shared_ptr< Property > & tThermalExpProp =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::THERMAL_EXPANSION ) );

            const std::shared_ptr< Property > & tRefTempProp =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::REF_TEMP ) );

            const std::shared_ptr< Property > & tInvPermeabProp   =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::INV_PERMEABILITY ) );

            const std::shared_ptr< Property > & tMassSourceProp   =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::MASS_SOURCE ) );

            // get the body load property
            const std::shared_ptr< Property > & tLoadProp =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::BODY_LOAD ) );

            // get the incompressible fluid constitutive model
            const std::shared_ptr< Constitutive_Model > & tIncFluidCM =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::INCOMPRESSIBLE_FLUID ) );

            // get the density property from CM
            const std::shared_ptr< Property > & tDensityProp =
                    tIncFluidCM->get_property( "Density" );

            // get the density value
            real tDensity = tDensityProp->val()( 0 );

            if( mIQITypeIndex == 0 )
            {
                // compute the residual strong form of momentum equation
                aQI = tDensity * trans( tVelocityFI->gradt( 1 ) ) +
                        tDensity * trans( tVelocityFI->gradx( 1 ) ) * tVelocityFI->val() -
                        tIncFluidCM->divflux();

                // if body load
                if ( tLoadProp != nullptr )
                {
                    // add contribution of body load term to momentum residual
                    aQI -= tLoadProp->val();
                }

                // if permeability
                if ( tInvPermeabProp != nullptr )
                {
                    // add Brinkman term to residual strong form
                    aQI += tInvPermeabProp->val()( 0 ) * tVelocityFI->val();
                }

                // if gravity
                if ( tGravityProp != nullptr )
                {
                    // add gravity to residual strong form
                    aQI += tDensity * tGravityProp->val();

                    // if thermal expansion and reference temperature
                    if( tThermalExpProp != nullptr && tRefTempProp != nullptr )
                    {
                        // get the temperature field interpolator
                        // FIXME protect FI
                        Field_Interpolator * tTempFI =
                                mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::TEMP );

                        // add contribution to residual
                        aQI -= tDensity * tGravityProp->val() * tThermalExpProp->val() *
                                ( tTempFI->val() - tRefTempProp->val() );
                    }
                }
            }
            else if( mIQITypeIndex == 1 )
            {
                // compute the residual strong form of continuity equation
                aQI = tVelocityFI->div();

                // if mass source
                if ( tMassSourceProp != nullptr )
                {
                    // add mass source to continuity residual
                    aQI -= tMassSourceProp->val()( 0 ) / tDensity;
                }
            }
            else
            {
                aQI = MORIS_REAL_MAX;
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Strong_Residual_Incompressible_NS::compute_QI( real aWStar )
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

        void IQI_Strong_Residual_Incompressible_NS::compute_dQIdu( real aWStar )
        {
            MORIS_ERROR(false,
                    "IQI_Strong_Residual_Incompressible_NS::compute_dQIdu - not implemented\n.");
        }

        //------------------------------------------------------------------------------

        void IQI_Strong_Residual_Incompressible_NS::compute_dQIdu(
                moris::Cell< MSI::Dof_Type > & aDofType,
                Matrix< DDRMat >             & adQIdu )
        {
            MORIS_ERROR(false,
                    "IQI_Strong_Residual_Incompressible_NS::compute_dQIdu - not implemented\n.");
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */
