/*
 * cl_FEM_IQI_Thermal_Energy.cpp
 *
 *  Created on: Sep 27, 2020
 *      Author: noel
 */
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Thermal_Energy.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Thermal_Energy::IQI_Thermal_Energy()
        {
            // set fem IQI type
            mFEMIQIType = fem::IQI_Type::THERMAL_ENERGY;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Fluid" ] = IQI_Constitutive_Type::FLUID;
        }

        //------------------------------------------------------------------------------

        void IQI_Thermal_Energy::set_constitutive_model(
                std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                std::string                           aConstitutiveString,
                mtk::Master_Slave                     aIsMaster )
        {
            // check that aPropertyString makes sense
            std::string tErrMsg =
                    std::string( "IQI_Thermal_Energy::set_constitutive_model - Unknown aConstitutiveModel: " ) +
                    aConstitutiveString;
            MORIS_ERROR( mConstitutiveMap.find( aConstitutiveString ) != mConstitutiveMap.end(), tErrMsg.c_str() );

            // check no slave allowed
            MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                    "IQI_Thermal_Energy::set_constitutive_model - No slave allowed." );

            // set the property in the property cell
            this->get_constitutive_models( aIsMaster )( static_cast< uint >( mConstitutiveMap[ aConstitutiveString ] ) ) = aConstitutiveModel;
        }

        //------------------------------------------------------------------------------

        void IQI_Thermal_Energy::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get the fluid CM
            std::shared_ptr< Constitutive_Model > tCMFluid =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::FLUID ) );

            // get density from CM
            std::shared_ptr< Property > tPropDensity = tCMFluid->get_property( "Density" );

            // FIXME protect dof type
            // get velocity field interpolator
            Field_Interpolator * tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // FIXME protect dof type
            // get temperature field interpolator
            Field_Interpolator * tFITemp =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::TEMP );

            // evaluate the QI
            aQI = tPropDensity->val() * tFITemp->val() * trans( tFIVelocity->val() ) * mNormal;
        }

        //------------------------------------------------------------------------------

        void IQI_Thermal_Energy::compute_dQIdu(
                moris::Cell< MSI::Dof_Type > & aDofType,
                Matrix< DDRMat >             & adQIdu )
        {
            // get the fluid CM
            std::shared_ptr< Constitutive_Model > tCMFluid =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::FLUID ) );

            // get density from CM
            std::shared_ptr< Property > tPropDensity = tCMFluid->get_property( "Density" );

            // FIXME protect dof type
            // get velocity field interpolator
            Field_Interpolator * tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // FIXME protect dof type
            // get temperature field interpolator
            Field_Interpolator * tFITemp =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::TEMP );

            // if dof type is velocity
            // FIXME protect dof type
            if ( aDofType( 0 ) == MSI::Dof_Type::VX )
            {
                adQIdu +=
                        tPropDensity->val()( 0 ) * tFITemp->val()( 0 ) *
                        trans( tFIVelocity->N() ) * mNormal;
            }

            // if dof type is temperature
            // FIXME protect dof type
            if ( aDofType( 0 ) == MSI::Dof_Type::TEMP )
            {
                adQIdu +=
                        tPropDensity->val()( 0 ) *
                        dot( tFIVelocity->val(), mNormal ) *
                        trans( tFITemp->N() );
            }

            // if density depends on dof type
            if ( tPropDensity->check_dof_dependency( aDofType ) )
            {
                // compute dQIdu
                adQIdu +=
                        tFITemp->val()( 0 ) *
                        trans( tFIVelocity->val() ) * mNormal *
                        tPropDensity->dPropdDOF( aDofType );
            }
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */



