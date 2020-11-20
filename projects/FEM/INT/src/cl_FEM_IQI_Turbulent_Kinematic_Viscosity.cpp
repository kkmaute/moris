/*
 * cl_FEM_IQI_Turbulent_Kinematic_Viscosity.cpp
 *
 *  Created on: Jul 20, 2020
 *      Author: noel
 */
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Turbulent_Kinematic_Viscosity.hpp"
#include "fn_norm.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Turbulent_Kinematic_Viscosity::IQI_Turbulent_Kinematic_Viscosity()
        {
            // set fem IQI type
            mFEMIQIType = fem::IQI_Type::TURBULENT_KINEMATIC_VISCOSITY;

            // set the property pointer cell size
            mMasterProp.resize( static_cast< uint >( Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "Density" ]          = static_cast< uint >( Property_Type::DENSITY );
            mPropertyMap[ "DynamicViscosity" ] = static_cast< uint >( Property_Type::DYNAMIC_VISCOSITY );
        }

        //------------------------------------------------------------------------------

        void IQI_Turbulent_Kinematic_Viscosity::set_dof_type_list(
                moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                moris::Cell< std::string >                  & aDofStrings,
                mtk::Master_Slave                             aIsMaster )
        {
            // switch on master slave
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER :
                {
                    // set dof type list
                    mMasterDofTypes = aDofTypes;

                    // loop on dof type
                    for( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
                    {
                        // get dof string
                        std::string tDofString = aDofStrings( iDof );

                        // get dof type
                        MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

                        // if viscosity
                        if( tDofString == "Viscosity" )
                        {
                            mMasterDofViscosity = tDofType;
                        }
                        else
                        {
                            // create error message
                            std::string tErrMsg =
                                    std::string( "IQI_Turbulent_Kinematic_Viscosity::set_dof_type_list - Unknown aDofString : ") +
                                    tDofString;
                            MORIS_ERROR( false , tErrMsg.c_str() );
                        }
                    }
                    break;
                }

                case mtk::Master_Slave::SLAVE :
                {
                    // set dof type list
                    mSlaveDofTypes = aDofTypes;
                    break;
                }

                default:
                    MORIS_ERROR( false, "IQI_Turbulent_Kinematic_Viscosity::set_dof_type_list - unknown master slave type." );
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Turbulent_Kinematic_Viscosity::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get field interpolator for modified viscosity
            Field_Interpolator * tFIModViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );

            // get the density property
            std::shared_ptr< Property > & tPropDensity =
                    mMasterProp( static_cast< uint >( Property_Type::DENSITY ) );

            // compute fv1
            real tFv1 = this->compute_fv1();

            // compute turbulent kinematic viscosity
            aQI = tPropDensity->val() * tFIModViscosity->val() * tFv1;
        }

        //------------------------------------------------------------------------------

        void IQI_Turbulent_Kinematic_Viscosity::compute_dQIdu(
                moris::Cell< MSI::Dof_Type > & aDofType,
                Matrix< DDRMat >             & adQIdu )
        {
            MORIS_ERROR( false, "compute_dQIdu() not implemented for turbulent kinematic viscosity IQI." );
        }

        //------------------------------------------------------------------------------

        real IQI_Turbulent_Kinematic_Viscosity::compute_fv1()
        {
            // compute chi
            real tChi = this->compute_chi();

            // compute chi³
            real tChi3 = std::pow( tChi, 3.0 );

            // compute cv1³
            real tCv13 = std::pow( mCv1, 3.0 );

            // compute fv1
            return tChi3 / ( tChi3 + tCv13 );
        }

        //------------------------------------------------------------------------------

        void IQI_Turbulent_Kinematic_Viscosity::compute_dfv1du(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adfv1du )
        {
            // compute chi
            real tChi = this->compute_chi();

            // compute chi²
            real tChi2 = std::pow( tChi, 2.0 );

            // compute chi³
            real tChi3 = std::pow( tChi, 3.0 );

            // compute cv1³
            real tCv13 = std::pow( mCv1, 3.0 );

            // compute dchidu
            Matrix< DDRMat > tdchidu;
            this->compute_dchidu( aDofTypes, tdchidu );

            // compute adfv1du
            adfv1du = 3.0 * tCv13 * tChi2 * tdchidu / std::pow( tChi3 + tCv13, 2.0 );
        }

        //------------------------------------------------------------------------------

        real IQI_Turbulent_Kinematic_Viscosity::compute_chi()
        {
            // get the modified viscosity dof FI
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );

            // get the density and gravity properties
            std::shared_ptr< Property > & tPropDynViscosity =
                    mMasterProp( static_cast< uint >( Property_Type::DYNAMIC_VISCOSITY ) );

            // compute chi
            return tFIViscosity->val()( 0 ) / tPropDynViscosity->val()( 0 );
        }

        //------------------------------------------------------------------------------
        void IQI_Turbulent_Kinematic_Viscosity::compute_dchidu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adchidu )
        {
            // get the derivative dof FIs
            Field_Interpolator * tDerFI =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init adchidu
            adchidu.set_size( 1, tDerFI->get_number_of_space_time_coefficients(), 0.0 );

            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIModViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );

            // get the density and gravity properties
            std::shared_ptr< Property > & tPropDynViscosity =
                    mMasterProp( static_cast< uint >( Property_Type::DYNAMIC_VISCOSITY ) );

            // if dof type is viscosity
            if( aDofTypes( 0 ) == mMasterDofViscosity )
            {
                adchidu += tDerFI->N() / tPropDynViscosity->val()( 0 );
            }

            // if viscosity property depends on dof type
            if( tPropDynViscosity->check_dof_dependency( aDofTypes ) )
            {
                adchidu -=
                        tFIModViscosity->val() * tPropDynViscosity->dPropdDOF( aDofTypes ) /
                        std::pow( tPropDynViscosity->val()( 0 ), 2 );
            }
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */



