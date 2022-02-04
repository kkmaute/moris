/*
 * cl_FEM_IQI_Power_Dissipation.cpp
 *
 *  Created on: May 7, 2021
 *      Author: noel
 */
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Power_Dissipation.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Power_Dissipation::IQI_Power_Dissipation()
        {
            // set fem IQI type
            mFEMIQIType = fem::IQI_Type::POWER_DISSIPATION;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Fluid" ] = static_cast< uint >( IQI_Constitutive_Type::FLUID );
        }
        
        //------------------------------------------------------------------------------

        void IQI_Power_Dissipation::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get the fluid CM
            const std::shared_ptr< Constitutive_Model > & tCMFluid =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::FLUID ) );

            // get density from CM
            const std::shared_ptr< Property > & tPropDensity =
                    tCMFluid->get_property( "Density" );

            // FIXME protect dof type
            // get velocity field interpolator
            Field_Interpolator * tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // evaluate the QI
            aQI = //trans( tCMFluid->traction( mNormal ) ) * tFIVelocity->val() -
                    - 0.5 * tPropDensity->val() * ( tFIVelocity->val_trans() * mNormal ) *
                    ( tFIVelocity->val_trans() * tFIVelocity->val() );
        }

        //------------------------------------------------------------------------------

        void IQI_Power_Dissipation::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the fluid CM
            const std::shared_ptr< Constitutive_Model > & tCMFluid =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::FLUID ) );

            // get density from CM
            const std::shared_ptr< Property > & tPropDensity =
                    tCMFluid->get_property( "Density" );

            // FIXME protect dof type
            // get velocity field interpolator
            Field_Interpolator * tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * (
                    //trans( tCMFluid->traction( mNormal ) ) * tFIVelocity->val() -
                    - 0.5 * tPropDensity->val() * ( tFIVelocity->val_trans() * mNormal ) *
                    ( tFIVelocity->val_trans() * tFIVelocity->val() ) );
        }

        //------------------------------------------------------------------------------

        void IQI_Power_Dissipation::compute_dQIdu( real aWStar )
        {
            // get the column index to assemble in residual
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the fluid CM
            const std::shared_ptr< Constitutive_Model > & tCMFluid =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::FLUID ) );

            // get density from CM
            const std::shared_ptr< Property > & tPropDensity =
                    tCMFluid->get_property( "Density" );

            // FIXME protect dof type
            // get velocity field interpolator
            Field_Interpolator * tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the number of master dof type dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // compute dQIdu for indirect dof dependencies
            for( uint iDof = 0; iDof < tNumDofDependencies; iDof++ )
            {
                // get the treated dof type
                const Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDof );

                // get master index for residual dof type, indices for assembly
                uint tMasterDofIndex      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
                uint tMasterDepStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

//                // if fluid CM depends on dof type
//                if ( tCMFluid->check_dof_dependency( tDofType ) )
//                {
//                    mSet->get_residual()( tQIIndex )(
//                            { tMasterDepStartIndex, tMasterDepStopIndex },
//                            { 0, 0 } ) += aWStar * ( trans( tCMFluid->dTractiondDOF( tDofType, mNormal ) ) * tFIVelocity->val() );
//                }

                // if dof type is velocity
                // FIXME protect dof type
                if ( tDofType( 0 ) == MSI::Dof_Type::VX )
                {
                    mSet->get_residual()( tQIIndex )(
                            { tMasterDepStartIndex, tMasterDepStopIndex },
                            { 0, 0 } ) += aWStar * (
                                    //trans( tFIVelocity->N() ) * tCMFluid->traction( mNormal ) -
                                    - 0.5 * tPropDensity->val()( 0 ) * ( tFIVelocity->N_trans() * mNormal ) *
                                    ( tFIVelocity->val_trans() * tFIVelocity->val() )
                                    - 1.0 * tPropDensity->val()( 0 ) * ( tFIVelocity->N_trans() * tFIVelocity->val() ) *
                                    ( tFIVelocity->val_trans() * mNormal ) );
                }

//                // if density depends on dof type
//                if ( tPropDensity->check_dof_dependency( tDofType ) )
//                {
//                    // compute dQIdu
//                    mSet->get_residual()( tQIIndex )(
//                            { tMasterDepStartIndex, tMasterDepStopIndex },
//                            { 0, 0 } ) -= aWStar * ( 0.5 * tFIVelocity->val_trans() * tFIVelocity->val() * tFIVelocity->val_trans() * mNormal *
//                                    tPropDensity->dPropdDOF( tDofType ) );
//                }
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Power_Dissipation::compute_dQIdu(
                moris::Cell< MSI::Dof_Type > & aDofType,
                Matrix< DDRMat >             & adQIdu )
        {
            // get the fluid CM
            const std::shared_ptr< Constitutive_Model > & tCMFluid =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::FLUID ) );

            // get density from CM
            const std::shared_ptr< Property > & tPropDensity =
                    tCMFluid->get_property( "Density" );

            // FIXME protect dof type
            // get velocity field interpolator
            Field_Interpolator * tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // if dof type is velocity
            // FIXME protect dof type
            if ( aDofType( 0 ) == MSI::Dof_Type::VX )
            {
                adQIdu += //trans( tCMFluid->traction( mNormal ) ) * tFIVelocity->N() -
                        - 0.5 * tPropDensity->val()( 0 ) * ( tFIVelocity->N_trans() * mNormal ) *
                        ( tFIVelocity->val_trans() * tFIVelocity->val() )
                        - tPropDensity->val()( 0 ) * ( tFIVelocity->N_trans() * tFIVelocity->val() ) *
                        ( tFIVelocity->val_trans() * mNormal );
                        //- 1.5 * tPropDensity->val()( 0 ) * tFIVelocity->N_trans() * tFIVelocity->val() * tFIVelocity->val_trans() * mNormal;
            }

//            // if density depends on dof type
//            if ( tPropDensity->check_dof_dependency( aDofType ) )
//            {
//                // compute dQIdu
//                adQIdu -= 0.5 * tFIVelocity->val_trans() * tFIVelocity->val() * tFIVelocity->val_trans() * mNormal *
//                        tPropDensity->dPropdDOF( aDofType );
//            }
//
//            // if fluid CM depends on dof type
//            if ( tCMFluid->check_dof_dependency( aDofType ) )
//            {
//                adQIdu +=  trans( tCMFluid->dTractiondDOF( aDofType, mNormal ) ) * tFIVelocity->val();
//            }
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */



