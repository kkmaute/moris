/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Thermal_Energy.cpp
 *
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
            mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Fluid" ] = static_cast< uint >( IQI_Constitutive_Type::FLUID );
        }

        //------------------------------------------------------------------------------

        void IQI_Thermal_Energy::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get the fluid CM
            const std::shared_ptr< Constitutive_Model > & tCMFluid =
                    mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::FLUID ) );

            // get density from CM
            const std::shared_ptr< Property > & tPropDensity =
                    tCMFluid->get_property( "Density" );

            // FIXME protect dof type
            // get velocity field interpolator
            Field_Interpolator * tFIVelocity =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // FIXME protect dof type
            // get temperature field interpolator
            Field_Interpolator * tFITemp =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::TEMP );

            // evaluate the QI
            aQI = tPropDensity->val() * tFITemp->val() * tFIVelocity->val_trans() * mNormal;
        }

        //------------------------------------------------------------------------------

        void IQI_Thermal_Energy::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the fluid CM
            const std::shared_ptr< Constitutive_Model > & tCMFluid =
                    mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::FLUID ) );

            // get density from CM
            const std::shared_ptr< Property > & tPropDensity =
                    tCMFluid->get_property( "Density" );

            // FIXME protect dof type
            // get velocity field interpolator
            Field_Interpolator * tFIVelocity =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // FIXME protect dof type
            // get temperature field interpolator
            Field_Interpolator * tFITemp =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::TEMP );

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * (
                    tPropDensity->val() * tFITemp->val() * tFIVelocity->val_trans() * mNormal );
        }

        //------------------------------------------------------------------------------

        void IQI_Thermal_Energy::compute_dQIdu( real aWStar )
        {
            // get the column index to assemble in residual
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the fluid CM
            const std::shared_ptr< Constitutive_Model > & tCMFluid =
                    mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::FLUID ) );

            // get density from CM
            const std::shared_ptr< Property > & tPropDensity =
                    tCMFluid->get_property( "Density" );

            // FIXME protect dof type
            // get velocity field interpolator
            Field_Interpolator * tFIVelocity =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // FIXME protect dof type
            // get temperature field interpolator
            Field_Interpolator * tFITemp =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::TEMP );

            // get the number of leader dof type dependencies
            uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            // compute dQIdu for indirect dof dependencies
            for( uint iDof = 0; iDof < tNumDofDependencies; iDof++ )
            {
                // get the treated dof type
                Cell< MSI::Dof_Type > & tDofType = mRequestedLeaderGlobalDofTypes( iDof );

                // get leader index for residual dof type, indices for assembly
                uint tLeaderDofIndex      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
                uint tLeaderDepStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

                // if dof type is velocity
                // FIXME protect dof type
                if ( tDofType( 0 ) == MSI::Dof_Type::VX )
                {
                    mSet->get_residual()( tQIIndex )(
                            { tLeaderDepStartIndex, tLeaderDepStopIndex },
                            { 0, 0 } ) += aWStar * (
                                    tPropDensity->val()( 0 ) * tFITemp->val()( 0 ) *
                                    tFIVelocity->N_trans() * mNormal );
                }
                // if dof type is temperature
                // FIXME protect dof type
                else if ( tDofType( 0 ) == MSI::Dof_Type::TEMP )
                {
                    mSet->get_residual()( tQIIndex )(
                            { tLeaderDepStartIndex, tLeaderDepStopIndex },
                            { 0, 0 } ) += aWStar * (
                                    tPropDensity->val()( 0 ) *
                                    dot( tFIVelocity->val(), mNormal ) *
                                    tFITemp->N_trans() );
                }

                // if density depends on dof type
                if ( tPropDensity->check_dof_dependency( tDofType ) )
                {
                    // compute dQIdu
                    mSet->get_residual()( tQIIndex )(
                            { tLeaderDepStartIndex, tLeaderDepStopIndex },
                            { 0, 0 } ) += aWStar * (
                                    tFITemp->val()( 0 ) *
                                    trans( tFIVelocity->val() ) * mNormal *
                                    tPropDensity->dPropdDOF( tDofType ) );
                }
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Thermal_Energy::compute_dQIdu(
                moris::Cell< MSI::Dof_Type > & aDofType,
                Matrix< DDRMat >             & adQIdu )
        {
            // get the fluid CM
            const std::shared_ptr< Constitutive_Model > & tCMFluid =
                    mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::FLUID ) );

            // get density from CM
            const std::shared_ptr< Property > & tPropDensity =
                    tCMFluid->get_property( "Density" );

            // FIXME protect dof type
            // get velocity field interpolator
            Field_Interpolator * tFIVelocity =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // FIXME protect dof type
            // get temperature field interpolator
            Field_Interpolator * tFITemp =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::TEMP );

            // if dof type is velocity
            // FIXME protect dof type
            if ( aDofType( 0 ) == MSI::Dof_Type::VX )
            {
                adQIdu +=
                        tPropDensity->val()( 0 ) * tFITemp->val()( 0 ) *
                        tFIVelocity->N_trans() * mNormal;
            }
            // if dof type is temperature
            // FIXME protect dof type
            else if ( aDofType( 0 ) == MSI::Dof_Type::TEMP )
            {
                adQIdu +=
                        tPropDensity->val()( 0 ) *
                        tFIVelocity->val_trans() * mNormal *
                        trans( tFITemp->N() );
            }

            // if density depends on dof type
            if ( tPropDensity->check_dof_dependency( aDofType ) )
            {
                // compute dQIdu
                adQIdu +=
                        tFITemp->val()( 0 ) *
                        tFIVelocity->val_trans() * mNormal *
                        tPropDensity->dPropdDOF( aDofType );
            }
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */

