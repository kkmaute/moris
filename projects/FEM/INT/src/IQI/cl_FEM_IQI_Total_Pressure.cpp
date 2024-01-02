/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Total_Pressure.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Total_Pressure.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Total_Pressure::IQI_Total_Pressure()
        {
            // set fem IQI type
            mFEMIQIType = fem::IQI_Type::TOTAL_PRESSURE;

            // set size for the constitutive model pointer cell
            mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Fluid" ] = static_cast< uint >( IQI_Constitutive_Type::FLUID );
        }

        //------------------------------------------------------------------------------

        void IQI_Total_Pressure::compute_QI( Matrix< DDRMat > & aQI )
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
            // get pressure field interpolator
            Field_Interpolator * tFIPressure =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::P );

            // evaluate the QI
            aQI = tFIPressure->val() +
                    0.5 * tPropDensity->val() * tFIVelocity->val_trans() * tFIVelocity->val();
        }

        //------------------------------------------------------------------------------

        void IQI_Total_Pressure::compute_QI( real aWStar )
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
            // get pressure field interpolator
            Field_Interpolator * tFIPressure =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::P );

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * ( tFIPressure->val() +
                    0.5 * tPropDensity->val() * tFIVelocity->val_trans() * tFIVelocity->val() );
        }

        //------------------------------------------------------------------------------

        void IQI_Total_Pressure::compute_dQIdu( real aWStar )
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
            // get pressure field interpolator
            Field_Interpolator * tFIPressure =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::P );

            // get the number of leader dof type dependencies
            uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            // compute dQIdu for indirect dof dependencies
            for( uint iDof = 0; iDof < tNumDofDependencies; iDof++ )
            {
                // get the treated dof type
                const Vector< MSI::Dof_Type > & tDofType = mRequestedLeaderGlobalDofTypes( iDof );

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
                                    tPropDensity->val()( 0 ) * tFIVelocity->N_trans() * tFIVelocity->val() );
                }
                // if dof type is pressure
                // FIXME protect dof type
                else if ( tDofType( 0 ) == MSI::Dof_Type::P )
                {
                    mSet->get_residual()( tQIIndex )(
                            { tLeaderDepStartIndex, tLeaderDepStopIndex },
                            { 0, 0 } ) += aWStar * ( tFIPressure->N_trans() );
                }

                // if density depends on dof type
                if ( tPropDensity->check_dof_dependency( tDofType ) )
                {
                    // compute dQIdu
                    mSet->get_residual()( tQIIndex )(
                            { tLeaderDepStartIndex, tLeaderDepStopIndex },
                            { 0, 0 } ) += aWStar * ( 0.5 * tFIVelocity->val_trans() * tFIVelocity->val() *
                                    tPropDensity->dPropdDOF( tDofType ) );
                }
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Total_Pressure::compute_dQIdu(
                Vector< MSI::Dof_Type > & aDofType,
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
            // get pressure field interpolator
            Field_Interpolator * tFIPressure =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::P );

            // if dof type is velocity
            // FIXME protect dof type
            if ( aDofType( 0 ) == MSI::Dof_Type::VX )
            {
                adQIdu += tPropDensity->val()( 0 ) * tFIVelocity->N_trans() * tFIVelocity->val();
            }
            // if dof type is pressure
            // FIXME protect dof type
            else if ( aDofType( 0 ) == MSI::Dof_Type::P )
            {
                adQIdu += tFIPressure->N_trans();
            }

            // if density depends on dof type
            if ( tPropDensity->check_dof_dependency( aDofType ) )
            {
                // compute dQIdu
                adQIdu += 0.5 * tFIVelocity->val_trans() * tFIVelocity->val() *
                        tPropDensity->dPropdDOF( aDofType );
            }
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */

