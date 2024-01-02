/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Thermal_Energy_Convective_Flux.cpp
 *
 */

#include "cl_FEM_IQI_Thermal_Energy_Convective_Flux.hpp"

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Thermal_Energy_Convective_Flux::IQI_Thermal_Energy_Convective_Flux()
        {
            // set fem IQI type
            mFEMIQIType = fem::IQI_Type::THERMAL_ENERGY_CONVECTIVE_FLUX;

            // set size for the constitutive model pointer cell
            mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Diffusion" ] = static_cast< uint >( IQI_Constitutive_Type::DIFFUSION );
        }

        //------------------------------------------------------------------------------

        void IQI_Thermal_Energy_Convective_Flux::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > & tCMDiffusion =
                    mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::DIFFUSION ) );

            // FIXME protect dof type
            // get velocity field interpolator
            Field_Interpolator * tFIVelocity =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // evaluate the QI
            aQI = tCMDiffusion->Energy()(0) * tFIVelocity->val_trans() * mNormal;
        }

        //------------------------------------------------------------------------------

        void IQI_Thermal_Energy_Convective_Flux::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > & tCMDiffusion =
                    mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::DIFFUSION ) );

            // FIXME protect dof type
            // get velocity field interpolator
            Field_Interpolator * tFIVelocity =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * (
                    tCMDiffusion->Energy()(0) * tFIVelocity->val_trans() * mNormal );
        }

        //------------------------------------------------------------------------------

        void IQI_Thermal_Energy_Convective_Flux::compute_dQIdu( real aWStar )
        {
            // get the column index to assemble in residual
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > & tCMDiffusion =
                    mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::DIFFUSION ) );

            // FIXME protect dof type
            // get velocity field interpolator
            Field_Interpolator * tFIVelocity =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the number of leader dof type dependencies
            uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            // compute dQIdu for indirect dof dependencies
            for( uint iDof = 0; iDof < tNumDofDependencies; iDof++ )
            {
                // get the treated dof type
                Vector< MSI::Dof_Type > & tDofType = mRequestedLeaderGlobalDofTypes( iDof );

                // get leader index for residual dof type, indices for assembly
                uint tLeaderDofIndex      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
                uint tLeaderDepStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

                // if dof type is velocity
                // FIXME protect dof type
                if ( tDofType( 0 ) == MSI::Dof_Type::VX )
                {
                    mSet->get_residual()( tQIIndex )(
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * (
                                    tCMDiffusion->Energy()(0) * tFIVelocity->N_trans() * mNormal );
                }

                // if density depends on dof type
                if ( tCMDiffusion->check_dof_dependency( tDofType ) )
                {
                    // compute dQIdu
                    mSet->get_residual()( tQIIndex )(
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * (
                                    dot( tFIVelocity->val(), mNormal ) *
                                    trans( tCMDiffusion->dEnergydDOF( tDofType ) ) );
                }
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Thermal_Energy_Convective_Flux::compute_dQIdu(
                Vector< MSI::Dof_Type > & aDofType,
                Matrix< DDRMat >             & adQIdu )
        {
            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > & tCMDiffusion =
                    mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::DIFFUSION ) );

            // FIXME protect dof type
            // get velocity field interpolator
            Field_Interpolator * tFIVelocity =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // if dof type is velocity
            // FIXME protect dof type
            if ( aDofType( 0 ) == MSI::Dof_Type::VX )
            {
                adQIdu += tCMDiffusion->Energy()(0) * tFIVelocity->N_trans() * mNormal;
            }

            // if density depends on dof type
            if ( tCMDiffusion->check_dof_dependency( aDofType ) )
            {
                // compute dQIdu
                adQIdu += dot( tFIVelocity->val(), mNormal ) * trans( tCMDiffusion->dEnergydDOF( aDofType ) );
            }
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */

