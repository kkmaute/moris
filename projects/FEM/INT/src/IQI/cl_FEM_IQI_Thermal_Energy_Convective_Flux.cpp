/*
 * cl_FEM_IQI_Thermal_Energy_Convective_Flux.cpp
 *
 *  Created on: Sep 27, 2020
 *      Author: noel
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
            mMasterCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Diffusion" ] = static_cast< uint >( IQI_Constitutive_Type::DIFFUSION );
        }
        
        //------------------------------------------------------------------------------

        void IQI_Thermal_Energy_Convective_Flux::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > & tCMDiffusion =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::DIFFUSION ) );

            // FIXME protect dof type
            // get velocity field interpolator
            Field_Interpolator * tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

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
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::DIFFUSION ) );

            // FIXME protect dof type
            // get velocity field interpolator
            Field_Interpolator * tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

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
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::DIFFUSION ) );

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
                Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDof );

                // get master index for residual dof type, indices for assembly
                uint tMasterDofIndex      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
                uint tMasterDepStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

                // if dof type is velocity
                // FIXME protect dof type
                if ( tDofType( 0 ) == MSI::Dof_Type::VX )
                {
                    mSet->get_residual()( tQIIndex )(
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    tCMDiffusion->Energy()(0) * tFIVelocity->N_trans() * mNormal );
                }

                // if density depends on dof type
                if ( tCMDiffusion->check_dof_dependency( tDofType ) )
                {
                    // compute dQIdu
                    mSet->get_residual()( tQIIndex )(
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    dot( tFIVelocity->val(), mNormal ) *
                                    trans( tCMDiffusion->dEnergydDOF( tDofType ) ) );
                }
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Thermal_Energy_Convective_Flux::compute_dQIdu(
                moris::Cell< MSI::Dof_Type > & aDofType,
                Matrix< DDRMat >             & adQIdu )
        {
            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > & tCMDiffusion =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::DIFFUSION ) );

            // FIXME protect dof type
            // get velocity field interpolator
            Field_Interpolator * tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

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

