/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Compressible_NS_Density_Bulk.cpp
 *
 */

#include "cl_FEM_IWG_Compressible_NS_Density_Bulk.hpp"

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Compressible_NS_Density_Bulk::IWG_Compressible_NS_Density_Bulk()
        {
            // set size for the constitutive model pointer cell
            mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Fluid" ] = static_cast< uint >( IWG_Constitutive_Type::FLUID );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Density_Bulk::compute_residual( real aWStar )
        {
            // check leader field interpolators
#ifdef MORIS_HAVE_DEBUG
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here velocity), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get the field interpolators
            Field_Interpolator * tDensityFI = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));
            Field_Interpolator * tVelocityFI = mLeaderFIManager->get_field_interpolators_for_type( mDofVelocity );

            // get the constitutive model
            std::shared_ptr< Constitutive_Model > tFluidCM =  mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::FLUID ) );

            // compute the residual weak form
            mSet->get_residual()( 0 )( { tLeaderResStartIndex, tLeaderResStopIndex }, { 0, 0 } ) += aWStar * (
                    trans( tDensityFI->N() ) * tDensityFI->gradt( 1 )
                    - trans( tDensityFI->dnNdxn( 1 ) ) * tDensityFI->val()( 0 ) * tVelocityFI->val() );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Compressible_NS_Density_Bulk::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------
        void IWG_Compressible_NS_Density_Bulk::compute_jacobian( real aWStar )
        {
            // check leader field interpolators
#ifdef MORIS_HAVE_DEBUG
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here pressure), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get the field interpolators
            Field_Interpolator * tDensityFI = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));
            Field_Interpolator * tVelocityFI = mLeaderFIManager->get_field_interpolators_for_type( mDofVelocity );

            // get the constitutive model
            std::shared_ptr< Constitutive_Model > tFluidCM =  mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::FLUID ) );

            // compute the jacobian for dof dependencies
            uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Vector< MSI::Dof_Type > & tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // if dof type is density, add diagonal term (density - density DoF types)
                if( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // add contibution
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * (
                                    trans( tDensityFI->N() ) * tDensityFI->dnNdtn( 1 )  -
                                    trans( tDensityFI->dnNdxn( 1 ) ) * tVelocityFI->val() * tDensityFI->N() );
                }

                // if dof type is velocity, add the mixed term (velocity - density DoF types)
                if( tDofType( 0 ) == mDofVelocity )
                {
                    // compute the jacobian
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * (
                                    - 1.0 * trans( tDensityFI->dnNdxn( 1 ) ) * tDensityFI->val()( 0 ) * tVelocityFI->N() );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Compressible_NS_Density_Bulk::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Density_Bulk::compute_jacobian_and_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Compressible_NS_Density_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Density_Bulk::compute_dRdp( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Compressible_NS_Density_Bulk::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Density_Bulk::compute_residual_strong_form(
                Matrix< DDRMat > & aRM,
                real             & aRC )
        {
            MORIS_ERROR( false, "IWG_Compressible_NS_Density_Bulk::compute_residual_strong_form - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Density_Bulk::compute_jacobian_strong_form(
                moris::Vector< MSI::Dof_Type >   aDofTypes,
                Matrix< DDRMat >             & aJM,
                Matrix< DDRMat >             & aJC )
        {
            MORIS_ERROR( false, "IWG_Compressible_NS_Density_Bulk::compute_jacobian_strong_form - Not implemented." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

