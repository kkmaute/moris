/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Linear_Interface_SLM_Constraint.cpp
 *
 */

#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Interface_SLM_Constraint.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IWG_Isotropic_Struc_Linear_Interface_SLM_Constraint::IWG_Isotropic_Struc_Linear_Interface_SLM_Constraint()
    {
        // set size for the property pointer cell
        mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );
        mFollowerCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

        // set size for the stabilization parameter pointer cell
        mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );
    }

    //------------------------------------------------------------------------------

    void
    IWG_Isotropic_Struc_Linear_Interface_SLM_Constraint::compute_residual( real aWStar )
    {
#ifdef MORIS_HAVE_DEBUG
            // check leader and follower field interpolators
            this->check_field_interpolators( mtk::Leader_Follower::LEADER );
            this->check_field_interpolators( mtk::Leader_Follower::FOLLOWER );
#endif

            // get leader index for residual dof type, indices for assembly
            const uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            const uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            const uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get follower index for residual dof type, indices for assembly
            const uint tFollowerDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );
            const uint tFollowerResStartIndex = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 0 );
            const uint tFollowerResStopIndex  = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 1 );

            Field_Interpolator* tFIDisplLeader =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get follower field interpolator for the residual dof type
            Field_Interpolator* tFIDisplFollower =
                    mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get leader field interpolator for the residual dof type
            Field_Interpolator* tFILambdaLeader =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get follower field interpolator for the residual dof type
            Field_Interpolator* tFILambdaFollower =
                    mFollowerFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            const auto tTestDisplLeader = tFIDisplLeader->N_trans();
            const auto tTestDisplFollower = tFIDisplFollower->N_trans();

            const auto tWeightedLeaderLambda = 0.5 * tFILambdaLeader->val();
            const auto tWeightedFollowerLambda = 0.5 * tFILambdaFollower->val();

            // compute leader residual
            mSet->get_residual()( 0 )( { tLeaderResStartIndex, tLeaderResStopIndex } ) -=
                    aWStar
                    * ( tTestDisplLeader * ( tWeightedLeaderLambda + tWeightedFollowerLambda ) );

            // compute follower residual
            mSet->get_residual()( 0 )( { tFollowerResStartIndex, tFollowerResStopIndex } ) +=
                    aWStar
                    * ( tTestDisplFollower * ( tWeightedLeaderLambda + tWeightedFollowerLambda ) );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Isotropic_Struc_Linear_Interface_SLM_Constraint::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
		IWG_Isotropic_Struc_Linear_Interface_SLM_Constraint::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader and follower field interpolators
            this->check_field_interpolators( mtk::Leader_Follower::LEADER );
            this->check_field_interpolators( mtk::Leader_Follower::FOLLOWER );
#endif

            // get leader index for residual dof type, indices for assembly
            const uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            const uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            const uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get follower index for residual dof type, indices for assembly
            const uint tFollowerDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );
            const uint tFollowerResStartIndex = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 0 );
            const uint tFollowerResStopIndex  = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 1 );

            Field_Interpolator* tFIDisplLeader =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get follower field interpolator for the residual dof type
            Field_Interpolator* tFIDisplFollower =
                    mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get leader field interpolator for the residual dof type
            Field_Interpolator* tFILambdaLeader =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get follower field interpolator for the residual dof type
            Field_Interpolator* tFILambdaFollower =
                    mFollowerFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // compute the jacobian for indirect dof dependencies through leader constitutive models
            uint tLeaderNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();
            for ( uint iDOF = 0; iDOF < tLeaderNumDofDependencies; iDOF++ )
            {
                // get the dof type
                const Vector< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                // get the index for the dof type
                const sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                const uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                const uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // extract sub-matrices
                auto tJacMM = mSet->get_jacobian()(
                        { tLeaderResStartIndex, tLeaderResStopIndex }, { tLeaderDepStartIndex, tLeaderDepStopIndex } );

                auto tJacSM = mSet->get_jacobian()(
                        { tFollowerResStartIndex, tFollowerResStopIndex }, { tLeaderDepStartIndex, tLeaderDepStopIndex } );

                // compute jacobian direct dependencies
                if ( tDofType( 0 ) == MSI::Dof_Type::VX )
                {
                    tJacMM -= aWStar
                            * ( 0.5 * tFIDisplLeader->N_trans() * tFILambdaLeader->N() );

                    tJacSM += aWStar
                            * ( 0.5 * tFIDisplFollower->N_trans() * tFILambdaLeader->N() );

                }
            }

            // compute the jacobian for indirect dof dependencies through follower constitutive models
            uint tFollowerNumDofDependencies = mRequestedFollowerGlobalDofTypes.size();
            for ( uint iDOF = 0; iDOF < tFollowerNumDofDependencies; iDOF++ )
            {
                // get dof type
                const Vector< MSI::Dof_Type >& tDofType = mRequestedFollowerGlobalDofTypes( iDOF );

                // get the index for the dof type
                const sint tDofDepIndex        = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::FOLLOWER );
                const uint tFollowerDepStartIndex = mSet->get_jac_dof_assembly_map()( tFollowerDofIndex )( tDofDepIndex, 0 );
                const uint tFollowerDepStopIndex  = mSet->get_jac_dof_assembly_map()( tFollowerDofIndex )( tDofDepIndex, 1 );

                // extract sub-matrices
                auto tJacMS = mSet->get_jacobian()(
                        { tLeaderResStartIndex, tLeaderResStopIndex }, { tFollowerDepStartIndex, tFollowerDepStopIndex } );

                auto tJacSS = mSet->get_jacobian()(
                        { tFollowerResStartIndex, tFollowerResStopIndex }, { tFollowerDepStartIndex, tFollowerDepStopIndex } );

                // if dof type is residual dof type
                if ( tDofType( 0 ) == MSI::Dof_Type::VX )
                {
                    tJacMS -= aWStar
                    		* ( 0.5 * tFIDisplLeader->N_trans() * tFILambdaFollower->N() );

                    tJacSS += aWStar
                            * ( 0.5 * tFIDisplFollower->N_trans() * tFILambdaFollower->N() );

                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Isotropic_Struc_Linear_Interface_SLM_Constraint::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
		IWG_Isotropic_Struc_Linear_Interface_SLM_Constraint::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false,
                    "IWG_Isotropic_Struc_Linear_Interface_SLM_Constraint::compute_jacobian_and_residual - This function does nothing." );
        }

        //------------------------------------------------------------------------------

        void
		IWG_Isotropic_Struc_Linear_Interface_SLM_Constraint::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Interface_SLM_Constraint::compute_dRdp - This function does nothing." );
        }

        //------------------------------------------------------------------------------
}    // namespace moris::fem
