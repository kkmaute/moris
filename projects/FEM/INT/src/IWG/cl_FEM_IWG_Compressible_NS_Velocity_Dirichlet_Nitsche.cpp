/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Compressible_NS_Velocity_Dirichlet_Nitsche.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Compressible_NS_Velocity_Dirichlet_Nitsche.hpp"
// LINALG/src
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Compressible_NS_Velocity_Dirichlet_Nitsche::IWG_Compressible_NS_Velocity_Dirichlet_Nitsche( sint aBeta )
        {
            // set sign for symmetric/unsymmetric Nitsche
            mBeta = aBeta;
            init_property( "PrescribedValue", IWG_Property_Type::PRESCRIBED_VALUE );
            init_property( "Select", IWG_Property_Type::SELECT );
            init_constitutive_model( "Fluid", IWG_Constitutive_Type::FLUID );
            init_stabilization_parameter( "NitschePenaltyParameter", IWG_Stabilization_Type::NITSCHE_PENALTY_PARAMETER );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Velocity_Dirichlet_Nitsche::compute_residual( real aWStar )
        {
            // check leader field interpolators
#ifdef MORIS_HAVE_DEBUG
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get the leader field interpolator for the residual dof type
            Field_Interpolator *tFIVelocity = get_leader_fi_manager()->get_field_interpolators_for_type( mDofVelocity );

            // get the selection matrix property
            std::shared_ptr< Property > const &tPropSelect = get_leader_property( IWG_Property_Type::SELECT );

            // set a default selection matrix if needed
            Matrix< DDRMat > tM;
            if ( tPropSelect == nullptr )
            {
                // get spatial dimension of velocity field
                uint tSpaceDim = tFIVelocity->get_number_of_fields();

                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = tPropSelect->val();
            }

            // get the imposed velocity property
            std::shared_ptr< Property > const &tPropVelocity = get_leader_property( IWG_Property_Type::PRESCRIBED_VALUE );

            // get the fluid constitutive model
            std::shared_ptr< Constitutive_Model > const &tCMFluid = get_leader_constitutive_model( IWG_Constitutive_Type::FLUID );

            // get the Nitsche stabilization parameter
            std::shared_ptr< Stabilization_Parameter > const &tSPNitsche = get_stabilization_parameter( IWG_Stabilization_Type::NITSCHE_PENALTY_PARAMETER );

            // compute the jump
            Matrix< DDRMat > tVelocityJump = tFIVelocity->val() - tPropVelocity->val();

            // compute leader residual
            mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex },
                    { 0, 0 } ) -= aWStar * ( mBeta * trans( tCMFluid->testTraction( get_normal(), mResidualDofType( 0 ), CM_Function_Type::MECHANICAL ) ) * tM * tVelocityJump );

            // if residual dof type is velocity
            if ( mResidualDofType( 0 )( 0 ) == mDofVelocity )
            {
                // compute leader residual
                mSet->get_residual()( 0 )(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { 0, 0 } ) -= aWStar * ( tFIVelocity->N_trans() * tM * tCMFluid->traction( get_normal(), CM_Function_Type::MECHANICAL ) + tSPNitsche->val()( 0 ) * tFIVelocity->N_trans() * tM * tVelocityJump );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Compressible_NS_Velocity_Dirichlet_Nitsche::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Velocity_Dirichlet_Nitsche::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get the leader field interpolator for residual dof type
            Field_Interpolator *tFIVelocity = get_leader_fi_manager()->get_field_interpolators_for_type( mDofVelocity );

            // get the selection matrix property
           std::shared_ptr< Property > const &tPropSelect = get_leader_property( IWG_Property_Type::SELECT );

            // set a default selection matrix if needed
            Matrix< DDRMat > tM;
            if ( tPropSelect == nullptr )
            {
                // get spatial dimension of velocity field
                uint tSpaceDim = tFIVelocity->get_dof_type().size();

                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = tPropSelect->val();
            }

            // get the imposed velocity property
            std::shared_ptr< Property > const &tPropVelocity = get_leader_property( IWG_Property_Type::PRESCRIBED_VALUE );

            // get the fluid constitutive model
            std::shared_ptr< Constitutive_Model > const &tCMFluid = get_leader_constitutive_model( IWG_Constitutive_Type::FLUID );

            // get the Nitsche stabilization parameter
            std::shared_ptr< Stabilization_Parameter > const &tSPNitsche = get_stabilization_parameter( IWG_Stabilization_Type::NITSCHE_PENALTY_PARAMETER );

            // compute the jump
            Matrix< DDRMat > tVelocityJump = tFIVelocity->val() - tPropVelocity->val();

            // get number of leader dependencies
            uint tLeaderNumDofDependencies = get_requested_leader_dof_types().size();

            // compute the jacobian for indirect dof dependencies through leader
            for ( uint iDOF = 0; iDOF < tLeaderNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Vector< MSI::Dof_Type > const &tDofType = get_requested_leader_dof_types()( iDOF );

                // get the index for the dof type
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // test-traction contribution
                mSet->get_jacobian()(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -= aWStar * ( mBeta * tCMFluid->dTestTractiondDOF( tDofType, get_normal(), tM * tVelocityJump, mResidualDofType( 0 ), CM_Function_Type::MECHANICAL ) );

                //---------------------------------------------------------------------
                // if residual dof type is velocity
                if ( mResidualDofType( 0 )( 0 ) == mDofVelocity )
                {
                    // compute leader residual
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -= aWStar * ( tFIVelocity->N_trans() * tM * tCMFluid->dTractiondDOF( tDofType, get_normal(), CM_Function_Type::MECHANICAL ) );

                    // if dof type is velocity
                    if ( tDofType( 0 ) == mDofVelocity )
                    {
                        // compute jacobian direct dependencies
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -= aWStar * ( mBeta * trans( tCMFluid->testTraction( get_normal(), mResidualDofType( 0 ), CM_Function_Type::MECHANICAL ) ) * tM * tFIVelocity->N() + tSPNitsche->val()( 0 ) * tFIVelocity->N_trans() * tM * tFIVelocity->N() );
                    }

                    // if imposed velocity depends on dof type
                    if ( tPropVelocity->check_dof_dependency( tDofType ) )
                    {
                        // add contribution from property to jacobian
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -= aWStar * ( mBeta * trans( tCMFluid->testTraction( get_normal(), mResidualDofType( 0 ), CM_Function_Type::MECHANICAL ) ) * tM * tPropVelocity->dPropdDOF( tDofType ) + tSPNitsche->val()( 0 ) * tFIVelocity->N_trans() * tM * tPropVelocity->dPropdDOF( tDofType ) );
                    }

                    // if stabilization parameter depends on the dof type
                    if ( tSPNitsche->check_dof_dependency( tDofType ) )
                    {
                        // add contribution of SP to jacobian
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -= aWStar * ( tFIVelocity->N_trans() * tM * tVelocityJump * tSPNitsche->dSPdLeaderDOF( tDofType ) );
                    }
                }
                //---------------------------------------------------------------------
                else    // residual dof type is density or temperature
                {
                    // if dof type is velocity
                    if ( tDofType( 0 ) == mDofVelocity )
                    {
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -= aWStar * ( mBeta * trans( tCMFluid->testTraction( get_normal(), mResidualDofType( 0 ), CM_Function_Type::MECHANICAL ) ) * tM * tFIVelocity->N() );
                    }

                    // if imposed velocity depends on dof type
                    if ( tPropVelocity->check_dof_dependency( tDofType ) )
                    {
                        // add contribution from property to jacobian
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -= aWStar * ( mBeta * trans( tCMFluid->testTraction( get_normal(), mResidualDofType( 0 ), CM_Function_Type::MECHANICAL ) ) * tM * tPropVelocity->dPropdDOF( tDofType ) );
                    }
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Compressible_NS_Velocity_Dirichlet_Nitsche::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Velocity_Dirichlet_Nitsche::compute_jacobian_and_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Compressible_NS_Velocity_Dirichlet_Nitsche::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Velocity_Dirichlet_Nitsche::compute_dRdp( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Compressible_NS_Velocity_Dirichlet_Nitsche::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
