/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Spalart_Allmaras_Turbulence_Interface.cpp
 *
 */

#include "cl_FEM_IWG_Spalart_Allmaras_Turbulence_Interface.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
//LINALG/src
#include "fn_trans.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Spalart_Allmaras_Turbulence_Interface::IWG_Spalart_Allmaras_Turbulence_Interface( sint aBeta )
        {
            // set beta for symmetric/unsymmetric Nitsche formulation
            mBeta = aBeta;
            init_constitutive_model("SpalartAllmarasTurbulence", IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE);
            init_stabilization_parameter("NitscheInterface", IWG_Stabilization_Type::NITSCHE_INTERFACE);
        }

        //------------------------------------------------------------------------------

        void IWG_Spalart_Allmaras_Turbulence_Interface::compute_residual( real aWStar )
        {
            // check leader field interpolators
#ifdef MORIS_HAVE_DEBUG
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get follower index for residual dof type, indices for assembly
            uint tFollowerDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );
            uint tFollowerResStartIndex = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 0 );
            uint tFollowerResStopIndex  = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 1 );

            // get leader field interpolator for the residual dof type
            Field_Interpolator * tFILeader = get_leader_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get follower field interpolator for the residual dof type
            Field_Interpolator * tFIFollower  = get_follower_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the SA turbulence CM
            const std::shared_ptr< Constitutive_Model > & tCMLeaderSATurbulence = get_leader_constitutive_model(IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE);
            const std::shared_ptr< Constitutive_Model > &tCMFollowerSATurbulence = get_follower_constitutive_model(IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE);

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPNitsche = get_stabilization_parameter(IWG_Stabilization_Type::NITSCHE_INTERFACE);
            real tNitsche      = tSPNitsche->val()( 0 );
            real tLeaderWeight = tSPNitsche->val()( 1 );
            real tFollowerWeight  = tSPNitsche->val()( 2 );

            // evaluate average traction
            Matrix< DDRMat > tTraction =
                    tLeaderWeight * tCMLeaderSATurbulence->traction( get_normal() ) +
                    tFollowerWeight  * tCMFollowerSATurbulence->traction( get_normal() );

            // evaluate temperature jump
            Matrix< DDRMat > tJumpViscosity = tFILeader->val() - tFIFollower->val();

            // compute leader residual
            mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex },
                    { 0, 0 } ) += aWStar * (
                            - tFILeader->N_trans() * tTraction
                            - mBeta * tLeaderWeight * tCMLeaderSATurbulence->testTraction( get_normal(), mResidualDofType( 0 ) ) * tJumpViscosity
                            + tNitsche * tFILeader->N_trans() * tJumpViscosity ) ;

            // compute follower residual
            mSet->get_residual()( 0 )(
                    { tFollowerResStartIndex, tFollowerResStopIndex },
                    { 0, 0 } ) += aWStar * (
                            + tFIFollower->N_trans() * tTraction
                            - mBeta * tFollowerWeight * tCMFollowerSATurbulence->testTraction( get_normal(), mResidualDofType( 0 ) ) * tJumpViscosity
                            - tNitsche * tFIFollower->N_trans() * tJumpViscosity );

            // check for nan, infinity
            MORIS_ASSERT( isfinite(  mSet->get_residual()( 0 ) ),
                    "IWG_Spalart_Allmaras_Turbulence_Interface::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Interface::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get follower index for residual dof type, indices for assembly
            uint tFollowerDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );
            uint tFollowerResStartIndex = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 0 );
            uint tFollowerResStopIndex  = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 1 );

            // get leader field interpolator for the residual dof type
            Field_Interpolator * tFILeader = get_leader_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get follower field interpolator for the residual dof type
            Field_Interpolator * tFIFollower  = get_follower_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the SA turbulence CM
            const std::shared_ptr< Constitutive_Model > & tCMLeaderSATurbulence = get_leader_constitutive_model(IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE);
            const std::shared_ptr< Constitutive_Model > &tCMFollowerSATurbulence = get_follower_constitutive_model(IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE);

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPNitsche = get_stabilization_parameter(IWG_Stabilization_Type::NITSCHE_INTERFACE);
            real tNitsche      = tSPNitsche->val()( 0 );
            real tLeaderWeight = tSPNitsche->val()( 1 );
            real tFollowerWeight  = tSPNitsche->val()( 2 );

            // evaluate average traction
            Matrix< DDRMat > tTraction =
                    tLeaderWeight * tCMLeaderSATurbulence->traction( get_normal() ) +
                    tFollowerWeight  * tCMFollowerSATurbulence->traction( get_normal() );

            // evaluate temperature jump
            Matrix< DDRMat > tJumpViscosity = tFILeader->val() - tFIFollower->val();

            // get number of leader dof dependencies
            uint tLeaderNumDofDependencies = get_requested_leader_dof_types().size();

            // compute the jacobian for indirect dof dependencies through leader constitutive models
            for( uint iDOF = 0; iDOF < tLeaderNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Vector< MSI::Dof_Type > const &tDofType = get_requested_leader_dof_types()( iDOF );

                // get the index for the dof type
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // compute jacobian direct dependencies
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * (
                                    - mBeta * tLeaderWeight * tCMLeaderSATurbulence->testTraction( get_normal(), mResidualDofType( 0 ) ) * tFILeader->N()
                                    + tNitsche * tFILeader->N_trans() * tFILeader->N() );

                    mSet->get_jacobian()(
                            { tFollowerResStartIndex,  tFollowerResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * (
                                    - mBeta * tFollowerWeight * tCMFollowerSATurbulence->testTraction( get_normal(), mResidualDofType( 0 ) ) * tFILeader->N()
                                    - tNitsche * tFIFollower->N_trans() * tFILeader->N() );
                }

                // if dependency of turbulence CM on the dof type
                if ( tCMLeaderSATurbulence->check_dof_dependency( tDofType ) )
                {
                    // add contribution from the derivative of the traction to jacobian
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * (
                                    - tFILeader->N_trans() * tLeaderWeight * tCMLeaderSATurbulence->dTractiondDOF( tDofType, get_normal() ) )
                                    - mBeta * tLeaderWeight * tCMLeaderSATurbulence->dTestTractiondDOF( tDofType, get_normal(), mResidualDofType( 0 ) ) * tJumpViscosity( 0 );

                    mSet->get_jacobian()(
                            { tFollowerResStartIndex,  tFollowerResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * (
                                    + tFIFollower->N_trans() * tLeaderWeight * tCMLeaderSATurbulence->dTractiondDOF( tDofType, get_normal() ) );
                }

                // if dependency of stabilization parameters on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Leader_Follower::LEADER ) )
                {
                    // get the derivatives of the SPs
                    Matrix< DDRMat > tNitscheDer      = tSPNitsche->dSPdLeaderDOF( tDofType ).get_row( 0 );
                    Matrix< DDRMat > tLeaderWeightDer = tSPNitsche->dSPdLeaderDOF( tDofType ).get_row( 1 );
                    Matrix< DDRMat > tFollowerWeightDer  = tSPNitsche->dSPdLeaderDOF( tDofType ).get_row( 2 );

                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * (
                                      tFILeader->N_trans() * tJumpViscosity * tNitscheDer
                                    - tFILeader->N_trans() * tCMLeaderSATurbulence->traction( get_normal() ) * tLeaderWeightDer
                                    - mBeta * tCMLeaderSATurbulence->testTraction( get_normal(), mResidualDofType( 0 ) ) * tJumpViscosity * tLeaderWeightDer
                                    - tFILeader->N_trans() * tCMFollowerSATurbulence->traction( get_normal() ) * tFollowerWeightDer );

                    mSet->get_jacobian()(
                            { tFollowerResStartIndex,  tFollowerResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * (
                                    - tFIFollower->N_trans() * tJumpViscosity * tNitscheDer
                                    + tFIFollower->N_trans() * tCMLeaderSATurbulence->traction( get_normal() ) * tLeaderWeightDer
                                    + tFIFollower->N_trans() * tCMFollowerSATurbulence->traction( get_normal() ) * tFollowerWeightDer
                                    - mBeta * tCMFollowerSATurbulence->testTraction( get_normal(), mResidualDofType( 0 ) ) * tJumpViscosity * tFollowerWeightDer );
                }
            }

            // compute the jacobian for indirect dof dependencies through follower constitutive models
            uint tFollowerNumDofDependencies = get_requested_follower_dof_types().size();
            for( uint iDOF = 0; iDOF < tFollowerNumDofDependencies; iDOF++ )
            {
                // get dof type
                Vector< MSI::Dof_Type > tDofType = get_requested_follower_dof_types()( iDOF );

                // get the index for the dof type
                sint tDofDepIndex        = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::FOLLOWER );
                uint tFollowerDepStartIndex = mSet->get_jac_dof_assembly_map()( tFollowerDofIndex )( tDofDepIndex, 0 );
                uint tFollowerDepStopIndex  = mSet->get_jac_dof_assembly_map()( tFollowerDofIndex )( tDofDepIndex, 1 );

                // if dof type is residual dof type
                if( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tFollowerDepStartIndex,  tFollowerDepStopIndex  } ) += aWStar * (
                                    + mBeta * tLeaderWeight * tCMLeaderSATurbulence->testTraction( get_normal(), mResidualDofType( 0 ) ) * tFIFollower->N()
                                    - tNitsche * tFILeader->N_trans() * tFIFollower->N() );

                    mSet->get_jacobian()(
                            { tFollowerResStartIndex, tFollowerResStopIndex },
                            { tFollowerDepStartIndex, tFollowerDepStopIndex } ) += aWStar * (
                                    + mBeta * tFollowerWeight * tCMFollowerSATurbulence->testTraction( get_normal(), mResidualDofType( 0 ) ) * tFIFollower->N()
                                    + tNitsche * tFIFollower->N_trans() * tFIFollower->N() );
                }

                // if dependency of turbulence CM on the dof type
                if ( tCMFollowerSATurbulence->check_dof_dependency( tDofType ) )
                {
                    // add contribution from the derivative of the traction to jacobian
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tFollowerDepStartIndex,  tFollowerDepStopIndex  } ) += aWStar * (
                                    - tFILeader->N_trans() * tFollowerWeight * tCMFollowerSATurbulence->dTractiondDOF( tDofType, get_normal() ) );

                    mSet->get_jacobian()(
                            { tFollowerResStartIndex, tFollowerResStopIndex },
                            { tFollowerDepStartIndex, tFollowerDepStopIndex } ) += aWStar * (
                                    + tFIFollower->N_trans() * tFollowerWeight * tCMFollowerSATurbulence->dTractiondDOF( tDofType, get_normal() ) )
                                    - mBeta * tFollowerWeight * tCMFollowerSATurbulence->dTestTractiondDOF( tDofType, get_normal(), mResidualDofType( 0 ) ) * tJumpViscosity( 0 );
                }

                // if dependency of stabilization parameters on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Leader_Follower::FOLLOWER ) )
                {
                    // get derivatives of SP
                    Matrix< DDRMat > tNitscheDer      = tSPNitsche->dSPdFollowerDOF( tDofType ).get_row( 0 );
                    Matrix< DDRMat > tLeaderWeightDer = tSPNitsche->dSPdFollowerDOF( tDofType ).get_row( 1 );
                    Matrix< DDRMat > tFollowerWeightDer  = tSPNitsche->dSPdFollowerDOF( tDofType ).get_row( 2 );

                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tFollowerDepStartIndex,  tFollowerDepStopIndex  } ) += aWStar * (
                                    + tFILeader->N_trans() * tJumpViscosity * tNitscheDer
                                    - tFILeader->N_trans() * tCMLeaderSATurbulence->traction( get_normal() ) * tLeaderWeightDer
                                    - mBeta * tCMLeaderSATurbulence->testTraction( get_normal(), mResidualDofType( 0 ) ) * tJumpViscosity * tLeaderWeightDer
                                    - tFILeader->N_trans() * tCMFollowerSATurbulence->traction( get_normal() ) * tFollowerWeightDer );

                    mSet->get_jacobian()(
                            { tFollowerResStartIndex, tFollowerResStopIndex },
                            { tFollowerDepStartIndex, tFollowerDepStopIndex } ) += aWStar * (
                                    - tFIFollower->N_trans() * tJumpViscosity * tNitscheDer
                                    + tFIFollower->N_trans() * tCMLeaderSATurbulence->traction( get_normal() ) * tLeaderWeightDer
                                    + tFIFollower->N_trans() * tCMFollowerSATurbulence->traction( get_normal() ) * tFollowerWeightDer
                                    - mBeta * tCMFollowerSATurbulence->testTraction( get_normal(), mResidualDofType( 0 ) ) * tJumpViscosity * tFollowerWeightDer );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Spalart_Allmaras_Turbulence_Interface::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Spalart_Allmaras_Turbulence_Interface::compute_jacobian_and_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Spalart_Allmaras_Turbulence_Interface::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Spalart_Allmaras_Turbulence_Interface::compute_dRdp( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Spalart_Allmaras_Turbulence_Interface::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

