/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Incompressible_NS_Velocity_Interface.cpp
 *
 */

#include "cl_FEM_IWG_Incompressible_NS_Velocity_Interface.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_trans.hpp"
#include "fn_isfinite.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IWG_Incompressible_NS_Velocity_Interface::IWG_Incompressible_NS_Velocity_Interface( sint aBeta )
    {
        // set mBeta for symmetric/unsymmetric Nitsche
        mBeta = aBeta;

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );
        mFollowerCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "IncompressibleFluid" ] = static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE );

        // set size for the stabilization parameter pointer cell
        mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

        // populate the stabilization map
        mStabilizationMap[ "NitscheInterface" ] = static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE );
    }

    //------------------------------------------------------------------------------

    void IWG_Incompressible_NS_Velocity_Interface::compute_residual( real aWStar )
    {
#ifdef MORIS_HAVE_DEBUG
            // check leader and follower field interpolators
            this->check_field_interpolators( mtk::Leader_Follower::LEADER );
            this->check_field_interpolators( mtk::Leader_Follower::FOLLOWER );
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
            Field_Interpolator * tFILeader =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get follower field interpolator for the residual dof type
            Field_Interpolator * tFIFollower  =
                    mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the fluid constitutive model
            const std::shared_ptr< Constitutive_Model > & tCMLeaderFluid =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );
            const std::shared_ptr< Constitutive_Model > & tCMFollowerFluid =
                    mFollowerCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );
            real tNitsche      = tSPNitsche->val()( 0 );
            real tLeaderWeight = tSPNitsche->val()( 1 );
            real tFollowerWeight  = tSPNitsche->val()( 2 );

            // evaluate average traction
            Matrix< DDRMat > tTractionFluid =
                    tLeaderWeight * tCMLeaderFluid->traction( mNormal ) +
                    tFollowerWeight  * tCMFollowerFluid->traction( mNormal );

            // evaluate temperature jump
            Matrix< DDRMat > tJumpVelocity = tFILeader->val() - tFIFollower->val();

            // compute leader residual
            mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex },
                    { 0, 0 } ) += aWStar * (
                            - tFILeader->N_trans() * tTractionFluid
                            - mBeta * tLeaderWeight * trans( tCMLeaderFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tJumpVelocity
                            + tNitsche * tFILeader->N_trans() * tJumpVelocity ) ;

            // compute follower residual
            mSet->get_residual()( 0 )(
                    { tFollowerResStartIndex, tFollowerResStopIndex },
                    { 0, 0 } ) += aWStar * (
                            + tFIFollower->N_trans() * tTractionFluid
                            - mBeta * tFollowerWeight * trans( tCMFollowerFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tJumpVelocity
                            - tNitsche * tFIFollower->N_trans() * tJumpVelocity );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Incompressible_NS_Velocity_Interface::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Velocity_Interface::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader and follower field interpolators
            this->check_field_interpolators( mtk::Leader_Follower::LEADER );
            this->check_field_interpolators( mtk::Leader_Follower::FOLLOWER );
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
            Field_Interpolator * tFILeader =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get follower field interpolator for the residual dof type
            Field_Interpolator * tFIFollower  =
                    mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the fluid constitutive model
            const std::shared_ptr< Constitutive_Model > & tCMLeaderFluid =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );
            const std::shared_ptr< Constitutive_Model > & tCMFollowerFluid =
                    mFollowerCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );
            real tNitsche      = tSPNitsche->val()( 0 );
            real tLeaderWeight = tSPNitsche->val()( 1 );
            real tFollowerWeight  = tSPNitsche->val()( 2 );

            // evaluate average traction
            Matrix< DDRMat > tTractionFluid =
                    tLeaderWeight * tCMLeaderFluid->traction( mNormal ) +
                    tFollowerWeight  * tCMFollowerFluid->traction( mNormal );

            // evaluate temperature jump
            Matrix< DDRMat > tJumpVelocity = tFILeader->val() - tFIFollower->val();

            // get number of leader dof dependencies
            uint tLeaderNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            // compute the jacobian for indirect dof dependencies through leader constitutive models
            for( uint iDOF = 0; iDOF < tLeaderNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Vector< MSI::Dof_Type > & tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

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
                                    - mBeta * tLeaderWeight * trans( tCMLeaderFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tFILeader->N()
                                    + tNitsche * tFILeader->N_trans() * tFILeader->N() );

                    mSet->get_jacobian()(
                            { tFollowerResStartIndex,  tFollowerResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * (
                                    - mBeta * tFollowerWeight * trans( tCMFollowerFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tFILeader->N()
                                    - tNitsche * tFIFollower->N_trans() * tFILeader->N() );
                }

                // if dependency on the dof type
                if ( tCMLeaderFluid->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * (
                                    - tFILeader->N_trans() * tLeaderWeight * tCMLeaderFluid->dTractiondDOF( tDofType, mNormal )
                                    - mBeta * tLeaderWeight * tCMLeaderFluid->dTestTractiondDOF(
                                            tDofType, mNormal, tJumpVelocity, mResidualDofType( 0 ) ) );

                    mSet->get_jacobian()(
                            { tFollowerResStartIndex,  tFollowerResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * (
                                    + tFIFollower->N_trans() * tLeaderWeight * tCMLeaderFluid->dTractiondDOF( tDofType, mNormal ) );
                }

                // if dependency of stabilization parameters on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Leader_Follower::LEADER ) )
                {
                    // get the derivatives of the SPs
                    Matrix< DDRMat > tNitscheDer      = tSPNitsche->dSPdLeaderDOF( tDofType ).get_row( 0 );
                    Matrix< DDRMat > tLeaderWeightDer = tSPNitsche->dSPdLeaderDOF( tDofType ).get_row( 1 );
                    Matrix< DDRMat > tFollowerWeightDer  = tSPNitsche->dSPdLeaderDOF( tDofType ).get_row( 2 );

                    // get traction derivative
                    Matrix< DDRMat > tTractionDer =
                            tCMLeaderFluid->traction( mNormal ) * tLeaderWeightDer +
                            tCMFollowerFluid->traction( mNormal )  * tFollowerWeightDer;

                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * (
                                    - tFILeader->N_trans() * tTractionDer
                                    - mBeta * trans( tCMLeaderFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tJumpVelocity * tLeaderWeightDer
                                    + tFILeader->N_trans() * tJumpVelocity * tNitscheDer );

                    mSet->get_jacobian()(
                            { tFollowerResStartIndex,  tFollowerResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * (
                                    + tFIFollower->N_trans() * tTractionDer
                                    - mBeta * trans( tCMFollowerFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tJumpVelocity * tFollowerWeightDer
                                    - tFIFollower->N_trans() * tJumpVelocity * tNitscheDer );
                }
            }

            // compute the jacobian for indirect dof dependencies through follower constitutive models
            uint tFollowerNumDofDependencies = mRequestedFollowerGlobalDofTypes.size();

            for( uint iDOF = 0; iDOF < tFollowerNumDofDependencies; iDOF++ )
            {
                // get dof type
                Vector< MSI::Dof_Type > tDofType = mRequestedFollowerGlobalDofTypes( iDOF );

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
                                    + mBeta * tLeaderWeight * trans( tCMLeaderFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tFIFollower->N()
                                    - tNitsche * tFILeader->N_trans() * tFIFollower->N() );

                    mSet->get_jacobian()(
                            { tFollowerResStartIndex, tFollowerResStopIndex },
                            { tFollowerDepStartIndex, tFollowerDepStopIndex } ) += aWStar * (
                                    + mBeta * tFollowerWeight * trans( tCMFollowerFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tFIFollower->N()
                                    + tNitsche * tFIFollower->N_trans() * tFIFollower->N() );
                }

                // if dependency on the dof type
                if ( tCMFollowerFluid->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tFollowerDepStartIndex,  tFollowerDepStopIndex  } ) += aWStar * (
                                    - tFILeader->N_trans() * tFollowerWeight * tCMFollowerFluid->dTractiondDOF( tDofType, mNormal ) );

                    mSet->get_jacobian()(
                            { tFollowerResStartIndex, tFollowerResStopIndex },
                            { tFollowerDepStartIndex, tFollowerDepStopIndex } ) += aWStar * (
                                    + tFIFollower->N_trans() * tFollowerWeight * tCMFollowerFluid->dTractiondDOF( tDofType, mNormal )
                                    - mBeta * tFollowerWeight * tCMFollowerFluid->dTestTractiondDOF(
                                            tDofType, mNormal, tJumpVelocity, mResidualDofType( 0 ) ) );
                }

                // if dependency of stabilization parameters on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Leader_Follower::FOLLOWER ) )
                {
                    // get the derivatives of the SPs
                    Matrix< DDRMat > tNitscheDer      = tSPNitsche->dSPdFollowerDOF( tDofType ).get_row( 0 );
                    Matrix< DDRMat > tLeaderWeightDer = tSPNitsche->dSPdFollowerDOF( tDofType ).get_row( 1 );
                    Matrix< DDRMat > tFollowerWeightDer  = tSPNitsche->dSPdFollowerDOF( tDofType ).get_row( 2 );

                    // get traction derivative
                    Matrix< DDRMat > tTractionDer =
                            tCMLeaderFluid->traction( mNormal ) * tLeaderWeightDer +
                            tCMFollowerFluid->traction( mNormal )  * tFollowerWeightDer;

                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tFollowerDepStartIndex,  tFollowerDepStopIndex  } ) += aWStar * (
                                    - tFILeader->N_trans() * tTractionDer
                                    - mBeta * trans( tCMLeaderFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tJumpVelocity * tLeaderWeightDer
                                    + tFILeader->N_trans() * tJumpVelocity * tNitscheDer );

                    mSet->get_jacobian()(
                            { tFollowerResStartIndex, tFollowerResStopIndex },
                            { tFollowerDepStartIndex, tFollowerDepStopIndex } ) += aWStar * (
                                    + tFIFollower->N_trans() * tTractionDer
                                    - mBeta * trans( tCMFollowerFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tJumpVelocity * tFollowerWeightDer
                                    - tFIFollower->N_trans() * tJumpVelocity * tNitscheDer );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Incompressible_NS_Velocity_Interface::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Velocity_Interface::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_Interface::compute_jacobian_and_residual - This function does nothing.");
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Velocity_Interface::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_Interface::compute_dRdp - This function does nothing.");
        }

        //------------------------------------------------------------------------------
}    // namespace moris::fem
