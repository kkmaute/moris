/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Incompressible_NS_Pressure_Interface.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Incompressible_NS_Pressure_Interface.hpp"
//LINALG/src
#include "fn_trans.hpp"

namespace moris::fem
{

    //------------------------------------------------------------------------------

    IWG_Incompressible_NS_Pressure_Interface::IWG_Incompressible_NS_Pressure_Interface( sint aBeta )
    {
        // set mBeta for symmetric/skew symmetric Nitsche
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

    void IWG_Incompressible_NS_Pressure_Interface::compute_residual( real aWStar )
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

            // get the velocity dof type FI
            // FIXME protect dof type
            Field_Interpolator * tFILeader =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );
            Field_Interpolator * tFIFollower =
                    mFollowerFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get leader/follower fluid constitutive model
            const std::shared_ptr< Constitutive_Model > & tCMLeaderFluid =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );
            const std::shared_ptr< Constitutive_Model > & tCMFollowerFluid =
                    mFollowerCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );
            real tLeaderWeight = tSPNitsche->val()( 1 );
            real tFollowerWeight  = tSPNitsche->val()( 2 );

            // compute the jump
            Matrix< DDRMat > tVelocityJump = tFILeader->val() - tFIFollower->val();

            // compute leader residual
            mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex },
                    { 0, 0 } ) -= aWStar * (
                            mBeta * tLeaderWeight * trans( tCMLeaderFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tVelocityJump );

            // compute follower residual
            mSet->get_residual()( 0 )(
                    { tFollowerResStartIndex, tFollowerResStopIndex },
                    { 0, 0 } ) -= aWStar * (
                            mBeta * tFollowerWeight * trans( tCMFollowerFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tVelocityJump );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Incompressible_NS_Pressure_Interface::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Pressure_Interface::compute_jacobian( real aWStar )
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

            // get the velocity dof type FI
            // FIXME protect dof type
            Field_Interpolator * tFILeader =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );
            Field_Interpolator * tFIFollower =
                    mFollowerFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get leader/follower fluid constitutive model
            const std::shared_ptr< Constitutive_Model > & tCMLeaderFluid =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );
            const std::shared_ptr< Constitutive_Model > & tCMFollowerFluid =
                    mFollowerCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );
            real tLeaderWeight = tSPNitsche->val()( 1 );
            real tFollowerWeight  = tSPNitsche->val()( 2 );

            // compute the jump
            Matrix< DDRMat > tVelocityJump = tFILeader->val() - tFIFollower->val();

            // get number of leader dependencies
            uint tLeaderNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            // compute the jacobian for indirect dof dependencies through leader
            for( uint iDOF = 0; iDOF < tLeaderNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Vector< MSI::Dof_Type > & tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // if dof type is velocity
                // FIXME protect dof type
                if ( tDofType( 0 ) == MSI::Dof_Type::VX )
                {
                    // add contribution to leader jacobian
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -= aWStar * (
                                    mBeta * tLeaderWeight * trans( tCMLeaderFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tFILeader->N() );

                    // add contribution to follower jacobian
                    mSet->get_jacobian()(
                            { tFollowerResStartIndex, tFollowerResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -= aWStar * (
                                    mBeta * tFollowerWeight * trans( tCMFollowerFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tFILeader->N() );
                }

                // if fluid constitutive model depends on dof type
                if ( tCMLeaderFluid->check_dof_dependency( tDofType ) )
                {
                    // add contribution to leader jacobian
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -= aWStar * (
                                    mBeta * tLeaderWeight * tCMLeaderFluid->dTestTractiondDOF(
                                            tDofType, mNormal, tVelocityJump, mResidualDofType( 0 ) ) );
                }

                // if leader SP depends on dof type
                if( tSPNitsche->check_dof_dependency( tDofType, mtk::Leader_Follower::LEADER ) )
                {
                    // get derivatives of SP
                    Matrix< DDRMat > tLeaderWeightDer = tSPNitsche->dSPdFollowerDOF( tDofType ).get_row( 1 );
                    Matrix< DDRMat > tFollowerWeightDer  = tSPNitsche->dSPdFollowerDOF( tDofType ).get_row( 2 );

                    // add contribution to leader jacobian
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -= aWStar * (
                                    mBeta * trans( tCMLeaderFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tVelocityJump * tLeaderWeightDer
                                    + mBeta * trans( tCMFollowerFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tVelocityJump * tFollowerWeightDer );
                }
            }

            // get number of follower dependencies
            uint tFollowerNumDofDependencies = mRequestedFollowerGlobalDofTypes.size();

            // compute the jacobian for indirect dof dependencies through leader
            for( uint iDOF = 0; iDOF < tFollowerNumDofDependencies; iDOF++ )
            {
                // get dof type
                Vector< MSI::Dof_Type > tDofType = mRequestedFollowerGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex        = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::FOLLOWER );
                uint tFollowerDepStartIndex = mSet->get_jac_dof_assembly_map()( tFollowerDofIndex )( tDofDepIndex, 0 );
                uint tFollowerDepStopIndex  = mSet->get_jac_dof_assembly_map()( tFollowerDofIndex )( tDofDepIndex, 1 );

                // if dof type is velocity
                // FIXME protect dof type
                if ( tDofType( 0 ) == MSI::Dof_Type::VX )
                {
                    // add contribution to leader jacobian
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tFollowerDepStartIndex, tFollowerDepStopIndex } ) += aWStar * (
                                    mBeta * tLeaderWeight * trans( tCMLeaderFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tFIFollower->N() );

                    // add contribution to follower jacobian
                    mSet->get_jacobian()(
                            { tFollowerResStartIndex, tFollowerResStopIndex },
                            { tFollowerDepStartIndex, tFollowerDepStopIndex } ) += aWStar * (
                                    mBeta * tFollowerWeight * trans( tCMFollowerFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tFIFollower->N() );
                }

                // if fluid constitutive model depends on dof type
                if ( tCMFollowerFluid->check_dof_dependency( tDofType ) )
                {
                    // add contribution to follower jacobian
                    mSet->get_jacobian()(
                            { tFollowerResStartIndex, tFollowerResStopIndex },
                            { tFollowerDepStartIndex, tFollowerDepStopIndex } ) -= aWStar * (
                                    mBeta * tFollowerWeight * tCMFollowerFluid->dTestTractiondDOF(
                                            tDofType, mNormal, tVelocityJump, mResidualDofType( 0 ) ) );
                }

                // if leader SP depends on dof type
                if( tSPNitsche->check_dof_dependency( tDofType, mtk::Leader_Follower::FOLLOWER ) )
                {
                    // get derivatives of SP
                    Matrix< DDRMat > tLeaderWeightDer = tSPNitsche->dSPdFollowerDOF( tDofType ).get_row( 1 );
                    Matrix< DDRMat > tFollowerWeightDer  = tSPNitsche->dSPdFollowerDOF( tDofType ).get_row( 2 );

                    // add contribution to leader jacobian
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tFollowerDepStartIndex, tFollowerDepStopIndex } ) -= aWStar * (
                                    mBeta * trans( tCMLeaderFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tVelocityJump * tLeaderWeightDer
                                    + mBeta * trans( tCMFollowerFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tVelocityJump * tFollowerWeightDer );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Incompressible_NS_Pressure_Interface::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Pressure_Interface::compute_jacobian_and_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Pressure_Interface::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Pressure_Interface::compute_dRdp( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Pressure_Interface::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------
}    // namespace moris::fem
