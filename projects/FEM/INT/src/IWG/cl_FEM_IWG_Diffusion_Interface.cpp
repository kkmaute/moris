/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Diffusion_Interface.cpp
 *
 */

#include "cl_FEM_IWG_Diffusion_Interface.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IWG_Diffusion_Interface::IWG_Diffusion_Interface( sint aBeta )
    {
        // set sint for symmetric/unsymmetric Nitsche
        mBeta = aBeta;

        // set size for the constitutive model pointer cell
        mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );
        mFollowerProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

        mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );
        mFollowerCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "Thickness" ] = static_cast< uint >( IWG_Property_Type::THICKNESS );
        mPropertyMap[ "Select" ]    = static_cast< uint >( IWG_Property_Type::SELECT );

        // populate the constitutive map
        mConstitutiveMap[ "Diffusion" ] = static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO );

        // set size for the stabilization parameter pointer cell
        mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

        // populate the stabilization map
        mStabilizationMap[ "NitscheInterface" ] = static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE );
    }

    //------------------------------------------------------------------------------

    void
    IWG_Diffusion_Interface::compute_residual( real aWStar )
    {
#ifdef MORIS_HAVE_DEBUG
            // check leader and follower field interpolators
            this->check_field_interpolators( mtk::Leader_Follower::LEADER );
            this->check_field_interpolators( mtk::Leader_Follower::FOLLOWER );
#endif

            // Check that leader and follower CMs are different
            MORIS_ASSERT( &mLeaderCM( 0 ) != &mFollowerCM( 0 ), "Leader and Follower constitutive model are the same. This will cause problems " );

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get follower index for residual dof type, indices for assembly
            uint tFollowerDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );
            uint tFollowerResStartIndex = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 0 );
            uint tFollowerResStopIndex  = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 1 );

            // get the leader field interpolator for the residual dof type
            Field_Interpolator* tFILeader =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the follower field interpolator for the residual dof type
            Field_Interpolator* tFIFollower =
                    mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the elasticity constitutive model
            const std::shared_ptr< Constitutive_Model >& tCMLeaderDiffusion =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) );

            const std::shared_ptr< Constitutive_Model >& tCMFollowerDiffusion =
                    mFollowerCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );

            real tNitsche      = tSPNitsche->val()( 0 );
            real tLeaderWeight = tSPNitsche->val()( 1 );
            real tFollowerWeight  = tSPNitsche->val()( 2 );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // get the selection matrix property
            const std::shared_ptr< Property >& tPropSelectLeader =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            const std::shared_ptr< Property >& tPropSelectFollower =
                    mFollowerProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            // get selection property from leader and follower
            real tMLeader = 1.0;
            if ( tPropSelectLeader != nullptr )
            {
                tMLeader = tPropSelectLeader->val()( 0 );
            }

            real tMFollower = 1.0;
            if ( tPropSelectFollower != nullptr )
            {
                tMFollower = tPropSelectFollower->val()( 0 );
            }

            // compute effective selection property
            real tM = tMLeader * tMFollower;

            // skip computing residual if selection property is zero
            if ( std::abs( tM ) < MORIS_REAL_EPS ) return;

            // evaluate average traction
            real tTraction =
                    tLeaderWeight * tCMLeaderDiffusion->traction( mNormal )( 0 ) +    //
                    tFollowerWeight * tCMFollowerDiffusion->traction( mNormal )( 0 );

            // evaluate temperature jump
            real tJump = tFILeader->val()( 0 ) - tFIFollower->val()( 0 );

            // compute leader residual
            mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex } ) +=
                    aWStar * (                                                                                                           //
                            -tFILeader->N_trans() * tM * tTraction                                                                       //
                            + mBeta * tLeaderWeight * tCMLeaderDiffusion->testTraction( mNormal, mResidualDofType( 0 ) ) * tM * tJump    //
                            + tNitsche * tFILeader->N_trans() * tM * tJump );

            // compute follower residual
            mSet->get_residual()( 0 )(
                    { tFollowerResStartIndex, tFollowerResStopIndex } ) +=
                    aWStar * (                                                                                                         //
                            +tFIFollower->N_trans() * tM * tTraction                                                                      //
                            + mBeta * tFollowerWeight * tCMFollowerDiffusion->testTraction( mNormal, mResidualDofType( 0 ) ) * tM * tJump    //
                            - tNitsche * tFIFollower->N_trans() * tM * tJump );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Diffusion_Interface::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Interface::compute_jacobian( real aWStar )
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
            Field_Interpolator* tFILeader =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get follower field interpolator for the residual dof type
            Field_Interpolator* tFIFollower =
                    mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the elasticity constitutive model
            const std::shared_ptr< Constitutive_Model >& tCMLeaderDiffusion =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) );
            const std::shared_ptr< Constitutive_Model >& tCMFollowerDiffusion =
                    mFollowerCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );

            real tNitsche      = tSPNitsche->val()( 0 );
            real tLeaderWeight = tSPNitsche->val()( 1 );
            real tFollowerWeight  = tSPNitsche->val()( 2 );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // get the selection matrix property
            const std::shared_ptr< Property >& tPropSelect =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            // set a default selection matrix if needed
            real tM = 1.0;
            if ( tPropSelect != nullptr )
            {
                tM = tPropSelect->val()( 0 );

                // skip computing residual if projection matrix is zero
                if ( std::abs( tM ) < MORIS_REAL_EPS ) return;
            }

            // evaluate temperature jump
            real tJump = tFILeader->val()( 0 ) - tFIFollower->val()( 0 );

            // get number of leader dof dependencies
            uint tLeaderNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            // loop over leader dof dependencies
            for ( uint iDOF = 0; iDOF < tLeaderNumDofDependencies; iDOF++ )
            {
                // get the dof type
                const Vector< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // get sub-matrices
                auto tJacLeader = mSet->get_jacobian()(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } );

                auto tJacFollower = mSet->get_jacobian()(
                        { tFollowerResStartIndex, tFollowerResStopIndex },
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } );

                // compute Jacobian direct dependencies
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    tJacLeader +=                                                                                                                        //
                            aWStar * (                                                                                                                   //
                                    +mBeta * tLeaderWeight * tCMLeaderDiffusion->testTraction( mNormal, mResidualDofType( 0 ) ) * tM * tFILeader->N()    //
                                    + tNitsche * tFILeader->N_trans() * tM * tFILeader->N() );

                    tJacFollower +=                                                                                                                       //
                            aWStar * (                                                                                                                 //
                                    +mBeta * tFollowerWeight * tCMFollowerDiffusion->testTraction( mNormal, mResidualDofType( 0 ) ) * tM * tFILeader->N()    //
                                    - tNitsche * tFIFollower->N_trans() * tM * tFILeader->N() );
                }

                // if dependency of constitutive models on the dof type
                if ( tCMLeaderDiffusion->check_dof_dependency( tDofType ) )
                {
                    // add contribution to Jacobian
                    tJacLeader +=                                                                                                          //
                            aWStar * (                                                                                                     //
                                    -tFILeader->N_trans() * tLeaderWeight * tM * tCMLeaderDiffusion->dTractiondDOF( tDofType, mNormal )    //
                                    + mBeta * tLeaderWeight * tCMLeaderDiffusion->dTestTractiondDOF( tDofType, mNormal, mResidualDofType( 0 ) ) * tM * tJump );

                    tJacFollower +=          //
                            aWStar * (    //
                                    +tFIFollower->N_trans() * tLeaderWeight * tM * tCMLeaderDiffusion->dTractiondDOF( tDofType, mNormal ) );
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
                            tCMLeaderDiffusion->traction( mNormal )( 0 ) * tLeaderWeightDer    //
                            + tCMFollowerDiffusion->traction( mNormal )( 0 ) * tFollowerWeightDer;

                    // add contribution to Jacobian
                    tJacLeader +=                                                                                                                   //
                            aWStar * (                                                                                                              //
                                    -tFILeader->N_trans() * tTractionDer                                                                            //
                                    + mBeta * tCMLeaderDiffusion->testTraction( mNormal, mResidualDofType( 0 ) ) * tM * tJump * tLeaderWeightDer    //
                                    + tFILeader->N_trans() * tM * tJump * tNitscheDer );

                    tJacFollower +=                                                                                                                  //
                            aWStar * (                                                                                                            //
                                    +tFIFollower->N_trans() * tTractionDer                                                                           //
                                    + mBeta * tCMFollowerDiffusion->testTraction( mNormal, mResidualDofType( 0 ) ) * tM * tJump * tFollowerWeightDer    //
                                    - tFIFollower->N_trans() * tM * tJump * tNitscheDer );
                }
            }

            // get number of follower dof dependencies
            uint tFollowerNumDofDependencies = mRequestedFollowerGlobalDofTypes.size();

            // loop over follower dof dependencies
            for ( uint iDOF = 0; iDOF < tFollowerNumDofDependencies; iDOF++ )
            {
                // get dof type
                const Vector< MSI::Dof_Type > tDofType = mRequestedFollowerGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex        = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::FOLLOWER );
                uint tFollowerDepStartIndex = mSet->get_jac_dof_assembly_map()( tFollowerDofIndex )( tDofDepIndex, 0 );
                uint tFollowerDepStopIndex  = mSet->get_jac_dof_assembly_map()( tFollowerDofIndex )( tDofDepIndex, 1 );

                // get sub-matrices
                auto tJacLeader = mSet->get_jacobian()(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { tFollowerDepStartIndex, tFollowerDepStopIndex } );

                auto tJacFollower = mSet->get_jacobian()(
                        { tFollowerResStartIndex, tFollowerResStopIndex },
                        { tFollowerDepStartIndex, tFollowerDepStopIndex } );

                // compute Jacobian direct dependencies
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    tJacLeader +=                                                                                                                       //
                            aWStar * (                                                                                                                  //
                                    -mBeta * tLeaderWeight * tCMLeaderDiffusion->testTraction( mNormal, mResidualDofType( 0 ) ) * tM * tFIFollower->N()    //
                                    - tNitsche * tFILeader->N_trans() * tM * tFIFollower->N() );

                    tJacFollower +=                                                                                                                      //
                            aWStar * (                                                                                                                //
                                    -mBeta * tFollowerWeight * tCMFollowerDiffusion->testTraction( mNormal, mResidualDofType( 0 ) ) * tM * tFIFollower->N()    //
                                    + tNitsche * tFIFollower->N_trans() * tM * tFIFollower->N() );
                }

                // if dependency on the dof type
                if ( tCMFollowerDiffusion->check_dof_dependency( tDofType ) )
                {
                    // add contribution to Jacobian
                    tJacLeader +=         //
                            aWStar * (    //
                                    -tFILeader->N_trans() * tFollowerWeight * tM * tCMFollowerDiffusion->dTractiondDOF( tDofType, mNormal ) );

                    tJacFollower +=                                                                                                        //
                            aWStar * (                                                                                                  //
                                    +tFIFollower->N_trans() * tFollowerWeight * tM * tCMFollowerDiffusion->dTractiondDOF( tDofType, mNormal )    //
                                    + mBeta * tFollowerWeight * tCMFollowerDiffusion->dTestTractiondDOF( tDofType, mNormal, mResidualDofType( 0 ) ) * tM * tJump );
                }

                // if dependency of stabilization parameters on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Leader_Follower::FOLLOWER ) )
                {
                    // get the derivatives of the SPs
                    Matrix< DDRMat > tNitscheDer      = tSPNitsche->dSPdFollowerDOF( tDofType ).get_row( 0 );
                    Matrix< DDRMat > tLeaderWeightDer = tSPNitsche->dSPdFollowerDOF( tDofType ).get_row( 1 );
                    Matrix< DDRMat > tFollowerWeightDer  = tSPNitsche->dSPdFollowerDOF( tDofType ).get_row( 2 );

                    // get traction derivative
                    Matrix< DDRMat > tTractionDer =                                            //
                            tCMLeaderDiffusion->traction( mNormal )( 0 ) * tLeaderWeightDer    //
                            + tCMFollowerDiffusion->traction( mNormal )( 0 ) * tFollowerWeightDer;

                    // add contribution to Jacobian
                    tJacLeader +=                                                                                                                   //
                            aWStar * (                                                                                                              //
                                    -tFILeader->N_trans() * tTractionDer                                                                            //
                                    + mBeta * tCMLeaderDiffusion->testTraction( mNormal, mResidualDofType( 0 ) ) * tM * tJump * tLeaderWeightDer    //
                                    + tFILeader->N_trans() * tM * tJump * tNitscheDer );

                    tJacFollower +=                                                                                                                  //
                            aWStar * (                                                                                                            //
                                    +tFIFollower->N_trans() * tTractionDer                                                                           //
                                    + mBeta * tCMFollowerDiffusion->testTraction( mNormal, mResidualDofType( 0 ) ) * tM * tJump * tFollowerWeightDer    //
                                    - tFIFollower->N_trans() * tM * tJump * tNitscheDer );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Diffusion_Interface::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Interface::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Diffusion_Interface::compute_jacobian_and_residual - This function does nothing." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Interface::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Diffusion_Interface::compute_dRdp - This function does nothing." );
        }

        //------------------------------------------------------------------------------
}    // namespace moris::fem
