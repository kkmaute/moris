/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Nonlocal_Interface.cpp
 *
 */

#include "cl_FEM_IWG_Nonlocal_Interface.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IWG_Nonlocal_Interface::IWG_Nonlocal_Interface( sint aBeta )
        {
            // set sint for symmetric/unsymmetric Nitsche
            mBeta = aBeta;

            // set size for the property pointer cell
            mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );
            mFollowerProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Weight" ] = static_cast< uint >( IWG_Property_Type::WEIGHT );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "NitscheInterface" ] = static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Nonlocal_Interface::set_parameters( const Vector< Matrix< DDRMat > >& aParameters )
        {
            // set parameters
            mParameters = aParameters;

            // check characteristic length provided
            MORIS_ERROR( mParameters.size() >= 1 && mParameters.size() < 3,
                    "IWG_Nonlocal_Interface::set_parameters - requires a characteristic length.\n" );

            // set a characteristic length
            mCharacteristicLength = aParameters( 0 )( 0 );

            // FIXME add higher order
            //            // set a order if provided by default 1
            //            if ( mParameters.size() > 1 )
            //            {
            //                mOrder = static_cast< uint >( aParameters( 1 )( 0 ) );
            //            }

            //            // FIXME for now only 1st order model supported
            //            MORIS_ERROR( mOrder < 2.0,
            //                    "IWG_Nonlocal_Interface::set_parameters - higher order model not supported.\n" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Nonlocal_Interface::compute_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check Leader and follower field interpolators
            this->check_field_interpolators( mtk::Leader_Follower::LEADER );
            this->check_field_interpolators( mtk::Leader_Follower::FOLLOWER );
#endif

            // get Leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get follower index for residual dof type, indices for assembly
            uint tFollowerDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );
            uint tFollowerResStartIndex = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 0 );
            uint tFollowerResStopIndex  = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 1 );

            // get the Leader field interpolator for the residual dof type
            Field_Interpolator* tFILeader =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the follower field interpolator for the residual dof type
            Field_Interpolator* tFIFollower =
                    mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get weight coefficient
            const std::shared_ptr< Property >& tPropWeight =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::WEIGHT ) );

            // if weight property is provided
            real tWeight = 1.0;
            if( tPropWeight )
            {
                tWeight = tPropWeight->val()( 0 );
            }

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );

            real tNitsche      = tSPNitsche->val()( 0 );
            real tLeaderWeight = tSPNitsche->val()( 1 );
            real tFollowerWeight = tSPNitsche->val()( 2 );

            // evaluate average traction
            Matrix< DDRMat > tTraction =
                    tLeaderWeight * std::pow( mCharacteristicLength, 2.0 * mOrder ) * trans( tFILeader->gradx( 1 ) ) * mNormal / mOrderCoeff( mOrder )
                    + tFollowerWeight * std::pow( mCharacteristicLength, 2.0 * mOrder ) * trans( tFIFollower->gradx( 1 ) ) * mNormal / mOrderCoeff( mOrder );

            // evaluate jump
            Matrix< DDRMat > tJump = tFILeader->val() - tFIFollower->val();

            // compute leader residual
            mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex } ) +=
                    aWStar * tWeight * ( -tFILeader->N_trans() * tTraction                                                                                                                        //
                                         + mBeta * tLeaderWeight * std::pow( mCharacteristicLength, 2.0 * mOrder ) * trans( tFILeader->dnNdxn( 1 ) ) * mNormal * tJump / mOrderCoeff( mOrder )    //
                                         + tNitsche * tFILeader->N_trans() * tJump );                                                                                                             //

            // compute follower residual
            mSet->get_residual()( 0 )(
                    { tFollowerResStartIndex, tFollowerResStopIndex } ) +=
                    aWStar * tWeight * ( +tFIFollower->N_trans() * tTraction                                                                                                                          //
                                         + mBeta * tFollowerWeight * std::pow( mCharacteristicLength, 2.0 * mOrder ) * trans( tFIFollower->dnNdxn( 1 ) ) * mNormal * tJump / mOrderCoeff( mOrder )    //
                                         - tNitsche * tFIFollower->N_trans() * tJump );                                                                                                               //

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Nonlocal_Interface::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Nonlocal_Interface::compute_jacobian( real aWStar )
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

            // get weight coefficient property
            const std::shared_ptr< Property >& tPropWeight =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::WEIGHT ) );

            // if weight property is provided
            real tWeight = 1.0;
            if( tPropWeight )
            {
                tWeight = tPropWeight->val()( 0 );
            }

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );

            real tNitsche      = tSPNitsche->val()( 0 );
            real tLeaderWeight = tSPNitsche->val()( 1 );
            real tFollowerWeight = tSPNitsche->val()( 2 );

            // evaluate temperature jump
            Matrix< DDRMat > tJump = tFILeader->val() - tFIFollower->val();

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
                    // compute derivative of traction
                    Matrix< DDRMat > tTractionDer =    //
                            tLeaderWeight * std::pow( mCharacteristicLength, 2.0 * mOrder ) * trans( mNormal ) * tFILeader->dnNdxn( 1 ) / mOrderCoeff( mOrder );

                    // add contribution to leader
                    tJacLeader +=
                            aWStar * tWeight * ( -tFILeader->N_trans() * tTractionDer                                                                                                                              //
                                                 + mBeta * tLeaderWeight * std::pow( mCharacteristicLength, 2.0 * mOrder ) * trans( tFILeader->dnNdxn( 1 ) ) * mNormal * tFILeader->N() / mOrderCoeff( mOrder )    //
                                                 + tNitsche * tFILeader->N_trans() * tFILeader->N() );

                    // add contribution to follower
                    tJacFollower +=
                            aWStar * tWeight * ( +tFIFollower->N_trans() * tTractionDer                                                                                                                                //
                                                 + mBeta * tFollowerWeight * std::pow( mCharacteristicLength, 2.0 * mOrder ) * trans( tFIFollower->dnNdxn( 1 ) ) * mNormal * tFILeader->N() / mOrderCoeff( mOrder )    //
                                                 - tNitsche * tFIFollower->N_trans() * tFILeader->N() );
                }

                // if dependency of constitutive models on the dof type
                if( tPropWeight )
                {
                    if ( tPropWeight->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to Jacobian
                        MORIS_ASSERT( false,
                                "IWG_Nonlocal_Interface::compute_jacobian - tPropWeight dependency on dof not accounted for." );
                    }
                }

                // if dependency of stabilization parameters on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Leader_Follower::LEADER ) )
                {
                    // get the derivatives of the SPs
                    Matrix< DDRMat > tNitscheDer        = tSPNitsche->dSPdLeaderDOF( tDofType ).get_row( 0 );
                    Matrix< DDRMat > tLeaderWeightDer   = tSPNitsche->dSPdLeaderDOF( tDofType ).get_row( 1 );
                    Matrix< DDRMat > tFollowerWeightDer = tSPNitsche->dSPdLeaderDOF( tDofType ).get_row( 2 );

                    // get traction derivative
                    Matrix< DDRMat > tTractionDer =
                            std::pow( mCharacteristicLength, 2.0 * mOrder ) * trans( mNormal ) * tFILeader->gradx( 1 ) * tLeaderWeightDer / mOrderCoeff( mOrder )    //
                            + std::pow( mCharacteristicLength, 2.0 * mOrder ) * trans( mNormal ) * tFIFollower->gradx( 1 ) * tFollowerWeightDer / mOrderCoeff( mOrder );

                    // add contribution to Jacobian
                    tJacLeader += aWStar * tWeight * ( -tFILeader->N_trans() * tTractionDer                                                                                                                        //
                                                       + mBeta * std::pow( mCharacteristicLength, 2.0 * mOrder ) * trans( tFILeader->dnNdxn( 1 ) ) * mNormal * tJump * tLeaderWeightDer / mOrderCoeff( mOrder )    //
                                                       + tFILeader->N_trans() * tJump * tNitscheDer );

                    tJacFollower += aWStar * tWeight * ( +tFIFollower->N_trans() * tTractionDer                                                                                                                          //
                                                         + mBeta * std::pow( mCharacteristicLength, 2.0 * mOrder ) * trans( tFIFollower->dnNdxn( 1 ) ) * mNormal * tJump * tFollowerWeightDer / mOrderCoeff( mOrder )    //
                                                         - tFIFollower->N_trans() * tJump * tNitscheDer );
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
                sint tDofDepIndex           = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::FOLLOWER );
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

                    // compute derivative of traction
                    Matrix< DDRMat > tTractionDer =    //
                            tFollowerWeight * std::pow( mCharacteristicLength, 2.0 * mOrder ) * trans( mNormal ) * tFIFollower->dnNdxn( 1 ) / mOrderCoeff( mOrder );

                    // add contribution to leader
                    tJacLeader += aWStar * tWeight * ( -tFILeader->N_trans() * tTractionDer                                                                                                                                //
                                                       - mBeta * tLeaderWeight * std::pow( mCharacteristicLength, 2.0 * mOrder ) * trans( tFILeader->dnNdxn( 1 ) ) * mNormal * tFIFollower->N() / mOrderCoeff( mOrder )    //
                                                       - tNitsche * tFILeader->N_trans() * tFIFollower->N() );

                    tJacFollower += aWStar * tWeight * ( +tFIFollower->N_trans() * tTractionDer                                                                                                                                  //
                                                         - mBeta * tFollowerWeight * std::pow( mCharacteristicLength, 2.0 * mOrder ) * trans( tFIFollower->dnNdxn( 1 ) ) * mNormal * tFIFollower->N() / mOrderCoeff( mOrder )    //
                                                         + tNitsche * tFIFollower->N_trans() * tFIFollower->N() );
                }

                // if dependency on the dof type
                if( tPropWeight )
                {
                    if ( tPropWeight->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to Jacobian
                        MORIS_ASSERT( false,
                                "IWG_Nonlocal_Interface::compute_jacobian - tPropWeight dependency on dof not accounted for." );
                    }
                }

                // if dependency of stabilization parameters on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Leader_Follower::FOLLOWER ) )
                {
                    // get the derivatives of the SPs
                    Matrix< DDRMat > tNitscheDer        = tSPNitsche->dSPdFollowerDOF( tDofType ).get_row( 0 );
                    Matrix< DDRMat > tLeaderWeightDer   = tSPNitsche->dSPdFollowerDOF( tDofType ).get_row( 1 );
                    Matrix< DDRMat > tFollowerWeightDer = tSPNitsche->dSPdFollowerDOF( tDofType ).get_row( 2 );

                    // get traction derivative
                    Matrix< DDRMat > tTractionDer =
                            std::pow( mCharacteristicLength, 2.0 * mOrder ) * trans( mNormal ) * tFILeader->gradx( 1 ) * tLeaderWeightDer / mOrderCoeff( mOrder )    //
                            + std::pow( mCharacteristicLength, 2.0 * mOrder ) * trans( mNormal ) * tFIFollower->gradx( 1 ) * tFollowerWeightDer / mOrderCoeff( mOrder );

                    // add contribution to Jacobian
                    tJacLeader += aWStar * tWeight * ( -tFILeader->N_trans() * tTractionDer                                                                                                                        //
                                                       + mBeta * std::pow( mCharacteristicLength, 2.0 * mOrder ) * trans( tFILeader->dnNdxn( 1 ) ) * mNormal * tJump * tLeaderWeightDer / mOrderCoeff( mOrder )    //
                                                       + tFILeader->N_trans() * tJump * tNitscheDer );

                    tJacFollower += aWStar * tWeight * ( +tFIFollower->N_trans() * tTractionDer                                                                                                                          //
                                                         + mBeta * std::pow( mCharacteristicLength, 2.0 * mOrder ) * trans( tFIFollower->dnNdxn( 1 ) ) * mNormal * tJump * tFollowerWeightDer / mOrderCoeff( mOrder )    //
                                                         - tFIFollower->N_trans() * tJump * tNitscheDer );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Nonlocal_Interface::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Nonlocal_Interface::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Nonlocal_Interface::compute_jacobian_and_residual - This function does nothing." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Nonlocal_Interface::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Nonlocal_Interface::compute_dRdp - This function does nothing." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

