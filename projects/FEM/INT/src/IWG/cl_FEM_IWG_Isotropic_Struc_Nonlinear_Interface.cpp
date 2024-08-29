/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Nonlinear_Interface.cpp
 *
 */

#include "cl_FEM_IWG_Isotropic_Struc_Nonlinear_Interface.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IWG_Isotropic_Struc_Nonlinear_Interface::IWG_Isotropic_Struc_Nonlinear_Interface(
            enum CM_Function_Type aStressType,
            enum CM_Function_Type aStrainType,
            sint                  aBeta )
    {
        // set sint for symmetric/unsymmetric Nitsche
        mBeta = aBeta;

        // assign stress and strain type to evaluate the IWG
        mStressType = aStressType;
        mStrainType = aStrainType;

        // set size for the property pointer cell
        mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "Thickness" ] = static_cast< uint >( IWG_Property_Type::THICKNESS );

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );
        mFollowerCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "ElastLinIso" ] = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );

        // set size for the stabilization parameter pointer cell
        mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

        // populate the stabilization map
        mStabilizationMap[ "NitscheInterface" ] = static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE );
    }

    //------------------------------------------------------------------------------

    void
    IWG_Isotropic_Struc_Nonlinear_Interface::compute_residual( real aWStar )
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

            // get leader field interpolator for the residual dof type
            Field_Interpolator* tFILeader =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get follower field interpolator for the residual dof type
            Field_Interpolator* tFIFollower =
                    mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the elasticity constitutive model
            const std::shared_ptr< Constitutive_Model >& tCMLeaderElasticity =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );
            const std::shared_ptr< Constitutive_Model >& tCMFollowerElasticity =
                    mFollowerCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );

            MORIS_ASSERT( tSPNitsche != nullptr,
                    "IWG_Isotropic_Struc_Nonlinear_Interface::compute_residual - Nitsche parameter missing." );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            const real tNitsche      = tSPNitsche->val()( 0 );
            const real tLeaderWeight = tSPNitsche->val()( 1 );
            const real tFollowerWeight  = tSPNitsche->val()( 2 );

            // evaluate average traction
            const Matrix< DDRMat > tTraction =
                    tLeaderWeight * tCMLeaderElasticity->traction( mNormal, mStressType )
                    + tFollowerWeight * tCMFollowerElasticity->traction( mNormal, mStressType );

            // evaluate temperature jump
            const auto tJump = tFILeader->val() - tFIFollower->val();

            // compute leader residual
            mSet->get_residual()( 0 )( { tLeaderResStartIndex, tLeaderResStopIndex } ) +=
                    aWStar
                    * ( -tFILeader->N_trans() * tTraction
                            + mBeta * tLeaderWeight
                                      * tCMLeaderElasticity->testTraction_trans( mNormal, mResidualDofType( 0 ), mStressType ) * tJump
                            + tNitsche * tFILeader->N_trans() * tJump );

            // compute follower residual
            mSet->get_residual()( 0 )( { tFollowerResStartIndex, tFollowerResStopIndex } ) +=
                    aWStar
                    * ( +tFIFollower->N_trans() * tTraction
                            + mBeta * tFollowerWeight * tCMFollowerElasticity->testTraction_trans( mNormal, mResidualDofType( 0 ), mStressType )
                                      * tJump
                            - tNitsche * tFIFollower->N_trans() * tJump );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Isotropic_Struc_Nonlinear_Interface::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Isotropic_Struc_Nonlinear_Interface::compute_jacobian( real aWStar )
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

            // get leader field interpolator for the residual dof type
            Field_Interpolator* tFILeader =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get follower field interpolator for the residual dof type
            Field_Interpolator* tFIFollower =
                    mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the elasticity constitutive model
            const std::shared_ptr< Constitutive_Model >& tCMLeaderElasticity =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );
            const std::shared_ptr< Constitutive_Model >& tCMFollowerElasticity =
                    mFollowerCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            const real tNitsche      = tSPNitsche->val()( 0 );
            const real tLeaderWeight = tSPNitsche->val()( 1 );
            const real tFollowerWeight  = tSPNitsche->val()( 2 );

            // get number of leader dof dependencies
            const uint tLeaderNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            // evaluate displacement jump
            Matrix< DDRMat > tJump = tFILeader->val() - tFIFollower->val();

            // compute the jacobian for indirect dof dependencies through leader constitutive models
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
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    tJacMM += aWStar
                            * ( +mBeta * tLeaderWeight
                                            * tCMLeaderElasticity->testTraction_trans( mNormal, mResidualDofType( 0 ), mStressType )
                                            * tFILeader->N()
                                    + tNitsche * tFILeader->N_trans() * tFILeader->N() );

                    tJacSM += aWStar
                            * ( +mBeta * tFollowerWeight
                                            * tCMFollowerElasticity->testTraction_trans( mNormal, mResidualDofType( 0 ), mStressType )
                                            * tFILeader->N()
                                    - tNitsche * tFIFollower->N_trans() * tFILeader->N() );
                }

                // if dependency on the dof type
                if ( tCMLeaderElasticity->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    tJacMM += aWStar
                            * ( -tFILeader->N_trans() * tLeaderWeight
                                            * tCMLeaderElasticity->dTractiondDOF( tDofType, mNormal, mStressType )
                                    + mBeta * tLeaderWeight
                                              * tCMLeaderElasticity->dTestTractiondDOF(
                                                      tDofType, mNormal, tJump, mResidualDofType( 0 ), mStressType ) );

                    tJacSM += aWStar
                            * ( tFIFollower->N_trans() * tLeaderWeight
                                    * tCMLeaderElasticity->dTractiondDOF( tDofType, mNormal, mStressType ) );
                }

                // if dependency of stabilization parameters on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Leader_Follower::LEADER ) )
                {
                    // get the derivatives of the SPs
                    const Matrix< DDRMat > tNitscheDer      = tSPNitsche->dSPdLeaderDOF( tDofType ).get_row( 0 );
                    const Matrix< DDRMat > tLeaderWeightDer = tSPNitsche->dSPdLeaderDOF( tDofType ).get_row( 1 );
                    const Matrix< DDRMat > tFollowerWeightDer  = tSPNitsche->dSPdLeaderDOF( tDofType ).get_row( 2 );

                    // get traction derivative
                    const Matrix< DDRMat > tTractionDer =
                            tCMLeaderElasticity->traction( mNormal, mStressType ) * tLeaderWeightDer
                            + tCMFollowerElasticity->traction( mNormal, mStressType ) * tFollowerWeightDer;

                    // add contribution to jacobian
                    tJacMM += aWStar
                            * ( -tFILeader->N_trans() * tTractionDer
                                    + mBeta * tCMLeaderElasticity->testTraction_trans( mNormal, mResidualDofType( 0 ), mStressType )
                                              * tJump * tLeaderWeightDer
                                    + tFILeader->N_trans() * tJump * tNitscheDer );

                    tJacSM += aWStar
                            * ( +tFIFollower->N_trans() * tTractionDer
                                    + mBeta * tCMFollowerElasticity->testTraction_trans( mNormal, mResidualDofType( 0 ), mStressType ) * tJump
                                              * tFollowerWeightDer
                                    - tFIFollower->N_trans() * tJump * tNitscheDer );
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
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    tJacMS += aWStar
                            * ( -mBeta * tLeaderWeight
                                            * tCMLeaderElasticity->testTraction_trans( mNormal, mResidualDofType( 0 ), mStressType )
                                            * tFIFollower->N()
                                    - tNitsche * tFILeader->N_trans() * tFIFollower->N() );

                    tJacSS += aWStar
                            * ( -mBeta * tFollowerWeight
                                            * tCMFollowerElasticity->testTraction_trans( mNormal, mResidualDofType( 0 ), mStressType )
                                            * tFIFollower->N()
                                    + tNitsche * tFIFollower->N_trans() * tFIFollower->N() );
                }

                // if dependency on the dof type
                if ( tCMFollowerElasticity->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    tJacMS += aWStar
                            * ( -tFILeader->N_trans() * tFollowerWeight
                                    * tCMFollowerElasticity->dTractiondDOF( tDofType, mNormal, mStressType ) );

                    tJacSS += aWStar
                            * ( +tFIFollower->N_trans() * tFollowerWeight
                                            * tCMFollowerElasticity->dTractiondDOF( tDofType, mNormal, mStressType )
                                    + mBeta * tFollowerWeight
                                              * tCMFollowerElasticity->dTestTractiondDOF(
                                                      tDofType, mNormal, tJump, mResidualDofType( 0 ), mStressType ) );
                }

                // if dependency of stabilization parameters on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Leader_Follower::FOLLOWER ) )
                {
                    // get the derivatives of the SPs
                    const Matrix< DDRMat > tNitscheDer      = tSPNitsche->dSPdFollowerDOF( tDofType ).get_row( 0 );
                    const Matrix< DDRMat > tLeaderWeightDer = tSPNitsche->dSPdFollowerDOF( tDofType ).get_row( 1 );
                    const Matrix< DDRMat > tFollowerWeightDer  = tSPNitsche->dSPdFollowerDOF( tDofType ).get_row( 2 );

                    // get traction derivative
                    const Matrix< DDRMat > tTractionDer =
                            tCMLeaderElasticity->traction( mNormal, mStressType ) * tLeaderWeightDer
                            + tCMFollowerElasticity->traction( mNormal, mStressType ) * tFollowerWeightDer;

                    // add contribution to jacobian
                    tJacMS += aWStar
                            * ( -tFILeader->N_trans() * tTractionDer
                                    + mBeta * tCMLeaderElasticity->testTraction_trans( mNormal, mResidualDofType( 0 ), mStressType )
                                              * tJump * tLeaderWeightDer
                                    + tFILeader->N_trans() * tJump * tNitscheDer );

                    tJacSS += aWStar
                            * ( +tFIFollower->N_trans() * tTractionDer
                                    + mBeta * tCMFollowerElasticity->testTraction_trans( mNormal, mResidualDofType( 0 ), mStressType ) * tJump
                                              * tFollowerWeightDer
                                    - tFIFollower->N_trans() * tJump * tNitscheDer );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Isotropic_Struc_Nonlinear_Interface::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Isotropic_Struc_Nonlinear_Interface::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false,
                    "IWG_Isotropic_Struc_Nonlinear_Interface::compute_jacobian_and_residual - This function does nothing." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Isotropic_Struc_Nonlinear_Interface::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Nonlinear_Interface::compute_dRdp - This function does nothing." );
        }

        //------------------------------------------------------------------------------
}    // namespace moris::fem
