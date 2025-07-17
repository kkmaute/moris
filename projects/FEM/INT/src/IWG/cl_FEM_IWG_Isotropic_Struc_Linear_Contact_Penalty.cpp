/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Penalty.cpp
 *
 */

#include "cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Penalty.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_isfinite.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IWG_Isotropic_Struc_Linear_Contact_Penalty::IWG_Isotropic_Struc_Linear_Contact_Penalty()
    {
        // set size for the property pointer cell
        mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "Thickness" ] = static_cast< uint >( IWG_Property_Type::THICKNESS );

        // set size for the constitutive model pointer cell
        // .resize: gives aValue:(The value to initialize the new elements with) and aCount:(new size of the Cell)
        mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );
        mFollowerCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "ElastLinIso" ] = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );

        // set size for the stabilization parameter pointer cell
        mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

        // populate the stabilization map
        mStabilizationMap[ "PenaltyContact" ]     = static_cast< uint >( IWG_Stabilization_Type::PENALTY_CONTACT );
        mStabilizationMap[ "StabPenaltyContact" ] = static_cast< uint >( IWG_Stabilization_Type::STAB_PENALTY_CONTACT );

        //            mStabilizationMap[ "LeaderWeightInterface" ] = IWG_Stabilization_Type::LEADER_WEIGHT_INTERFACE;
        //            mStabilizationMap[ "FollowerWeightInterface" ]  = IWG_Stabilization_Type::FOLLOWER_WEIGHT_INTERFACE;
    }

    //------------------------------------------------------------------------------

    void IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_residual( real aWStar )
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
            Field_Interpolator *tFILeader = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get follower field interpolator for the residual dof type
            Field_Interpolator *tFIFollower = mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the elasticity constitutive model
            std::shared_ptr< Constitutive_Model > &tCMLeaderElasticity =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );
            std::shared_ptr< Constitutive_Model > &tCMFollowerElasticity =
                    mFollowerCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // get the Nitsche stabilization parameter
            std::shared_ptr< Stabilization_Parameter > &tSPPenalty = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::PENALTY_CONTACT ) );

            // get the Nitsche stabilization parameter
            std::shared_ptr< Stabilization_Parameter > &tSPStabPen = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::STAB_PENALTY_CONTACT ) );

            // get thickness property
            const std::shared_ptr< Property > &tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // evaluate average traction
            Matrix< DDRMat > tJumpTraction = tCMLeaderElasticity->traction( mNormal ) - tCMFollowerElasticity->traction( mNormal );
            Matrix< DDRMat > tJumpPressure = trans( mNormal ) * tJumpTraction;

            // moris::print( tJumpTraction, "tJumpTraction" );

            // evaluate gap
            Matrix< DDRMat > tGap = trans( mNormal ) * ( tFIFollower->val() - tFILeader->val() );    // mNormal is normal on Leaderside

            if ( tGap( 0, 0 ) > 0 )
            {
                tGap( 0, 0 ) = 0;
            }

            // compute contact residual on follower side
            mSet->get_residual()( 0 )( { tFollowerResStartIndex, tFollowerResStopIndex }, { 0, 0 } ) +=
                    tSPPenalty->val()( 0 ) * trans( tFIFollower->N() ) * mNormal * tGap( 0, 0 ) * aWStar;

            // compute contact residual on leader side
            mSet->get_residual()( 0 )( { tLeaderResStartIndex, tLeaderResStopIndex }, { 0, 0 } ) +=
                    ( -1 ) * tSPPenalty->val()( 0 ) * trans( tFILeader->N() ) * mNormal * tGap( 0, 0 ) * aWStar;

            // if stabilized Penalty
            if ( tSPStabPen != nullptr )
            {
                // compute contact residual on follower side
                mSet->get_residual()( 0 )(
                        { tFollowerResStartIndex, tFollowerResStopIndex },
                        { 0, 0 } ) += aWStar * ( tSPStabPen->val()( 0 ) * tCMFollowerElasticity->testTraction_trans( mNormal, mResidualDofType( 0 ) ) * mNormal * tJumpPressure );

                // compute contact residual on leader side
                mSet->get_residual()( 0 )(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { 0, 0 } ) += aWStar * ( tSPStabPen->val()( 0 ) * tCMLeaderElasticity->testTraction_trans( mNormal, mResidualDofType( 0 ) ) * mNormal * tJumpPressure );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader and follower field interpolators
            this->check_field_interpolators( mtk::Leader_Follower::LEADER );
            this->check_field_interpolators( mtk::Leader_Follower::FOLLOWER );
#endif

            //             // get leader index for residual dof type
            //             uint tDofIndexLeader = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            //
            //             // get follower index for residual dof type
            //             uint tDofIndexFollower  = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );
            //
            //             // get leader field interpolator for the residual dof type
            //             Field_Interpolator * tFILeader = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));
            //
            //             // get follower field interpolator for the residual dof type
            //             Field_Interpolator * tFIFollower  = mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));
            //
            //             // get indices for SP, CM and properties
            //             uint tElastLinIsoIndex  = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );
            //             uint tPenIndex          = static_cast< uint >( IWG_Stabilization_Type::PENALTY_CONTACT );
            //             uint tStabPenIndex      = static_cast< uint >( IWG_Stabilization_Type::STAB_PENALTY_CONTACT );

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get follower index for residual dof type, indices for assembly
            uint tFollowerDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );
            uint tFollowerResStartIndex = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 0 );
            uint tFollowerResStopIndex  = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 1 );

            // get leader field interpolator for the residual dof type
            Field_Interpolator *tFILeader = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get follower field interpolator for the residual dof type
            Field_Interpolator *tFIFollower = mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the elasticity constitutive model
            std::shared_ptr< Constitutive_Model > &tCMLeaderElasticity =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );
            std::shared_ptr< Constitutive_Model > &tCMFollowerElasticity =
                    mFollowerCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // get the Nitsche stabilization parameter
            std::shared_ptr< Stabilization_Parameter > &tSPPenalty =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::PENALTY_CONTACT ) );

            // get the Nitsche stabilization parameter
            std::shared_ptr< Stabilization_Parameter > &tSPStabPen =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::STAB_PENALTY_CONTACT ) );

            // get thickness property
            const std::shared_ptr< Property > &tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // evaluate average traction
            Matrix< DDRMat > tJumpTraction =
                    tCMLeaderElasticity->traction( mNormal ) - tCMFollowerElasticity->traction( mNormal );

            // get number of leader dof dependencies
            uint tLeaderNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            // create quadrature of normal
            Matrix< DDRMat > tNormQ = mNormal * trans( mNormal );

            // evaluate gap
            Matrix< DDRMat > tGap = trans( mNormal ) * ( tFIFollower->val() - tFILeader->val() );    // mNormal is normal on Leader side

            if ( tGap( 0, 0 ) <= 0 )
            {
                // hier
                // compute the jacobian through leader constitutive models
                for ( uint iDOF = 0; iDOF < tLeaderNumDofDependencies; iDOF++ )
                {
                    // get the dof type
                    Vector< MSI::Dof_Type > &tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                    // get the index for the dof type
                    sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                    uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                    uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                    // compute jacobian direct dependencies
                    if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                    {
                        mSet->get_jacobian()( { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) +=
                                tSPPenalty->val()( 0 ) * trans( tFILeader->N() ) * tNormQ * tFILeader->N() * aWStar;

                        mSet->get_jacobian()( { tFollowerResStartIndex, tFollowerResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) +=
                                ( -1 ) * tSPPenalty->val()( 0 ) * trans( tFILeader->N() ) * tNormQ * tFIFollower->N() * aWStar;
                    }
                }

                // compute the jacobian through follower constitutive models
                for ( uint iDOF = 0; iDOF < tLeaderNumDofDependencies; iDOF++ )
                {
                    // get the dof type
                    Vector< MSI::Dof_Type > &tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                    // get the index for the dof type
                    sint tDofDepIndex           = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::FOLLOWER );
                    uint tFollowerDepStartIndex = mSet->get_jac_dof_assembly_map()( tFollowerDofIndex )( tDofDepIndex, 0 );
                    uint tFollowerDepStopIndex  = mSet->get_jac_dof_assembly_map()( tFollowerDofIndex )( tDofDepIndex, 1 );

                    // compute jacobian direct dependencies
                    if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                    {
                        mSet->get_jacobian()( { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tFollowerDepStartIndex, tFollowerDepStopIndex } ) +=
                                ( -1 ) * tSPPenalty->val()( 0 ) * trans( tFIFollower->N() ) * tNormQ * tFILeader->N() * aWStar;

                        mSet->get_jacobian()( { tFollowerResStartIndex, tFollowerResStopIndex },
                                { tFollowerDepStartIndex, tFollowerDepStopIndex } ) +=
                                tSPPenalty->val()( 0 ) * trans( tFIFollower->N() ) * tNormQ * tFIFollower->N() * aWStar;
                    }
                }

                // hier

                //                 // compute the jacobian for direct dof dependencies
                //                 if ( mResidualDofTypeRequested )
                //                 {
                //                 // compute jacobian fon follower side
                //                 // follower follower
                //                 mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexFollower )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexFollower )( 0, 1 ) },
                //                         { mSet->get_jac_dof_assembly_map()( tDofIndexFollower )( tDofIndexFollower, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexFollower )( tDofIndexFollower, 1 ) } )
                //                         += tSPPenalty->val()( 0 ) * trans( tFIFollower -> N()) * tNormQ * tFIFollower -> N() * aWStar;
                //
                //                 // follower leader
                //                 mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexFollower )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexFollower )( 0, 1 ) },
                //                         { mSet->get_jac_dof_assembly_map()( tDofIndexFollower )( tDofIndexLeader, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexFollower )( tDofIndexLeader, 1 ) } )
                //                         += (-1) * tSPPenalty->val()( 0 ) * trans( tFIFollower -> N()) * tNormQ * tFILeader -> N() * aWStar;
                //
                //
                //                 // compute the jacobian on leader side
                //                 //leader leader
                //                 mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexLeader )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexLeader )( 0, 1 ) },
                //                         { mSet->get_jac_dof_assembly_map()( tDofIndexLeader )( tDofIndexLeader, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexLeader )( tDofIndexLeader, 1 ) } )
                //                         += mStabilizationParam( tPenIndex )->val()( 0 ) * trans( tFILeader -> N()) * tNormQ * tFILeader -> N() * tWStar;
                //
                //                 //leader follower
                //                 mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexLeader )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexLeader )( 0, 1 ) },
                //                         { mSet->get_jac_dof_assembly_map()( tDofIndexLeader  )( tDofIndexFollower, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexLeader )( tDofIndexFollower, 1 ) } )
                //                         += (-1) * mStabilizationParam( tPenIndex )->val()( 0 ) * trans( tFILeader -> N()) * tNormQ * tFIFollower -> N() * tWStar;
                //                 }
            }

            // if stabilized Penalty
            if ( tSPStabPen != nullptr )
            {
                // hier

                // compute the jacobian for indirect dof dependencies through leader constitutive models
                for ( uint iDOF = 0; iDOF < tLeaderNumDofDependencies; iDOF++ )
                {
                    // get the dof type
                    Vector< MSI::Dof_Type > &tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                    // get the index for the dof type
                    sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                    uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                    uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                    // if dependency on the dof type
                    if ( tCMLeaderElasticity->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()( { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) +=
                                tSPStabPen->val()( 0 ) * tCMLeaderElasticity->testTraction_trans( mNormal, mResidualDofType( 0 ) )
                                * tNormQ * tCMLeaderElasticity->dTractiondDOF( tDofType, mNormal ) * aWStar;

                        mSet->get_jacobian()( { tFollowerResStartIndex, tFollowerResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) +=
                                tSPStabPen->val()( 0 ) * tCMLeaderElasticity->testTraction_trans( mNormal, mResidualDofType( 0 ) )
                                * tNormQ * tCMFollowerElasticity->dTractiondDOF( tDofType, mNormal ) * aWStar;
                    }
                }

                // compute the jacobian for indirect dof dependencies through follower constitutive models
                uint tFollowerNumDofDependencies = mRequestedFollowerGlobalDofTypes.size();
                for ( uint iDOF = 0; iDOF < tFollowerNumDofDependencies; iDOF++ )
                {
                    // get the dof type
                    Vector< MSI::Dof_Type > &tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                    // get the index for the dof type
                    sint tDofDepIndex           = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::FOLLOWER );
                    uint tFollowerDepStartIndex = mSet->get_jac_dof_assembly_map()( tFollowerDofIndex )( tDofDepIndex, 0 );
                    uint tFollowerDepStopIndex  = mSet->get_jac_dof_assembly_map()( tFollowerDofIndex )( tDofDepIndex, 1 );

                    // if dependency on the dof type
                    if ( tCMFollowerElasticity->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()( { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tFollowerDepStartIndex, tFollowerDepStopIndex } ) -=
                                tSPStabPen->val()( 0 ) * tCMFollowerElasticity->testTraction_trans( mNormal, mResidualDofType( 0 ) )
                                * tNormQ * tCMLeaderElasticity->dTractiondDOF( tDofType, mNormal ) * aWStar;

                        mSet->get_jacobian()( { tFollowerResStartIndex, tFollowerResStopIndex },
                                { tFollowerDepStartIndex, tFollowerDepStopIndex } ) -=
                                tSPStabPen->val()( 0 ) * tCMFollowerElasticity->testTraction_trans( mNormal, mResidualDofType( 0 ) )
                                * tNormQ * tCMFollowerElasticity->dTractiondDOF( tDofType, mNormal ) * aWStar;
                    }
                }
                // hier ende

                //                 // compute the jacobian for indirect dof dependencies through leader constitutive models
                //                 for( uint iDOF = 0; iDOF < tLeaderNumDofDependencies; iDOF++ )
                //                 {
                //                     // get the dof type
                //                     Vector< MSI::Dof_Type > & tDofType = mRequestedLeaderGlobalDofTypes( iDOF );
                //
                //                     // get index for the dof type
                //                     sint tIndexDep = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                //
                //                     // if dependency on the dof type
                //
                //                     if ( mLeaderCM( tElastLinIsoIndex )->check_dof_dependency( tDofType ) )
                //                     {
                //                         // add contribution to jacobian
                //
                //                         // leader leader
                //                         mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexLeader )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexLeader )( 0, 1 ) },
                //                                               { mSet->get_jac_dof_assembly_map()( tDofIndexLeader )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexLeader )( tIndexDep, 1 ) } )
                //                         += mStabilizationParam( tStabPenIndex )->val()( 0 ) * mLeaderCM( tElastLinIsoIndex )->testTraction( mNormal )
                //                             * mLeaderCM( tElastLinIsoIndex )->dTractiondDOF( tDofType, mNormal ) * tWStar;
                //
                //                         // leader follower
                //                         mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexFollower )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexFollower )( 0, 1 ) },
                //                                               { mSet->get_jac_dof_assembly_map()( tDofIndexFollower )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexFollower )( tIndexDep, 1 ) } )
                //                         += mStabilizationParam( tStabPenIndex )->val()( 0 ) * mLeaderCM( tElastLinIsoIndex )->testTraction( mNormal )
                //                             * mFollowerCM( tElastLinIsoIndex )->dTractiondDOF( tDofType, mNormal ) * tWStar;
                //                     }
                //                 }
                //
                //                 // compute the jacobian for indirect dof dependencies through follower constitutive models
                //                 uint tFollowerNumDofDependencies = mRequestedFollowerGlobalDofTypes.size();
                //                 for( uint iDOF = 0; iDOF < tFollowerNumDofDependencies; iDOF++ )
                //                 {
                //                     // get dof type
                //                     Vector< MSI::Dof_Type > tDofType = mRequestedFollowerGlobalDofTypes( iDOF );
                //
                //                     // get index for dof type
                //                     sint tIndexDep = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::FOLLOWER );
                //
                //                     // if dependency on the dof type
                //                     if ( mFollowerCM( tElastLinIsoIndex )->check_dof_dependency( tDofType ) )
                //                     {
                //                         // add contribution to jacobian
                //
                //                         // follower leader
                //                         mSet->get_jacobian()( { tLeaderResStartIndex, tLeaderResStopIndex },
                //                                               { tFollowerDepStartIndex,  tFollowerDepStopIndex  } )
                //                         += mStabilizationParam( tStabPenIndex )->val()( 0 ) * mFollowerCM( tElastLinIsoIndex )->testTraction( mNormal )
                //                             * mLeaderCM( tElastLinIsoIndex )->dTractiondDOF( tDofType, mNormal ) * tWStar;
                //
                //                         //follower follower
                //                         mSet->get_jacobian()( { tFollowerResStartIndex, tFollowerResStopIndex },
                //                                               { tFollowerDepStartIndex, tFollowerDepStopIndex } )
                //                         += mStabilizationParam( tStabPenIndex )->val()( 0 ) * mFollowerCM( tElastLinIsoIndex )->testTraction( mNormal )
                //                             * mFollowerCM( tElastLinIsoIndex )->dTractiondDOF( tDofType, mNormal ) * tWStar;
                //                     }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_jacobian_and_residual - This function does nothing." );
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_dRdp - This function does nothing." );
        }

        //------------------------------------------------------------------------------
}    // namespace moris::fem
