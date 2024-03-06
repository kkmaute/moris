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
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IWG_Isotropic_Struc_Linear_Contact_Penalty::IWG_Isotropic_Struc_Linear_Contact_Penalty()
        {
            init_property( "Thickness", IWG_Property_Type::THICKNESS );
            init_constitutive_model( "ElastLinIso", IWG_Constitutive_Type::ELAST_LIN_ISO );
            init_stabilization_parameter( "PenaltyContact", IWG_Stabilization_Type::PENALTY_CONTACT );
            init_stabilization_parameter( "StabPenaltyContact", IWG_Stabilization_Type::STAB_PENALTY_CONTACT );

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
            Field_Interpolator *tFILeader = get_leader_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get follower field interpolator for the residual dof type
            Field_Interpolator *tFIFollower = get_follower_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the elasticity constitutive model
            std::shared_ptr< Constitutive_Model > const &tCMLeaderElasticity   = get_leader_constitutive_model( IWG_Constitutive_Type::ELAST_LIN_ISO );
            std::shared_ptr< Constitutive_Model > const &tCMFollowerElasticity = get_follower_constitutive_model( IWG_Constitutive_Type::ELAST_LIN_ISO );

            // get the Nitsche stabilization parameter
            std::shared_ptr< Stabilization_Parameter > const &tSPPenalty = get_stabilization_parameter( IWG_Stabilization_Type::PENALTY_CONTACT );

            // get the Nitsche stabilization parameter
            std::shared_ptr< Stabilization_Parameter > const &tSPStabPen = get_stabilization_parameter( IWG_Stabilization_Type::STAB_PENALTY_CONTACT );

            // get thickness property
            const std::shared_ptr< Property > &tPropThickness = get_leader_property( IWG_Property_Type::THICKNESS );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // evaluate average traction
            Matrix< DDRMat > tJumpTraction = tCMLeaderElasticity->traction( get_normal() ) - tCMFollowerElasticity->traction( get_normal() );
            Matrix< DDRMat > tJumpPressure = trans( get_normal() ) * tJumpTraction;

            // moris::print( tJumpTraction, "tJumpTraction" );

            // evaluate gap
            Matrix< DDRMat > tGap = trans( get_normal() ) * ( tFIFollower->val() - tFILeader->val() );    // get_normal() is normal on Leaderside

            if ( tGap( 0, 0 ) > 0 )
            {
                tGap( 0, 0 ) = 0;
            }

            // compute contact residual on follower side
            mSet->get_residual()( 0 )( { tFollowerResStartIndex, tFollowerResStopIndex }, { 0, 0 } ) += tSPPenalty->val()( 0 ) * trans( tFIFollower->N() ) * get_normal() * tGap( 0, 0 ) * aWStar;

            // compute contact residual on leader side
            mSet->get_residual()( 0 )( { tLeaderResStartIndex, tLeaderResStopIndex }, { 0, 0 } ) += ( -1 ) * tSPPenalty->val()( 0 ) * trans( tFILeader->N() ) * get_normal() * tGap( 0, 0 ) * aWStar;

            // if stabilized Penalty
            if ( tSPStabPen != nullptr )
            {
                // compute contact residual on follower side
                mSet->get_residual()( 0 )(
                        { tFollowerResStartIndex, tFollowerResStopIndex },
                        { 0, 0 } ) += aWStar * ( tSPStabPen->val()( 0 ) * tCMFollowerElasticity->testTraction_trans( get_normal(), mResidualDofType( 0 ) ) * get_normal() * tJumpPressure );

                // compute contact residual on leader side
                mSet->get_residual()( 0 )(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { 0, 0 } ) += aWStar * ( tSPStabPen->val()( 0 ) * tCMLeaderElasticity->testTraction_trans( get_normal(), mResidualDofType( 0 ) ) * get_normal() * tJumpPressure );
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
            //             Field_Interpolator * tFILeader = get_leader_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));
            //
            //             // get follower field interpolator for the residual dof type
            //             Field_Interpolator * tFIFollower  = get_follower_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));
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
            Field_Interpolator *tFILeader = get_leader_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get follower field interpolator for the residual dof type
            Field_Interpolator *tFIFollower = get_follower_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the elasticity constitutive model
            std::shared_ptr< Constitutive_Model > const &tCMLeaderElasticity   = get_leader_constitutive_model( IWG_Constitutive_Type::ELAST_LIN_ISO );
            std::shared_ptr< Constitutive_Model > const &tCMFollowerElasticity = get_follower_constitutive_model( IWG_Constitutive_Type::ELAST_LIN_ISO );

            // get the Nitsche stabilization parameter
            std::shared_ptr< Stabilization_Parameter > const &tSPPenalty = get_stabilization_parameter( IWG_Stabilization_Type::PENALTY_CONTACT );

            // get the Nitsche stabilization parameter
            std::shared_ptr< Stabilization_Parameter > const &tSPStabPen = get_stabilization_parameter( IWG_Stabilization_Type::STAB_PENALTY_CONTACT );

            // get thickness property
            const std::shared_ptr< Property > &tPropThickness = get_leader_property( IWG_Property_Type::THICKNESS );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // evaluate average traction
            Matrix< DDRMat > tJumpTraction = tCMLeaderElasticity->traction( get_normal() ) - tCMFollowerElasticity->traction( get_normal() );

            // get number of leader dof dependencies
            uint tLeaderNumDofDependencies = get_requested_leader_dof_types().size();

            // create quadrature of normal
            Matrix< DDRMat > tNormQ = get_normal() * trans( get_normal() );

            // evaluate gap
            Matrix< DDRMat > tGap = trans( get_normal() ) * ( tFIFollower->val() - tFILeader->val() );    // get_normal() is normal on Leader side

            if ( tGap( 0, 0 ) <= 0 )
            {
                // hier
                // compute the jacobian through leader constitutive models
                for ( uint iDOF = 0; iDOF < tLeaderNumDofDependencies; iDOF++ )
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
                        mSet->get_jacobian()( { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += tSPPenalty->val()( 0 ) * trans( tFILeader->N() ) * tNormQ * tFILeader->N() * aWStar;

                        mSet->get_jacobian()( { tFollowerResStartIndex, tFollowerResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += ( -1 ) * tSPPenalty->val()( 0 ) * trans( tFILeader->N() ) * tNormQ * tFIFollower->N() * aWStar;
                    }
                }

                // compute the jacobian through follower constitutive models
                for ( uint iDOF = 0; iDOF < tLeaderNumDofDependencies; iDOF++ )
                {
                    // get the dof type
                    Vector< MSI::Dof_Type > const &tDofType = get_requested_leader_dof_types()( iDOF );

                    // get the index for the dof type
                    sint tDofDepIndex           = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::FOLLOWER );
                    uint tFollowerDepStartIndex = mSet->get_jac_dof_assembly_map()( tFollowerDofIndex )( tDofDepIndex, 0 );
                    uint tFollowerDepStopIndex  = mSet->get_jac_dof_assembly_map()( tFollowerDofIndex )( tDofDepIndex, 1 );

                    // compute jacobian direct dependencies
                    if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                    {
                        mSet->get_jacobian()( { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tFollowerDepStartIndex, tFollowerDepStopIndex } ) += ( -1 ) * tSPPenalty->val()( 0 ) * trans( tFIFollower->N() ) * tNormQ * tFILeader->N() * aWStar;

                        mSet->get_jacobian()( { tFollowerResStartIndex, tFollowerResStopIndex },
                                { tFollowerDepStartIndex, tFollowerDepStopIndex } ) += tSPPenalty->val()( 0 ) * trans( tFIFollower->N() ) * tNormQ * tFIFollower->N() * aWStar;
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
                    Vector< MSI::Dof_Type > const &tDofType = get_requested_leader_dof_types()( iDOF );

                    // get the index for the dof type
                    sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                    uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                    uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                    // if dependency on the dof type
                    if ( tCMLeaderElasticity->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()( { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += tSPStabPen->val()( 0 ) * tCMLeaderElasticity->testTraction_trans( get_normal(), mResidualDofType( 0 ) )
                                                                                 * tNormQ * tCMLeaderElasticity->dTractiondDOF( tDofType, get_normal() ) * aWStar;

                        mSet->get_jacobian()( { tFollowerResStartIndex, tFollowerResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += tSPStabPen->val()( 0 ) * tCMLeaderElasticity->testTraction_trans( get_normal(), mResidualDofType( 0 ) )
                                                                                 * tNormQ * tCMFollowerElasticity->dTractiondDOF( tDofType, get_normal() ) * aWStar;
                    }
                }

                // compute the jacobian for indirect dof dependencies through follower constitutive models
                uint tFollowerNumDofDependencies = get_requested_follower_dof_types().size();
                for ( uint iDOF = 0; iDOF < tFollowerNumDofDependencies; iDOF++ )
                {
                    // get the dof type
                    Vector< MSI::Dof_Type > const &tDofType = get_requested_leader_dof_types()( iDOF );

                    // get the index for the dof type
                    sint tDofDepIndex           = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::FOLLOWER );
                    uint tFollowerDepStartIndex = mSet->get_jac_dof_assembly_map()( tFollowerDofIndex )( tDofDepIndex, 0 );
                    uint tFollowerDepStopIndex  = mSet->get_jac_dof_assembly_map()( tFollowerDofIndex )( tDofDepIndex, 1 );

                    // if dependency on the dof type
                    if ( tCMFollowerElasticity->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()( { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tFollowerDepStartIndex, tFollowerDepStopIndex } ) -= tSPStabPen->val()( 0 ) * tCMFollowerElasticity->testTraction_trans( get_normal(), mResidualDofType( 0 ) )
                                                                                     * tNormQ * tCMLeaderElasticity->dTractiondDOF( tDofType, get_normal() ) * aWStar;

                        mSet->get_jacobian()( { tFollowerResStartIndex, tFollowerResStopIndex },
                                { tFollowerDepStartIndex, tFollowerDepStopIndex } ) -= tSPStabPen->val()( 0 ) * tCMFollowerElasticity->testTraction_trans( get_normal(), mResidualDofType( 0 ) )
                                                                                     * tNormQ * tCMFollowerElasticity->dTractiondDOF( tDofType, get_normal() ) * aWStar;
                    }
                }
                // hier ende

                //                 // compute the jacobian for indirect dof dependencies through leader constitutive models
                //                 for( uint iDOF = 0; iDOF < tLeaderNumDofDependencies; iDOF++ )
                //                 {
                //                     // get the dof type
                //                     Vector< MSI::Dof_Type > const &tDofType = get_requested_leader_dof_types()( iDOF );
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
                //                         += mStabilizationParam( tStabPenIndex )->val()( 0 ) * mLeaderCM( tElastLinIsoIndex )->testTraction( get_normal() )
                //                             * mLeaderCM( tElastLinIsoIndex )->dTractiondDOF( tDofType, get_normal() ) * tWStar;
                //
                //                         // leader follower
                //                         mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexFollower )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexFollower )( 0, 1 ) },
                //                                               { mSet->get_jac_dof_assembly_map()( tDofIndexFollower )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexFollower )( tIndexDep, 1 ) } )
                //                         += mStabilizationParam( tStabPenIndex )->val()( 0 ) * mLeaderCM( tElastLinIsoIndex )->testTraction( get_normal() )
                //                             * mFollowerCM( tElastLinIsoIndex )->dTractiondDOF( tDofType, get_normal() ) * tWStar;
                //                     }
                //                 }
                //
                //                 // compute the jacobian for indirect dof dependencies through follower constitutive models
                //                 uint tFollowerNumDofDependencies = get_requested_follower_dof_types().size();
                //                 for( uint iDOF = 0; iDOF < tFollowerNumDofDependencies; iDOF++ )
                //                 {
                //                     // get dof type
                //                     Vector< MSI::Dof_Type > tDofType = get_requested_follower_dof_types()( iDOF );
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
                //                         += mStabilizationParam( tStabPenIndex )->val()( 0 ) * mFollowerCM( tElastLinIsoIndex )->testTraction( get_normal() )
                //                             * mLeaderCM( tElastLinIsoIndex )->dTractiondDOF( tDofType, get_normal() ) * tWStar;
                //
                //                         //follower follower
                //                         mSet->get_jacobian()( { tFollowerResStartIndex, tFollowerResStopIndex },
                //                                               { tFollowerDepStartIndex, tFollowerDepStopIndex } )
                //                         += mStabilizationParam( tStabPenIndex )->val()( 0 ) * mFollowerCM( tElastLinIsoIndex )->testTraction( get_normal() )
                //                             * mFollowerCM( tElastLinIsoIndex )->dTractiondDOF( tDofType, get_normal() ) * tWStar;
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
    } /* namespace fem */
} /* namespace moris */
