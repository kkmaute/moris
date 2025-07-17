/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Gap_Nitsche.cpp
 *
 */

#include "cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Gap_Nitsche.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_dot.hpp"
#include "fn_isfinite.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IWG_Isotropic_Struc_Linear_Contact_Gap_Nitsche::IWG_Isotropic_Struc_Linear_Contact_Gap_Nitsche( sint aBeta )
    {
        // sign for symmetric/unsymmetric Nitsche
        mBeta = aBeta;

        // set size for the property pointer cell
        mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "Thickness" ] = static_cast< uint >( IWG_Property_Type::THICKNESS );
        mPropertyMap[ "Gap" ]       = static_cast< uint >( IWG_Property_Type::GAP );

        // set size for the constitutive model pointer cell
        // .resize: gives aValue:(The value to initialize the new elements with) and aCount:(new size of the Cell)
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
    IWG_Isotropic_Struc_Linear_Contact_Gap_Nitsche::compute_residual( real aWStar )
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
            Field_Interpolator* tFILeader = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get follower field interpolator for the residual dof type
            Field_Interpolator* tFIFollower = mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

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

            // normals on leader and follower side (needs to be generalized)
            const Matrix< DDRMat > tNormalLeader = mNormal;
            const Matrix< DDRMat > tNormalFollower  = -1.0 * mNormal;

            // projections from follower to leader and leader to follower (needs to be generalized)
            const real tPrjFollowerToLeader      = 1.0;    // should be matrix
            const real tPrjLeaderToFollower      = 1.0;    // should be matrix
            const real tPrjFollowerToLeaderTrans = 1.0;    // should be matrix
            const real tPrjLeaderToFollowerTrans = 1.0;    // should be matrix

            // normal projection operator
            const Matrix< DDRMat > tNormalProjectorLeader = tNormalLeader * trans( tNormalLeader );
            const Matrix< DDRMat > tNormalProjectorFollower  = tNormalFollower * trans( tNormalFollower );

            // compute the jump
            const Matrix< DDRMat > tJumpLeader = tFILeader->val() - tPrjFollowerToLeader * tFIFollower->val();
            const Matrix< DDRMat > tJumpFollower  = tFIFollower->val() - tPrjLeaderToFollower * tFILeader->val();

            // compute projection of displacement jump onto normal
            const real tNormalJumpLeader = dot( tJumpLeader, tNormalLeader );
            const real tNormalJumpFollower  = dot( tJumpFollower, tNormalFollower );

            // evaluate tractions
            const Matrix< DDRMat > tTractionLeader = tCMLeaderElasticity->traction( tNormalLeader );
            const Matrix< DDRMat > tTractionFollower  = tCMFollowerElasticity->traction( tNormalFollower );

            // compute contact pressure
            const real tIfcPressureLeader = dot( tTractionLeader, tNormalLeader );
            const real tIfcPressureFollower  = dot( tTractionFollower, tNormalFollower );

            // check for contact on leader side
            if ( tIfcPressureLeader - tNitsche * tNormalJumpLeader < 0 )
            {
                // compute leader residual
                mSet->get_residual()( 0 )(
                        { tLeaderResStartIndex, tLeaderResStopIndex } ) +=                                                                                          //
                        aWStar * tLeaderWeight * (                                                                                                                  //
                                -tFILeader->N_trans() * tNormalProjectorLeader * tTractionLeader                                                                    //
                                + mBeta * tCMLeaderElasticity->testTraction_trans( tNormalLeader, mResidualDofType( 0 ) ) * tNormalProjectorLeader * tJumpLeader    //
                                + tNitsche * tFILeader->N_trans() * tNormalProjectorLeader * tJumpLeader );

                mSet->get_residual()( 0 )(
                        { tFollowerResStartIndex, tFollowerResStopIndex } ) +=                                                    //
                        aWStar * tLeaderWeight * (                                                                          //
                                +tFIFollower->N_trans() * tPrjFollowerToLeaderTrans * tNormalProjectorLeader * tTractionLeader    //
                                - tNitsche * tFIFollower->N_trans() * tPrjFollowerToLeaderTrans * tNormalProjectorLeader * tJumpLeader );
            }
            else
            {
                mSet->get_residual()( 0 )(
                        { tLeaderResStartIndex, tLeaderResStopIndex } ) +=    //
                        aWStar * tLeaderWeight * (                            //
                                -mBeta / tNitsche * tCMLeaderElasticity->testTraction_trans( tNormalLeader, mResidualDofType( 0 ) ) * tNormalProjectorLeader * tTractionLeader );
            }

            // check for contact on follower side
            if ( tIfcPressureFollower - tNitsche * tNormalJumpFollower < 0 )
            {
                // compute follower residual
                mSet->get_residual()( 0 )(
                        { tLeaderResStartIndex, tLeaderResStopIndex } ) +=                                                 //
                        aWStar * tFollowerWeight * (                                                                          //
                                +tFILeader->N_trans() * tPrjLeaderToFollowerTrans * tNormalProjectorFollower * tTractionFollower    //
                                - tNitsche * tFILeader->N_trans() * tPrjLeaderToFollowerTrans * tNormalProjectorFollower * tJumpFollower );

                mSet->get_residual()( 0 )(
                        { tFollowerResStartIndex, tFollowerResStopIndex } ) +=                                                                                        //
                        aWStar * tFollowerWeight * (                                                                                                               //
                                -tFIFollower->N_trans() * tNormalProjectorFollower * tTractionFollower                                                                   //
                                + mBeta * tCMFollowerElasticity->testTraction_trans( tNormalFollower, mResidualDofType( 0 ) ) * tNormalProjectorFollower * tJumpFollower    //
                                + tNitsche * tFIFollower->N_trans() * tNormalProjectorFollower * tJumpFollower );
            }
            else
            {
                mSet->get_residual()( 0 )(
                        { tFollowerResStartIndex, tFollowerResStopIndex } ) +=    //
                        aWStar * tFollowerWeight * (                           //
                                -mBeta / tNitsche * tCMFollowerElasticity->testTraction_trans( tNormalFollower, mResidualDofType( 0 ) ) * tNormalProjectorFollower * tTractionFollower );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Isotropic_Struc_Linear_Contact_Gap::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Isotropic_Struc_Linear_Contact_Gap_Nitsche::compute_jacobian( real aWStar )
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
            Field_Interpolator* tFILeader = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get follower field interpolator for the residual dof type
            Field_Interpolator* tFIFollower = mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

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

            // normals on leader and follower side (needs to be generalized)
            const Matrix< DDRMat > tNormalLeader = mNormal;
            const Matrix< DDRMat > tNormalFollower  = -1.0 * mNormal;

            // projections from follower to leader and leader to follower (needs to be generalized)
            const real tPrjFollowerToLeader      = 1.0;    // should be matrix
            const real tPrjLeaderToFollower      = 1.0;    // should be matrix
            const real tPrjFollowerToLeaderTrans = 1.0;    // should be matrix
            const real tPrjLeaderToFollowerTrans = 1.0;    // should be matrix

            // normal projection operator
            const Matrix< DDRMat > tNormalProjectorLeader = tNormalLeader * trans( tNormalLeader );
            const Matrix< DDRMat > tNormalProjectorFollower  = tNormalFollower * trans( tNormalFollower );

            // compute the jump
            const Matrix< DDRMat > tJumpLeader = tFILeader->val() - tPrjFollowerToLeader * tFIFollower->val();
            const Matrix< DDRMat > tJumpFollower  = tFIFollower->val() - tPrjLeaderToFollower * tFILeader->val();

            // compute projection of displacement jump onto normal
            const real tNormalJumpLeader = dot( tJumpLeader, tNormalLeader );
            const real tNormalJumpFollower  = dot( tJumpFollower, tNormalFollower );

            // evaluate tractions
            const Matrix< DDRMat > tTractionLeader = tCMLeaderElasticity->traction( tNormalLeader );
            const Matrix< DDRMat > tTractionFollower  = tCMFollowerElasticity->traction( tNormalFollower );

            // compute contact pressure
            const real tIfcPressureLeader = dot( tTractionLeader, tNormalLeader );
            const real tIfcPressureFollower  = dot( tTractionFollower, tNormalFollower );

            // get number of leader dof dependencies
            const uint tLeaderNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            // compute the Jacobian for indirect dof dependencies through leader constitutive models
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
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } );

                auto tJacSM = mSet->get_jacobian()(
                        { tFollowerResStartIndex, tFollowerResStopIndex },
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } );

                // compute Jacobian direct dependencies
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    if ( tIfcPressureLeader - tNitsche * tNormalJumpLeader < 0 )
                    {
                        tJacMM += aWStar * tLeaderWeight * (                                                                                                                    //
                                          +mBeta * tCMLeaderElasticity->testTraction_trans( tNormalLeader, mResidualDofType( 0 ) ) * tNormalProjectorLeader * tFILeader->N()    //
                                          + tNitsche * tFILeader->N_trans() * tNormalProjectorLeader * tFILeader->N() );

                        tJacSM += aWStar * tLeaderWeight * (    //
                                          -tNitsche * tFIFollower->N_trans() * tPrjFollowerToLeaderTrans * tNormalProjectorLeader * tFILeader->N() );
                    }

                    if ( tIfcPressureFollower - tNitsche * tNormalJumpFollower < 0 )
                    {
                        tJacMM += aWStar * tFollowerWeight * (    //
                                          +tNitsche * tFILeader->N_trans() * tPrjLeaderToFollowerTrans * tNormalProjectorFollower * tPrjLeaderToFollower * tFILeader->N() );

                        tJacSM += aWStar * tFollowerWeight * (                                                                                                                                      //
                                          -mBeta * tCMFollowerElasticity->testTraction_trans( tNormalFollower, mResidualDofType( 0 ) ) * tNormalProjectorFollower * tPrjLeaderToFollower * tFILeader->N()    //
                                          - tNitsche * tFIFollower->N_trans() * tNormalProjectorFollower * tPrjLeaderToFollower * tFILeader->N() );
                    }
                }

                // if dependency on the dof type
                if ( tCMLeaderElasticity->check_dof_dependency( tDofType ) )
                {
                    if ( tIfcPressureLeader - tNitsche * tNormalJumpLeader < 0 )
                    {
                        // add contribution to Jacobian
                        tJacMM += aWStar * tLeaderWeight * (                                                                                                //
                                          -tFILeader->N_trans() * tNormalProjectorLeader * tCMLeaderElasticity->dTractiondDOF( tDofType, tNormalLeader )    //
                                          + mBeta * tCMLeaderElasticity->dTestTractiondDOF( tDofType, tNormalLeader, tNormalProjectorLeader * tJumpLeader, mResidualDofType( 0 ) ) );

                        tJacSM += aWStar * tLeaderWeight * (    //
                                          +tFIFollower->N_trans() * tPrjFollowerToLeaderTrans * tNormalProjectorLeader * tCMLeaderElasticity->dTractiondDOF( tDofType, tNormalLeader ) );
                    }
                    else
                    {
                        tJacMM += aWStar * tLeaderWeight * -mBeta / tNitsche * (                                                                                                                                      //
                                          tCMLeaderElasticity->testTraction_trans( tNormalLeader, mResidualDofType( 0 ) ) * tNormalProjectorLeader * tCMLeaderElasticity->dTractiondDOF( tDofType, tNormalLeader )    //
                                          + tCMLeaderElasticity->dTestTractiondDOF( tDofType, tNormalLeader, tNormalProjectorLeader * tTractionLeader, mResidualDofType( 0 ) ) );
                    }
                }

                // if dependency of stabilization parameters on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Leader_Follower::LEADER ) )
                {
                    // get the derivatives of the SPs
                    const Matrix< DDRMat > tNitscheDer      = tSPNitsche->dSPdLeaderDOF( tDofType ).get_row( 0 );
                    const Matrix< DDRMat > tLeaderWeightDer = tSPNitsche->dSPdLeaderDOF( tDofType ).get_row( 1 );
                    const Matrix< DDRMat > tFollowerWeightDer  = tSPNitsche->dSPdLeaderDOF( tDofType ).get_row( 2 );

                    if ( tIfcPressureLeader - tNitsche * tNormalJumpLeader < 0 )
                    {
                        // add contribution to Jacobian
                        tJacMM += aWStar * (                                                                                                                                          //
                                          ( -tFILeader->N_trans() * tNormalProjectorLeader * tTractionLeader                                                                          //
                                                  + mBeta * tCMLeaderElasticity->testTraction_trans( tNormalLeader, mResidualDofType( 0 ) ) * tNormalProjectorLeader * tJumpLeader    //
                                                  + tNitsche * tFILeader->N_trans() * tNormalProjectorLeader * tJumpLeader )                                                          //
                                                  * tLeaderWeightDer                                                                                                                  //
                                          + tLeaderWeight * tFILeader->N_trans() * tNormalProjectorLeader * tJumpLeader * tNitscheDer );

                        tJacSM += aWStar * (                                                                                                            //
                                          ( +tFIFollower->N_trans() * tPrjFollowerToLeaderTrans * tNormalProjectorLeader * tTractionLeader                    //
                                                  - tNitsche * tFIFollower->N_trans() * tPrjFollowerToLeaderTrans * tNormalProjectorLeader * tJumpLeader )    //
                                                  * tLeaderWeightDer                                                                                    //
                                          - tLeaderWeight * tFIFollower->N_trans() * tPrjFollowerToLeader * tNormalProjectorLeader * tJumpLeader * tNitscheDer );
                    }
                    else
                    {
                        tJacMM += aWStar * (    //
                                          -mBeta * tCMLeaderElasticity->testTraction_trans( tNormalLeader, mResidualDofType( 0 ) ) * tNormalProjectorLeader * tTractionLeader )
                                * ( 1.0 / tNitsche * tLeaderWeightDer - tLeaderWeight / tNitsche / tNitsche * tNitscheDer );
                    }

                    if ( tIfcPressureFollower - tNitsche * tNormalJumpFollower < 0 )
                    {
                        tJacMM += aWStar * (                                                                                                           //
                                          ( +tFILeader->N_trans() * tPrjLeaderToFollowerTrans * tNormalProjectorFollower * tTractionFollower                    //
                                                  - tNitsche * tFILeader->N_trans() * tPrjLeaderToFollowerTrans * tNormalProjectorFollower * tJumpFollower )    //
                                                  * tFollowerWeightDer                                                                                    //
                                          - tFollowerWeight * tFILeader->N_trans() * tPrjLeaderToFollowerTrans * tNormalProjectorFollower * tJumpFollower * tNitscheDer );

                        tJacSM += aWStar * (                                                                                                                                      //
                                          ( -tFIFollower->N_trans() * tNormalProjectorFollower * tTractionFollower                                                                         //
                                                  + mBeta * tCMFollowerElasticity->testTraction_trans( tNormalFollower, mResidualDofType( 0 ) ) * tNormalProjectorFollower * tJumpFollower    //
                                                  + tNitsche * tFIFollower->N_trans() * tNormalProjectorFollower * tJumpFollower )                                                         //
                                                  * tFollowerWeightDer                                                                                                               //
                                          + tFollowerWeight * tFIFollower->N_trans() * tNormalProjectorFollower * tJumpFollower * tNitscheDer );
                    }
                    else
                    {
                        tJacSM += aWStar * (    //
                                          -mBeta * tCMFollowerElasticity->testTraction_trans( tNormalFollower, mResidualDofType( 0 ) ) * tNormalProjectorFollower * tTractionFollower )
                                * ( 1.0 / tNitsche * tFollowerWeightDer - tFollowerWeight / tNitsche / tNitsche * tNitscheDer );
                    }
                }
            }

            // compute the Jacobian for indirect dof dependencies through follower constitutive models
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
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { tFollowerDepStartIndex, tFollowerDepStopIndex } );

                auto tJacSS = mSet->get_jacobian()(
                        { tFollowerResStartIndex, tFollowerResStopIndex },
                        { tFollowerDepStartIndex, tFollowerDepStopIndex } );

                // if dof type is residual dof type
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    if ( tIfcPressureLeader - tNitsche * tNormalJumpLeader < 0 )
                    {
                        tJacMS += aWStar * tLeaderWeight * (                                                                                                                                       //
                                          -mBeta * tCMLeaderElasticity->testTraction_trans( tNormalLeader, mResidualDofType( 0 ) ) * tNormalProjectorLeader * tPrjFollowerToLeader * tFIFollower->N()    //
                                          - tNitsche * tFILeader->N_trans() * tNormalProjectorLeader * tPrjFollowerToLeader * tFIFollower->N() );

                        tJacSS += aWStar * tLeaderWeight * (    //
                                          +tNitsche * tFIFollower->N_trans() * tPrjFollowerToLeaderTrans * tNormalProjectorLeader * tPrjFollowerToLeader * tFIFollower->N() );
                    }

                    if ( tIfcPressureFollower - tNitsche * tNormalJumpFollower < 0 )
                    {
                        tJacMS += aWStar * tFollowerWeight * (    //
                                          -tNitsche * tFILeader->N_trans() * tPrjLeaderToFollowerTrans * tNormalProjectorFollower * tFIFollower->N() );

                        tJacSS += aWStar * tFollowerWeight * (                                                                                                                 //
                                          +mBeta * tCMFollowerElasticity->testTraction_trans( tNormalFollower, mResidualDofType( 0 ) ) * tNormalProjectorFollower * tFIFollower->N()    //
                                          + tNitsche * tFIFollower->N_trans() * tNormalProjectorFollower * tFIFollower->N() );
                    }
                }

                // if dependency on the dof type
                if ( tCMFollowerElasticity->check_dof_dependency( tDofType ) )
                {
                    if ( tIfcPressureFollower - tNitsche * tNormalJumpFollower < 0 )
                    {
                        // add contribution to Jacobian
                        tJacMS += aWStar * tFollowerWeight * (    //
                                          +tFILeader->N_trans() * tPrjLeaderToFollowerTrans * tNormalProjectorFollower * tCMFollowerElasticity->dTractiondDOF( tDofType, tNormalFollower ) );

                        tJacSS += aWStar * tFollowerWeight * (                                                                                             //
                                          -tFIFollower->N_trans() * tNormalProjectorFollower * tCMFollowerElasticity->dTractiondDOF( tDofType, tNormalFollower )    //
                                          + mBeta * tCMFollowerElasticity->dTestTractiondDOF( tDofType, tNormalFollower, tNormalProjectorFollower * tJumpFollower, mResidualDofType( 0 ) ) );
                    }
                    else
                    {
                        tJacSS += aWStar * tFollowerWeight * -mBeta / tNitsche * (                                                                                                                                  //
                                          tCMFollowerElasticity->testTraction_trans( tNormalFollower, mResidualDofType( 0 ) ) * tNormalProjectorFollower * tCMFollowerElasticity->dTractiondDOF( tDofType, tNormalFollower )    //
                                          + tCMFollowerElasticity->dTestTractiondDOF( tDofType, tNormalFollower, tNormalProjectorFollower * tTractionFollower, mResidualDofType( 0 ) ) );
                    }
                }

                // if dependency of stabilization parameters on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Leader_Follower::FOLLOWER ) )
                {
                    // get the derivatives of the SPs
                    const Matrix< DDRMat > tNitscheDer      = tSPNitsche->dSPdFollowerDOF( tDofType ).get_row( 0 );
                    const Matrix< DDRMat > tLeaderWeightDer = tSPNitsche->dSPdFollowerDOF( tDofType ).get_row( 1 );
                    const Matrix< DDRMat > tFollowerWeightDer  = tSPNitsche->dSPdFollowerDOF( tDofType ).get_row( 2 );

                    if ( tIfcPressureLeader - tNitsche * tNormalJumpLeader < 0 )
                    {
                        // add contribution to Jacobian
                        tJacMS += aWStar * (                                                                                                                                          //
                                          ( -tFILeader->N_trans() * tNormalProjectorLeader * tTractionLeader                                                                          //
                                                  + mBeta * tCMLeaderElasticity->testTraction_trans( tNormalLeader, mResidualDofType( 0 ) ) * tNormalProjectorLeader * tJumpLeader    //
                                                  + tNitsche * tFILeader->N_trans() * tNormalProjectorLeader * tJumpLeader )
                                                  * tLeaderWeightDer    //
                                          + tLeaderWeight * tFILeader->N_trans() * tNormalProjectorLeader * tJumpLeader * tNitscheDer );

                        tJacSS += aWStar * (                                                                                                            //
                                          ( +tFIFollower->N_trans() * tPrjFollowerToLeaderTrans * tNormalProjectorLeader * tTractionLeader                    //
                                                  - tNitsche * tFIFollower->N_trans() * tPrjFollowerToLeaderTrans * tNormalProjectorLeader * tJumpLeader )    //
                                                  * tLeaderWeightDer                                                                                    //
                                          - tLeaderWeight * tFIFollower->N_trans() * tPrjFollowerToLeaderTrans * tNormalProjectorLeader * tJumpLeader * tNitscheDer );
                    }
                    else
                    {
                        tJacMS += aWStar * (                                                                                                                                       //
                                          -mBeta * tCMLeaderElasticity->testTraction_trans( tNormalLeader, mResidualDofType( 0 ) ) * tNormalProjectorLeader * tTractionLeader )    //
                                * ( 1.0 / tNitsche * tLeaderWeightDer - tLeaderWeight / tNitsche / tNitsche * tNitscheDer );
                    }

                    if ( tIfcPressureFollower - tNitsche * tNormalJumpFollower < 0 )
                    {
                        tJacMS += aWStar * (                                                                                                           //
                                          ( +tFILeader->N_trans() * tPrjLeaderToFollowerTrans * tNormalProjectorFollower * tTractionFollower                    //
                                                  - tNitsche * tFILeader->N_trans() * tPrjLeaderToFollowerTrans * tNormalProjectorFollower * tJumpFollower )    //
                                                  * tFollowerWeightDer                                                                                    //
                                          - tFollowerWeight * tFILeader->N_trans() * tPrjLeaderToFollowerTrans * tNormalProjectorFollower * tJumpFollower * tNitscheDer );

                        tJacSS += aWStar * (                                                                                                                                      //
                                          ( -tFIFollower->N_trans() * tNormalProjectorFollower * tTractionFollower                                                                         //
                                                  + mBeta * tCMFollowerElasticity->testTraction_trans( tNormalFollower, mResidualDofType( 0 ) ) * tNormalProjectorFollower * tJumpFollower    //
                                                  + tNitsche * tFIFollower->N_trans() * tNormalProjectorFollower * tJumpFollower )                                                         //
                                                  * tFollowerWeightDer                                                                                                               //
                                          + tFIFollower->N_trans() * tNormalProjectorFollower * tJumpFollower * tNitscheDer );
                    }
                    else
                    {
                        tJacSS += aWStar * (    //
                                          -mBeta * tCMFollowerElasticity->testTraction_trans( tNormalFollower, mResidualDofType( 0 ) ) * tNormalProjectorFollower * tTractionFollower )
                                * ( 1.0 / tNitsche * tFollowerWeightDer - tFollowerWeight / tNitsche / tNitsche * tNitscheDer );
                    }
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Isotropic_Struc_Linear_Interface::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Isotropic_Struc_Linear_Contact_Gap_Nitsche::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Contact_Gap_Nitsche::compute_jacobian_and_residual - This function does nothing." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Isotropic_Struc_Linear_Contact_Gap_Nitsche::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Contact_Gap_Nitsche::compute_dRdp - This function does nothing." );
        }

        //------------------------------------------------------------------------------
}    // namespace moris::fem
