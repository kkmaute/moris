/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Nitsche.cpp
 *
 */

#include "cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Nitsche.hpp"
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

        IWG_Isotropic_Struc_Linear_Contact_Nitsche::IWG_Isotropic_Struc_Linear_Contact_Nitsche( sint aBeta )
        {
            // sign for symmetric/unsymmetric Nitsche
            mBeta = aBeta;
            init_property("Thickness", IWG_Property_Type::THICKNESS);
            init_constitutive_model("ElastLinIso", IWG_Constitutive_Type::ELAST_LIN_ISO);
            init_stabilization_parameter("NitscheInterface", IWG_Stabilization_Type::NITSCHE_INTERFACE);
        }

        //------------------------------------------------------------------------------

        void
        IWG_Isotropic_Struc_Linear_Contact_Nitsche::compute_residual( real aWStar )
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
            Field_Interpolator* tFILeader = get_leader_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get follower field interpolator for the residual dof type
            Field_Interpolator* tFIFollower = get_follower_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the elasticity constitutive model
            const std::shared_ptr< Constitutive_Model >& tCMLeaderElasticity = get_leader_constitutive_model(IWG_Constitutive_Type::ELAST_LIN_ISO);

            const std::shared_ptr< Constitutive_Model > &tCMFollowerElasticity = get_follower_constitutive_model(IWG_Constitutive_Type::ELAST_LIN_ISO);

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche = get_stabilization_parameter(IWG_Stabilization_Type::NITSCHE_INTERFACE);

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness = get_leader_property(IWG_Property_Type::THICKNESS);

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            const real tNitsche      = tSPNitsche->val()( 0 );
            const real tLeaderWeight = tSPNitsche->val()( 1 );
            const real tFollowerWeight  = tSPNitsche->val()( 2 );

            // normal projection operator
            const Matrix< DDRMat > tNormalProjector = get_normal() * trans( get_normal() );

            // compute the jump
            const Matrix< DDRMat > tJump = tFILeader->val() - tFIFollower->val();

            // compute projection of displacement jump onto normal
            const real tNormalJump = dot( tJump, get_normal() );

            // evaluate average traction
            const Matrix< DDRMat > tTraction =
                    tLeaderWeight * tCMLeaderElasticity->traction( get_normal() )    //
                    + tFollowerWeight * tCMFollowerElasticity->traction( get_normal() );

            // compute contact pressure
            const real tIfcPressure = dot( tTraction, get_normal() );

            // check for contact
            if ( tIfcPressure - tNitsche * tNormalJump < 0 )
            {
                // compute leader residual
                mSet->get_residual()( 0 )(
                        { tLeaderResStartIndex, tLeaderResStopIndex } ) +=                                                                                        //
                        aWStar * (                                                                                                                                //
                                -tFILeader->N_trans() * tNormalProjector * tTraction                                                                              //
                                + mBeta * tLeaderWeight * tCMLeaderElasticity->testTraction_trans( get_normal(), mResidualDofType( 0 ) ) * tNormalProjector * tJump    //
                                + tNitsche * tFILeader->N_trans() * tNormalProjector * tJump );

                // compute follower residual
                mSet->get_residual()( 0 )(
                        { tFollowerResStartIndex, tFollowerResStopIndex } ) +=                                                                                        //
                        aWStar * (                                                                                                                              //
                                +tFIFollower->N_trans() * tNormalProjector * tTraction                                                                             //
                                + mBeta * tFollowerWeight * tCMFollowerElasticity->testTraction_trans( get_normal(), mResidualDofType( 0 ) ) * tNormalProjector * tJump    //
                                - tNitsche * tFIFollower->N_trans() * tNormalProjector * tJump );
            }
            else
            {
                mSet->get_residual()( 0 )(
                        { tLeaderResStartIndex, tLeaderResStopIndex } ) +=    //
                        aWStar * (                                            //
                                -mBeta / tNitsche * tLeaderWeight * tCMLeaderElasticity->testTraction_trans( get_normal(), mResidualDofType( 0 ) ) * tNormalProjector * tTraction );

                mSet->get_residual()( 0 )(
                        { tFollowerResStartIndex, tFollowerResStopIndex } ) +=    //
                        aWStar * (                                          //
                                -mBeta / tNitsche * tFollowerWeight * tCMFollowerElasticity->testTraction_trans( get_normal(), mResidualDofType( 0 ) ) * tNormalProjector * tTraction );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Isotropic_Struc_Linear_Contact_Nitsche::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Isotropic_Struc_Linear_Contact_Nitsche::compute_jacobian( real aWStar )
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
            Field_Interpolator* tFILeader = get_leader_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get follower field interpolator for the residual dof type
            Field_Interpolator* tFIFollower = get_follower_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the elasticity constitutive model
            const std::shared_ptr< Constitutive_Model >& tCMLeaderElasticity = get_leader_constitutive_model(IWG_Constitutive_Type::ELAST_LIN_ISO);

            const std::shared_ptr< Constitutive_Model > &tCMFollowerElasticity = get_follower_constitutive_model(IWG_Constitutive_Type::ELAST_LIN_ISO);

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche = get_stabilization_parameter(IWG_Stabilization_Type::NITSCHE_INTERFACE);

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness = get_leader_property(IWG_Property_Type::THICKNESS);

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            const real tNitsche      = tSPNitsche->val()( 0 );
            const real tLeaderWeight = tSPNitsche->val()( 1 );
            const real tFollowerWeight  = tSPNitsche->val()( 2 );

            // normal projection operator
            const Matrix< DDRMat > tNormalProjector = get_normal() * trans( get_normal() );

            // compute the jump
            const Matrix< DDRMat > tJump = tFILeader->val() - tFIFollower->val();

            // compute projection of displacement jump onto normal
            const real tNormalJump = dot( tJump, get_normal() );

            // evaluate average traction
            const Matrix< DDRMat > tTraction =
                    tLeaderWeight * tCMLeaderElasticity->traction( get_normal() ) + tFollowerWeight * tCMFollowerElasticity->traction( get_normal() );

            // compute contact pressure
            const real tIfcPressure = dot( tTraction, get_normal() );

            // get number of leader dof dependencies
            const uint tLeaderNumDofDependencies = get_requested_leader_dof_types().size();

            // compute the Jacobian for indirect dof dependencies through leader constitutive models
            for ( uint iDOF = 0; iDOF < tLeaderNumDofDependencies; iDOF++ )
            {
                // get the dof type
                const Vector< MSI::Dof_Type >& tDofType = get_requested_leader_dof_types()( iDOF );

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
                    // check for contact
                    if ( tIfcPressure - tNitsche * tNormalJump < 0 )
                    {
                        tJacMM += aWStar * (                                                                                                                                        //
                                          +mBeta * tLeaderWeight * tCMLeaderElasticity->testTraction_trans( get_normal(), mResidualDofType( 0 ) ) * tNormalProjector * tFILeader->N()    //
                                          + tNitsche * tFILeader->N_trans() * tNormalProjector * tFILeader->N() );

                        tJacSM += aWStar * (                                                                                                                                      //
                                          +mBeta * tFollowerWeight * tCMFollowerElasticity->testTraction_trans( get_normal(), mResidualDofType( 0 ) ) * tNormalProjector * tFILeader->N()    //
                                          - tNitsche * tFIFollower->N_trans() * tNormalProjector * tFILeader->N() );
                    }
                }

                // if dependency on the dof type
                if ( tCMLeaderElasticity->check_dof_dependency( tDofType ) )
                {
                    if ( tIfcPressure - tNitsche * tNormalJump < 0 )
                    {
                        // add contribution to Jacobian
                        tJacMM += aWStar * (                                                                                                                    //
                                          -tLeaderWeight * tFILeader->N_trans() * tNormalProjector * tCMLeaderElasticity->dTractiondDOF( tDofType, get_normal() )    //
                                          + mBeta * tLeaderWeight * tCMLeaderElasticity->dTestTractiondDOF( tDofType, get_normal(), tNormalProjector * tJump, mResidualDofType( 0 ) ) );

                        tJacSM += aWStar * ( tLeaderWeight * tFIFollower->N_trans() * tNormalProjector * tCMLeaderElasticity->dTractiondDOF( tDofType, get_normal() ) );
                    }
                    else
                    {
                        tJacMM += aWStar * -mBeta / tNitsche * tLeaderWeight * (                                                                                                                                    //
                                          tLeaderWeight * tCMLeaderElasticity->testTraction_trans( get_normal(), mResidualDofType( 0 ) ) * tNormalProjector * tCMLeaderElasticity->dTractiondDOF( tDofType, get_normal() )    //
                                          + tCMLeaderElasticity->dTestTractiondDOF( tDofType, get_normal(), tNormalProjector * tTraction, mResidualDofType( 0 ) ) );

                        tJacSM += aWStar * -mBeta / tNitsche * tFollowerWeight * (    //
                                          tLeaderWeight * tCMFollowerElasticity->testTraction_trans( get_normal(), mResidualDofType( 0 ) ) * tNormalProjector * tCMLeaderElasticity->dTractiondDOF( tDofType, get_normal() ) );
                    }
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
                            tCMLeaderElasticity->traction( get_normal() ) * tLeaderWeightDer + tCMFollowerElasticity->traction( get_normal() ) * tFollowerWeightDer;

                    if ( tIfcPressure - tNitsche * tNormalJump < 0 )
                    {
                        // add contribution to Jacobian
                        tJacMM += aWStar * (                                                                                                                                   //
                                          -tFILeader->N_trans() * tNormalProjector * tTractionDer                                                                              //
                                          + mBeta * tCMLeaderElasticity->testTraction_trans( get_normal(), mResidualDofType( 0 ) ) * tNormalProjector * tJump * tLeaderWeightDer    //
                                          + tFILeader->N_trans() * tNormalProjector * tJump * tNitscheDer );

                        tJacSM += aWStar * (                                                                                                                                 //
                                          +tFIFollower->N_trans() * tNormalProjector * tTractionDer                                                                             //
                                          + mBeta * tCMFollowerElasticity->testTraction_trans( get_normal(), mResidualDofType( 0 ) ) * tNormalProjector * tJump * tFollowerWeightDer    //
                                          - tFIFollower->N_trans() * tNormalProjector * tJump * tNitscheDer );
                    }
                    else
                    {
                        MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Contact_Nitsche::compute_jacobian - state dependent stabilization not implemented." );
                    }
                }
            }

            // compute the Jacobian for indirect dof dependencies through follower constitutive models
            uint tFollowerNumDofDependencies = get_requested_follower_dof_types().size();
            for ( uint iDOF = 0; iDOF < tFollowerNumDofDependencies; iDOF++ )
            {
                // get dof type
                const Vector< MSI::Dof_Type >& tDofType = get_requested_follower_dof_types()( iDOF );

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
                    if ( tIfcPressure - tNitsche * tNormalJump < 0 )
                    {
                        tJacMS += aWStar * (                                                                                                                                       //
                                          -mBeta * tLeaderWeight * tCMLeaderElasticity->testTraction_trans( get_normal(), mResidualDofType( 0 ) ) * tNormalProjector * tFIFollower->N()    //
                                          - tNitsche * tFILeader->N_trans() * tNormalProjector * tFIFollower->N() );

                        tJacSS += aWStar * (                                                                                                                                     //
                                          -mBeta * tFollowerWeight * tCMFollowerElasticity->testTraction_trans( get_normal(), mResidualDofType( 0 ) ) * tNormalProjector * tFIFollower->N()    //
                                          + tNitsche * tFIFollower->N_trans() * tNormalProjector * tFIFollower->N() );
                    }
                }

                // if dependency on the dof type
                if ( tCMFollowerElasticity->check_dof_dependency( tDofType ) )
                {
                    if ( tIfcPressure - tNitsche * tNormalJump < 0 )
                    {
                        // add contribution to Jacobian
                        tJacMS += aWStar * ( -tFollowerWeight * tFILeader->N_trans() * tNormalProjector * tCMFollowerElasticity->dTractiondDOF( tDofType, get_normal() ) );

                        tJacSS += aWStar * (                                                                                                                 //
                                          +tFollowerWeight * tFIFollower->N_trans() * tNormalProjector * tCMFollowerElasticity->dTractiondDOF( tDofType, get_normal() )    //
                                          + mBeta * tFollowerWeight * tCMFollowerElasticity->dTestTractiondDOF( tDofType, get_normal(), tNormalProjector * tJump, mResidualDofType( 0 ) ) );
                    }
                    else
                    {
                        tJacMS += aWStar * -mBeta / tNitsche * tLeaderWeight * (    //
                                          tFollowerWeight * tCMLeaderElasticity->testTraction_trans( get_normal(), mResidualDofType( 0 ) ) * tNormalProjector * tCMFollowerElasticity->dTractiondDOF( tDofType, get_normal() ) );

                        tJacSS += aWStar * -mBeta / tNitsche * tFollowerWeight * (                                                                                                                                  //
                                          tFollowerWeight * tCMFollowerElasticity->testTraction_trans( get_normal(), mResidualDofType( 0 ) ) * tNormalProjector * tCMFollowerElasticity->dTractiondDOF( tDofType, get_normal() )    //
                                          + tCMFollowerElasticity->dTestTractiondDOF( tDofType, get_normal(), tNormalProjector * tTraction, mResidualDofType( 0 ) ) );
                    }
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
                            tCMLeaderElasticity->traction( get_normal() ) * tLeaderWeightDer + tCMFollowerElasticity->traction( get_normal() ) * tFollowerWeightDer;

                    if ( tIfcPressure - tNitsche * tNormalJump < 0 )
                    {
                        // add contribution to Jacobian
                        tJacMS += aWStar * (                                                                                                                                   //
                                          -tFILeader->N_trans() * tNormalProjector * tTractionDer                                                                              //
                                          + mBeta * tCMLeaderElasticity->testTraction_trans( get_normal(), mResidualDofType( 0 ) ) * tNormalProjector * tJump * tLeaderWeightDer    //
                                          + tFILeader->N_trans() * tNormalProjector * tJump * tNitscheDer );

                        tJacSS += aWStar * (                                                                                                                                 //
                                          tFIFollower->N_trans() * tTractionDer                                                                                                 //
                                          + mBeta * tCMFollowerElasticity->testTraction_trans( get_normal(), mResidualDofType( 0 ) ) * tNormalProjector * tJump * tFollowerWeightDer    //
                                          - tFIFollower->N_trans() * tNormalProjector * tJump * tNitscheDer );
                    }
                    else
                    {
                        MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Contact_Nitsche::compute_jacobian - state dependent stabilization not implemented." );
                    }
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Isotropic_Struc_Linear_Contact_Nitsche::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Isotropic_Struc_Linear_Contact_Nitsche::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Contact_Nitsche::compute_jacobian_and_residual - This function does nothing." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Isotropic_Struc_Linear_Contact_Nitsche::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Contact_Nitsche::compute_dRdp - This function does nothing." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

