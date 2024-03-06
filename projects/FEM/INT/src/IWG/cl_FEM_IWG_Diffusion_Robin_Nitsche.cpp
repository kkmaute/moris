/*
 * Copyright (c) 2022 University of Colorado
 *Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Diffusion_Robin_Nitsche.cpp
 *
 */

// FEM/INT/src
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Diffusion_Robin_Nitsche.hpp"
// LINALG/src
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IWG_Diffusion_Robin_Nitsche::IWG_Diffusion_Robin_Nitsche( sint aBeta )
        {
            // set sign for symmetric/unsymmetric Nitsche
            mBeta = aBeta;
            init_property("Dirichlet", IWG_Property_Type::DIRICHLET);
            init_property("NeumannPenalty", IWG_Property_Type::NEUMANN_PENALTY);
            init_property("Traction", IWG_Property_Type::TRACTION);
            init_constitutive_model("Diffusion", IWG_Constitutive_Type::DIFFUSION);
            init_stabilization_parameter("RobinNitsche", IWG_Stabilization_Type::ROBIN_NITSCHE);
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Robin_Nitsche::compute_residual( real aWStar )
        {
            // check leader field interpolators
#ifdef MORIS_HAVE_DEBUG
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get the leader field interpolator for the residual dof type
            Field_Interpolator* tFITemp = get_leader_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the imposed velocity property
            const std::shared_ptr< Property >& tPropDirichlet = get_leader_property(IWG_Property_Type::DIRICHLET);

            // get the slip length property
            const std::shared_ptr< Property >& tPropNeumannPen = get_leader_property(IWG_Property_Type::NEUMANN_PENALTY);

            // get the traction property
            const std::shared_ptr< Property >& tPropTraction = get_leader_property(IWG_Property_Type::TRACTION);

            // get the fluid constitutive model
            const std::shared_ptr< Constitutive_Model >& tCMDiffusion = get_leader_constitutive_model(IWG_Constitutive_Type::DIFFUSION);

            // get the dynamic viscosity property
            const std::shared_ptr< Property >& tPropConductivity = tCMDiffusion->get_property( "Conductivity" );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche = get_stabilization_parameter(IWG_Stabilization_Type::ROBIN_NITSCHE);

            // check that slip length is defined
            MORIS_ASSERT( tPropNeumannPen and tPropConductivity,
                    "IWG_Diffusion_Robin_Nitsche::compute_residual - Slip length not defined.\n" );

            // get the slip length
            const real tNeumannPen = tPropNeumannPen->val()( 0 );

            // compute the dirichlet jump
            Matrix< DDRMat > tJump = tPropConductivity->val() * tFITemp->val();
            if ( tPropDirichlet )
            {
                // subtract the prescribed dirichlet , by default is zero
                tJump -= tPropConductivity->val() * tPropDirichlet->val();
            }

            // compute the traction jump
            tJump += tNeumannPen * tCMDiffusion->traction( get_normal() );
            if ( tPropTraction )
            {
                // subtract the prescribed traction , by default is zero
                tJump -= tNeumannPen * tPropTraction->val();
            }

            // penalty parameters
            const real tStabilityPenalty = tSPNitsche->val()( 0 );
            const real tAdjointPenalty   = tSPNitsche->val()( 1 );

            // get sub-matrix
            auto tRes = mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex } );

            // compute leader residual
            tRes += aWStar * (                                                                            //
                            tFITemp->N_trans() * (                                                        //
                                    -tCMDiffusion->traction( get_normal() )                                    //                                     //
                                    + tStabilityPenalty * tJump )                                         //
                            - mBeta * tCMDiffusion->testTraction( get_normal(), mResidualDofType( 0 ) ) * (    //                                                 //
                                      tAdjointPenalty * tJump ) );                                        //

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Diffusion_Robin_Nitsche::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Robin_Nitsche::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get the leader field interpolator for the residual dof type
            Field_Interpolator* tFITemp = get_leader_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the imposed velocity property
            const std::shared_ptr< Property >& tPropDirichlet = get_leader_property(IWG_Property_Type::DIRICHLET);

            // get the slip length property
            const std::shared_ptr< Property >& tPropNeumannPen = get_leader_property(IWG_Property_Type::NEUMANN_PENALTY);

            // get the traction property
            const std::shared_ptr< Property >& tPropTraction = get_leader_property(IWG_Property_Type::TRACTION);

            // get the fluid constitutive model
            const std::shared_ptr< Constitutive_Model >& tCMDiffusion = get_leader_constitutive_model(IWG_Constitutive_Type::DIFFUSION);

            // get the dynamic viscosity property
            const std::shared_ptr< Property >& tPropConductivity = tCMDiffusion->get_property( "Conductivity" );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche = get_stabilization_parameter(IWG_Stabilization_Type::ROBIN_NITSCHE);

            // check that slip length is defined
            MORIS_ASSERT( tPropNeumannPen and tPropConductivity,
                    "IWG_Diffusion_Robin_Nitsche::compute_residual - Slip length not defined.\n" );

            // get the slip length
            const real tNeumannPen = tPropNeumannPen->val()( 0 );

            // compute the dirichlet jump
            Matrix< DDRMat > tJump = tPropConductivity->val() * tFITemp->val();
            if ( tPropDirichlet )
            {
                // subtract the prescribed dirichlet , by default is zero
                tJump -= tPropConductivity->val() * tPropDirichlet->val();
            }

            // compute the traction jump
            tJump += tNeumannPen * tCMDiffusion->traction( get_normal() );
            if ( tPropTraction )
            {
                // subtract the prescribed traction , by default is zero
                tJump -= tNeumannPen * tPropTraction->val();
            }

            // penalty parameters
            const real tStabilityPenalty = tSPNitsche->val()( 0 );
            const real tAdjointPenalty   = tSPNitsche->val()( 1 );

            // get number of leader dependencies
            const uint tLeaderNumDofDependencies = get_requested_leader_dof_types().size();

            // compute the Jacobian for indirect dof dependencies through leader
            for ( uint iDOF = 0; iDOF < tLeaderNumDofDependencies; iDOF++ )
            {
                // get the dof type
                const Vector< MSI::Dof_Type >& tDofType = get_requested_leader_dof_types()( iDOF );

                // get the index for the dof type
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // get sub-matrix
                auto tJac = mSet->get_jacobian()(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } );

                // if dof type is residual dof type
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // compute Jacobian direct dependencies
                    tJac += aWStar * (                                                                            //
                                    tFITemp->N_trans() * (                                                        //
                                            tStabilityPenalty * tPropConductivity->val() * tFITemp->N() )         //
                                    - mBeta * tCMDiffusion->testTraction( get_normal(), mResidualDofType( 0 ) ) * (    //
                                              tAdjointPenalty * tPropConductivity->val() * tFITemp->N() ) );      //
                }

                // if fluid constitutive model depends on dof type
                if ( tCMDiffusion->check_dof_dependency( tDofType ) )
                {
                    // compute Jacobian direct dependencies
                    tJac += aWStar * (                                                                                           //
                                    tFITemp->N_trans() * (                                                                       //
                                            -tCMDiffusion->dTractiondDOF( tDofType, get_normal() ) )                                  //
                                    - mBeta * tCMDiffusion->dTestTractiondDOF( tDofType, get_normal(), mResidualDofType( 0 ) ) * (    //
                                              tAdjointPenalty * tJump( 0 ) ) );

                    // compute the dependencies of the jacobian on the jump term which has traction in it
                    tJac += aWStar * (                                                                                              //
                                    tFITemp->N_trans() * (                                                                          //
                                            tStabilityPenalty * tNeumannPen * tCMDiffusion->dTractiondDOF( tDofType, get_normal() ) )    //
                                    - mBeta * tCMDiffusion->testTraction( get_normal(), mResidualDofType( 0 ) ) * (                      //
                                              tAdjointPenalty * tNeumannPen * tCMDiffusion->dTractiondDOF( tDofType, get_normal() ) ) );
                }

                // if prescribed traction depends on the dof type
                if ( tPropDirichlet )
                {
                    if ( tPropDirichlet->check_dof_dependency( tDofType ) )
                    {
                        MORIS_ERROR( false, "IWG_Diffusion_Robin_Nitsche::compute_jacobian - %s.\n", "Dof dependency of prescribed traction not implemented yet" );
                    }
                }

                // if prescribed traction depends on the dof type
                if ( tPropTraction )
                {
                    if ( tPropTraction->check_dof_dependency( tDofType ) )
                    {
                        MORIS_ERROR( false, "IWG_Diffusion_Robin_Nitsche::compute_jacobian - %s.\n", "Dof dependency of prescribed traction not implemented yet" );
                    }
                }

                if ( tPropNeumannPen->check_dof_dependency( tDofType ) )
                {
                    MORIS_ERROR( false, "IWG_Diffusion_Robin_Nitsche::compute_jacobian - %s.\n", "Dof dependency of prescribed traction not implemented yet" );
                }

                // if stabilization parameter depends on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType ) )
                {
                    MORIS_ERROR( false, "IWG_Diffusion_Robin_Nitsche::compute_jacobian - %s.\n", "Dof dependency of prescribed traction not implemented yet" );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Diffusion_Robin_Nitsche::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Robin_Nitsche::compute_jacobian_and_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Diffusion_Robin_Nitsche::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Robin_Nitsche::compute_dRdp( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Diffusion_Robin_Nitsche::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
