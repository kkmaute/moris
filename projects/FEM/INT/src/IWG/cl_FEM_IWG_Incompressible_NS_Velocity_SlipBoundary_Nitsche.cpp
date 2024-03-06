/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche.hpp"
//LINALG/src
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche::IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche( sint aBeta )
                {
            // set sign for symmetric/unsymmetric Nitsche
            mBeta = aBeta;
            init_property("Dirichlet", IWG_Property_Type::DIRICHLET);
            init_property("SlipLength", IWG_Property_Type::SLIPLENGTH);
            init_property("Traction", IWG_Property_Type::TRACTION);
            init_constitutive_model("IncompressibleFluid", IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE);
            init_stabilization_parameter("DirichletNitsche", IWG_Stabilization_Type::VELOCITY_SPLIPLENGTH_NITSCHE);
                }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche::compute_residual( real aWStar )
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
            Field_Interpolator * tFIVelocity = get_leader_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the imposed velocity property
            const std::shared_ptr< Property > & tPropVelocity = get_leader_property(IWG_Property_Type::DIRICHLET);

            // get the slip length property
            const std::shared_ptr< Property > & tPropSlipLength = get_leader_property(IWG_Property_Type::SLIPLENGTH);

            // get the traction property
            const std::shared_ptr< Property > & tPropTraction = get_leader_property(IWG_Property_Type::TRACTION);

            // get the fluid constitutive model
            const std::shared_ptr< Constitutive_Model > & tCMFluid = get_leader_constitutive_model(IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE);

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPNitsche = get_stabilization_parameter(IWG_Stabilization_Type::VELOCITY_SPLIPLENGTH_NITSCHE);

            // get the dynamic viscosity property
            const std::shared_ptr< Property > & tPropViscosity = tCMFluid->get_property( "Viscosity" );

            // check that slip length is defined
            MORIS_ASSERT( tPropSlipLength,
                    "IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche::compute_residual - Slip length not defined.\n");

            // check that prescribed velocity parameter is defined
            MORIS_ASSERT( tPropVelocity,
                    "IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::compute_residual - Prescribed velocity not defined.\n");

            // get the dynamic viscosity value
            const real tViscosity = tPropViscosity->val()( 0 );

            // get the slip length
            const real tSplipLength = tPropSlipLength->val()( 0 );

            // get spatial dimension
            uint tSpaceDim = tFIVelocity->get_number_of_fields();

            // build projector
            const Matrix<DDRMat> tNormalProjector  = get_normal() * trans( get_normal() );
            const Matrix<DDRMat> tTangentProjector = moris::eye(tSpaceDim,tSpaceDim) - tNormalProjector;

            // compute velocity jump in normal direction
            const Matrix<DDRMat> tNormalVelocityJump = tNormalProjector * ( tFIVelocity->val() - tPropVelocity->val() );

            // slip condition violation
            Matrix<DDRMat> tSlipVelocityJump = tTangentProjector * ( tSplipLength * tCMFluid->traction( get_normal() )
                    + tViscosity * ( tFIVelocity->val() - tPropVelocity->val() ) );

            // add contribution of prescribed traction to slip condition violation
            if ( tPropTraction )
            {
                tSlipVelocityJump -= tSplipLength * tTangentProjector * tPropTraction->val();
            }

            // penalty parameters
            const real tNormalPenalty   = tSPNitsche->val()( 0 );
            const real tTangentPenalty1 = tSPNitsche->val()( 1 );
            const real tTangentPenalty2 = tSPNitsche->val()( 2 );

            // compute leader residual
            mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex } ) += aWStar * (
                            + tFIVelocity->N_trans() * (
                                    - tCMFluid->traction( get_normal() )
                                    + tNormalPenalty   * tNormalVelocityJump
                                    + tTangentPenalty1 * tSlipVelocityJump )
                            - mBeta * trans( tCMFluid->testTraction( get_normal(), mResidualDofType( 0 ) ) ) * (
                                    + tNormalVelocityJump
                                    + tTangentPenalty2 * tSlipVelocityJump ) );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get the leader field interpolator for residual dof type
            Field_Interpolator * tFIVelocity = get_leader_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the imposed velocity property
            const std::shared_ptr< Property > & tPropVelocity = get_leader_property(IWG_Property_Type::DIRICHLET);

            // get the slip length property
            const std::shared_ptr< Property > & tPropSlipLength = get_leader_property(IWG_Property_Type::SLIPLENGTH);

            // get the traction property
            const std::shared_ptr< Property > & tPropTraction = get_leader_property(IWG_Property_Type::TRACTION);

            // get the fluid constitutive model
            const std::shared_ptr< Constitutive_Model > & tCMFluid = get_leader_constitutive_model(IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE);

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPNitsche = get_stabilization_parameter(IWG_Stabilization_Type::VELOCITY_SPLIPLENGTH_NITSCHE);

            // get the dynamic viscosity property
            const std::shared_ptr< Property > & tPropViscosity = tCMFluid->get_property( "Viscosity" );

            // get the dynamic viscosity value
            const real tViscosity = tPropViscosity->val()( 0 );

            // get the slip length
            const real tSplipLength = tPropSlipLength->val()( 0 );

            // get spatial dimension
            uint tSpaceDim = tFIVelocity->get_number_of_fields();

            // build projector
            const Matrix<DDRMat> tNormalProjector  = get_normal() * trans( get_normal() );
            const Matrix<DDRMat> tTangentProjector = moris::eye(tSpaceDim,tSpaceDim) - tNormalProjector;

            // compute velocity jump in normal direction
            const Matrix<DDRMat> tNormalVelocityJump = tNormalProjector * ( tFIVelocity->val() - tPropVelocity->val() );

            // slip condition violation
            Matrix<DDRMat> tSlipVelocityJump = tTangentProjector * ( tSplipLength * tCMFluid->traction( get_normal() )
                    + tViscosity * ( tFIVelocity->val() - tPropVelocity->val() ) );

            // add contribution of prescribed traction to slip condition violation
            if ( tPropTraction )
            {
                tSlipVelocityJump -= tSplipLength * tTangentProjector * tPropTraction->val();
            }

            // penalty parameters
            const real tNormalPenalty   = tSPNitsche->val()( 0 );
            const real tTangentPenalty1 = tSPNitsche->val()( 1 );
            const real tTangentPenalty2 = tSPNitsche->val()( 2 );

            // get number of leader dependencies
            const uint tLeaderNumDofDependencies = get_requested_leader_dof_types().size();

            // compute the Jacobian for indirect dof dependencies through leader
            for( uint iDOF = 0; iDOF < tLeaderNumDofDependencies; iDOF++ )
            {
                // get the dof type
                const Vector< MSI::Dof_Type > & tDofType = get_requested_leader_dof_types()( iDOF );

                // get the index for the dof type
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // if dof type is residual dof type
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // compute Jacobian direct dependencies
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * (
                                    + tFIVelocity->N_trans() * (
                                            + tNormalPenalty   * tNormalProjector
                                            + tTangentPenalty1 * tTangentProjector * tViscosity ) * tFIVelocity->N()
                                    - mBeta * trans( tCMFluid->testTraction( get_normal(), mResidualDofType( 0 ) ) ) * (
                                            + tNormalProjector
                                            + tTangentPenalty2 * tTangentProjector * tViscosity  ) * tFIVelocity->N() );
                }

                // if imposed velocity depends on dof type
                if ( tPropVelocity->check_dof_dependency( tDofType ) )
                {
                    // add contribution from property to Jacobian
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -= aWStar * (
                                    + tFIVelocity->N_trans() * (
                                            + tNormalPenalty   * tNormalProjector
                                            + tTangentPenalty1 * tTangentProjector * tViscosity ) * tPropVelocity->dPropdDOF( tDofType )
                                    - mBeta * trans( tCMFluid->testTraction( get_normal(), mResidualDofType( 0 ) ) ) * (
                                            + tNormalProjector
                                            + tTangentPenalty2 * tTangentProjector * tViscosity  ) * tPropVelocity->dPropdDOF( tDofType ) );
                }

                // if fluid constitutive model depends on dof type
                if ( tCMFluid->check_dof_dependency( tDofType ) )
                {
                    // add contribution of CM to Jacobian
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * (
                                    - tFIVelocity->N_trans() * tCMFluid->dTractiondDOF( tDofType, get_normal() )
                                    - mBeta * tCMFluid->dTestTractiondDOF(
                                            tDofType,
                                            get_normal(),
                                            tNormalVelocityJump + tTangentPenalty2 * tSlipVelocityJump,
                                            mResidualDofType( 0 ) ) );

                    // add contribution due to dependency of SlipVelocityJump on velocity and pressure
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * tSplipLength * ( (
                                    + tTangentPenalty1 * tFIVelocity->N_trans()
                                    - tTangentPenalty2 * mBeta * trans( tCMFluid->testTraction( get_normal(), mResidualDofType( 0 ) ) ) ) *
                                    tTangentProjector * tCMFluid->dTractiondDOF( tDofType, get_normal() ) );
                }

                // if prescribed traction depends on the dof type
                if ( tPropTraction )
                {
                    if ( tPropTraction->check_dof_dependency( tDofType ) )
                    {
                        MORIS_ERROR(false,"IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche::compute_jacobian - %s.\n",
                                "Dof dependency of prescribed traction not implemented yet");
                    }
                }

                // if viscosity depends on the dof type
                if ( tPropViscosity->check_dof_dependency( tDofType ) )
                {
                    MORIS_ERROR(false,"IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche::compute_jacobian - %s.\n",
                            "Dof dependency of viscosity not implemented yet");
                }

                // if stabilization parameter depends on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType ) )
                {
                    // get derivative of penalty parameter
                    const Matrix<DDRMat> & tDSPNitsche = tSPNitsche->dSPdLeaderDOF( tDofType );

                    // add contribution of SP to Jacobian
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * (
                                    + tFIVelocity->N_trans() * (
                                            + tNormalVelocityJump * tDSPNitsche.get_row( 0 ) ) );

                    // add the following lines if Nitsche penalty for tangential direction depends on dofs
                    //                        + tSlipVelocityJump   * tDSPNitsche.get_row( 1 ) )
                    //                - mBeta * trans( tCMFluid->testTraction( get_normal(), mResidualDofType( 0 ) ) ) *
                    //                        tSlipVelocityJump * tDSPNitsche.get_row( 2 ) );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche::compute_jacobian_and_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche::compute_dRdp( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

