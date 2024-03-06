/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Spalart_Allmaras_Turbulence_Dirichlet.cpp
 *
 */

#include "cl_FEM_IWG_Spalart_Allmaras_Turbulence_Dirichlet.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
// LINALG/src
#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_norm.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Spalart_Allmaras_Turbulence_Dirichlet::IWG_Spalart_Allmaras_Turbulence_Dirichlet( sint aBeta )
        {
            // set mBeta for symmetric/unsymmetric Nitsche
            mBeta = aBeta;
            init_property("Dirichlet", IWG_Property_Type::DIRICHLET);
            init_property("Select", IWG_Property_Type::SELECT);
            init_property("Upwind", IWG_Property_Type::UPWIND);
            init_constitutive_model("SpalartAllmarasTurbulence", IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE);
            init_stabilization_parameter("Nitsche", IWG_Stabilization_Type::NITSCHE);
        }

        //------------------------------------------------------------------------------

        void
        IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_residual( real aWStar )
        {
            // check leader field interpolators
#ifdef MORIS_HAVE_DEBUG
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here velocity), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get the residual viscosity FI
            Field_Interpolator* tFIViscosity = get_leader_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the velocity dof FI
            Field_Interpolator* tFIVelocity = get_leader_fi_manager()->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the imposed viscosity property
            const std::shared_ptr< Property >& tPropDirichlet = get_leader_property(IWG_Property_Type::DIRICHLET);

            // get the SA turbulence CM
            const std::shared_ptr< Constitutive_Model >& tCMSATurbulence = get_leader_constitutive_model(IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE);

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche = get_stabilization_parameter(IWG_Stabilization_Type::NITSCHE);

            // get the selection matrix property
            const std::shared_ptr< Property >& tPropSelect = get_leader_property(IWG_Property_Type::SELECT);

            // set a default selection matrix if needed
            Matrix< DDRMat > tM;
            if ( tPropSelect == nullptr )
            {
                // get number of fields which should be one
                const uint tNumFields = tFIViscosity->get_number_of_fields();

                // set selection matrix as identity
                eye( tNumFields, tNumFields, tM );
            }
            else
            {
                tM = tPropSelect->val();
            }

            // compute the jump
            Matrix< DDRMat > tJump = tFIViscosity->val() - tPropDirichlet->val();

            // compute the residual weak form
            mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex }, { 0, 0 } ) +=                                      //
                    aWStar * (                                                                                        //
                            -tFIViscosity->N_trans() * tM * tCMSATurbulence->traction( get_normal() )                      //
                            - mBeta * tCMSATurbulence->testTraction( get_normal(), mResidualDofType( 0 ) ) * tM * tJump    //
                            + tSPNitsche->val()( 0 ) * tFIViscosity->N_trans() * tM * tJump );

            // get the upwind property
            const std::shared_ptr< Property >& tPropUpwind = get_leader_property(IWG_Property_Type::UPWIND);

            // upwind term
            if ( tPropUpwind )
            {
                // compute modified velocity
                Matrix< DDRMat > tModVelocity =
                        tFIVelocity->val() - mCb2 * tFIViscosity->gradx( 1 ) / mSigma;

                // add upwind contribution to residual
                mSet->get_residual()( 0 )(
                        { tLeaderResStartIndex, tLeaderResStopIndex }, { 0, 0 } ) -=    //
                        aWStar * (                                                      //
                                tPropUpwind->val()( 0 ) * tFIViscosity->N_trans() * dot( tModVelocity, get_normal() ) * tM * tJump );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here velocity), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get the residual dof FI (here viscosity)
            Field_Interpolator* tFIViscosity = get_leader_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the velocity dof FI
            Field_Interpolator* tFIVelocity = get_leader_fi_manager()->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the dirichlet property
            const std::shared_ptr< Property >& tPropDirichlet = get_leader_property(IWG_Property_Type::DIRICHLET);

            // get the SA turbulence CM
            const std::shared_ptr< Constitutive_Model >& tCMSATurbulence = get_leader_constitutive_model(IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE);

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche = get_stabilization_parameter(IWG_Stabilization_Type::NITSCHE);

            // get the selection matrix property
            const std::shared_ptr< Property >& tPropSelect = get_leader_property(IWG_Property_Type::SELECT);

            // set a default selection matrix if needed
            Matrix< DDRMat > tM;
            if ( tPropSelect == nullptr )
            {
                // get number of fields which should be one
                const uint tSpaceDim = tFIViscosity->get_number_of_fields();

                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = tPropSelect->val();
            }

            // compute the jump
            Matrix< DDRMat > tJump = tFIViscosity->val() - tPropDirichlet->val();

            // get number of dof dependencies
            uint tNumDofDependencies = get_requested_leader_dof_types().size();

            // loop over the dof dependencies
            for ( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Vector< MSI::Dof_Type > const &tDofType = get_requested_leader_dof_types()( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // if residual dof type (here viscosity)
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) +=                                                           //
                            aWStar * (                                                                                                   //
                                    -mBeta * tCMSATurbulence->testTraction( get_normal(), mResidualDofType( 0 ) ) * tM * tFIViscosity->N()    //
                                    + tSPNitsche->val()( 0 ) * tFIViscosity->N_trans() * tM * tFIViscosity->N() );
                }

                // if imposed viscosity depends on dof type
                if ( tPropDirichlet->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) +=                                                                               //
                            aWStar * (                                                                                                                       //
                                    +mBeta * tCMSATurbulence->testTraction( get_normal(), mResidualDofType( 0 ) ) * tM * tPropDirichlet->dPropdDOF( tDofType )    //
                                    - tSPNitsche->val()( 0 ) * tFIViscosity->N_trans() * tM * tPropDirichlet->dPropdDOF( tDofType ) );
                }

                // if Nitsche SP depends on dof type
                if ( tSPNitsche->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) +=    //
                            aWStar * (                                            //
                                    tFIViscosity->N_trans() * tM * tJump * tSPNitsche->dSPdLeaderDOF( tDofType ) );
                }

                // if turbulence CM depends on dof type
                if ( tCMSATurbulence->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) +=                                             //
                            aWStar * (                                                                                     //
                                    -tFIViscosity->N_trans() * tM * tCMSATurbulence->dTractiondDOF( tDofType, get_normal() )    //
                                    - mBeta * tM( 0 ) * tJump( 0 ) * tCMSATurbulence->dTestTractiondDOF( tDofType, get_normal(), mResidualDofType( 0 ) ) );
                }

                // get the upwind property
                const std::shared_ptr< Property >& tPropUpwind = get_leader_property(IWG_Property_Type::UPWIND);

                // upwind term
                if ( tPropUpwind )
                {
                    // compute modified velocity
                    Matrix< DDRMat > tModVelocity =
                            tFIVelocity->val() - mCb2 * tFIViscosity->gradx( 1 ) / mSigma;

                    // if dof type is residual dof type
                    if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                    {
                        // compute dModVelocitydModViscosity
                        Matrix< DDRMat > tModVelocityDer = -mCb2 * tFIViscosity->dnNdxn( 1 ) / mSigma;

                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -=                                                                   //
                                aWStar * (                                                                                                           //
                                        tPropUpwind->val()( 0 ) * tFIViscosity->N_trans() * dot( tModVelocity, get_normal() ) * tM * tFIViscosity->N()    //
                                        + tPropUpwind->val()( 0 ) * tFIViscosity->N_trans() * tM * tJump * trans( get_normal() ) * tModVelocityDer );
                    }

                    // if dof type is residual dof type
                    if ( tDofType( 0 ) == MSI::Dof_Type::VX )
                    {
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -=    //
                                aWStar * (                                            //
                                        tPropUpwind->val()( 0 ) * tFIViscosity->N_trans() * tM * tJump * trans( get_normal() ) * tFIVelocity->N() );
                    }

                    // if imposed velocity depends on dof type
                    if ( tPropDirichlet->check_dof_dependency( tDofType ) )
                    {
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) +=    //
                                aWStar * (                                            //
                                        tPropUpwind->val()( 0 ) * tFIViscosity->N_trans() * dot( tModVelocity, get_normal() ) * tM * tPropDirichlet->dPropdDOF( tDofType ) );
                    }

                    // if upwind parameter depends on the dof type
                    if ( tPropUpwind->check_dof_dependency( tDofType ) )
                    {
                        // add contribution of SP to jacobian
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -=    //
                                aWStar * (                                            //
                                        tFIViscosity->N_trans() * dot( tModVelocity, get_normal() ) * tM * tJump * tPropUpwind->dPropdDOF( tDofType ) );
                    }
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_jacobian_and_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_dRdp( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Spalart_Allmaras_Turbulence_Dirichlet::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

