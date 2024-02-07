/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche.hpp"
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

        IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche( sint aBeta )
        {
            // set sign for symmetric/unsymmetric Nitsche
            mBeta = aBeta;

            // set size for the property pointer cell
            mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Dirichlet" ] = static_cast< uint >( IWG_Property_Type::DIRICHLET );
            mPropertyMap[ "Select" ]    = static_cast< uint >( IWG_Property_Type::SELECT );
            mPropertyMap[ "Upwind" ]    = static_cast< uint >( IWG_Property_Type::UPWIND );

            // set size for the constitutive model pointer cell
            mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "IncompressibleFluid" ] = static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "DirichletNitsche" ] = static_cast< uint >( IWG_Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::compute_residual( real aWStar )
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
            Field_Interpolator* tFIVelocity =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the selection matrix property
            const std::shared_ptr< Property >& tPropSelect =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            // set a default selection matrix if needed
            Matrix< DDRMat > tM;
            if ( tPropSelect == nullptr )
            {
                // get number of fields which should equal spatial dimension
                const uint tSpaceDim = tFIVelocity->get_number_of_fields();

                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = tPropSelect->val();

                // skip computing residual if projection matrix is zero
                if ( norm( tM ) < MORIS_REAL_EPS ) return;
            }

            // get the imposed velocity property
            const std::shared_ptr< Property >& tPropVelocity =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) );

            // get the upwind property
            const std::shared_ptr< Property >& tPropUpwind =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::UPWIND ) );

            // get the fluid constitutive model
            const std::shared_ptr< Constitutive_Model >& tCMFluid =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE ) );

            // check that prescribed velocity parameter is defined
            MORIS_ASSERT( tPropVelocity,
                    "IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::compute_residual - Prescribed velocity not defined.\n" );

            // compute the jump
            const auto tVelocityJump = tFIVelocity->val() - tPropVelocity->val();

            // compute leader residual
            mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex }, { 0, 0 } ) +=                                                //
                    aWStar * (                                                                                                  //
                            -tFIVelocity->N_trans() * tM * tCMFluid->traction( mNormal )                                        //
                            - mBeta * trans( tCMFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tM * tVelocityJump    //
                            + tSPNitsche->val()( 0 ) * tFIVelocity->N_trans() * tM * tVelocityJump );

            // if upwind
            if ( tPropUpwind )
            {
                // get the density property
                const std::shared_ptr< Property >& tDensityProp = tCMFluid->get_property( "Density" );

                // get the density value
                const real tDensity = tDensityProp->val()( 0 );

                mSet->get_residual()( 0 )(
                        { tLeaderResStartIndex, tLeaderResStopIndex }, { 0, 0 } ) -=    //
                        aWStar * (                                                      //
                                tPropUpwind->val()( 0 ) * tDensity * tFIVelocity->N_trans() * dot( tFIVelocity->val(), mNormal ) * tM * tVelocityJump );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::compute_jacobian( real aWStar )
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
            Field_Interpolator* tFIVelocity =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the selection matrix property
            const std::shared_ptr< Property >& tPropSelect =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            // set a default selection matrix if needed
            Matrix< DDRMat > tM;
            if ( tPropSelect == nullptr )
            {
                // get number of fields which should equal spatial dimension
                const uint tSpaceDim = tFIVelocity->get_number_of_fields();

                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = tPropSelect->val();

                // skip computing residual if projection matrix is zero
                if ( norm( tM ) < MORIS_REAL_EPS ) return;
            }

            // get the imposed velocity property
            const std::shared_ptr< Property >& tPropVelocity =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) );

            // get the upwind property
            const std::shared_ptr< Property >& tPropUpwind =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::UPWIND ) );

            // get the fluid constitutive model
            const std::shared_ptr< Constitutive_Model >& tCMFluid =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE ) );

            // compute the jump
            const auto tVelocityJump = tFIVelocity->val() - tPropVelocity->val();

            // get number of leader dependencies
            uint tLeaderNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            // compute the jacobian for indirect dof dependencies through leader
            for ( uint iDOF = 0; iDOF < tLeaderNumDofDependencies; iDOF++ )
            {
                // get the dof type
                const Vector< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // if dof type is residual dof type
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // compute jacobian direct dependencies
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) +=                                                            //
                            aWStar * (                                                                                                    //
                                    -mBeta * trans( tCMFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tM * tFIVelocity->N()    //
                                    + tSPNitsche->val()( 0 ) * tFIVelocity->N_trans() * tM * tFIVelocity->N() );
                }

                // if imposed velocity depends on dof type
                if ( tPropVelocity->check_dof_dependency( tDofType ) )
                {
                    // add contribution from property to jacobian
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -=                                                                                //
                            aWStar * (                                                                                                                        //
                                    -mBeta * trans( tCMFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tM * tPropVelocity->dPropdDOF( tDofType )    //
                                    + tSPNitsche->val()( 0 ) * tFIVelocity->N_trans() * tM * tPropVelocity->dPropdDOF( tDofType ) );
                }

                // if fluid constitutive model depends on dof type
                if ( tCMFluid->check_dof_dependency( tDofType ) )
                {
                    // add contribution of CM to jacobian
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) +=                                     //
                            aWStar * (                                                                             //
                                    -tFIVelocity->N_trans() * tM * tCMFluid->dTractiondDOF( tDofType, mNormal )    //
                                    - mBeta * tCMFluid->dTestTractiondDOF( tDofType, mNormal, tM * tVelocityJump, mResidualDofType( 0 ) ) );
                }

                // if stabilization parameter depends on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType ) )
                {
                    // add contribution of SP to jacobian
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) +=    //
                            aWStar * (                                            //
                                    tFIVelocity->N_trans() * tM * tVelocityJump * tSPNitsche->dSPdLeaderDOF( tDofType ) );
                }

                // if upwind
                if ( tPropUpwind )
                {
                    // get the density property
                    const std::shared_ptr< Property >& tDensityProp = tCMFluid->get_property( "Density" );

                    // get the density value
                    const real tDensity = tDensityProp->val()( 0 );

                    // if dof type is residual dof type
                    if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                    {
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -=                            //
                                aWStar * (                                                                    //
                                        tPropUpwind->val()( 0 ) * tDensity * tFIVelocity->N_trans() * (       //
                                                dot( tFIVelocity->val(), mNormal ) * tM * tFIVelocity->N()    //
                                                + tM * tVelocityJump * ( trans( mNormal ) * tFIVelocity->N() ) ) );
                    }

                    // if imposed velocity depends on dof type
                    if ( tPropVelocity->check_dof_dependency( tDofType ) )
                    {
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) +=                       //
                                aWStar * (                                                               //
                                        tPropUpwind->val()( 0 ) * tDensity * tFIVelocity->N_trans() *    //
                                        dot( tFIVelocity->val(), mNormal ) * tM * tPropVelocity->dPropdDOF( tDofType ) );
                    }

                    // if upwind parameter depends on the dof type
                    if ( tPropUpwind->check_dof_dependency( tDofType ) )
                    {
                        // add contribution of SP to jacobian
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -=                                  //
                                aWStar * (                                                                          //
                                        tFIVelocity->N_trans() * tDensity * dot( tFIVelocity->val(), mNormal ) *    //
                                        tM * tVelocityJump * tPropUpwind->dPropdDOF( tDofType ) );
                    }

                    // if density depends on dof type
                    if ( tDensityProp->check_dof_dependency( tDofType ) )
                    {
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -=                                                 //
                                aWStar * (                                                                                         //
                                        tFIVelocity->N_trans() * tPropUpwind->val()( 0 ) * dot( tFIVelocity->val(), mNormal ) *    //
                                        tM * tVelocityJump * tDensityProp->dPropdDOF( tDofType ) );
                    }
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::compute_jacobian_and_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::compute_dRdp( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
