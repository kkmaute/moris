/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Advection_Bulk.cpp
 *
 */

#include "cl_FEM_IWG_Advection_Bulk.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_dot.hpp"
#include "fn_eye.hpp"
#include "fn_FEM_IWG_Crosswind_Stabilization_Tools.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IWG_Advection_Bulk::IWG_Advection_Bulk()
    {
        // set size for the property pointer cell
        mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "Load" ] = static_cast< uint >( IWG_Property_Type::BODY_LOAD );

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "Diffusion" ] = static_cast< uint >( IWG_Constitutive_Type::DIFFUSION );

        // set size for the stabilization parameter pointer cell
        mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

        // populate the stabilization map
        mStabilizationMap[ "SUPG" ]               = static_cast< uint >( IWG_Stabilization_Type::SUPG );
        mStabilizationMap[ "YZBeta" ]             = static_cast< uint >( IWG_Stabilization_Type::YZBETA );
        mStabilizationMap[ "DiffusionCrosswind" ] = static_cast< uint >( IWG_Stabilization_Type::DIFFUSION_CROSSWIND );
        mStabilizationMap[ "DiffusionIsotropic" ] = static_cast< uint >( IWG_Stabilization_Type::DIFFUSION_ISOTROPIC );
    }

    //------------------------------------------------------------------------------

    void IWG_Advection_Bulk::compute_residual( real aWStar )
    {
        // check leader field interpolators
#ifdef MORIS_HAVE_DEBUG
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type, indices for assembly
            const uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            const uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            const uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get the residual dof (here temperature) field interpolator
            Field_Interpolator* tFITemp =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the velocity dof field interpolator
            // FIXME protect dof type
            Field_Interpolator* tFIVelocity =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > & tCMDiffusion =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::DIFFUSION ) );

            // get the SUPG stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPSUPG =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SUPG ) );

            // get the YZBeta stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPYZBeta =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::YZBETA ) );

            // get the crosswind diffusion stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPCrosswind =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIFFUSION_CROSSWIND ) );

            // get the isotropic diffusion stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPIsotropic =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIFFUSION_ISOTROPIC ) );

            // compute the residual strong form if either SUPG or YZBeta is used
            Matrix< DDRMat > tRT;
            if ( tSPSUPG || tSPYZBeta || tSPCrosswind || tSPIsotropic )
            {
                this->compute_residual_strong_form( tRT );
            }

            // get sub-matrix of residual
            auto tRes = mSet->get_residual()( 0 )( { tLeaderResStartIndex, tLeaderResStopIndex } );

            // compute the residual
            tRes += aWStar *
                    tFITemp->N_trans() * tFIVelocity->val_trans() * tCMDiffusion->gradEnergy();

            // compute SUPG contribution to residual
            if ( tSPSUPG )
            {
                tRes += aWStar *
                        tSPSUPG->val()( 0 ) * trans( tFITemp->dnNdxn( 1 ) ) * tFIVelocity->val() * tRT( 0, 0 );
            }

            // compute YZBeta contribution to residual
            if ( tSPYZBeta )
            {
                tRes += aWStar *
                        tSPYZBeta->val()( 0 ) * std::abs(tRT( 0, 0 )) * trans( tFITemp->dnNdxn( 1 ) ) * tCMDiffusion->gradEnergy();
            }

            // if crosswind stabilization
            if( tSPCrosswind || tSPIsotropic )
            {
                // bool for crosswind or isotropic
                bool tIsCrosswind = tSPCrosswind != nullptr;

                // get conductivity property from CM
                const std::shared_ptr< Property >& tPropConductivity = //
                        tCMDiffusion->get_property( "Conductivity" );

                // get zero tolerance for isotropic or crosswind
                real tEpsilon;
                real tSPValue;
                if( tIsCrosswind )
                {
                    // get zero tolerance for crosswind
                    tEpsilon = tSPCrosswind->val()( 1 );

                    // get the value of SP
                    tSPValue = tSPCrosswind->val()( 0 );
                }
                else
                {
                    // get zero tolerance for isotropic
                    tEpsilon = tSPIsotropic->val()( 1 );

                    // get the value of SP
                    tSPValue = tSPIsotropic->val()( 0 );
                }

                // compute the norm of the viscosity gradient
                real tNorm = std::max( norm( tFITemp->gradx( 1 ) ), tEpsilon );

                // compute the abs of the strong form of the residual
                real tRAbs = std::max( std::abs( tRT( 0 ) ), tEpsilon );

                // compute full crosswind stabilization parameter value
                real tCrosswind = std::max( tSPValue * tRAbs / tNorm - tPropConductivity->val()( 0 ), 0.0 );

                // id crosswind stabilization parameter is greater than zero
                if( tCrosswind > 0.0 )
                {
                    // get the velocity FI
                    // FIXME protect dof type
                    Field_Interpolator* tFIVelocity = //
                            mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

                    // get space dimension
                    uint tSpaceDim = tFIVelocity->get_number_of_fields();

                    // compute crosswind projection of temperature gradient
                    // FIXME protect velocity dof type
                    Matrix< DDRMat > tcgradxt;
                    compute_cgradxw( { MSI::Dof_Type::VX }, mResidualDofType( 0 ), //
                            mLeaderFIManager, tSpaceDim, tEpsilon, tIsCrosswind,   //
                            tcgradxt );

                    // add contribution to residual
                    tRes += aWStar * trans( tFITemp->dnNdxn( 1 ) ) * tCrosswind * tcgradxt;
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Advection_Bulk::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Advection_Bulk::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get the residual dof (here temperature) field interpolator
            Field_Interpolator* tFITemp = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the velocity dof field interpolator
            // FIXME protect dof type
            Field_Interpolator* tFIVelocity =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > & tCMDiffusion =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::DIFFUSION ) );

            // get the SUPG stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPSUPG =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SUPG ) );

            // get the YZBeta stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPYZBeta =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::YZBETA ) );

            // get the crosswind diffusion stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPCrosswind =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIFFUSION_CROSSWIND ) );

            // get the isotropic diffusion stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPIsotropic =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIFFUSION_ISOTROPIC ) );

            // compute the residual strong form if either SUPG or YZBeta is used
            Matrix< DDRMat > tRT;
            if ( tSPSUPG || tSPYZBeta || tSPCrosswind || tSPIsotropic )
            {
                this->compute_residual_strong_form( tRT );
            }

            // get number of leader dof dependencies
            const  uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            // loop over leader dof dependencies
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                const Vector< MSI::Dof_Type > & tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                const sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                const uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                const uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // get sub-matrix
                auto tJac = mSet->get_jacobian()(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } );

                // if velocity dof type
                // FIXME protect dof type
                if( tDofType( 0 ) == MSI::Dof_Type::VX )
                {
                    // add contribution to Jacobian
                    tJac += aWStar *
                            tFITemp->N_trans() * trans( tCMDiffusion->gradEnergy() ) * tFIVelocity->N();

                    // consider SUPG term
                    if ( tSPSUPG )
                    {
                        tJac += aWStar *
                                tSPSUPG->val()( 0 ) * trans( tFITemp->dnNdxn( 1 ) ) * tFIVelocity->N() * tRT( 0, 0 );
                    }
                }

                // if diffusion CM depends on dof type
                if( tCMDiffusion->check_dof_dependency( tDofType ) )
                {
                    // add contribution to Jacobian
                    tJac += aWStar *
                            ( tFITemp->N_trans() * tFIVelocity->val_trans() * tCMDiffusion->dGradEnergydDOF( tDofType ) );
                }

                // compute the Jacobian strong form

                Matrix< DDRMat > tJT;
                if ( tSPSUPG || tSPYZBeta || tSPCrosswind || tSPIsotropic )
                {
                    this->compute_jacobian_strong_form( tDofType, tJT );
                }

                // compute the SUPG contribution
                if ( tSPSUPG )
                {
                    // contribution due to the dof dependence of strong form
                    tJac += aWStar * (
                            tSPSUPG->val()( 0 ) * trans( tFITemp->dnNdxn( 1 ) ) * tFIVelocity->val() * tJT );

                    // if SUPG stabilization parameter depends on dof type
                    if( tSPSUPG->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to Jacobian
                        tJac += aWStar * (
                                trans( tFITemp->dnNdxn( 1 ) ) * tFIVelocity->val() * tRT( 0, 0 ) * tSPSUPG->dSPdLeaderDOF( tDofType ) );
                    }
                }

                // compute the YZBeta contribution
                if ( tSPYZBeta )
                {
                    // contribution from strong form of residual
                    const real tSign = tRT( 0, 0 ) < 0 ? -1.0 : 1.0;

                    tJac += tSign * aWStar * (
                            tSPYZBeta->val()( 0 ) * trans( tFITemp->dnNdxn( 1 ) ) * tCMDiffusion->gradEnergy() * tJT );

                    // contribution from spatial gradient of energy
                    if( tCMDiffusion->check_dof_dependency( tDofType ) )
                    {
                        tJac += aWStar * std::abs( tRT( 0, 0 ) ) * (
                                tSPYZBeta->val()( 0 ) * trans( tFITemp->dnNdxn( 1 ) ) * tCMDiffusion->dGradEnergydDOF( tDofType ) );
                    }

                    // if YZBeta stabilization parameter depends on dof type
                    if( tSPYZBeta->check_dof_dependency( tDofType ) )
                    {
                        tJac += aWStar * std::abs( tRT( 0, 0 ) ) * (
                                trans( tFITemp->dnNdxn( 1 ) ) * tCMDiffusion->gradEnergy() * tSPYZBeta->dSPdLeaderDOF( tDofType ) );
                    }
                }

                // if isotropic or crosswind diffusion stabilization
                if( tSPCrosswind || tSPIsotropic )
                {
                    // bool for crosswind or isotropic
                    bool tIsCrosswind = tSPCrosswind != nullptr;

                    // get conductivity property from CM
                    const std::shared_ptr< Property >& tPropConductivity = //
                            tCMDiffusion->get_property( "Conductivity" );

                    // get zero tolerance for isotropic or crosswind
                    real tEpsilon;
                    real tSPValue;
                    if( tIsCrosswind )
                    {
                        // get zero tolerance for crosswind
                        tEpsilon = tSPCrosswind->val()( 1 );

                        // get the value of SP
                        tSPValue = tSPCrosswind->val()( 0 );
                    }
                    else
                    {
                        // get zero tolerance for isotropic
                        tEpsilon = tSPIsotropic->val()( 1 );

                        // get the value of SP
                        tSPValue = tSPIsotropic->val()( 0 );
                    }

                    // compute the norm of the gradient of the viscosity
                    real tNorm = std::max( norm( tFITemp->gradx( 1 ) ), tEpsilon );

                    // compute the norm of the gradient of the strong form of the residual
                    real tRAbs = std::max( std::abs( tRT( 0 ) ), tEpsilon );

                    // compute full crosswind stabilization parameter value
                    real tCrosswind = std::max( tSPValue * tRAbs / tNorm - tPropConductivity->val()( 0 ), 0.0 );

                    // if crosswind stabilization is greater than zero
                    if ( tCrosswind > 0.0 )
                    {
                        // get the velocity FI
                        // FIXME protect dof type
                        Field_Interpolator* tFIVelocity = //
                                mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

                        // get space dimension
                        uint tSpaceDim = tFIVelocity->get_number_of_fields();

                        // compute crosswind projection of velocity gradient
                        // FIXME protect dof type
                        Matrix< DDRMat > tcgradxt;
                        compute_cgradxw( { MSI::Dof_Type::VX }, mResidualDofType( 0 ), //
                                mLeaderFIManager, tSpaceDim, tEpsilon, tIsCrosswind,   //
                                tcgradxt );

                        // compute derivative of crosswind projection of velocity gradient
                        // FIXME protect dof type
                        Matrix< DDRMat > tdcgradtdu;
                        compute_dcgradxwdu( { MSI::Dof_Type::VX }, mResidualDofType( 0 ),      //
                                mLeaderFIManager, tDofType, tSpaceDim, tEpsilon, tIsCrosswind, //
                                tdcgradtdu );

                        // add contribution to jacobian per direction
                        tJac += aWStar * trans( tFITemp->dnNdxn( 1 ) ) * tCrosswind * tdcgradtdu;

                        // if derivative dof is velocity dof
                        if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) && tNorm > tEpsilon )
                        {
                            // check deno value not too small
                            // FIXME introduce inconsistent derivative
                            const real tNormDeno = std::max( std::pow( tNorm, 3.0 ), tEpsilon );

                            // add contribution of derivative of crosswind stabilization to jacobian
                            tJac -= aWStar * trans( tFITemp->dnNdxn( 1 ) ) * tcgradxt * //
                                    tSPValue * tRAbs * trans( tFITemp->gradx( 1 ) ) * tFITemp->dnNdxn( 1 ) / tNormDeno;
                        }

                        // if absolute value of strong form residual greater than zero
                        if( tRAbs > tEpsilon )
                        {
                            // add contribution of derivative of crosswind stabilization to jacobian
                            tJac += aWStar * trans( tFITemp->dnNdxn( 1 ) ) * tcgradxt * //
                                    tSPValue * tRT( 0 ) * tJT / ( tNorm * tRAbs );
                        }

                        // if conductivity depends on dof
                        if ( tPropConductivity->check_dof_dependency( tDofType ) )
                        {
                            // add contribution of derivative of crosswind stabilization to jacobian
                            tJac -= aWStar * trans( tFITemp->dnNdxn( 1 ) ) * tcgradxt *//
                                    tPropConductivity->dPropdDOF( tDofType );
                        }

                        // if crosswind diffusion and SP depends on dof
                        if ( tIsCrosswind && tSPCrosswind->check_dof_dependency( tDofType ) )
                        {
                            // add contribution of derivative of crosswind stabilization to jacobian
                            tJac += aWStar * trans( tFITemp->dnNdxn( 1 ) ) * tcgradxt * //
                                    tRAbs * tSPCrosswind->dSPdLeaderDOF( tDofType ) / tNorm ;
                        }

                        // if isotropic diffusion and SP depends on dof
                        if ( !tIsCrosswind && tSPIsotropic->check_dof_dependency( tDofType ) )
                        {
                            // add contribution of derivative of crosswind stabilization to jacobian
                            tJac += aWStar * trans( tFITemp->dnNdxn( 1 ) ) * tcgradxt * //
                                    tRAbs * tSPIsotropic->dSPdLeaderDOF( tDofType ) / tNorm ;
                        }
                    }
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Advection_Bulk::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Advection_Bulk::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Advection_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Advection_Bulk::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Advection_Bulk::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Advection_Bulk::compute_residual_strong_form( Matrix< DDRMat > & aRT )
        {
            // get the velocity dof field interpolator
            // FIXME protect dof type
            Field_Interpolator* tFIVelocity =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > & tCMDiffusion =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::DIFFUSION ) );

            // get body load property
            const std::shared_ptr< Property > & tPropLoad =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::BODY_LOAD ) );

            aRT = tCMDiffusion->EnergyDot() +
                    tFIVelocity->val_trans() * tCMDiffusion->gradEnergy() -
                    tCMDiffusion->divflux();

            // if body load exists
            if ( tPropLoad != nullptr )
            {
                // add energy source term
                aRT -= tPropLoad->val();
            }
        }

        //------------------------------------------------------------------------------

        void IWG_Advection_Bulk::compute_jacobian_strong_form (
                const Vector< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & aJT )
        {
            // get the res dof and the derivative dof FIs
            Field_Interpolator * tFIDer =
                    mLeaderFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get the velocity dof field interpolator
            // FIXME protect dof type
            Field_Interpolator* tFIVelocity =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // initialize aJT
            aJT.set_size( 1, tFIDer->get_number_of_space_time_coefficients());

            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > & tCMDiffusion =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::DIFFUSION ) );

            // get body load property
            const std::shared_ptr< Property > & tPropLoad =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::BODY_LOAD ) );

            // if CM diffusion depends on dof type
            if( tCMDiffusion->check_dof_dependency( aDofTypes ) )
            {
                // compute contribution to Jacobian strong form
                aJT =   tCMDiffusion->dEnergyDotdDOF( aDofTypes ) +
                        tFIVelocity->val_trans() * tCMDiffusion->dGradEnergydDOF( aDofTypes ) -
                        tCMDiffusion->ddivfluxdu( aDofTypes );
            }
            else
            {
                aJT.fill(0.0);
            }

            // if derivative wrt to velocity dof type
            if( aDofTypes( 0 ) == MSI::Dof_Type::VX )
            {
                aJT += trans( tCMDiffusion->gradEnergy() ) * tFIVelocity->N() ;
            }

            // if body load exists and depends on DOFs
            if ( tPropLoad != nullptr )
            {
                if ( tPropLoad->check_dof_dependency( aDofTypes ) )
                {
                    aJT -= tPropLoad->dPropdDOF( aDofTypes );
                }
            }
        }

        //------------------------------------------------------------------------------
}    // namespace moris::fem
