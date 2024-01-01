/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Spalart_Allmaras_Turbulence_Bulk.cpp
 *
 */

#include "cl_FEM_IWG_Spalart_Allmaras_Turbulence_Bulk.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "fn_FEM_IWG_Crosswind_Stabilization_Tools.hpp"
#include "cl_FEM_CM_Spalart_Allmaras_Turbulence.hpp"

// LINALG/src
#include "fn_trans.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Spalart_Allmaras_Turbulence_Bulk::IWG_Spalart_Allmaras_Turbulence_Bulk()
        {
            // set size for the constitutive model pointer cell
            mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "SpalartAllmarasTurbulence" ] = static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "SUPG" ] = static_cast< uint >( IWG_Stabilization_Type::SUPG );
            mStabilizationMap[ "DiffusionCrosswind" ] = static_cast< uint >( IWG_Stabilization_Type::DIFFUSION_CROSSWIND );
            mStabilizationMap[ "DiffusionIsotropic" ] = static_cast< uint >( IWG_Stabilization_Type::DIFFUSION_ISOTROPIC );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Spalart_Allmaras_Turbulence_Bulk::compute_residual( real aWStar )
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
            Field_Interpolator* tFIViscosity = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the SA turbulence CM
            const std::shared_ptr< Constitutive_Model >& tCMSATurbulence =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );

            // get the SUPG stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPSUPG =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SUPG ) );

            // get the crosswind diffusion stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPCrosswind =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIFFUSION_CROSSWIND ) );

            // get the isotropic diffusion  stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPIsotropic =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIFFUSION_ISOTROPIC ) );

            // compute residual of the strong form
            Matrix< DDRMat > tR;
            this->compute_residual_strong_form( tR );

            // get sub-matrix of residual
            auto tRes = mSet->get_residual()( 0 )( { tLeaderResStartIndex, tLeaderResStopIndex } );

            // compute the residual weak form
            tRes += aWStar
                  * ( tFIViscosity->N_trans()
                                  * ( tFIViscosity->gradt( 1 ) + trans( tCMSATurbulence->modified_velocity() ) * tFIViscosity->gradx( 1 )
                                          - tCMSATurbulence->production_term() + tCMSATurbulence->wall_destruction_term() )
                          + trans( tFIViscosity->dnNdxn( 1 ) )
                                    * ( tCMSATurbulence->diffusion_coefficient()( 0 ) * tFIViscosity->gradx( 1 )
                                            + tCMSATurbulence->modified_velocity() * tSPSUPG->val()( 0 ) * tR( 0 ) ) );

            // if crosswind stabilization
            if( tSPCrosswind || tSPIsotropic )
            {
                // bool for crosswind or isotropic
                bool tIsCrosswind = tSPCrosswind != nullptr;

                // cast constitutive model base class pointer to SA constitutive model
                CM_Spalart_Allmaras_Turbulence* tCMSATurbulencePtr =
                        dynamic_cast< CM_Spalart_Allmaras_Turbulence* >( tCMSATurbulence.get() );

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
                real tNorm = std::max( norm( tFIViscosity->gradx( 1 ) ), tEpsilon );

                // compute the abs of the strong form of the residual
                real tRAbs = std::max( std::abs( tR( 0, 0 ) ), tEpsilon );

                // compute full crosswind stabilization parameter value
                real tCrosswind = std::max( tSPValue * tRAbs / tNorm - tCMSATurbulencePtr->diffusion_coefficient()( 0 ), 0.0) ;

                // id crosswind stabilization parameter is greater than zero
                if( tCrosswind > 0.0 )
                {
                    // get the velocity FI
                    // FIXME protect dof type
                    Field_Interpolator* tFIVelocity = //
                            mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

                    // get space dimension
                    uint tSpaceDim = tFIVelocity->get_number_of_fields();

                    // compute crosswind projection of modified viscosity gradient
                    // FIXME protect velocity dof type
                    Matrix< DDRMat > tcgradxnu;
                    compute_cgradxw( { MSI::Dof_Type::VX }, mResidualDofType( 0 ), //
                            mLeaderFIManager, tSpaceDim, tEpsilon, tIsCrosswind,   //
                            tcgradxnu );

                    // add contribution to residual
                    tRes += aWStar * trans( tFIViscosity->dnNdxn( 1 ) ) * tCrosswind * tcgradxnu;
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Spalart_Allmaras_Turbulence_Bulk::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Spalart_Allmaras_Turbulence_Bulk::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here velocity), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get the residual viscosity FI
            Field_Interpolator* tFIViscosity = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the SA turbulence CM
            const std::shared_ptr< Constitutive_Model >& tCMSATurbulence =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );

            // get the SUPG stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPSUPG =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SUPG ) );

            // get the crosswind diffusion stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPCrosswind =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIFFUSION_CROSSWIND ) );

            // get the isotropic diffusion stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPIsotropic =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIFFUSION_ISOTROPIC ) );

            // compute residual of the strng form
            Matrix< DDRMat > tR;
            this->compute_residual_strong_form( tR );

            // get number of dof dependencies
            uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            // loop over the dof dependencies
            for ( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                const Vector< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                const sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                const uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                const uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // extract sub-matrix
                auto tJac = mSet->get_jacobian()(
                        { tLeaderResStartIndex, tLeaderResStopIndex }, { tLeaderDepStartIndex, tLeaderDepStopIndex } );

                // if residual dof type (here viscosity)
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // compute the jacobian
                    tJac += aWStar
                            * ( tFIViscosity->N_trans() * tFIViscosity->dnNdtn( 1 )
                                    + tFIViscosity->N_trans() * trans( tCMSATurbulence->modified_velocity() ) * tFIViscosity->dnNdxn( 1 )
                                    + trans( tFIViscosity->dnNdxn( 1 ) ) * tCMSATurbulence->diffusion_coefficient()( 0 ) * tFIViscosity->dnNdxn( 1 ) );
                }

                // if turbulence model depends on dof
                if ( tCMSATurbulence->check_dof_dependency( tDofType ) )
                {
                    // compute the jacobian
                    tJac += aWStar
                            * ( tFIViscosity->N_trans() * ( tCMSATurbulence->dwalldestructiontermdu( tDofType ) - tCMSATurbulence->dproductiontermdu( tDofType ) )
                                    + trans( tFIViscosity->dnNdxn( 1 ) ) * tFIViscosity->gradx( 1 ) * tCMSATurbulence->ddiffusioncoeffdu( tDofType )
                                    + tFIViscosity->N_trans() * trans( tFIViscosity->gradx( 1 ) ) * tCMSATurbulence->dmodvelocitydu( tDofType )
                                    + trans( tFIViscosity->dnNdxn( 1 ) ) * tCMSATurbulence->dmodvelocitydu( tDofType ) * tSPSUPG->val()( 0 ) * tR( 0 ) );
                }

                // if SP SUPG depends on dof
                if ( tSPSUPG->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    tJac += aWStar * ( trans( tFIViscosity->dnNdxn( 1 ) ) * tCMSATurbulence->modified_velocity() * tR( 0 ) * tSPSUPG->dSPdLeaderDOF( tDofType ) );
                }

                // compute jacobian of the strong form
                Matrix< DDRMat > tJ;
                this->compute_jacobian_strong_form( tDofType, tJ );

                // compute the jacobian
                tJac += aWStar * ( trans( tFIViscosity->dnNdxn( 1 ) ) * tCMSATurbulence->modified_velocity() * tSPSUPG->val()( 0 ) * tJ );

                // if crosswind stabilization
                if( tSPCrosswind || tSPIsotropic )
                {
                    // bool for crosswind or isotropic
                    bool tIsCrosswind = tSPCrosswind != nullptr;

                    // cast constitutive model base class pointer to SA constitutive model
                    CM_Spalart_Allmaras_Turbulence* tCMSATurbulencePtr =
                            dynamic_cast< CM_Spalart_Allmaras_Turbulence* >( tCMSATurbulence.get() );

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
                    real tNorm = std::max( norm( tFIViscosity->gradx( 1 ) ), tEpsilon );

                    // compute the norm of the gradient of the strong form of the residual
                    real tRAbs = std::max( std::abs( tR( 0, 0 ) ), tEpsilon );

                    // compute full crosswind stabilization parameter value
                    real tCrosswind = std::max( tSPValue * tRAbs / tNorm - tCMSATurbulencePtr->diffusion_coefficient()( 0 ), 0.0 );

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
                        Matrix< DDRMat > tcgradxnu;
                        compute_cgradxw( { MSI::Dof_Type::VX }, mResidualDofType( 0 ), //
                                mLeaderFIManager, tSpaceDim, tEpsilon, tIsCrosswind,   //
                                tcgradxnu );

                        // compute derivative of crosswind projection of velocity gradient
                        // FIXME protect dof type
                        Matrix< DDRMat > tdcgradnudu;
                        compute_dcgradxwdu( { MSI::Dof_Type::VX }, mResidualDofType( 0 ),      //
                                mLeaderFIManager, tDofType, tSpaceDim, tEpsilon, tIsCrosswind, //
                                tdcgradnudu );

                        // add contribution to jacobian per direction
                        tJac += aWStar * trans( tFIViscosity->dnNdxn( 1 ) ) * tCrosswind * tdcgradnudu;

                        // if derivative dof is velocity dof
                        if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) && tNorm > tEpsilon )
                        {
                            // check deno value not too small
                            // FIXME introduce inconsistent derivative
                            real tNormDeno = std::max( std::pow( tNorm, 3.0 ), tEpsilon );

                            // add contribution of derivative of the norm to jacobian
                            tJac -= aWStar * trans( tFIViscosity->dnNdxn( 1 ) ) * tcgradxnu * //
                                    tSPValue * tRAbs * trans( tFIViscosity->gradx( 1 ) ) * tFIViscosity->dnNdxn( 1 ) / tNormDeno;
                        }

                        // if absolute value of strong form of residual is greater than zero
                        if( tRAbs > tEpsilon )
                        {
                            // add contribution of the derivative of the strong of the residual to jacobian
                            tJac += aWStar * trans( tFIViscosity->dnNdxn( 1 ) ) * tcgradxnu * //
                                    tSPValue * tR( 0, 0 ) * tJ / ( tNorm * tRAbs );
                        }

                        // if dynamic viscosity depends on dof
                        if ( tCMSATurbulencePtr->check_dof_dependency( tDofType ) )
                        {
                            // add contribution of derivative of dynamic viscosity to jacobian
                            tJac -= aWStar * trans( tFIViscosity->dnNdxn( 1 ) ) * tcgradxnu * //
                                    tCMSATurbulencePtr->ddiffusioncoeffdu( tDofType );
                        }

                        // if crosswind diffusion and SP depends on dof
                        if ( tIsCrosswind && tSPCrosswind->check_dof_dependency( tDofType ) )
                        {
                            // add contribution of derivative of crosswind stabilization parameter
                            tJac += aWStar * trans( tFIViscosity->dnNdxn( 1 ) ) * tcgradxnu * //
                                    tRAbs * tSPCrosswind->dSPdLeaderDOF( tDofType ) / tNorm ;
                        }

                        // if isotropic diffusion and SP depends on dof
                        if ( !tIsCrosswind && tSPIsotropic->check_dof_dependency( tDofType ) )
                        {
                            // add contribution of derivative of crosswind stabilization parameter
                            tJac += aWStar * trans( tFIViscosity->dnNdxn( 1 ) ) * tcgradxnu * //
                                    tRAbs * tSPIsotropic->dSPdLeaderDOF( tDofType ) / tNorm ;
                        }
                    }
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Spalart_Allmaras_Turbulence_Bulk::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------
        void
        IWG_Spalart_Allmaras_Turbulence_Bulk::compute_jacobian_and_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------
        void
        IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dRdp( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_Bulk::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Spalart_Allmaras_Turbulence_Bulk::compute_residual_strong_form( Matrix< DDRMat >& aR )
        {
            // get the residual viscosity FI
            Field_Interpolator* tFIViscosity = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the SA turbulence CM
            const std::shared_ptr< Constitutive_Model >& tCMSATurbulence =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );

            // compute strong form of residual
            aR = tFIViscosity->gradt( 1 )
               + trans( tCMSATurbulence->modified_velocity() ) * tFIViscosity->gradx( 1 )
               - tCMSATurbulence->production_term()
               + tCMSATurbulence->wall_destruction_term()
               - tCMSATurbulence->divflux();

            MORIS_ASSERT( isfinite( aR ),
                    "IWG_Spalart_Allmaras_Turbulence_Bulk::compute_residual_strong_form - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Spalart_Allmaras_Turbulence_Bulk::compute_jacobian_strong_form(
                const moris::Vector< MSI::Dof_Type >& aDofTypes,
                Matrix< DDRMat >&                   aJ )
        {
            // get the residual viscosity FI
            Field_Interpolator* tFIViscosity = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the SA turbulence CM
            const std::shared_ptr< Constitutive_Model >& tCMSATurbulence =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );

            // get the der FI
            Field_Interpolator* tFIDer = mLeaderFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init aJ
            aJ.set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // if dof type is residual dof type (here viscosity)
            if ( aDofTypes( 0 ) == mResidualDofType( 0 )( 0 ) )
            {
                aJ = tFIViscosity->dnNdtn( 1 ) + trans( tCMSATurbulence->modified_velocity() ) * tFIViscosity->dnNdxn( 1 );
            }
            else
            {
                aJ.fill( 0.0 );
            }

            // if turbulence model depends on dof
            if ( tCMSATurbulence->check_dof_dependency( aDofTypes ) )
            {
                // add contribution to jacobian
                aJ += - tCMSATurbulence->dproductiontermdu( aDofTypes )
                      + tCMSATurbulence->dwalldestructiontermdu( aDofTypes )
                      - tCMSATurbulence->ddivfluxdu( aDofTypes )
                      + trans( tFIViscosity->gradx( 1 ) ) * tCMSATurbulence->dmodvelocitydu( aDofTypes );
            }

            MORIS_ASSERT( isfinite( aJ ),
                    "IWG_Spalart_Allmaras_Turbulence_Bulk::compute_jacobian_strong_form - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

