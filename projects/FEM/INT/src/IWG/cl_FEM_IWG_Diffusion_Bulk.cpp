/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Diffusion_Bulk.cpp
 *
 */

#include "cl_FEM_IWG_Diffusion_Bulk.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
// LINALG/src
#include "fn_trans.hpp"

namespace moris::fem
{

    //------------------------------------------------------------------------------

    IWG_Diffusion_Bulk::IWG_Diffusion_Bulk()
    {
        // set size for the property pointer cell
        mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "Load" ]       = static_cast< uint >( IWG_Property_Type::BODY_LOAD );
        mPropertyMap[ "Thickness" ]  = static_cast< uint >( IWG_Property_Type::THICKNESS );
        mPropertyMap[ "H2Penalty" ]  = static_cast< uint >( IWG_Property_Type::H2_PENALTY );
        mPropertyMap[ "H3Penalty" ]  = static_cast< uint >( IWG_Property_Type::H3_PENALTY );
        mPropertyMap[ "PhaseField" ] = static_cast< uint >( IWG_Property_Type::Phase_Field );
        mPropertyMap[ "Select" ]     = static_cast< uint >( IWG_Property_Type::SELECT );

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "Diffusion" ] = static_cast< uint >( IWG_Constitutive_Type::DIFFUSION );

        // set size for the stabilization parameter pointer cell
        mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

        // populate the stabilization map
        mStabilizationMap[ "GGLSParam" ] = static_cast< uint >( IWG_Stabilization_Type::GGLS_DIFFUSION );
    }

    //------------------------------------------------------------------------------

    void
    IWG_Diffusion_Bulk::compute_residual( real aWStar )
    {
        // check leader field interpolators
#ifdef MORIS_HAVE_DEBUG
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get residual dof type field interpolator (here temperature)
            Field_Interpolator* tFITemp = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get body load property
            const std::shared_ptr< Property >& tPropLoad =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::BODY_LOAD ) );

            // get select property
            const std::shared_ptr< Property >& tPropSelect =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            real tSelectValue = 1.0;
            if ( tPropSelect != nullptr )
            {
                tSelectValue = tPropSelect->val()( 0 );
            }

            // get the elasticity CM
            const std::shared_ptr< Constitutive_Model >& tCMDiffusion =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::DIFFUSION ) );

            // get the Stabilization Parameter
            const std::shared_ptr< Stabilization_Parameter >& tGGLSParam =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::GGLS_DIFFUSION ) );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // get H2 penalty property
            const std::shared_ptr< Property >& tPropH2Pen =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::H2_PENALTY ) );

            // get H3 penalty property
            const std::shared_ptr< Property >& tPropH3Pen =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::H3_PENALTY ) );

            // get phase field property
            const std::shared_ptr< Property >& tPropPhaseField =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::Phase_Field ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // get sub-matrix
            auto tRes = mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex } );

            // compute the residual
            tRes += aWStar * (                                                           //
                            tCMDiffusion->testStrain_trans() * tCMDiffusion->flux() +    //
                            tFITemp->N_trans() * tCMDiffusion->EnergyDot() );

            // if body load
            if ( tPropLoad != nullptr )
            {
                // compute contribution of body load to residual
                tRes -= aWStar * tSelectValue * ( tFITemp->N_trans() * tPropLoad->val()( 0 ) );
            }

            // if H2 penalty
            if ( tPropH2Pen != nullptr )
            {
                // compute contribution of body load to residual
                tRes += aWStar * tPropH2Pen->val()( 0 ) * ( trans( tFITemp->dnNdxn( 2 ) ) * tFITemp->gradx( 2 ) );
            }

            // if H3 penalty
            if ( tPropH3Pen != nullptr )
            {
                // compute contribution of body load to residual
                tRes += aWStar * tPropH3Pen->val()( 0 ) * ( trans( tFITemp->dnNdxn( 3 ) ) * tFITemp->gradx( 3 ) );
            }

            // if phase field term
            if ( tPropPhaseField )
            {
                tRes += aWStar * tPropPhaseField->val()( 0 ) * tFITemp->N_trans()    //
                      * ( tFITemp->val()( 0 ) * ( 1.0 - tFITemp->val()( 0 ) ) * ( 1.0 - 2.0 * tFITemp->val()( 0 ) ) );
            }

            // if stabilization parameter is defined
            if ( tGGLSParam != nullptr )
            {
                // compute the residual from bulk diffusion term
                tRes += aWStar * (                                                             //
                                tCMDiffusion->testStrain_trans() * tGGLSParam->val()( 0 ) *    //
                                ( tCMDiffusion->gradEnergyDot() - tCMDiffusion->graddivflux() ) );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Diffusion_Bulk::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Bulk::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get field interpolator for a given dof type
            Field_Interpolator* tFITemp = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get body load property
            const std::shared_ptr< Property >& tPropLoad =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::BODY_LOAD ) );

            // get select property
            const std::shared_ptr< Property >& tPropSelect =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            real tSelectValue = 1.0;
            if ( tPropSelect != nullptr )
            {
                tSelectValue = tPropSelect->val()( 0 );
            }

            // get the elasticity CM
            const std::shared_ptr< Constitutive_Model >& tCMDiffusion =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::DIFFUSION ) );

            // get the Stabilization Parameter
            const std::shared_ptr< Stabilization_Parameter >& tGGLSParam =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::GGLS_DIFFUSION ) );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // get H2 penalty property
            const std::shared_ptr< Property >& tPropH2Pen =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::H2_PENALTY ) );

            // get H3 penalty property
            const std::shared_ptr< Property >& tPropH3Pen =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::H3_PENALTY ) );

            // get phase field property
            const std::shared_ptr< Property >& tPropPhaseField =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::Phase_Field ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // get the number of leader dof type dependencies
            uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            // loop over leader dof type dependencies
            for ( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                const Vector< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // get sub-matrix
                auto tJac = mSet->get_jacobian()(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } );

                // if body load
                if ( tPropLoad != nullptr )
                {
                    // if body load property has dependency on the dof type
                    if ( tPropLoad->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to Jacobian
                        tJac -= aWStar * tSelectValue * ( tFITemp->N_trans() * tPropLoad->dPropdDOF( tDofType ) );
                    }
                }

                // for explicit dependence on residual dof type
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // if H2 penalty
                    if ( tPropH2Pen != nullptr )
                    {
                        // compute contribution of body load to residual
                        tJac += aWStar * tPropH2Pen->val()( 0 ) * ( trans( tFITemp->dnNdxn( 2 ) ) * tFITemp->dnNdxn( 2 ) );
                    }

                    // if H2 penalty
                    if ( tPropH3Pen != nullptr )
                    {
                        // compute contribution of body load to residual
                        tJac += aWStar * tPropH3Pen->val()( 0 ) * ( trans( tFITemp->dnNdxn( 3 ) ) * tFITemp->dnNdxn( 3 ) );
                    }

                    // if phase field term
                    if ( tPropPhaseField )
                    {
                        tJac += aWStar * tPropPhaseField->val()( 0 ) * tFITemp->N_trans() * (                    //
                                        ( ( 1.0 - tFITemp->val()( 0 ) ) * ( 1.0 - 2.0 * tFITemp->val()( 0 ) )    //
                                                - tFITemp->val()( 0 ) * ( 1.0 - 2.0 * tFITemp->val()( 0 ) )      //
                                                - 2.0 * tFITemp->val()( 0 ) * ( 1.0 - tFITemp->val()( 0 ) ) )    //
                                        * tFITemp->N() );
                    }
                }

                // if constitutive model has dependency on the dof type
                if ( tCMDiffusion->check_dof_dependency( tDofType ) )
                {
                    // compute the Jacobian
                    tJac += aWStar * (                                                                          //
                                    tCMDiffusion->testStrain_trans() * tCMDiffusion->dFluxdDOF( tDofType ) +    //
                                    tFITemp->N_trans() * tCMDiffusion->dEnergyDotdDOF( tDofType ) );
                    // FIXME add derivative of the test strain
                }

                // if stabilization parameter is defined
                if ( tGGLSParam != nullptr )
                {
                    // FIXME: spatial derivative of body load property needed
                    // // if body load
                    // if ( tPropLoad != nullptr )
                    // {
                    //     // if body load property has dependency on the dof type
                    //     if ( tPropLoad->check_dof_dependency( tDofType ) )
                    //     {
                    //         // compute contribution of body load to Jacobian
                    //         mSet->get_jacobian()( { tLeaderResStartIndex, tLeaderResStopIndex },
                    //                 { tLeaderDepStartIndex, tLeaderDepStopIndex } )
                    //                 -= aWStar * ( trans( tFITemp->N() ) * tGGLSParam->val()( 0 ) * tPropLoad->dPropdDOF( tDofType ) );
                    //     }
                    // }

                    // if constitutive model has dependency on the dof type
                    if ( tCMDiffusion->check_dof_dependency( tDofType ) )
                    {
                        // compute the Jacobian
                        tJac += aWStar * (                                                             //
                                        tGGLSParam->val()( 0 ) * tCMDiffusion->testStrain_trans() *    //
                                        ( tCMDiffusion->dGradEnergyDotdDOF( tDofType ) - tCMDiffusion->dGradDivFluxdDOF( tDofType ) ) );
                    }

                    // add contribution from stabilization parameter
                    if ( tGGLSParam->check_dof_dependency( tDofType ) )
                    {
                        tJac += aWStar * (                                                                                                      //
                                        tCMDiffusion->testStrain_trans() * ( tCMDiffusion->gradEnergyDot() - tCMDiffusion->graddivflux() ) *    //
                                        tGGLSParam->dSPdLeaderDOF( tDofType ) );
                    }
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Diffusion_Bulk::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Bulk::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Diffusion_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Bulk::compute_dRdp( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Bulk::compute_dRdp - Not implemented." );

            //            // get leader index for residual dof type, indices for assembly
            //            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            //            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            //            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );
            //
            //            // get residual dof type field interpolator (here temperature)
            //            Field_Interpolator * tFITemp = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));
            //
            //            // get body load property
            //            std::shared_ptr< Property > tPropLoad
            //            = mLeaderProp( static_cast< uint >( IWG_Property_Type::BODY_LOAD ) );
            //
            //            // get density property
            //            std::shared_ptr< Property > tPropDensity
            //            = mLeaderProp( static_cast< uint >( IWG_Property_Type::DENSITY ) );
            //
            //            // get heat capacity property
            //            std::shared_ptr< Property > tPropHeatCapacity
            //            = mLeaderProp( static_cast< uint >( IWG_Property_Type::HEAT_CAPACITY ) );
            //
            //            // get the elasticity CM
            //            std::shared_ptr< Constitutive_Model > tCMDiffusion
            //            = mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::DIFFUSION ) );
            //
            //            // get number of leader dv dependencies
            //            uint tNumDvDependencies = mLeaderGlobalDvTypes.size();
            //            // loop over leader dv dependencies
            //            for( uint iDv = 0; iDv < tNumDvDependencies; iDv++ )
            //            {
            //                // get the treated dv type
            //                Vector< gen::PDV_Type > tDvType = mLeaderGlobalDvTypes( iDv );
            //
            //                // get the index for dof type, indices for assembly
            //                sint tDvDepIndex          = ;
            //                uint tLeaderDepStartIndex = ;
            //                uint tLeaderDepStopIndex  = ;
            //
            //                // get index for the treated dv type
            //                uint tIndexDep = mSet->get_dv_index_for_type( tDvType( 0 ), mtk::Leader_Follower::LEADER );
            //
            //                // if body load
            //                if ( tPropLoad != nullptr )
            //                {
            //                    // if load property has dependency on the dv type
            //                    if ( tPropLoad->check_dv_dependency( tDvType ) )
            //                    {
            //                        // compute drdpdv
            //                        mSet->get_drdpdv()( { tLeaderResStartIndex, tLeaderResStopIndex },
            //                                            { tLeaderDepStartIndex, tLeaderDepStopIndex } )
            //                        -= aWStar * ( tFI->N_trans() * tPropLoad->dPropdDV( tDvType ) );
            //                    }
            //                }
            //
            //                // if diffusion constitutive model has dependency on the dv type
            //                if ( tCMDiffusion->check_dv_dependency( tDvType ) )
            //                {
            //                    // compute the Jacobian
            //                    mSet->get_jacobian()( { tLeaderResStartIndex, tLeaderResStopIndex },
            //                                          { tLeaderDepStartIndex, tLeaderDepStopIndex } )
            //                    += aWStar * ( tCMDiffusion->testStrain_trans() * tCMDiffusion->dFluxdDV( tDvType ) );
            //                }
            //            }
        }

        //------------------------------------------------------------------------------
}    // namespace moris::fem
