/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Compressible_NS_Temperature_Dirichlet_Nitsche.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Compressible_NS_Temperature_Dirichlet_Nitsche.hpp"
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

        IWG_Compressible_NS_Temperature_Dirichlet_Nitsche::IWG_Compressible_NS_Temperature_Dirichlet_Nitsche( sint aBeta )
        {
            // set sign for symmetric/unsymmetric Nitsche
            mBeta = aBeta;

            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "PrescribedValue" ] = static_cast< uint >( IWG_Property_Type::PRESCRIBED_VALUE );

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Fluid" ] = static_cast< uint >( IWG_Constitutive_Type::FLUID );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "NitschePenaltyParameter" ] = static_cast< uint >( IWG_Stabilization_Type::NITSCHE_PENALTY_PARAMETER );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Temperature_Dirichlet_Nitsche::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the master field interpolator for the temperature field
            Field_Interpolator * tFITemp =
                    mMasterFIManager->get_field_interpolators_for_type( mDofTemperature );

            // get the imposed velocity property
            std::shared_ptr< Property > & tPropTemperature =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::PRESCRIBED_VALUE ) );

            // get the fluid constitutive model
            std::shared_ptr< Constitutive_Model > & tCMFluid =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID ) );

            // get the Nitsche stabilization parameter
            std::shared_ptr< Stabilization_Parameter > & tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_PENALTY_PARAMETER ) );

            // compute the jump
            Matrix< DDRMat > tTemperatureJump = tFITemp->val() - tPropTemperature->val();

            // compute master residual
            mSet->get_residual()( 0 )(
                    { tMasterResStartIndex, tMasterResStopIndex },
                    { 0, 0 } ) += aWStar * (
                            mBeta * trans( tCMFluid->testTraction( mNormal, mResidualDofType( 0 ), CM_Function_Type::THERMAL ) ) * tTemperatureJump );

            // if residual dof type is velocity
            if ( mResidualDofType( 0 )( 0 ) == mDofTemperature )
            {
                // compute master residual
                mSet->get_residual()( 0 )(
                        { tMasterResStartIndex, tMasterResStopIndex },
                        { 0, 0 } ) += aWStar * (
                                tFITemp->N_trans() * tCMFluid->traction( mNormal, CM_Function_Type::THERMAL )
                                + tSPNitsche->val()( 0 ) * tFITemp->N_trans() * tTemperatureJump );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Compressible_NS_Temperature_Dirichlet_Nitsche::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Temperature_Dirichlet_Nitsche::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the master field interpolator for temperature dof type
            Field_Interpolator * tFITemp =
                    mMasterFIManager->get_field_interpolators_for_type( mDofTemperature );

            // get the imposed velocity property
            std::shared_ptr< Property > & tPropTemperature =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::PRESCRIBED_VALUE ) );

            // get the fluid constitutive model
            std::shared_ptr< Constitutive_Model > & tCMFluid =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID ) );

            // get the Nitsche stabilization parameter
            std::shared_ptr< Stabilization_Parameter > & tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_PENALTY_PARAMETER ) );

            // compute the jump
            Matrix< DDRMat > tTemperatureJump = tFITemp->val() - tPropTemperature->val();

            // get number of master dependencies
            uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // compute the jacobian for indirect dof dependencies through master
            for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // test-traction contribution
                mSet->get_jacobian()(
                        { tMasterResStartIndex, tMasterResStopIndex },
                        { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                mBeta * tCMFluid->dTestTractiondDOF(
                                        tDofType, mNormal, tTemperatureJump, mResidualDofType( 0 ), CM_Function_Type::THERMAL ) );

                //---------------------------------------------------------------------
                // if residual dof type is velocity
                if ( mResidualDofType( 0 )( 0 ) == mDofTemperature )
                {
                    // compute master residual
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    tFITemp->N_trans() * tCMFluid->dTractiondDOF( tDofType, mNormal, CM_Function_Type::THERMAL ) );

                    // if dof type is velocity
                    if ( tDofType( 0 ) == mDofTemperature )
                    {
                        // compute jacobian direct dependencies
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                        mBeta * trans( tCMFluid->testTraction( mNormal, mResidualDofType( 0 ), CM_Function_Type::THERMAL ) ) * tFITemp->N()
                                        + tSPNitsche->val()( 0 ) * tFITemp->N_trans() * tFITemp->N() );
                    }

                    // if imposed velocity depends on dof type
                    if ( tPropTemperature->check_dof_dependency( tDofType ) )
                    {
                        // add contribution from property to jacobian
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                        mBeta * trans( tCMFluid->testTraction( mNormal, mResidualDofType( 0 ), CM_Function_Type::THERMAL ) ) * tPropTemperature->dPropdDOF( tDofType )
                                        + tSPNitsche->val()( 0 ) * tFITemp->N_trans() * tPropTemperature->dPropdDOF( tDofType ) );
                    }

                    // if stabilization parameter depends on the dof type
                    if ( tSPNitsche->check_dof_dependency( tDofType ) )
                    {
                        // add contribution of SP to jacobian
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                        tFITemp->N_trans() * tTemperatureJump * tSPNitsche->dSPdMasterDOF( tDofType ) );
                    }
                }
                //---------------------------------------------------------------------
                else // residual dof type is density or temperature
                {
                    // if dof type is velocity
                    if ( tDofType( 0 ) == mDofTemperature )
                    {
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                        mBeta * trans( tCMFluid->testTraction( mNormal, mResidualDofType( 0 ), CM_Function_Type::THERMAL ) ) * tFITemp->N() );
                    }

                    // if imposed velocity depends on dof type
                    if ( tPropTemperature->check_dof_dependency( tDofType ) )
                    {
                        // add contribution from property to jacobian
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                        mBeta * trans( tCMFluid->testTraction( mNormal, mResidualDofType( 0 ), CM_Function_Type::THERMAL ) ) * tPropTemperature->dPropdDOF( tDofType ) );
                    }
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Compressible_NS_Temperature_Dirichlet_Nitsche::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Temperature_Dirichlet_Nitsche::compute_jacobian_and_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Compressible_NS_Temperature_Dirichlet_Nitsche::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Temperature_Dirichlet_Nitsche::compute_dRdp( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Compressible_NS_Temperature_Dirichlet_Nitsche::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

