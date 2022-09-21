/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Compressible_NS_Velocity_Dirichlet_Nitsche.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Compressible_NS_Velocity_Dirichlet_Nitsche.hpp"
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

        IWG_Compressible_NS_Velocity_Dirichlet_Nitsche::IWG_Compressible_NS_Velocity_Dirichlet_Nitsche( sint aBeta )
        {
            // set sign for symmetric/unsymmetric Nitsche
            mBeta = aBeta;

            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "PrescribedValue" ] = static_cast< uint >( IWG_Property_Type::PRESCRIBED_VALUE );
            mPropertyMap[ "Select" ]    = static_cast< uint >( IWG_Property_Type::SELECT );

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

        void IWG_Compressible_NS_Velocity_Dirichlet_Nitsche::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef MORIS_HAVE_DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the master field interpolator for the residual dof type
            Field_Interpolator * tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );

            // get the selection matrix property
            std::shared_ptr< Property > & tPropSelect =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            // set a default selection matrix if needed
            Matrix< DDRMat > tM;
            if ( tPropSelect == nullptr )
            {
                // get spatial dimension of velocity field
                uint tSpaceDim = tFIVelocity->get_number_of_fields();

                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = tPropSelect->val();
            }

            // get the imposed velocity property
            std::shared_ptr< Property > & tPropVelocity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::PRESCRIBED_VALUE ) );

            // get the fluid constitutive model
            std::shared_ptr< Constitutive_Model > & tCMFluid =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID ) );

            // get the Nitsche stabilization parameter
            std::shared_ptr< Stabilization_Parameter > & tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_PENALTY_PARAMETER ) );

            // compute the jump
            Matrix< DDRMat > tVelocityJump = tFIVelocity->val() - tPropVelocity->val();

            // compute master residual
            mSet->get_residual()( 0 )(
                    { tMasterResStartIndex, tMasterResStopIndex },
                    { 0, 0 } ) -= aWStar * (
                            mBeta * trans( tCMFluid->testTraction( mNormal, mResidualDofType( 0 ), CM_Function_Type::MECHANICAL ) ) * tM * tVelocityJump );

            // if residual dof type is velocity
            if ( mResidualDofType( 0 )( 0 ) == mDofVelocity )
            {
                // compute master residual
                mSet->get_residual()( 0 )(
                        { tMasterResStartIndex, tMasterResStopIndex },
                        { 0, 0 } ) -= aWStar * (
                                tFIVelocity->N_trans() * tM * tCMFluid->traction( mNormal, CM_Function_Type::MECHANICAL) +
                                tSPNitsche->val()( 0 ) * tFIVelocity->N_trans() * tM * tVelocityJump );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Compressible_NS_Velocity_Dirichlet_Nitsche::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Velocity_Dirichlet_Nitsche::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the master field interpolator for residual dof type
            Field_Interpolator * tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );

            // get the selection matrix property
            std::shared_ptr< Property > & tPropSelect =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            // set a default selection matrix if needed
            Matrix< DDRMat > tM;
            if ( tPropSelect == nullptr )
            {
                // get spatial dimension of velocity field
                uint tSpaceDim = tFIVelocity->get_dof_type().size();

                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = tPropSelect->val();
            }

            // get the imposed velocity property
            std::shared_ptr< Property > & tPropVelocity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::PRESCRIBED_VALUE ) );

            // get the fluid constitutive model
            std::shared_ptr< Constitutive_Model > & tCMFluid =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID ) );

            // get the Nitsche stabilization parameter
            std::shared_ptr< Stabilization_Parameter > & tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_PENALTY_PARAMETER ) );

            // compute the jump
            Matrix< DDRMat > tVelocityJump = tFIVelocity->val() - tPropVelocity->val();

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
                        { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                mBeta * tCMFluid->dTestTractiondDOF(
                                        tDofType, mNormal, tM * tVelocityJump, mResidualDofType( 0 ), CM_Function_Type::MECHANICAL ) );

                //---------------------------------------------------------------------
                // if residual dof type is velocity
                if ( mResidualDofType( 0 )( 0 ) == mDofVelocity )
                {
                    // compute master residual
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                    tFIVelocity->N_trans() * tM * tCMFluid->dTractiondDOF( tDofType, mNormal, CM_Function_Type::MECHANICAL ) );

                    // if dof type is velocity
                    if ( tDofType( 0 ) == mDofVelocity )
                    {
                        // compute jacobian direct dependencies
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                        mBeta * trans( tCMFluid->testTraction( mNormal, mResidualDofType( 0 ), CM_Function_Type::MECHANICAL ) ) * tM * tFIVelocity->N()
                                        + tSPNitsche->val()( 0 )  * tFIVelocity->N_trans() * tM * tFIVelocity->N() );
                    }

                    // if imposed velocity depends on dof type
                    if ( tPropVelocity->check_dof_dependency( tDofType ) )
                    {
                        // add contribution from property to jacobian
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                        mBeta * trans( tCMFluid->testTraction( mNormal, mResidualDofType( 0 ), CM_Function_Type::MECHANICAL ) ) * tM * tPropVelocity->dPropdDOF( tDofType )
                                        + tSPNitsche->val()( 0 ) * tFIVelocity->N_trans() * tM * tPropVelocity->dPropdDOF( tDofType ) );
                    }

                    // if stabilization parameter depends on the dof type
                    if ( tSPNitsche->check_dof_dependency( tDofType ) )
                    {
                        // add contribution of SP to jacobian
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                        tFIVelocity->N_trans() * tM * tVelocityJump * tSPNitsche->dSPdMasterDOF( tDofType ) );
                    }
                }
                //---------------------------------------------------------------------
                else // residual dof type is density or temperature
                {
                    // if dof type is velocity
                    if ( tDofType( 0 ) == mDofVelocity )
                    {
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                        mBeta * trans( tCMFluid->testTraction( mNormal, mResidualDofType( 0 ), CM_Function_Type::MECHANICAL ) ) * tM * tFIVelocity->N() );
                    }

                    // if imposed velocity depends on dof type
                    if ( tPropVelocity->check_dof_dependency( tDofType ) )
                    {
                        // add contribution from property to jacobian
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                        mBeta * trans( tCMFluid->testTraction( mNormal, mResidualDofType( 0 ), CM_Function_Type::MECHANICAL ) ) * tM * tPropVelocity->dPropdDOF( tDofType ) );
                    }
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Compressible_NS_Velocity_Dirichlet_Nitsche::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Velocity_Dirichlet_Nitsche::compute_jacobian_and_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Compressible_NS_Velocity_Dirichlet_Nitsche::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Velocity_Dirichlet_Nitsche::compute_dRdp( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Compressible_NS_Velocity_Dirichlet_Nitsche::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

