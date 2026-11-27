/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Nonlinear_Contact_Seitz.cpp
 *
 */

#include "cl_FEM_Contact_Gap_Data.hpp"
#include "cl_FEM_Contact_Geometry.hpp"
#include "cl_FEM_IWG_Isotropic_Struc_Nonlinear_Contact_Seitz.hpp"
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Cluster.hpp"

#include "cl_FEM_Stabilization_Parameter.hpp"
#include "cl_MSI_Dof_Type_Enums.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_Vector.hpp"
#include "fn_inv.hpp"
#include "fn_isfinite.hpp"
#include "fn_assert.hpp"
#include "cl_Matrix.hpp"
#include "fn_trans.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"

#include "fn_dot.hpp"
#include "fn_eye.hpp"

#include <string>
#include <memory>
#include <utility>

namespace moris::fem
{
    IWG_Isotropic_Struc_Nonlinear_Contact_Seitz::IWG_Isotropic_Struc_Nonlinear_Contact_Seitz(
            sint             aTheta,
            CM_Function_Type aCMFunctionType )
            : mNitscheTheta( aTheta )
            , mTheta( aTheta )    // switch for symmetric/unsymmetric/neutral Nitsche
            , mCMFunctionType( aCMFunctionType )
    {
        // set time continuity flag
        mTimeContinuity = true;

        // set IWG type -- needs to be consistent with IWG::factory
        if ( aTheta == 0 )
        {
            mIWGType = fem::IWG_Type::STRUC_NONLINEAR_CONTACT_SEITZ_UNBIASED_NEUTRAL;
        }
        else if ( aTheta == 1 )
        {
            mIWGType = fem::IWG_Type::STRUC_NONLINEAR_CONTACT_SEITZ_UNBIASED_SYMMETRIC;
        }
        else if ( aTheta == -1 )
        {
            mIWGType = fem::IWG_Type::STRUC_NONLINEAR_CONTACT_SEITZ_UNBIASED_UNSYMMETRIC;
        }
        else
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Nonlinear_Contact_Seitz::IWG_Isotropic_Struc_Nonlinear_Contact_Seitz - theta must be 0, 1, or -1." );
        }

        // set size for the property pointer cell
        mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );
        mFollowerProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "Thickness" ] = static_cast< uint >( IWG_Property_Type::THICKNESS );
        mPropertyMap[ "Select" ]    = static_cast< uint >( IWG_Property_Type::SELECT );

        // set size for the constitutive model pointer cell
        // .resize: gives aValue:(The value to initialize the new elements with) and aCount:(new size of the Cell)
        mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );
        mFollowerCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "ElastLinIso" ] = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );

        // set size for the stabilization parameter pointer cell
        mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

        // populate the stabilization map
        mStabilizationMap[ "NitscheInterface" ] = static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE );

        // set flag to use displacement for gap
        mUseDeformedGeometryForGap = true;

        // scaling factor to enable pure penalty formulation
        mTractionScaling = 1.0;
    }

    //------------------------------------------------------------------------------

    void
    IWG_Isotropic_Struc_Nonlinear_Contact_Seitz::set_parameters( const Vector< Matrix< DDRMat > >& aParameters )
    {
        // set parameters
        mParameters = aParameters;

        // Read friction, penalty-only and debug settings.
        if ( mParameters.size() > 0 )
        {
            MORIS_ERROR( mParameters.size() == 4 &&                  //
                                 mParameters( 0 ).numel() == 1 &&    //
                                 mParameters( 1 ).numel() == 1 &&    //
                                 mParameters( 2 ).numel() == 1 &&    //
                                 mParameters( 3 ).numel() == 1,
                    "IWG_Isotropic_Struc_Nonlinear_Contact_Seitz::IWG_Isotropic_Struc_Nonlinear_Contact_Seitz - Expected four scalar parameters: friction coefficient, tangential stiffness, penalty-only flag, and debug flag." );

            // parameter 0: friction coefficient
            mFrictionCoefficient = mParameters( 0 )( 0 );

            // parameter 1: frictional stiffness
            mGammaTau = mParameters( 1 )( 0 );

            // Recompute both effective settings on every complete update so
            // disabling penalty-only restores the constructor's Nitsche variant.
            const bool tPenaltyOnly = mParameters( 2 )( 0 ) > 0.5;
            mTractionScaling        = tPenaltyOnly ? 0.0 : 1.0;
            mTheta                  = tPenaltyOnly ? 0 : mNitscheTheta;

            // parameter 3: for debug purposes, to trigger debug output
            mDebugFlag = mParameters( 3 )( 0 );
        }
    }

    //------------------------------------------------------------------------------

    /**
     * Frictional contact based on Seitz, Wall and Popp (2019),
     * Eqs. (33)-(35). MORIS integrates on the leader side and retains its
     * reference-traction and geometry conventions.
     */
    void
    IWG_Isotropic_Struc_Nonlinear_Contact_Seitz::compute_residual( real aWStar )
    {
#ifdef MORIS_HAVE_DEBUG
        // check leader and follower field interpolators
        this->check_field_interpolators( mtk::Leader_Follower::LEADER );
        this->check_field_interpolators( mtk::Leader_Follower::FOLLOWER );
#endif
        // check whether contact should be enforced
        const std::shared_ptr< Property >& tSelectLeader =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

        const std::shared_ptr< Property >& tSelectFollower =
                mFollowerProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

        const real tImposeContactLeader   = ( tSelectLeader != nullptr ) ? tSelectLeader->val()( 0 ) : 1.0;
        const real tImposeContactFollower = ( tSelectFollower != nullptr ) ? tSelectFollower->val()( 0 ) : 1.0;

        if ( tImposeContactLeader < MORIS_REAL_EPS || tImposeContactFollower < MORIS_REAL_EPS )
        {
            return;
        }

        // get leader index for residual dof type, indices for assembly
        Vector< MSI::Dof_Type > const tDisplDofTypes = mResidualDofType( 0 );

        // get leader index for residual dof type, indices for assembly
        const uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
        const uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
        const uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

        // get follower index for residual dof type, indices for assembly
        const uint tFollowerDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );
        const uint tFollowerResStartIndex = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 0 );
        const uint tFollowerResStopIndex  = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 1 );

        // build gap data and remap follower coordinates
        const Matrix< DDRMat > tRemappedFollowerCoords = Contact_Geometry::remap_nonconformal_rays(
                mUseDeformedGeometryForGap,
                tDisplDofTypes,
                mLeaderFIManager,
                mFollowerFIManager,
                mLeaderPreviousFIManager,
                mFollowerPreviousFIManager,
                mGapData );

        // check whether the remapping is successful
        if ( std::abs( tRemappedFollowerCoords( 0 ) ) > 1 )
        {
            fprintf( stderr, "IWG_Isotropic_Struc_Nonlinear_Contact_Seitz::compute_residual - remapping of follower coordinates failed.\n" );
            return;    // exit if remapping was not successful
        }

        // set integration point for follower side
        mFollowerFIManager->set_space_time_from_local_IG_point( tRemappedFollowerCoords );
        mFollowerPreviousFIManager->set_space_time_from_local_IG_point( tRemappedFollowerCoords );

        // get leader field interpolator for the residual dof type
        Field_Interpolator* tFILeader =
                mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

        // get follower field interpolator for the residual dof type
        Field_Interpolator* tFIFollower =
                mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

        // get the elasticity constitutive models
        const std::shared_ptr< Constitutive_Model >& tConstitutiveModelLeader =
                mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

        // get the Nitsche stabilization parameter
        const std::shared_ptr< Stabilization_Parameter >& tSPNitsche =
                mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );

        // get thickness property
        const std::shared_ptr< Property >& tPropThickness =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

        // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
        aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

        // Nitsche parameter gamma
        const real tNitscheParam = tSPNitsche->val()( 0 );

        // Nitsche parameter for tangential contributions
        const real tNitscheTangentParam = tNitscheParam;

        // Tangential penalty parameter for frictional contributions, if not set by user, use the Nitsche parameter
        const real tGammaTau = ( mGammaTau > 0.0 ) ? mGammaTau : tNitscheTangentParam;

        // get Piola traction
        const Matrix< DDRMat > tTraction =
                mTractionScaling * tConstitutiveModelLeader->traction( mGapData->mLeaderRefNormal, mCMFunctionType );

        // get test traction
        const Matrix< DDRMat > tTestTraction =
                tConstitutiveModelLeader->testTraction_trans( mGapData->mLeaderRefNormal, tDisplDofTypes, mCMFunctionType );

        // compute contact pressure using current normal
        const real tContactPressure = mTractionScaling * dot( tTraction, mGapData->mLeaderNormal );

        // compute augmented Lagrangian term
        const real tAugLagrTerm = tContactPressure + tNitscheParam * mGapData->mGap;

        // check whether contact is enforced, i.e. tAugLagrTerm >= 0
        if ( !mGapData->is_contact( tAugLagrTerm <= 0 ) )
        {
            if ( mDebugFlag > 0 )
            {
                fprintf( stderr, "IWG_Isotropic_Struc_Nonlinear_Contact_Seitz::compute_residual - no contact enforced.\n" );
            }

            // no contact: stress penalty including both normal and tangential contributions
            mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex } ) +=
                    -0.5 * aWStar * mTheta / tNitscheParam * tTestTraction * tTraction;
        }
        else
        {
            mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex } ) +=                        //
                    0.5 * aWStar * (                                                          //
                            -tAugLagrTerm * tFILeader->N_trans() * mGapData->mLeaderNormal    //
                            + mTheta * mGapData->mGap * tTestTraction * mGapData->mLeaderNormal );

            mSet->get_residual()( 0 )(
                    { tFollowerResStartIndex, tFollowerResStopIndex } ) +=                         //
                    0.5 * aWStar * (                                                               //
                            +tAugLagrTerm * tFIFollower->N_trans() * mGapData->mLeaderNormal );    //

            // compute current tangential traction
            const Matrix< DDRMat > tTangentTraction = mGapData->mLeaderTangentProjection * tTraction;

            // tangential contributions
            if ( mFrictionCoefficient <= 0.0 )
            {
                if ( mDebugFlag > 0 )
                {
                    fprintf( stderr, "IWG_Isotropic_Struc_Nonlinear_Contact_Seitz::compute_residual - no contact enforced.\n" );
                }

                // add penalty on tangential stress in case of frictionless contact
                mSet->get_residual()( 0 )(
                        { tLeaderResStartIndex, tLeaderResStopIndex } ) +=
                        -0.5 * aWStar * mTheta / tNitscheParam * tTestTraction * tTangentTraction;
            }
            else
            {
                // make sure time continuity is set
                MORIS_ASSERT( get_time_continuity(),
                        "IWG_Isotropic_Struc_Nonlinear_Contact_Seitz::compute_residual - time continuity flag not set.\n" );

                // compute slip traction magnitude
                real tSplipTraction = -mFrictionCoefficient * tAugLagrTerm;

                // compute trial tangential traction
                Matrix< DDRMat > tTrailTangentTraction = tGammaTau * mGapData->mLeaderTangentVelocity;

                // Eq. (33) uses current bulk traction. Previous configurations enter only the slip increment.
                tTrailTangentTraction += tTangentTraction;

                // check whether the tangential traction is in the slip or stick regime and project accordingly
                Matrix< DDRMat > tProjectedTangentTraction;

                real tTrailTangentTractionNorm = norm( tTrailTangentTraction );

                // Select the stick-side derivative at equality, including q=0 and radius=0.
                // This permits an initially unloaded interface without norm division.
                if ( !mGapData->is_stick( tTrailTangentTractionNorm <= tSplipTraction ) )
                {
                    if ( mDebugFlag > 0 )
                    {
                        fprintf( stderr, "IWG_Isotropic_Struc_Nonlinear_Contact_Seitz::compute_residual - contact enforced - Slip.\n" );
                    }
                    // slip condition
                    tProjectedTangentTraction = tSplipTraction / tTrailTangentTractionNorm * tTrailTangentTraction;
                }
                else
                {
                    if ( mDebugFlag > 0 )
                    {
                        fprintf( stderr, "IWG_Isotropic_Struc_Nonlinear_Contact_Seitz::compute_residual - contact enforced - Stick.\n" );
                    }
                    // stick condition
                    tProjectedTangentTraction = tTrailTangentTraction;
                }

                mSet->get_residual()( 0 )(
                        { tLeaderResStartIndex, tLeaderResStopIndex } ) +=           //
                        0.5 * aWStar * (                                             //
                                -tFILeader->N_trans() * tProjectedTangentTraction    //
                                - mTheta / tNitscheTangentParam * tTestTraction * ( tTangentTraction - tProjectedTangentTraction ) );

                mSet->get_residual()( 0 )( { tFollowerResStartIndex, tFollowerResStopIndex } ) +=    //
                        0.5 * aWStar * (                                                             //
                                +tFIFollower->N_trans() * tProjectedTangentTraction );
            }
        }

        // call debug function
        if ( mDebugFlag > 0 )
        {
            this->debug_function();
        }

        // check for nan, infinity
        MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                "IWG_Isotropic_Struc_Nonlinear_Contact_Seitz::compute_residual - Residual contains NAN or INF, exiting!" );
    }

    //------------------------------------------------------------------------------

    void IWG_Isotropic_Struc_Nonlinear_Contact_Seitz::compute_jacobian( real aWStar )
    {
#ifdef MORIS_HAVE_DEBUG
        // check leader and follower field interpolators
        this->check_field_interpolators( mtk::Leader_Follower::LEADER );
        this->check_field_interpolators( mtk::Leader_Follower::FOLLOWER );
#endif
        // check whether contact should be enforced
        const std::shared_ptr< Property >& tSelectLeader =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

        const std::shared_ptr< Property >& tSelectFollower =
                mFollowerProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

        real tImposeContactLeader   = ( tSelectLeader != nullptr ) ? tSelectLeader->val()( 0 ) : 1.0;
        real tImposeContactFollower = ( tSelectFollower != nullptr ) ? tSelectFollower->val()( 0 ) : 1.0;

        if ( tImposeContactLeader < MORIS_REAL_EPS || tImposeContactFollower < MORIS_REAL_EPS )
        {
            return;
        }

        // get leader index for residual dof type, indices for assembly
        Vector< MSI::Dof_Type > const tDisplDofTypes = mResidualDofType( 0 );

        // get leader index for residual dof type, indices for assembly
        const uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
        const uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
        const uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

        // get follower index for residual dof type, indices for assembly
        const uint tFollowerDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );
        const uint tFollowerResStartIndex = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 0 );
        const uint tFollowerResStopIndex  = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 1 );

        // build gap data and remap follower coordinates
        const Matrix< DDRMat > tRemappedFollowerCoords = Contact_Geometry::remap_nonconformal_rays(
                mUseDeformedGeometryForGap,
                tDisplDofTypes,
                mLeaderFIManager,
                mFollowerFIManager,
                mLeaderPreviousFIManager,
                mFollowerPreviousFIManager,
                mGapData );

        // check whether the remapping is successful
        if ( std::abs( tRemappedFollowerCoords( 0 ) ) > 1 )
        {
            return;    // exit if remapping was not successful
        }

        // set integration point for follower side
        mFollowerFIManager->set_space_time_from_local_IG_point( tRemappedFollowerCoords );
        mFollowerPreviousFIManager->set_space_time_from_local_IG_point( tRemappedFollowerCoords );

        // get leader field interpolator for the residual dof type
        Field_Interpolator* tFILeader =
                mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

        // get follower field interpolator for the residual dof type
        Field_Interpolator* tFIFollower =
                mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

        // get the elasticity constitutive models
        const std::shared_ptr< Constitutive_Model >& tConstitutiveModelLeader =
                mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

        // get the Nitsche stabilization parameter
        const std::shared_ptr< Stabilization_Parameter >& tSPNitsche =
                mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );

        // get thickness property
        const std::shared_ptr< Property >& tPropThickness =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

        // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
        aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

        // Nitsche parameter gamma
        const real tNitscheParam = tSPNitsche->val()( 0 );

        // Nitsche parameter for tangential contributions
        const real tNitscheTangentParam = tNitscheParam;

        // Tangential penalty parameter for frictional contributions, if not set by user, use the Nitsche parameter
        const real tGammaTau = ( mGammaTau > 0.0 ) ? mGammaTau : tNitscheTangentParam;

        // get Piola traction
        const Matrix< DDRMat > tTraction =
                mTractionScaling * tConstitutiveModelLeader->traction( mGapData->mLeaderRefNormal, mCMFunctionType );

        // get test traction
        const Matrix< DDRMat > tTestTraction =
                tConstitutiveModelLeader->testTraction_trans( mGapData->mLeaderRefNormal, tDisplDofTypes, mCMFunctionType );

        // compute contact pressure using current normal
        const real tContactPressure = mTractionScaling * dot( tTraction, mGapData->mLeaderNormal );

        // compute augmented Lagrangian term
        const real tAugLagrTerm = tContactPressure + tNitscheParam * mGapData->mGap;

        // get number of leader dof dependencies
        const uint tLeaderNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

        // compute the Jacobian for indirect dof dependencies through leader constitutive models
        for ( uint iDOF = 0; iDOF < tLeaderNumDofDependencies; iDOF++ )
        {
            // get the dof type
            const Vector< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

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
            if ( tDofType( 0 ) == tDisplDofTypes( 0 ) )
            {
                // compute derivative of traction
                const Matrix< DDRMat > tdTractiondu =
                        tConstitutiveModelLeader->dTractiondDOF( tDofType, mGapData->mLeaderRefNormal, mCMFunctionType );

                // compute derivative of contact pressure with respect to the dof type
                const Matrix< DDRMat > tdContactPressuredu =
                        mTractionScaling * ( trans( mGapData->mLeaderNormal ) * tdTractiondu    //
                                             + trans( tTraction ) * mGapData->mLeaderdNormaldu );

                if ( tAugLagrTerm > 0 )
                {
                    tJacMM += -0.5 * aWStar * mTheta / tNitscheParam * (              //
                                      tConstitutiveModelLeader->dTestTractiondDOF(    //
                                              tDofType,                               //
                                              mGapData->mLeaderRefNormal,             //
                                              tTraction,                              //
                                              tDisplDofTypes,                         //
                                              mCMFunctionType )                       //
                                      + tTestTraction * tdTractiondu );
                }
                else
                {
                    tJacMM +=
                            0.5 * aWStar * (                                                                  //
                                    -tAugLagrTerm * tFILeader->N_trans() * mGapData->mLeaderdNormaldu         //
                                    - tFILeader->N_trans() * mGapData->mLeaderNormal * (                      //
                                              tdContactPressuredu + tNitscheParam * mGapData->mdGapdu )       //
                                    + mTheta * tConstitutiveModelLeader->dTestTractiondDOF(                   //
                                              tDofType,                                                       //
                                              mGapData->mLeaderRefNormal,                                     //
                                              mGapData->mGap * mGapData->mLeaderNormal,                       //
                                              tDisplDofTypes,                                                 //
                                              mCMFunctionType )                                               //
                                    + mTheta * mGapData->mGap * tTestTraction * mGapData->mLeaderdNormaldu    //
                                    + mTheta * tTestTraction * mGapData->mLeaderNormal * mGapData->mdGapdu );

                    tJacSM +=
                            0.5 * aWStar * (                                                                                   //
                                    +tAugLagrTerm * mGapData->multiply_follower_dshapefunctiondu( mGapData->mLeaderNormal )    //
                                    + tAugLagrTerm * tFIFollower->N_trans() * mGapData->mLeaderdNormaldu                       //
                                    + tFIFollower->N_trans() * mGapData->mLeaderNormal * (                                     //
                                              tdContactPressuredu + tNitscheParam * mGapData->mdGapdu ) );

                    // compute current tangential traction
                    const Matrix< DDRMat > tTangentTraction = mGapData->mLeaderTangentProjection * tTraction;

                    const Matrix< DDRMat > tdTangentTractiondu = mGapData->multiply_leader_dtangentprojectiondu( tTraction )    //
                                                               + mTractionScaling * mGapData->mLeaderTangentProjection * tdTractiondu;

                    if ( mFrictionCoefficient <= 0.0 )
                    {
                        tJacMM += -0.5 * aWStar * mTheta / tNitscheParam * (              //
                                          tConstitutiveModelLeader->dTestTractiondDOF(    //
                                                  tDofType,                               //
                                                  mGapData->mLeaderRefNormal,             //
                                                  tTangentTraction,                       //
                                                  tDisplDofTypes,                         //
                                                  mCMFunctionType )                       //
                                          + tTestTraction * tdTangentTractiondu );
                    }
                    else
                    {
                        // compute slip traction magnitude
                        real tSplipTraction = -mFrictionCoefficient * tAugLagrTerm;

                        // compute trial tangential traction
                        Matrix< DDRMat > tTrailTangentTraction = tGammaTau * mGapData->mLeaderTangentVelocity;

                        Matrix< DDRMat > tdTrailTangentTractiondu = tGammaTau * mGapData->mLeaderdTangentVelocitydu;

                        // Eq. (33) uses current bulk traction. Previous configurations enter only the slip increment.
                        tTrailTangentTraction += tTangentTraction;
                        tdTrailTangentTractiondu += tdTangentTractiondu;

                        Matrix< DDRMat > tProjectedTangentTraction;
                        Matrix< DDRMat > tdProjectedTangentTractiondu;

                        real tTrailTangentTractionNorm = norm( tTrailTangentTraction );

                        if ( tTrailTangentTractionNorm > tSplipTraction )
                        {
                            tProjectedTangentTraction = tSplipTraction / tTrailTangentTractionNorm * tTrailTangentTraction;

                            uint tSpaceDim = tTrailTangentTraction.n_rows();

                            tdProjectedTangentTractiondu =
                                    -mFrictionCoefficient / tTrailTangentTractionNorm * tTrailTangentTraction *                        //
                                            ( tdContactPressuredu + tNitscheParam * mGapData->mdGapdu )                                //
                                    + tSplipTraction / std::pow( tTrailTangentTractionNorm, 3.0 ) *                                    //
                                              ( tTrailTangentTractionNorm * tTrailTangentTractionNorm * eye( tSpaceDim, tSpaceDim )    //
                                                      - tTrailTangentTraction * trans( tTrailTangentTraction ) )
                                              * tdTrailTangentTractiondu;
                        }
                        else
                        {
                            tProjectedTangentTraction = tTrailTangentTraction;

                            tdProjectedTangentTractiondu = tdTrailTangentTractiondu;
                        }

                        tJacMM += 0.5 * aWStar * (                                                                                                    //
                                          -tFILeader->N_trans() * tdProjectedTangentTractiondu                                                        //
                                          - mTheta / tNitscheTangentParam * tTestTraction * ( tdTangentTractiondu - tdProjectedTangentTractiondu )    //
                                          - mTheta / tNitscheTangentParam * tConstitutiveModelLeader->dTestTractiondDOF(                              //
                                                    tDofType,                                                                                         //
                                                    mGapData->mLeaderRefNormal,                                                                       //
                                                    tTangentTraction - tProjectedTangentTraction,                                                     //
                                                    tDisplDofTypes,                                                                                   //
                                                    mCMFunctionType ) );

                        tJacSM += 0.5 * aWStar * (                                                                      //
                                          +mGapData->multiply_follower_dshapefunctiondu( tProjectedTangentTraction )    //
                                          + tFIFollower->N_trans() * tdProjectedTangentTractiondu );
                    }
                }
            }

            // if dependency on the dof type
            if ( tConstitutiveModelLeader->check_dof_dependency( tDofType ) )
            {
                if ( tAugLagrTerm > 0 )
                {
                    // move terms from above to here
                }
                else
                {
                    // move terms from above to here
                }
            }

            // if dependency of stabilization parameters on the dof type
            if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Leader_Follower::LEADER ) )
            {
                MORIS_ERROR( false, "IWG_Isotropic_Struc_Nonlinear_Contact_Seitz::compute_jacobian - state dependent stabilization not implemented." );
            }
        }

        // compute the Jacobian for indirect dof dependencies through follower constitutive models
        const uint tFollowerNumDofDependencies = mRequestedFollowerGlobalDofTypes.size();
        for ( uint iDOF = 0; iDOF < tFollowerNumDofDependencies; iDOF++ )
        {
            // get dof type
            const Vector< MSI::Dof_Type >& tDofType = mRequestedFollowerGlobalDofTypes( iDOF );

            // get the index for the dof type
            const sint tDofDepIndex           = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::FOLLOWER );
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
                if ( tAugLagrTerm > 0 )
                {
                    // nothing here
                }
                else
                {
                    tJacMS +=
                            0.5 * aWStar * (                                                                               //
                                    -tNitscheParam * tFILeader->N_trans() * mGapData->mLeaderNormal * mGapData->mdGapdv    //
                                    + mTheta * tTestTraction * mGapData->mLeaderNormal * mGapData->mdGapdv );

                    tJacSS +=
                            0.5 * aWStar * (                                                                                   //
                                    +tAugLagrTerm * mGapData->multiply_follower_dshapefunctiondv( mGapData->mLeaderNormal )    //
                                    + tNitscheParam * tFIFollower->N_trans() * mGapData->mLeaderNormal * mGapData->mdGapdv );

                    // compute current tangential traction
                    const Matrix< DDRMat > tTangentTraction = mGapData->mLeaderTangentProjection * tTraction;

                    if ( mFrictionCoefficient <= 0.0 )
                    {
                        // nothing here
                    }
                    else
                    {
                        // compute slip traction magnitude
                        real tSplipTraction = -mFrictionCoefficient * tAugLagrTerm;

                        // compute trial tangential traction
                        Matrix< DDRMat > tTrailTangentTraction = tGammaTau * mGapData->mLeaderTangentVelocity;

                        Matrix< DDRMat > tdTrailTangentTractiondv = tGammaTau * mGapData->mLeaderdTangentVelocitydv;

                        // Eq. (33) uses current bulk traction. Previous configurations enter only the slip increment.
                        tTrailTangentTraction += tTangentTraction;

                        Matrix< DDRMat > tProjectedTangentTraction;
                        Matrix< DDRMat > tdProjectedTangentTractiondv;

                        real tTrailTangentTractionNorm = norm( tTrailTangentTraction );

                        if ( tTrailTangentTractionNorm > tSplipTraction )
                        {
                            tProjectedTangentTraction = tSplipTraction / tTrailTangentTractionNorm * tTrailTangentTraction;

                            uint tSpaceDim = tTrailTangentTraction.n_rows();

                            tdProjectedTangentTractiondv =
                                    -mFrictionCoefficient * tNitscheParam / tTrailTangentTractionNorm * tTrailTangentTraction * mGapData->mdGapdv    //
                                    + tSplipTraction / std::pow( tTrailTangentTractionNorm, 3.0 ) *                                                  //
                                              ( tTrailTangentTractionNorm * tTrailTangentTractionNorm * eye( tSpaceDim, tSpaceDim )                  //
                                                      - tTrailTangentTraction * trans( tTrailTangentTraction ) )
                                              * tdTrailTangentTractiondv;
                        }
                        else
                        {
                            tProjectedTangentTraction = tTrailTangentTraction;

                            tdProjectedTangentTractiondv = tdTrailTangentTractiondv;
                        }

                        tJacMS += 0.5 * aWStar * (                                                //
                                          -tFILeader->N_trans() * tdProjectedTangentTractiondv    //
                                          + mTheta / tNitscheTangentParam * tTestTraction * tdProjectedTangentTractiondv );

                        tJacSS += 0.5 * aWStar * (                                                                      //
                                          +mGapData->multiply_follower_dshapefunctiondv( tProjectedTangentTraction )    //
                                          + tFIFollower->N_trans() * tdProjectedTangentTractiondv );
                    }
                }
            }

            // if dependency on the dof type: even if so, nothing here

            // if dependency of stabilization parameters on the dof type
            if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Leader_Follower::FOLLOWER ) )
            {
                MORIS_ERROR( false, "IWG_Isotropic_Struc_Nonlinear_Contact_Seitz::compute_jacobian - state dependent stabilization not implemented." );
            }
        }

        // check for nan, infinity
        MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                "IWG_Isotropic_Struc_Nonlinear_Contact_Seitz::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
    }

    //------------------------------------------------------------------------------

    void IWG_Isotropic_Struc_Nonlinear_Contact_Seitz::compute_jacobian_and_residual( real aWStar )
    {
        MORIS_ERROR( false, "IWG_Isotropic_Struc_Nonlinear_Contact_Seitz::compute_jacobian_and_residual - This function does nothing." );
    }

    //------------------------------------------------------------------------------

    void IWG_Isotropic_Struc_Nonlinear_Contact_Seitz::compute_dRdp( real aWStar )
    {
        MORIS_ERROR( false, "IWG_Isotropic_Struc_Nonlinear_Contact_Seitz::compute_dRdp - This function does nothing." );
    }

    //------------------------------------------------------------------------------

    void IWG_Isotropic_Struc_Nonlinear_Contact_Seitz::compute_jacobian_previous( real aWStar )
    {
#ifdef MORIS_HAVE_DEBUG
        // check leader and follower field interpolators
        this->check_field_interpolators( mtk::Leader_Follower::LEADER );
        this->check_field_interpolators( mtk::Leader_Follower::FOLLOWER );
#endif

        // check for friction
        if ( mFrictionCoefficient <= 0.0 )
        {
            return;
        }

        // check whether contact should be enforced
        const std::shared_ptr< Property >& tSelectLeader =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

        const std::shared_ptr< Property >& tSelectFollower =
                mFollowerProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

        real tImposeContactLeader   = ( tSelectLeader != nullptr ) ? tSelectLeader->val()( 0 ) : 1.0;
        real tImposeContactFollower = ( tSelectFollower != nullptr ) ? tSelectFollower->val()( 0 ) : 1.0;

        if ( tImposeContactLeader < MORIS_REAL_EPS || tImposeContactFollower < MORIS_REAL_EPS )
        {
            return;
        }

        // get leader index for residual dof type, indices for assembly
        Vector< MSI::Dof_Type > const tDisplDofTypes = mResidualDofType( 0 );

        // get leader index for residual dof type, indices for assembly
        const uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
        const uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
        const uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

        // get follower index for residual dof type, indices for assembly
        const uint tFollowerDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );
        const uint tFollowerResStartIndex = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 0 );
        const uint tFollowerResStopIndex  = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 1 );

        // build gap data and remap follower coordinates
        const Matrix< DDRMat > tRemappedFollowerCoords = Contact_Geometry::remap_nonconformal_rays(
                mUseDeformedGeometryForGap,
                tDisplDofTypes,
                mLeaderFIManager,
                mFollowerFIManager,
                mLeaderPreviousFIManager,
                mFollowerPreviousFIManager,
                mGapData );

        // check whether the remapping is successful
        if ( std::abs( tRemappedFollowerCoords( 0 ) ) > 1 )
        {
            return;    // exit if remapping was not successful
        }

        // set integration point for follower side
        mFollowerFIManager->set_space_time_from_local_IG_point( tRemappedFollowerCoords );
        mFollowerPreviousFIManager->set_space_time_from_local_IG_point( tRemappedFollowerCoords );

        // get leader field interpolator for the residual dof type
        Field_Interpolator* tFILeader =
                mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

        // get follower field interpolator for the residual dof type
        Field_Interpolator* tFIFollower =
                mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

        // get the elasticity constitutive models
        const std::shared_ptr< Constitutive_Model >& tConstitutiveModelLeader =
                mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

        // get the Nitsche stabilization parameter
        const std::shared_ptr< Stabilization_Parameter >& tSPNitsche =
                mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );

        // get thickness property
        const std::shared_ptr< Property >& tPropThickness =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

        // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
        aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

        // Nitsche parameter gamma
        const real tNitscheParam = tSPNitsche->val()( 0 );

        // Nitsche parameter for tangential contributions
        const real tNitscheTangentParam = tNitscheParam;

        // Tangential penalty parameter for frictional contributions, if not set by user, use the Nitsche parameter
        const real tGammaTau = ( mGammaTau > 0.0 ) ? mGammaTau : tNitscheTangentParam;

        // get Piola traction
        const Matrix< DDRMat > tTraction =
                mTractionScaling * tConstitutiveModelLeader->traction( mGapData->mLeaderRefNormal, mCMFunctionType );

        // get test traction
        const Matrix< DDRMat > tTestTraction =
                tConstitutiveModelLeader->testTraction_trans( mGapData->mLeaderRefNormal, tDisplDofTypes, mCMFunctionType );

        // compute contact pressure using current normal
        const real tContactPressure = mTractionScaling * dot( tTraction, mGapData->mLeaderNormal );

        // compute augmented Lagrangian term
        const real tAugLagrTerm = tContactPressure + tNitscheParam * mGapData->mGap;

        // get number of leader dof dependencies
        const uint tLeaderNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

        // compute the Jacobian for indirect dof dependencies through leader constitutive models
        for ( uint iDOF = 0; iDOF < tLeaderNumDofDependencies; iDOF++ )
        {
            // get the dof type
            const Vector< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

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
            if ( tDofType( 0 ) == tDisplDofTypes( 0 ) )
            {
                // Current traction is independent of the previous-slab DOFs.
                const Matrix< DDRMat > tTangentTraction = mGapData->mLeaderTangentProjection * tTraction;

                if ( tAugLagrTerm > 0 )
                {
                    // nothing here
                }
                else
                {
                    // compute slip traction magnitude
                    real tSplipTraction = -mFrictionCoefficient * tAugLagrTerm;

                    // compute trial tangential traction
                    Matrix< DDRMat > tTrailTangentTraction = tGammaTau * mGapData->mLeaderTangentVelocity;

                    Matrix< DDRMat > tdTrailTangentTractiondu = tGammaTau * mGapData->mLeaderdTangentVelocityduprev;

                    // Eq. (33) uses current bulk traction. Previous configurations enter only the slip increment.
                    tTrailTangentTraction += tTangentTraction;

                    Matrix< DDRMat > tdProjectedTangentTractiondu;

                    real tTrailTangentTractionNorm = norm( tTrailTangentTraction );

                    if ( tTrailTangentTractionNorm > tSplipTraction )
                    {

                        uint tSpaceDim = tTrailTangentTraction.n_rows();

                        tdProjectedTangentTractiondu =
                                +tSplipTraction / std::pow( tTrailTangentTractionNorm, 3.0 ) *                           //
                                ( tTrailTangentTractionNorm * tTrailTangentTractionNorm * eye( tSpaceDim, tSpaceDim )    //
                                        - tTrailTangentTraction * trans( tTrailTangentTraction ) )
                                * tdTrailTangentTractiondu;
                    }
                    else
                    {

                        tdProjectedTangentTractiondu = tdTrailTangentTractiondu;
                    }

                    tJacMM += 0.5 * aWStar * (                                                //
                                      -tFILeader->N_trans() * tdProjectedTangentTractiondu    //
                                      + mTheta / tNitscheTangentParam * tTestTraction * tdProjectedTangentTractiondu );

                    tJacSM += 0.5 * aWStar * (    //
                                      +tFIFollower->N_trans() * tdProjectedTangentTractiondu );
                }
            }

            // if dependency on the dof type
            if ( tConstitutiveModelLeader->check_dof_dependency( tDofType ) )
            {
                if ( tAugLagrTerm > 0 )
                {
                    // move terms from above to here
                }
                else
                {
                    // move terms from above to here
                }
            }

            // if dependency of stabilization parameters on the dof type
            if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Leader_Follower::LEADER ) )
            {
                MORIS_ERROR( false, "IWG_Isotropic_Struc_Nonlinear_Contact_Seitz::compute_jacobian - state dependent stabilization not implemented." );
            }
        }

        // compute the Jacobian for indirect dof dependencies through follower constitutive models
        const uint tFollowerNumDofDependencies = mRequestedFollowerGlobalDofTypes.size();

        for ( uint iDOF = 0; iDOF < tFollowerNumDofDependencies; iDOF++ )
        {
            // get dof type
            const Vector< MSI::Dof_Type >& tDofType = mRequestedFollowerGlobalDofTypes( iDOF );

            // get the index for the dof type
            const sint tDofDepIndex           = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::FOLLOWER );
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
                // compute current tangential traction
                const Matrix< DDRMat > tTangentTraction = mGapData->mLeaderTangentProjection * tTraction;

                if ( tAugLagrTerm > 0 )
                {
                    // nothing here
                }
                else
                {
                    // compute slip traction magnitude
                    real tSplipTraction = -mFrictionCoefficient * tAugLagrTerm;

                    // compute trial tangential traction
                    Matrix< DDRMat > tTrailTangentTraction = tGammaTau * mGapData->mLeaderTangentVelocity;

                    Matrix< DDRMat > tdTrailTangentTractiondv = tGammaTau * mGapData->mLeaderdTangentVelocitydvprev;

                    // Eq. (33) uses current bulk traction. Previous configurations enter only the slip increment.
                    tTrailTangentTraction += tTangentTraction;

                    Matrix< DDRMat > tProjectedTangentTraction;
                    Matrix< DDRMat > tdProjectedTangentTractiondv;

                    real tTrailTangentTractionNorm = norm( tTrailTangentTraction );

                    if ( tTrailTangentTractionNorm > tSplipTraction )
                    {
                        tProjectedTangentTraction = tSplipTraction / tTrailTangentTractionNorm * tTrailTangentTraction;

                        uint tSpaceDim = tTrailTangentTraction.n_rows();

                        tdProjectedTangentTractiondv =
                                +tSplipTraction / std::pow( tTrailTangentTractionNorm, 3.0 ) *                           //
                                ( tTrailTangentTractionNorm * tTrailTangentTractionNorm * eye( tSpaceDim, tSpaceDim )    //
                                        - tTrailTangentTraction * trans( tTrailTangentTraction ) )
                                * tdTrailTangentTractiondv;
                    }
                    else
                    {
                        tProjectedTangentTraction = tTrailTangentTraction;

                        tdProjectedTangentTractiondv = tdTrailTangentTractiondv;
                    }

                    tJacMS += 0.5 * aWStar * (                                                //
                                      -tFILeader->N_trans() * tdProjectedTangentTractiondv    //
                                      + mTheta / tNitscheTangentParam * tTestTraction * tdProjectedTangentTractiondv );

                    tJacSS += 0.5 * aWStar * (    //
                                      +tFIFollower->N_trans() * tdProjectedTangentTractiondv );
                }
            }

            // if dependency on the dof type: even if so, nothing here

            // if dependency of stabilization parameters on the dof type
            if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Leader_Follower::FOLLOWER ) )
            {
                MORIS_ERROR( false, "IWG_Isotropic_Struc_Nonlinear_Contact_Seitz::compute_jacobian - state dependent stabilization not implemented." );
            }
        }

        // check for nan, infinity
        MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                "IWG_Isotropic_Struc_Nonlinear_Contact_Seitz::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
    }

    //------------------------------------------------------------------------------

    void IWG_Isotropic_Struc_Nonlinear_Contact_Seitz::debug_function()
    {
        switch ( mDebugFlag )
        {
            case 1:
            {
                // get leader index for residual dof type, indices for assembly
                Vector< MSI::Dof_Type > const tDisplDofTypes = mResidualDofType( 0 );

                sint tNiter = (sint)gLogger.get_iteration( "NonLinearAlgorithm", "Newton", "Solve", true );

                Matrix< DDRMat > tRayCastPoint = mLeaderFIManager->get_field_interpolators_for_type( tDisplDofTypes( 0 ) )->val()    //
                                               + trans( mLeaderFIManager->get_IP_geometry_interpolator()->valx() );
                Matrix< DDRMat > tTgtPoint = mFollowerFIManager->get_field_interpolators_for_type( tDisplDofTypes( 0 ) )->val()    //
                                           + trans( mFollowerFIManager->get_IP_geometry_interpolator()->valx() );

                const std::shared_ptr< Constitutive_Model >& tConstitutiveModelLeader =
                        mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

                // get the Nitsche stabilization parameter
                const std::shared_ptr< Stabilization_Parameter >& tSPNitsche =
                        mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );

                // get Piola traction
                const Matrix< DDRMat > tTraction =
                        mTractionScaling * tConstitutiveModelLeader->traction( mGapData->mLeaderRefNormal, mCMFunctionType );

                // get Cauchy traction
                const Matrix< DDRMat > tCauchyTraction =
                        tConstitutiveModelLeader->traction( mGapData->mLeaderNormal, CM_Function_Type::CAUCHY );

                // compute contact pressure using current normal
                real tContactPressure = dot( tTraction, mGapData->mLeaderNormal );

                // Nitsche parameter gamma
                const real tNitscheParam = tSPNitsche->val()( 0 );

                // compute augmented Lagrangian term
                const real tAugLagrTerm = tContactPressure + tNitscheParam * mGapData->mGap;

                fprintf( stdout, "Niter = %d Seitz %e  %e  %e  %e  %e  %e  %e\n",    //
                        tNiter,
                        tRayCastPoint( 0 ),
                        tRayCastPoint( 1 ),
                        tTgtPoint( 0 ),
                        tTgtPoint( 1 ),
                        tContactPressure,
                        tAugLagrTerm,
                        mGapData->mGap );
                fprintf( stdout, "Niter = %d NormalSeitz %e  %e  %e  %e  %e  %e\n",    //
                        tNiter,
                        tRayCastPoint( 0 ),
                        tRayCastPoint( 1 ),
                        mGapData->mLeaderNormal( 0 ),
                        mGapData->mLeaderNormal( 1 ),
                        mGapData->mLeaderRefNormal( 0 ),
                        mGapData->mLeaderRefNormal( 1 ) );

                fprintf( stdout, "Niter = %d TractionSeitz %e  %e  %e  %e  %e  %e\n",    //
                        tNiter,
                        tRayCastPoint( 0 ),
                        tRayCastPoint( 1 ),
                        tTraction( 0 ),
                        tTraction( 1 ),
                        tCauchyTraction( 0 ),
                        tCauchyTraction( 1 ) );
                break;
            }
            case 2:
            {
                Matrix< DDRMat > tCheckPoint =
                        mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->valx();

                if ( std::abs( tCheckPoint( 1 ) - 1.0 ) < 0.005 && tCheckPoint( 0 ) > 2.253 && tCheckPoint( 0 ) < 2.374 )
                {
                    // check whether contact should be enforced
                    const std::shared_ptr< Property >& tSelectLeader =
                            mLeaderProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

                    const std::shared_ptr< Property >& tSelectFollower =
                            mFollowerProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

                    real tImposeContactLeader   = ( tSelectLeader != nullptr ) ? tSelectLeader->val()( 0 ) : 1.0;
                    real tImposeContactFollower = ( tSelectFollower != nullptr ) ? tSelectFollower->val()( 0 ) : 1.0;

                    real tvalueLeader   = mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::TEMP )->val()( 0 );
                    real tvalueFollower = mFollowerFIManager->get_field_interpolators_for_type( MSI::Dof_Type::TEMP )->val()( 0 );

                    fprintf( stdout, "valueLeader = %f  valueFollower = %f  tImposeContactLeader = %f  tImposeContactFollower = %f\n",    //
                            tvalueLeader,
                            tvalueFollower,
                            tImposeContactLeader,
                            tImposeContactFollower );
                }
                //
                //            const real tDelta = 1.0e-3;    // tolerance for checking the IG geometry interpolator
                //
                //            if ( std::abs( tCheckPoint( 0 ) - tRefPoint( 0 ) ) < tDelta && std::abs( tCheckPoint( 1 ) - tRefPoint( 1 ) ) < tDelta )
                //            {
                //                real tTimeWeightFactor = gLogger.get_action_data( "NonLinearAlgorithm", "Newton", "Solve", "LoadFactor" );
                //
                //
                //                fprintf( stdout, "tCheckPoint: %f  %f - %f - tLeaderDisp: %e  %e\n",    //
                //                        tCheckPoint( 0 ),
                //                        tCheckPoint( 1 ),
                //                        tTimeWeightFactor,
                //                        tLeaderDisp( 0 ),
                //                        tLeaderDisp( 1 ) );
                //                fprintf( stdout, "tCheckPoint: %f  %f - %f - tFollowerDisp: %e  %e\n",    //
                //                        tCheckPoint( 0 ),
                //                        tCheckPoint( 1 ),
                //                        tTimeWeightFactor,
                //                        tFollowerDisp( 0 ),
                //                        tFollowerDisp( 1 ) );
                //            }
                break;
            }
            default:
                MORIS_ERROR( false, "IWG_Isotropic_Struc_Nonlinear_Contact_Seitz::debug_function - Invalid debug flag." );
        }
    }

    //------------------------------------------------------------------------------

}    // namespace moris::fem
