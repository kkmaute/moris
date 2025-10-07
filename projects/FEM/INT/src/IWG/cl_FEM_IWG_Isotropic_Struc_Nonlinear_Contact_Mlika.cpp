/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Nonlinear_Contact_Mlika.cpp
 *
 */

#include "cl_FEM_IWG_Isotropic_Struc_Nonlinear_Contact_Mlika.hpp"
#include "armadillo"
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "cl_FEM_Stabilization_Parameter.hpp"
#include "cl_MSI_Dof_Type_Enums.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_Vector.hpp"
#include "fn_assert.hpp"
#include "cl_Matrix_Arma_Dynamic.hpp"
#include "fn_trans.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include "fn_dot.hpp"
#include <string>
#include <memory>
#include <utility>

namespace moris::fem
{
    IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::IWG_Isotropic_Struc_Nonlinear_Contact_Mlika(
            sint             aTheta,
            CM_Function_Type aCMFunctionType )
            : mTheta( aTheta )    // switch for symmetric/unsymmetric/neutral Nitsche
            , mCMFunctionType( aCMFunctionType )
    {
        // set time continuity flag
        mTimeContinuity = true;

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
        mUseDeformedGeometryForGap           = true;
        mUseConsistentDeformedGeometryForGap = true;

        // scaling factor to enable pure penalty formulation
        mTractionScaling = 1.0;
    }

    //------------------------------------------------------------------------------

    void
    IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::set_parameters( const Vector< Matrix< DDRMat > >& aParameters )
    {
        // set parameters
        mParameters = aParameters;

        // Read consistency and pure penalty flag from input parameters
        if ( mParameters.size() > 0 )
        {
            MORIS_ERROR( mParameters.size() == 3 &&                  //
                                 mParameters( 0 ).numel() == 1 &&    //
                                 mParameters( 1 ).numel() == 1 &&    //
                                 mParameters( 2 ).numel() == 1,
                    "IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::IWG_Isotropic_Struc_Nonlinear_Contact_Mlika - Only three parameters possible." );

            // parameter 0: consistency flag
            mUseConsistentDeformedGeometryForGap = mParameters( 0 )( 0 ) > 0.5 ? true : false;

            // parameter 1: pure penalty flag
            // create penalty only formulation
            if ( mParameters( 1 )( 0 ) > 0.5 )
            {
                mTractionScaling = 0.0;    // set scaling of all traction related terms to zero
                mTheta           = 0.0;    // theta overwritten to eliminate adjoint term for pure penalty formulation
            }

            mDebugFlag = mParameters( 2 )( 0 );
        }
    }

    //------------------------------------------------------------------------------

    /**
     * \brief Implementation of this IWG is mainly based on
     * 'Mlika, Rabii. “Nitsche Method for Frictional Contact and Self-Contact: Mathematical and Numerical Study.” PhD Thesis, Universite de Lyon, 2018. https://theses.hal.science/tel-02067118.'
     * @attention The implementation is leader-oriented (integration happens on leader side) while Mlika (and the literature in general) integrates on the follower side. Thus, the coordinates X will be on the leader side and Y will be on the follower side (in the literature, it is the other way around).
     * \param aWStar
     */
    void
    IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::compute_residual( real aWStar )
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
        const Matrix< DDRMat > tRemappedFollowerCoords = this->remap_nonconformal_rays(
                mUseDeformedGeometryForGap,
                mUseConsistentDeformedGeometryForGap,
                tDisplDofTypes,
                mLeaderFIManager,
                mFollowerFIManager,
                mGapData );

        // check whether the remapping is successful
        if ( std::abs( tRemappedFollowerCoords( 0 ) ) > 1 )
        {
            return;    // exit if remapping was not successful
        }

        // set integration point for follower side
        mFollowerFIManager->set_space_time_from_local_IG_point( tRemappedFollowerCoords );


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

        // get Piola traction
        const Matrix< DDRMat > tTraction =
                mTractionScaling * tConstitutiveModelLeader->traction( mGapData->mLeaderRefNormal, mCMFunctionType );

        // get test traction
        const Matrix< DDRMat >& tTestTraction =
                tConstitutiveModelLeader->testTraction_trans( mGapData->mLeaderRefNormal, tDisplDofTypes, mCMFunctionType );

        // compute contact pressure using current normal
        const real tContactPressure = mTractionScaling * dot( tTraction, mGapData->mLeaderNormal );

        // compute variation of contact pressure
        Matrix< DDRMat > tTestContactPressuredUleader =
                mTractionScaling * ( tTestTraction * mGapData->mLeaderNormal + trans( mGapData->mLeaderdNormaldu ) * tTraction );

        // compute augmented Lagrangian term
        const real tAugLagrTerm = tContactPressure + tNitscheParam * mGapData->mGap;

        if ( tAugLagrTerm > 0 )
        {
            mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex } ) +=
                    -0.5 * aWStar * mTheta / tNitscheParam * tTestContactPressuredUleader * tContactPressure;
        }
        else
        {
            mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex } ) +=                  //
                    0.5 * aWStar * (                                                    //
                            trans( mGapData->mdGapdu ) * tContactPressure               //
                            + mTheta * tTestContactPressuredUleader * mGapData->mGap    //
                            + tNitscheParam * trans( mGapData->mdGapdu ) * mGapData->mGap );

            mSet->get_residual()( 0 )(
                    { tFollowerResStartIndex, tFollowerResStopIndex } ) +=    //
                    0.5 * aWStar * (                                          //
                            trans( mGapData->mdGapdv ) * tContactPressure     //
                            + tNitscheParam * trans( mGapData->mdGapdv ) * mGapData->mGap );
        }

        // frictional mechanics -> tangential contributions
        const real tFrictionCoefficient = 0.0; // hardcoded for now
        if ( tFrictionCoefficient > 0.0 )
        {
            // make sure time continuity is set
            MORIS_ERROR( get_time_continuity(),
                    "IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::compute_residual - time continuity flag not set.\n" );

            // LEADER

            // get residual dof type field interpolator for current time step
            Field_Interpolator* tLeaderCurrentFI = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            MORIS_ERROR( tLeaderCurrentFI != nullptr,
                    "IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::compute_residual - previous time step field interpolator manager is not set.\n" );

            // get residual dof type field interpolator for previous time step
            Field_Interpolator* tLeaderPreviousFI = mLeaderPreviousFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            //  store field manager of previous time step with field manager of current time step (previous state might be used in property)
            mLeaderFIManager->set_field_interpolator_manager_previous( mLeaderPreviousFIManager );

            // FOLLOWER

            // get residual dof type field interpolator for current time step
            Field_Interpolator* tFollowerCurrentFI = mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            MORIS_ERROR( tFollowerCurrentFI != nullptr,
                    "IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::compute_residual - previous time step field interpolator manager is not set.\n" );

            //  get residual dof type field interpolator for previous time step
            //Field_Interpolator* tFollowerPreviousFI = mFollowerPreviousFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            //  store field manager of previous time step with field manager of current time step (previous state might be used in property)
            mFollowerFIManager->set_field_interpolator_manager_previous( mFollowerPreviousFIManager );

            // for the Leader Previous FI manager, the space time is set from the local IG point (see above)
            // it seems to be necessary for the Follower Previous FI manager as well
            mFollowerPreviousFIManager->set_space_time_from_local_IG_point( tRemappedFollowerCoords );

            // FIXME: set initial time
            real tInitTime = 0.0;

            // initialize jump
            Matrix< DDRMat > tJump = tLeaderCurrentFI->val();
            //std::cout << "current disp: " << tLeaderCurrentFI->val() << std::endl;

            if ( mLeaderFIManager->get_IP_geometry_interpolator()->valt()( 0 ) > tInitTime )
            {
                tJump -= tLeaderPreviousFI->val();
                //std::cout << "previous disp: " << tLeaderPreviousFI->val() << std::endl;
            }

            //std::cout << "Time: " << mLeaderFIManager->get_IP_geometry_interpolator()->valt() << "; Displacement jump: " << tJump << std::endl;
        }

        // call debug function
        if ( mDebugFlag > 0 )
        {
            this->debug_function();
        }

        // check for nan, infinity
        MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                "IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::compute_residual - Residual contains NAN or INF, exiting!" );
    }

    //------------------------------------------------------------------------------

    void IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::compute_jacobian( real aWStar )
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
        const Matrix< DDRMat > tRemappedFollowerCoords = this->remap_nonconformal_rays(
                mUseDeformedGeometryForGap,
                mUseConsistentDeformedGeometryForGap,
                tDisplDofTypes,
                mLeaderFIManager,
                mFollowerFIManager,
                mGapData );

        // check whether the remapping is successful
        if ( std::abs( tRemappedFollowerCoords( 0 ) ) > 1 )
        {
            return;    // exit if remapping was not successful
        }

        // set integration point for follower side
        mFollowerFIManager->set_space_time_from_local_IG_point( tRemappedFollowerCoords );

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

        // get Piola traction
        const Matrix< DDRMat > tTraction =
                mTractionScaling * tConstitutiveModelLeader->traction( mGapData->mLeaderRefNormal, mCMFunctionType );

        // get test traction
        const Matrix< DDRMat >& tTestTraction =
                tConstitutiveModelLeader->testTraction_trans( mGapData->mLeaderRefNormal, tDisplDofTypes, mCMFunctionType );

        // compute contact pressure using current normal
        const real tContactPressure = mTractionScaling * dot( tTraction, mGapData->mLeaderNormal );

        // compute variation of contact pressure
        const Matrix< DDRMat > tTestContactPressuredUleader =
                mTractionScaling * ( tTestTraction * mGapData->mLeaderNormal + trans( mGapData->mLeaderdNormaldu ) * tTraction );

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
                if ( tAugLagrTerm > 0 )
                {
                    tJacMM += -0.5 * aWStar * mTheta / tNitscheParam * (                                         //
                                      tContactPressure * tTestTraction * mGapData->mLeaderdNormaldu              //
                                      + tContactPressure * mGapData->multiply_leader_dnormal2du2( tTraction )    //
                                      + tTestContactPressuredUleader * trans( tTraction ) * mGapData->mLeaderdNormaldu );
                }
                else
                {
                    tJacMM +=
                            0.5 * aWStar * (                                                                          //
                                    tContactPressure * mGapData->mdGap2du2                                            //
                                    + trans( mGapData->mdGapdu ) * trans( tTraction ) * mGapData->mLeaderdNormaldu    //
                                    + mTheta * mGapData->mGap * mGapData->multiply_leader_dnormal2du2( tTraction )    //
                                    + mTheta * mGapData->mGap * tTestTraction * mGapData->mLeaderdNormaldu            //
                                    + mTheta * tTestContactPressuredUleader * mGapData->mdGapdu                       //
                                    + tNitscheParam * ( trans( mGapData->mdGapdu ) * mGapData->mdGapdu                //
                                                        + mGapData->mGap * mGapData->mdGap2du2 ) );

                    tJacSM +=
                            0.5 * aWStar * (                                                                          //
                                    tContactPressure * trans( mGapData->mdGap2duv )                                   //
                                    + trans( mGapData->mdGapdv ) * trans( tTraction ) * mGapData->mLeaderdNormaldu    //
                                    + tNitscheParam * ( trans( mGapData->mdGapdv ) * mGapData->mdGapdu                //
                                                        + mGapData->mGap * trans( mGapData->mdGap2duv ) ) );
                }
            }

            // if dependency on the dof type
            if ( tConstitutiveModelLeader->check_dof_dependency( tDofType ) )
            {
                if ( tAugLagrTerm > 0 )
                {
                    tJacMM += -0.5 * aWStar * mTheta / tNitscheParam * (                                                                                                                                  //
                                      tConstitutiveModelLeader->dTestTractiondDOF( tDofType, mGapData->mLeaderRefNormal, tContactPressure * mGapData->mLeaderNormal, tDisplDofTypes, mCMFunctionType )    //
                                      + tContactPressure * trans( mGapData->mLeaderdNormaldu ) * tConstitutiveModelLeader->dTractiondDOF( tDofType, mGapData->mLeaderRefNormal, mCMFunctionType )         //
                                      + tTestContactPressuredUleader * trans( mGapData->mLeaderNormal ) * tConstitutiveModelLeader->dTractiondDOF( tDofType, mGapData->mLeaderRefNormal, mCMFunctionType ) );
                }
                else
                {
                    tJacMM +=
                            0.5 * aWStar * (                                                                                                                                                                               //
                                    mTractionScaling * trans( mGapData->mdGapdu ) * trans( mGapData->mLeaderNormal ) * tConstitutiveModelLeader->dTractiondDOF( tDofType, mGapData->mLeaderRefNormal, mCMFunctionType )    //
                                    + mTheta * tConstitutiveModelLeader->dTestTractiondDOF( tDofType, mGapData->mLeaderRefNormal, mGapData->mGap * mGapData->mLeaderNormal, tDisplDofTypes, mCMFunctionType )              //
                                    + mTheta * mGapData->mGap * trans( mGapData->mLeaderdNormaldu ) * tConstitutiveModelLeader->dTractiondDOF( tDofType, mGapData->mLeaderRefNormal, mCMFunctionType ) );

                    tJacSM +=
                            0.5 * aWStar * (    //
                                    mTractionScaling * trans( mGapData->mdGapdv ) * trans( mGapData->mLeaderNormal ) * tConstitutiveModelLeader->dTractiondDOF( tDofType, mGapData->mLeaderRefNormal, mCMFunctionType ) );
                }
            }

            // if dependency of stabilization parameters on the dof type
            if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Leader_Follower::LEADER ) )
            {
                MORIS_ERROR( false, "IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::compute_jacobian - state dependent stabilization not implemented." );
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
                            0.5 * aWStar * (                                                              //
                                    tContactPressure * mGapData->mdGap2duv                                //
                                    + mTheta * tTestContactPressuredUleader * mGapData->mdGapdv           //
                                    + tNitscheParam * ( trans( mGapData->mdGapdu ) * mGapData->mdGapdv    //
                                                        + mGapData->mGap * mGapData->mdGap2duv ) );
                    tJacSS +=
                            0.5 * aWStar * (                                                              //
                                    tContactPressure * mGapData->mdGap2dv2                                //
                                    + tNitscheParam * ( trans( mGapData->mdGapdv ) * mGapData->mdGapdv    //
                                                        + mGapData->mGap * trans( mGapData->mdGap2dv2 ) ) );
                }
            }

            // if dependency on the dof type: even if so, nothing here

            // if dependency of stabilization parameters on the dof type
            if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Leader_Follower::FOLLOWER ) )
            {
                MORIS_ERROR( false, "IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::compute_jacobian - state dependent stabilization not implemented." );
            }
        }

        // check for nan, infinity
        MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                "IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
    }

    //------------------------------------------------------------------------------

    void IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::compute_jacobian_and_residual( real aWStar )
    {
        MORIS_ERROR( false, "IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::compute_jacobian_and_residual - This function does nothing." );
    }

    //------------------------------------------------------------------------------

    void IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::compute_dRdp( real aWStar )
    {
        MORIS_ERROR( false, "IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::compute_dRdp - This function does nothing." );
    }

    //------------------------------------------------------------------------------

    void
    IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::debug_function()
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
                const Matrix< DDRMat >& tCauchyTraction =
                        tConstitutiveModelLeader->traction( mGapData->mLeaderNormal, CM_Function_Type::CAUCHY );

                // compute contact pressure using current normal
                real tContactPressure = dot( tTraction, mGapData->mLeaderNormal );

                // Nitsche parameter gamma
                const real tNitscheParam = tSPNitsche->val()( 0 );

                // compute augmented Lagrangian term
                const real tAugLagrTerm = tContactPressure + tNitscheParam * mGapData->mGap;

                fprintf( stdout, "Niter = %d Mlika %e  %e  %e  %e  %e  %e  %e\n",    //
                        tNiter,
                        tRayCastPoint( 0 ),
                        tRayCastPoint( 1 ),
                        tTgtPoint( 0 ),
                        tTgtPoint( 1 ),
                        tContactPressure,
                        tAugLagrTerm,
                        mGapData->mGap );
                fprintf( stdout, "Niter = %d NormalMlika %e  %e  %e  %e  %e  %e\n",    //
                        tNiter,
                        tRayCastPoint( 0 ),
                        tRayCastPoint( 1 ),
                        mGapData->mLeaderNormal( 0 ),
                        mGapData->mLeaderNormal( 1 ),
                        mGapData->mLeaderRefNormal( 0 ),
                        mGapData->mLeaderRefNormal( 1 ) );

                fprintf( stdout, "Niter = %d TractionMlika %e  %e  %e  %e  %e  %e\n",    //
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
                //            const Matrix< DDRMat >& tRefPoint = { { 0.000012, -0.004000 } };
                //
                //            const real tDelta = 1.0e-3;    // tolerance for checking the IG geometry interpolator
                //
                //            if ( std::abs( tCheckPoint( 0 ) - tRefPoint( 0 ) ) < tDelta && std::abs( tCheckPoint( 1 ) - tRefPoint( 1 ) ) < tDelta )
                //            {
                //                real tTimeWeightFactor = gLogger.get_action_data( "NonLinearAlgorithm", "Newton", "Solve", "LoadFactor" );
                //
                //                Matrix< DDRMat > tLeaderDisp   = mLeaderFIManager->get_field_interpolators_for_type( tDisplDofTypes( 0 ) )->val();
                //                Matrix< DDRMat > tFollowerDisp = mFollowerFIManager->get_field_interpolators_for_type( tDisplDofTypes( 0 ) )->val();
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
                MORIS_ERROR( false, "IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::debug_function - Invalid debug flag." );
        }
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
