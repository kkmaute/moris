/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Normal_Nitsche_Unbiased_Unbiased.cpp
 *
 */

#include "cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Normal_Nitsche_Unbiased.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Stabilization_Parameter.hpp"
#include "fn_assert.hpp"
#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"
#include <memory>
#include <utility>

namespace moris::fem
{

    //------------------------------------------------------------------------------

    IWG_Isotropic_Struc_Linear_Contact_Normal_Nitsche_Unbiased::IWG_Isotropic_Struc_Linear_Contact_Normal_Nitsche_Unbiased( sint aBeta )
    {
        // sign for symmetric/unsymmetric Nitsche
        mBeta = aBeta;

        // set size for the property pointer cell
        mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "Thickness" ] = static_cast< uint >( IWG_Property_Type::THICKNESS );
        mPropertyMap[ "Gap" ]       = static_cast< uint >( IWG_Property_Type::GAP );

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
    }

    //------------------------------------------------------------------------------

    void
    IWG_Isotropic_Struc_Linear_Contact_Normal_Nitsche_Unbiased::compute_residual( real aWStar )
    {
#ifdef MORIS_HAVE_DEBUG
        // check leader and follower field interpolators
        this->check_field_interpolators( mtk::Leader_Follower::LEADER );
        this->check_field_interpolators( mtk::Leader_Follower::FOLLOWER );
#endif

        // get leader index for residual dof type, indices for assembly
        Vector< MSI::Dof_Type > const tDisplDofTypes = mResidualDofType( 0 );

        // get leader index for residual dof type, indices for assembly
        uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
        uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
        uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

        // get follower index for residual dof type, indices for assembly
        uint tFollowerDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );
        uint tFollowerResStartIndex = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 0 );
        uint tFollowerResStopIndex  = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 1 );

        // get field interpolator for the residual dof type
        Field_Interpolator*    tLeaderDofs     = mLeaderFIManager->get_field_interpolators_for_type( tDisplDofTypes( 0 ) );
        Geometry_Interpolator* tLeaderGeometry = mLeaderFIManager->get_IG_geometry_interpolator();

        // get follower field interpolator for the residual dof type
        Field_Interpolator*    tFollowerDofs     = mFollowerFIManager->get_field_interpolators_for_type( tDisplDofTypes( 0 ) );
        Geometry_Interpolator* tFollowerGeometry = mFollowerFIManager->get_IG_geometry_interpolator();

        // get leader field interpolator for the residual dof type
        Field_Interpolator* tFILeader = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

        // get follower field interpolator for the residual dof type
        Field_Interpolator* tFIFollower = mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

        // get user defined constitutive model, stabilization parameter and thickness property
        const std::shared_ptr< Constitutive_Model >& tConstitutiveModel =
                mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

        const std::shared_ptr< Stabilization_Parameter >& tStabilizationParameter =
                mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );

        const real tNitscheParam = tStabilizationParameter->val()( 0 );

        const std::shared_ptr< Property >& tThicknessProperty =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

        // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
        aWStar *= ( tThicknessProperty != nullptr ) ? tThicknessProperty->val()( 0 ) : 1;

        // initial gap
        const real tGap = dot( mNormal, tFollowerGeometry->valx() - tLeaderGeometry->valx() );

        // compute test functions in normal directions
        const Matrix< DDRMat > tTestULeader   = trans( mNormal ) * tFILeader->N();
        const Matrix< DDRMat > tTestUFollower = trans( mNormal ) * tFIFollower->N();

        // compute penetration, i.e. negative of gap in deformed configuration
        const real tPenetration = dot( mNormal, tLeaderDofs->val() - tFollowerDofs->val() ) - tGap;

        // evaluate contact pressure on leader side
        const real tPressure = dot( mNormal, tConstitutiveModel->traction( mNormal ) );

        // evaluate test contact pressure on leader side
        const Matrix< DDRMat > tTestPressure =
                tConstitutiveModel->testTraction_trans( mNormal, tDisplDofTypes ) * mNormal;

        // check for contact on leader side
        if ( tPressure - tNitscheParam * tPenetration <= 0 )
        {
            // compute leader residual
            mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex } ) +=    //
                    0.5 * aWStar * (                                      //
                            -trans( tTestULeader ) * tPressure            //
                            + mBeta * tTestPressure * tPenetration        //
                            + tNitscheParam * trans( tTestULeader ) * tPenetration );

            mSet->get_residual()( 0 )(
                    { tFollowerResStartIndex, tFollowerResStopIndex } ) +=    //
                    0.5 * aWStar * (                                          //
                            +trans( tTestUFollower ) * tPressure              //
                            - tNitscheParam * trans( tTestUFollower ) * tPenetration );
        }
        else
        {
            mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex } ) +=    //
                    aWStar * 0.5 * (                                      //
                            -mBeta / tNitscheParam * tTestPressure * tPressure );
        }

        // check for nan, infinity
        MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                "WG_Isotropic_Struc_Linear_Contact_Normal_Nitsche_Unbiased::compute_residual - Residual contains NAN or INF, exiting!" );
    }

    //------------------------------------------------------------------------------

    void
    IWG_Isotropic_Struc_Linear_Contact_Normal_Nitsche_Unbiased::compute_jacobian( real aWStar )
    {
#ifdef MORIS_HAVE_DEBUG
        // check leader and follower field interpolators
        this->check_field_interpolators( mtk::Leader_Follower::LEADER );
        this->check_field_interpolators( mtk::Leader_Follower::FOLLOWER );
#endif

        // get leader index for residual dof type, indices for assembly
        Vector< MSI::Dof_Type > const tDisplDofTypes = mResidualDofType( 0 );

        // get leader index for residual dof type, indices for assembly
        uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
        uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
        uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

        // get follower index for residual dof type, indices for assembly
        uint tFollowerDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );
        uint tFollowerResStartIndex = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 0 );
        uint tFollowerResStopIndex  = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 1 );

        // get field interpolator for the residual dof type
        Field_Interpolator*    tLeaderDofs     = mLeaderFIManager->get_field_interpolators_for_type( tDisplDofTypes( 0 ) );
        Geometry_Interpolator* tLeaderGeometry = mLeaderFIManager->get_IG_geometry_interpolator();

        // get follower field interpolator for the residual dof type
        Field_Interpolator*    tFollowerDofs     = mFollowerFIManager->get_field_interpolators_for_type( tDisplDofTypes( 0 ) );
        Geometry_Interpolator* tFollowerGeometry = mFollowerFIManager->get_IG_geometry_interpolator();

        // get leader field interpolator for the residual dof type
        Field_Interpolator* tFILeader = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

        // get follower field interpolator for the residual dof type
        Field_Interpolator* tFIFollower = mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

        // get user defined constitutive model, stabilization parameter and thickness property
        const std::shared_ptr< Constitutive_Model >& tConstitutiveModel =
                mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

        const std::shared_ptr< Stabilization_Parameter >& tStabilizationParameter =
                mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );

        const real tNitscheParam = tStabilizationParameter->val()( 0 );

        const std::shared_ptr< Property >& tThicknessProperty =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

        // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
        aWStar *= ( tThicknessProperty != nullptr ) ? tThicknessProperty->val()( 0 ) : 1;

        // initial gap
        const real tGap = dot( mNormal, tFollowerGeometry->valx() - tLeaderGeometry->valx() );

        // compute test functions in normal directions
        const Matrix< DDRMat > tTestULeader   = trans( mNormal ) * tFILeader->N();
        const Matrix< DDRMat > tTestUFollower = trans( mNormal ) * tFIFollower->N();

        // compute penetration, i.e. negative of gap in deformed configuration
        const real tPenetration = dot( mNormal, tLeaderDofs->val() - tFollowerDofs->val() ) - tGap;

        // evaluate contact pressure on leader side
        const real tPressure = dot( mNormal, tConstitutiveModel->traction( mNormal ) );

        // evaluate test contact pressure on leader side
        const Matrix< DDRMat > tTestPressure =
                tConstitutiveModel->testTraction_trans( mNormal, tDisplDofTypes ) * mNormal;

        // get number of leader dof dependencies
        const uint tLeaderNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

        // compute the Jacobian for indirect dof dependencies through leader constitutive models
        for ( uint iDOF = 0; iDOF < tLeaderNumDofDependencies; iDOF++ )
        {
            // get the dof type
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
            if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
            {
                if ( tPressure - tNitscheParam * tPenetration <= 0 )
                {
                    tJacMM += 0.5 * aWStar * (                                 //
                                      +mBeta * tTestPressure * tTestULeader    //
                                      + tNitscheParam * trans( tTestULeader ) * tTestULeader );

                    tJacSM += 0.5 * aWStar * (    //
                                      -tNitscheParam * trans( tTestUFollower ) * tTestULeader );
                }
            }
            // if dependency on the dof type
            if ( tConstitutiveModel->check_dof_dependency( tDofType ) )
            {
                if ( tPressure - tNitscheParam * tPenetration <= 0 )
                {
                    // add contribution to Jacobian
                    tJacMM += 0.5 * aWStar * (                                                                                              //
                                      -trans( tTestULeader ) * trans( mNormal ) * tConstitutiveModel->dTractiondDOF( tDofType, mNormal )    //
                                      + mBeta * tConstitutiveModel->dTestTractiondDOF( tDofType, mNormal, tPenetration * mNormal, mResidualDofType( 0 ) ) );

                    tJacSM += 0.5 * aWStar * (    //
                                      +trans( tTestUFollower ) * trans( mNormal ) * tConstitutiveModel->dTractiondDOF( tDofType, mNormal ) );
                }
                else
                {
                    tJacMM += 0.5 * aWStar * (                                                                                             //
                                      -mBeta / tNitscheParam * (                                                                           //
                                              tTestPressure * trans( mNormal ) * tConstitutiveModel->dTractiondDOF( tDofType, mNormal )    //
                                              + tConstitutiveModel->dTestTractiondDOF( tDofType, mNormal, tPressure * mNormal, mResidualDofType( 0 ) ) ) );
                }
            }

            // if dependency of stabilization parameters on the dof type
            if ( tStabilizationParameter->check_dof_dependency( tDofType, mtk::Leader_Follower::LEADER ) )
            {
                MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Contact_Normal_Nitsche_Unbiased::compute_jacobian - state dependent stabilization not implemented." );
            }
        }

        // compute the Jacobian for indirect dof dependencies through follower constitutive models
        uint tFollowerNumDofDependencies = mRequestedFollowerGlobalDofTypes.size();
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
                if ( tPressure - tNitscheParam * tPenetration <= 0 )
                {
                    tJacMS += 0.5 * aWStar * (                                   //
                                      -mBeta * tTestPressure * tTestUFollower    //
                                      - tNitscheParam * trans( tTestULeader ) * tTestUFollower );

                    tJacSS += 0.5 * aWStar * (    //
                                      +tNitscheParam * trans( tTestUFollower ) * tTestUFollower );
                }
            }

            // if dependency of stabilization parameters on the dof type
            if ( tStabilizationParameter->check_dof_dependency( tDofType, mtk::Leader_Follower::FOLLOWER ) )
            {
                MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Contact_Normal_Nitsche_Unbiased::compute_jacobian - state dependent stabilization not implemented." );
            }
        }

        // check for nan, infinity
        MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                "WG_Isotropic_Struc_Linear_Contact_Normal_Nitsche_Unbiased::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
    }

    //------------------------------------------------------------------------------

    void
    IWG_Isotropic_Struc_Linear_Contact_Normal_Nitsche_Unbiased::compute_jacobian_and_residual( real aWStar )
    {
        MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Contact_Normal_Nitsche_Unbiased::compute_jacobian_and_residual - This function does nothing." );
    }

    //------------------------------------------------------------------------------

    void
    IWG_Isotropic_Struc_Linear_Contact_Normal_Nitsche_Unbiased::compute_dRdp( real aWStar )
    {
        MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Contact_Normal_Nitsche_Unbiased::compute_dRdp - This function does nothing." );
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
