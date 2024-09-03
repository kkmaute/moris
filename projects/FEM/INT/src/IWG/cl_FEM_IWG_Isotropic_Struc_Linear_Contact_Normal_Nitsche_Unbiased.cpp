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

            uint const tDofIndex      = mSet->get_dof_index_for_type( tDisplDofTypes( 0 ), mtk::Leader_Follower::LEADER );
            uint const tResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint const tResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );
            auto const tResRange      = std::make_pair( tResStartIndex, tResStopIndex );

            // get field interpolator for the residual dof type
            Field_Interpolator*    tLeaderDofs     = mLeaderFIManager->get_field_interpolators_for_type( tDisplDofTypes( 0 ) );
            Geometry_Interpolator* tLeaderGeometry = mLeaderFIManager->get_IG_geometry_interpolator();

            // get follower field interpolator for the residual dof type
            Field_Interpolator*    tFollowerDofs     = mFollowerFIManager->get_field_interpolators_for_type( tDisplDofTypes( 0 ) );
            Geometry_Interpolator* tFollowerGeometry = mFollowerFIManager->get_IG_geometry_interpolator();

            // get user defined constitutive model, stabilization parameter and thickness property
            const std::shared_ptr< Constitutive_Model >&      tConstitutiveModel      = mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );
            const std::shared_ptr< Stabilization_Parameter >& tStabilizationParameter = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );
            const std::shared_ptr< Property >&                tThicknessProperty      = mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tThicknessProperty != nullptr ) ? tThicknessProperty->val()( 0 ) : 1;

            const real tNitscheParam = tStabilizationParameter->val()( 0 );    // stabilization parameter gamma

            const uint             tDim      = mNormal.numel();
            const Matrix< DDRMat > tIdentity = eye( tDim, tDim );

            // test functions of leader (L) and follower (F)
            const Matrix< DDRMat > tNL = tLeaderDofs->N();      // ( 2 x 8 )
            const Matrix< DDRMat > tNF = tFollowerDofs->N();    // ( 2 x 8 )

            // Determining the normal vector on the leader side is done by computing the difference between the deformed coordinates of the leader and follower.
            Matrix< DDRMat > const txL = tLeaderGeometry->valx_current( tLeaderDofs );
            Matrix< DDRMat > const txF = tFollowerGeometry->valx_current( tFollowerDofs );

            Matrix< DDRMat > tNormalLeader = tLeaderGeometry->get_normal_current( tLeaderDofs );

            // normal projection operator
            const Matrix< DDRMat > tNormalProjectorLeader = tNormalLeader * trans( tNormalLeader );

            // initial gap
            const real tGap = norm( tFollowerGeometry->valx() - tLeaderGeometry->valx() );

            // compute the jump
            const Matrix< DDRMat > tJumpLeader = tLeaderDofs->val() - tFollowerDofs->val() - tNormalLeader * tGap;

            // compute projection of displacement jump onto normal
            const real tNormalJumpLeader = dot( tJumpLeader, tNormalLeader );

            // evaluate traction and their pressure equivalent (normal component)
            const Matrix< DDRMat > tTraction       = tConstitutiveModel->traction( tNormalLeader );    // ( 2 x 1 )
            const real             tPressure       = dot( tTraction, tNormalLeader );
            const Matrix< DDRMat > tNormalTraction = tPressure * tNormalLeader;    // ( 2 x 1 ) normal component of traction

            // evaluate test traction and the corresponding pressure equivalent (normal component) by computing the dot product between
            // each column of tTestTraction and tNormalLeader. This is done by element-wise multiplication and then summing the columns.
            const Matrix< DDRMat > tTestTraction       = tConstitutiveModel->testTraction( tNormalLeader, tDisplDofTypes );    // ( 2 x 8 )
            const Matrix< DDRMat > tTestPressure       = trans( tNormalLeader ) * tTestTraction;                               // ( 1 x 8 )
            const Matrix< DDRMat > tNormalTestTraction = tNormalLeader * tTestPressure;                                        // ( 2 x 8 ) normal component of test traction

            // compute contact pressure
            const real tIfcPressureLeader = dot( tTraction, tNormalLeader );

            // check for contact on leader side
            if ( tIfcPressureLeader - tNitscheParam * tNormalJumpLeader < 0 )
            {
                // compute leader residual
                mSet->get_residual()( 0 )( tResRange ) +=
                        aWStar * 0.5 * (                                                                   //
                                -trans( tNL ) * tNormalProjectorLeader * tTraction                         //
                                + mBeta * trans( tTestTraction ) * tNormalProjectorLeader * tJumpLeader    //
                                + tNitscheParam * trans( tNL ) * tNormalProjectorLeader * tJumpLeader      //
                        );
            }
            else
            {
                mSet->get_residual()( 0 )( tResRange ) +=
                        aWStar * 0.5 * ( -mBeta / tNitscheParam * trans( tTestTraction ) * tNormalProjectorLeader * tTraction );
            }

            [[maybe_unused]] uint const tIteration = gLogger.get_iteration( "NonLinearAlgorithm", "Newton", "Solve" );
            const Matrix< DDRMat >      tDebugArea{
                     { 0.55, 0.60 },    // xmin, xmax
                     { 0.5, 0.8 }       // ymin, ymax
            };
            const Matrix< DDRMat > tXL = trans( tLeaderGeometry->valx() );      // ( 2 x 1 )
            const Matrix< DDRMat > tXF = trans( tFollowerGeometry->valx() );    // ( 2 x 1 )
            const Matrix< DDRMat > tUL = tLeaderDofs->val();                    // ( 2 x 1 )
            const Matrix< DDRMat > tUF = tFollowerDofs->val();                  // ( 2 x 1 )
            // check if point tXL is in the debug area
            [[maybe_unused]] const bool tIsDebugPoint = ( tXL( 0 ) > tDebugArea( 0, 0 )
                                                          && tXL( 0 ) < tDebugArea( 0, 1 )
                                                          && tXL( 1 ) > tDebugArea( 1, 0 )
                                                          && tXL( 1 ) < tDebugArea( 1, 1 ) );
            if ( tIsDebugPoint )
            {
                //                std::cout << "IWG:" << std::setprecision( 10 )
                //                          << tIteration << ","
                //                          << tGap << ","
                //                          << pmat1( tXL ) << ","
                //                          << pmat1( tXF ) << ","
                //                          << pmat1( tUL ) << ","
                //                          << pmat1( tUF ) << ","
                //                          << pmat1( txL ) << ","
                //                          << pmat1( txF ) << ","
                //                          << pmat1( tNormalLeader ) << ","
                //                          << std::endl;
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Isotropic_Struc_Linear_Contact_Gap::compute_residual - Residual contains NAN or INF, exiting!" );
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
            const uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            const uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            const uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get follower index for residual dof type, indices for assembly
            const uint tFollowerDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );

            // get leader field interpolator for the residual dof type
            Field_Interpolator*    tFILeader = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );
            Geometry_Interpolator* tGILeader = mLeaderFIManager->get_IG_geometry_interpolator();
            // get follower field interpolator for the residual dof type
            Field_Interpolator*    tFIFollower = mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );
            Geometry_Interpolator* tGIFollower = mFollowerFIManager->get_IG_geometry_interpolator();

            real tGap = norm( tGIFollower->valx() - tGILeader->valx() );

            // get the elasticity constitutive model
            const std::shared_ptr< Constitutive_Model >& tCMLeaderElasticity =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            const std::shared_ptr< Constitutive_Model >& tCMFollowerElasticity =
                    mFollowerCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            const real tNitsche      = tSPNitsche->val()( 0 );
            const real tLeaderWeight = tSPNitsche->val()( 1 );

            // normals on leader and follower side (needs to be generalized)
            const Matrix< DDRMat > tNormalLeader   = mNormal;
            const Matrix< DDRMat > tNormalFollower = -1.0 * mNormal;

            // normal projection operator
            const Matrix< DDRMat > tNormalProjectorLeader   = tNormalLeader * trans( tNormalLeader );
            const Matrix< DDRMat > tNormalProjectorFollower = tNormalFollower * trans( tNormalFollower );

            // compute the jump
            const Matrix< DDRMat > tJumpLeader   = tFILeader->val() - tFIFollower->val() - tNormalLeader * tGap;
            const Matrix< DDRMat > tJumpFollower = tFIFollower->val() - tFILeader->val() - tNormalFollower * tGap;

            // compute projection of displacement jump onto normal
            const real tNormalJumpLeader = dot( tJumpLeader, tNormalLeader );

            // evaluate tractions
            const Matrix< DDRMat > tTractionLeader   = tCMLeaderElasticity->traction( tNormalLeader );
            const Matrix< DDRMat > tTractionFollower = tCMFollowerElasticity->traction( tNormalFollower );

            // compute contact pressure
            const real tIfcPressureLeader = dot( tTractionLeader, tNormalLeader );

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

                // compute Jacobian direct dependencies
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    if ( tIfcPressureLeader - tNitsche * tNormalJumpLeader < 0 )
                    {
                        tJacMM += aWStar * tLeaderWeight * (                                                                                                                    //
                                          +mBeta * tCMLeaderElasticity->testTraction_trans( tNormalLeader, mResidualDofType( 0 ) ) * tNormalProjectorLeader * tFILeader->N()    //
                                          + tNitsche * tFILeader->N_trans() * tNormalProjectorLeader * tFILeader->N() );
                    }
                }

                // if dependency on the dof type
                if ( tCMLeaderElasticity->check_dof_dependency( tDofType ) )
                {
                    if ( tIfcPressureLeader - tNitsche * tNormalJumpLeader < 0 )
                    {
                        // add contribution to Jacobian
                        tJacMM += aWStar * tLeaderWeight * (                                                                                                //
                                          -tFILeader->N_trans() * tNormalProjectorLeader * tCMLeaderElasticity->dTractiondDOF( tDofType, tNormalLeader )    //
                                          + mBeta * tCMLeaderElasticity->dTestTractiondDOF( tDofType, tNormalLeader, tNormalProjectorLeader * tJumpLeader, mResidualDofType( 0 ) ) );
                    }
                    else
                    {
                        tJacMM += aWStar * tLeaderWeight * -mBeta / tNitsche * (                                                                                                                                      //
                                          tCMLeaderElasticity->testTraction_trans( tNormalLeader, mResidualDofType( 0 ) ) * tNormalProjectorLeader * tCMLeaderElasticity->dTractiondDOF( tDofType, tNormalLeader )    //
                                          + tCMLeaderElasticity->dTestTractiondDOF( tDofType, tNormalLeader, tNormalProjectorLeader * tTractionLeader, mResidualDofType( 0 ) ) );
                    }
                }

                // if dependency of stabilization parameters on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Leader_Follower::LEADER ) )
                {
                    // get the derivatives of the SPs
                    const Matrix< DDRMat > tNitscheDer        = tSPNitsche->dSPdLeaderDOF( tDofType ).get_row( 0 );
                    const Matrix< DDRMat > tLeaderWeightDer   = tSPNitsche->dSPdLeaderDOF( tDofType ).get_row( 1 );
                    const Matrix< DDRMat > tFollowerWeightDer = tSPNitsche->dSPdLeaderDOF( tDofType ).get_row( 2 );

                    if ( tIfcPressureLeader - tNitsche * tNormalJumpLeader < 0 )
                    {
                        // add contribution to Jacobian
                        tJacMM += aWStar * (                                                                                                                                          //
                                          ( -tFILeader->N_trans() * tNormalProjectorLeader * tTractionLeader                                                                          //
                                                  + mBeta * tCMLeaderElasticity->testTraction_trans( tNormalLeader, mResidualDofType( 0 ) ) * tNormalProjectorLeader * tJumpLeader    //
                                                  + tNitsche * tFILeader->N_trans() * tNormalProjectorLeader * tJumpLeader )                                                          //
                                                  * tLeaderWeightDer                                                                                                                  //
                                          + tLeaderWeight * tFILeader->N_trans() * tNormalProjectorLeader * tJumpLeader * tNitscheDer );
                    }
                    else
                    {
                        tJacMM += aWStar * (    //
                                          -mBeta * tCMLeaderElasticity->testTraction_trans( tNormalLeader, mResidualDofType( 0 ) ) * tNormalProjectorLeader * tTractionLeader )
                                * ( 1.0 / tNitsche * tLeaderWeightDer - tLeaderWeight / tNitsche / tNitsche * tNitscheDer );
                    }
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

                // if dof type is residual dof type
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    if ( tIfcPressureLeader - tNitsche * tNormalJumpLeader < 0 )
                    {
                        tJacMS += aWStar * tLeaderWeight * (                                                                                                                      //
                                          -mBeta * tCMLeaderElasticity->testTraction_trans( tNormalLeader, mResidualDofType( 0 ) ) * tNormalProjectorLeader * tFIFollower->N()    //
                                          - tNitsche * tFILeader->N_trans() * tNormalProjectorLeader * tFIFollower->N() );
                    }
                }

                // if dependency of stabilization parameters on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Leader_Follower::FOLLOWER ) )
                {
                    // get the derivatives of the SPs
                    const Matrix< DDRMat > tNitscheDer        = tSPNitsche->dSPdFollowerDOF( tDofType ).get_row( 0 );
                    const Matrix< DDRMat > tLeaderWeightDer   = tSPNitsche->dSPdFollowerDOF( tDofType ).get_row( 1 );
                    const Matrix< DDRMat > tFollowerWeightDer = tSPNitsche->dSPdFollowerDOF( tDofType ).get_row( 2 );

                    if ( tIfcPressureLeader - tNitsche * tNormalJumpLeader < 0 )
                    {
                        // add contribution to Jacobian
                        tJacMS += aWStar * ( ( -tFILeader->N_trans() * tNormalProjectorLeader * tTractionLeader + mBeta * tCMLeaderElasticity->testTraction_trans( tNormalLeader, mResidualDofType( 0 ) ) * tNormalProjectorLeader * tJumpLeader + tNitsche * tFILeader->N_trans() * tNormalProjectorLeader * tJumpLeader ) * tLeaderWeightDer + tLeaderWeight * tFILeader->N_trans() * tNormalProjectorLeader * tJumpLeader * tNitscheDer );
                    }
                    else
                    {
                        tJacMS += aWStar * ( -mBeta * tCMLeaderElasticity->testTraction_trans( tNormalLeader, mResidualDofType( 0 ) ) * tNormalProjectorLeader * tTractionLeader )
                                * ( 1.0 / tNitsche * tLeaderWeightDer - tLeaderWeight / tNitsche / tNitsche * tNitscheDer );
                    }
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Isotropic_Struc_Linear_Interface::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
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
