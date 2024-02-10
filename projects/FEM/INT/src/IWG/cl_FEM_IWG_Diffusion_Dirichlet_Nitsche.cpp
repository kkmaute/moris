/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Diffusion_Dirichlet_Nitsche.cpp
 *
 */

#include "cl_FEM_IWG_Diffusion_Dirichlet_Nitsche.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Diffusion_Dirichlet_Nitsche::IWG_Diffusion_Dirichlet_Nitsche( sint aBeta )
        {
            // set sign for symmetric/unsymmetric Nitsche
            mBeta = aBeta;

            // set size for the property pointer cell
            mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Dirichlet" ] = static_cast< uint >( IWG_Property_Type::DIRICHLET );
            mPropertyMap[ "Select" ]    = static_cast< uint >( IWG_Property_Type::SELECT );
            mPropertyMap[ "Thickness" ] = static_cast< uint >( IWG_Property_Type::THICKNESS );

            // set size for the constitutive model pointer cell
            mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Diffusion" ] = static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "DirichletNitsche" ] = static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Dirichlet_Nitsche::compute_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check field interpolators, properties, constitutive models
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get field interpolator for a given dof type
            Field_Interpolator* tFI =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the selection matrix property
            const std::shared_ptr< Property >& tPropSelect =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            // set a default selection matrix if needed
            real tM = 1.0;
            if ( tPropSelect != nullptr )
            {
                tM = tPropSelect->val()( 0 );

                // skip computing residual if projection matrix is zero
                if ( std::abs( tM ) < MORIS_REAL_EPS ) return;
            }

            // get imposed temperature property
            const std::shared_ptr< Property >& tPropDirichlet =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) );

            MORIS_ASSERT( tPropDirichlet != nullptr,
                    "IWG_Diffusion_Dirichlet_Nitsche::compute_residual - DBC property missing." );

            // get the elasticity CM
            const std::shared_ptr< Constitutive_Model >& tCMDiffusion =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) );

            MORIS_ASSERT( tCMDiffusion != nullptr,
                    "IWG_Diffusion_Dirichlet_Nitsche::compute_residual - Constitutive model missing." );

            // get the stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE ) );

            MORIS_ASSERT( tSPNitsche != nullptr,
                    "IWG_Diffusion_Dirichlet_Nitsche::compute_residual - Nitsche stabilization parameter missing." );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // compute jump
            real tJump = tFI->val()( 0 ) - tPropDirichlet->val()( 0 );

            // compute the residual
            mSet->get_residual()( 0 )( { tLeaderResStartIndex, tLeaderResStopIndex } ) +=
                    aWStar * ( -tFI->N_trans() * tM * tCMDiffusion->traction( mNormal )                               //
                               + mBeta * tCMDiffusion->testTraction( mNormal, mResidualDofType( 0 ) ) * tM * tJump    //
                               + tSPNitsche->val()( 0 ) * tFI->N_trans() * tM * tJump );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Diffusion_Dirichlet_Nitsche::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Dirichlet_Nitsche::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check field interpolators, properties, constitutive models
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get field interpolator for residual dof type
            Field_Interpolator* tFI =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the selection matrix property
            const std::shared_ptr< Property >& tPropSelect =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            // set a default selection matrix if needed
            real tM = 1.0;
            if ( tPropSelect != nullptr )
            {
                tM = tPropSelect->val()( 0 );

                // skip computing residual if projection matrix is zero
                if ( std::abs( tM ) < MORIS_REAL_EPS ) return;
            }

            // get imposed temperature property
            const std::shared_ptr< Property >& tPropDirichlet =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) );

            // get the elasticity CM
            const std::shared_ptr< Constitutive_Model >& tCMDiffusion =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) );

            // get the stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE ) );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // compute jump
            real tJump = tFI->val()( 0 ) - tPropDirichlet->val()( 0 );

            // compute the Jacobian for indirect dof dependencies through properties
            uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            for ( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the dof type
                const Vector< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // get sub-matrix
                auto tJac = mSet->get_jacobian()(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } );

                // if dof type is residual dof type
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    tJac += aWStar * ( mBeta * tCMDiffusion->testTraction( mNormal, mResidualDofType( 0 ) ) * tM * tFI->N() +    //
                                       tSPNitsche->val()( 0 ) * tFI->N_trans() * tM * tFI->N() );
                }

                // if dependency on the dof type
                if ( tPropDirichlet->check_dof_dependency( tDofType ) )
                {
                    // add contribution to Jacobian
                    tJac += aWStar * ( -mBeta * tCMDiffusion->testTraction( mNormal, mResidualDofType( 0 ) ) * tM * tPropDirichlet->dPropdDOF( tDofType )    //
                                       - tSPNitsche->val()( 0 ) * tFI->N_trans() * tM * tPropDirichlet->dPropdDOF( tDofType ) );
                }

                // if dependency on the dof type
                if ( tCMDiffusion->check_dof_dependency( tDofType ) )
                {
                    // add contribution to Jacobian
                    tJac += aWStar * ( -tFI->N_trans() * tM * tCMDiffusion->dTractiondDOF( tDofType, mNormal )    //
                                       + mBeta * tCMDiffusion->dTestTractiondDOF( tDofType, mNormal, mResidualDofType( 0 ) ) * tM * tJump );
                }

                // if dependency on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType ) )
                {
                    // add contribution to Jacobian
                    tJac += aWStar * ( tFI->N_trans() * tM * tJump * tSPNitsche->dSPdLeaderDOF( tDofType ) );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Diffusion_Dirichlet_Nitsche::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Dirichlet_Nitsche::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Diffusion_Dirichlet_Nitsche::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Dirichlet_Nitsche::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Diffusion_Dirichlet_Nitsche::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
