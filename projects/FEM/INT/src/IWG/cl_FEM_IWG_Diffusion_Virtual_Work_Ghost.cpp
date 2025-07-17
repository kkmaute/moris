/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Diffusion_Virtual_Work_Ghost.cpp
 *
 */

#include "cl_FEM_IWG_Diffusion_Virtual_Work_Ghost.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_isfinite.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IWG_Diffusion_Virtual_Work_Ghost::IWG_Diffusion_Virtual_Work_Ghost()
    {
        // set ghost flag
        mIsGhost = true;

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );
        mFollowerCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "Diffusion" ] = static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO );

        // set size for the stabilization parameter pointer cell
        mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

        // populate the stabilization map
        mStabilizationMap[ "GhostVWOrder" ] = static_cast< uint >( IWG_Stabilization_Type::GHOST_VW );
    }

    //------------------------------------------------------------------------------

    void IWG_Diffusion_Virtual_Work_Ghost::compute_residual( real aWStar )
    {
#ifdef MORIS_HAVE_DEBUG
            // check leader and follower field interpolators
            this->check_field_interpolators( mtk::Leader_Follower::LEADER );
            this->check_field_interpolators( mtk::Leader_Follower::FOLLOWER );
#endif
            // set interpolation order
            IWG::set_interpolation_order();

            // get leader index for residual dof type, indices for assembly
            uint tDofIndexLeader      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndexLeader )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndexLeader )( 0, 1 );

            // get follower index for residual dof type, indices for assembly
            uint tDofIndexFollower      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );
            uint tFollowerResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndexFollower )( 0, 0 );
            uint tFollowerResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndexFollower )( 0, 1 );

            // get leader field interpolator for the residual dof type
            Field_Interpolator * tLeaderFI =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get follower field interpolator for the residual dof type
            Field_Interpolator * tFollowerFI  =
                    mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the diffusion constitutive model for leader and follower
            const std::shared_ptr< Constitutive_Model > & tCMLeaderDiff =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) );
            const std::shared_ptr< Constitutive_Model > & tCMFollowerDiff =
                    mFollowerCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) );

            // get the stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPGhost =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::GHOST_VW ) );

            // FIXME get the conductivity for leader and follower
            real tLeaderConductivity = tCMLeaderDiff->constitutive()( 0, 0 );
            real tFollowerConductivity  = tCMFollowerDiff->constitutive()( 0, 0 );

            // get sub-matrices
            auto tResLeader = mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex } );

            auto tResFollower = mSet->get_residual()( 0 )(
                    { tFollowerResStartIndex, tFollowerResStopIndex } );

            // loop over the interpolation order
            for ( uint iOrder = 1; iOrder <= mOrder; iOrder++ )
            {
                // set the order for the stabilization parameter
                tSPGhost->set_interpolation_order( iOrder );

                // get flattened normal matrix
                Matrix< DDRMat > tNormalMatrix;
                this->get_flat_normal_matrix( tNormalMatrix, iOrder );

                // multiply common terms (note: needs to be stored as matrix
                Matrix< DDRMat > tPreMultiply =
                        tSPGhost->val()( 0 ) * trans( tNormalMatrix ) * tNormalMatrix *
                        ( tLeaderConductivity * tLeaderFI->gradx( iOrder ) - tFollowerConductivity * tFollowerFI->gradx( iOrder ) );

                // compute leader residual
                tResLeader += aWStar * (
                        trans( tLeaderFI->dnNdxn( iOrder ) ) * tPreMultiply );

                // compute follower residual
                tResFollower -= aWStar * (
                        trans( tFollowerFI->dnNdxn( iOrder ) ) * tPreMultiply );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Diffusion_Virtual_Work_Ghost::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------
        void IWG_Diffusion_Virtual_Work_Ghost::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader and follower field interpolators
            this->check_field_interpolators( mtk::Leader_Follower::LEADER );
            this->check_field_interpolators( mtk::Leader_Follower::FOLLOWER );
#endif

            // set interpolation order
            IWG::set_interpolation_order();

            // get leader index for residual dof type, indices for assembly
            uint tDofIndexLeader      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndexLeader )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndexLeader )( 0, 1 );

            // get follower index for residual dof type, indices for assembly
            uint tDofIndexFollower      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );
            uint tFollowerResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndexFollower )( 0, 0 );
            uint tFollowerResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndexFollower )( 0, 1 );

            // get leader field interpolator for the residual dof type
            Field_Interpolator * tLeaderFI =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get follower field interpolator for the residual dof type
            Field_Interpolator * tFollowerFI  =
                    mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the diffusion constitutive model for leader and follower
            const std::shared_ptr< Constitutive_Model > & tCMLeaderDiff =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) );
            const std::shared_ptr< Constitutive_Model > & tCMFollowerDiff =
                    mFollowerCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) );

            // get the stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPGhost =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::GHOST_VW ) );

            // FIXME get the conductivity for leader and follower
            real tLeaderConductivity = tCMLeaderDiff->constitutive()( 0, 0 );
            real tFollowerConductivity  = tCMFollowerDiff->constitutive()( 0, 0 );

            // get number of leader dof dependencies
            uint tLeaderNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            // loop over the order
            for ( uint iOrder = 1; iOrder <= mOrder; iOrder++ )
            {
                // set the order for the stabilization parameter
                tSPGhost->set_interpolation_order( iOrder );

                // get flattened normal matrix
                Matrix< DDRMat > tNormalMatrix;
                this->get_flat_normal_matrix( tNormalMatrix, iOrder );

                // compute the Jacobian for indirect dof dependencies through leader constitutive models
                for( uint iDOF = 0; iDOF < tLeaderNumDofDependencies; iDOF++ )
                {
                    // get the dof type
                    const Vector< MSI::Dof_Type > & tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                    // get the index for the dof type
                    sint tIndexDep      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                    uint tDepStartIndex = mSet->get_jac_dof_assembly_map()( tDofIndexLeader )( tIndexDep, 0 );
                    uint tDepStopIndex  = mSet->get_jac_dof_assembly_map()( tDofIndexLeader )( tIndexDep, 1 );

                    auto tJacLeader = mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tDepStartIndex,       tDepStopIndex } );

                    auto tJacFollower =  mSet->get_jacobian()(
                            { tFollowerResStartIndex, tFollowerResStopIndex },
                            { tDepStartIndex,      tDepStopIndex } );

                    // if dependency on the dof type
                    if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                    {
                        // add contribution to Jacobian
                        tJacLeader += aWStar * (
                                        tSPGhost->val()( 0 )
                                        * trans( tLeaderFI->dnNdxn( iOrder ) ) * trans( tNormalMatrix )
                                        * tNormalMatrix * tLeaderConductivity * tLeaderFI->dnNdxn( iOrder ) );

                       tJacFollower -= aWStar * (
                                        tSPGhost->val()( 0 )
                                        * trans( tFollowerFI->dnNdxn( iOrder ) ) * trans( tNormalMatrix )
                                        * tNormalMatrix * tLeaderConductivity * tLeaderFI->dnNdxn( iOrder ) );
                    }

                    // if diffusion CM depends on dof type
                    if( tCMLeaderDiff->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to Jacobian
                        tJacLeader -= aWStar * (
                                        tSPGhost->val()( 0 )
                                        * trans( tLeaderFI->dnNdxn( iOrder ) ) * trans( tNormalMatrix )
                                        * tNormalMatrix * tLeaderFI->gradx( iOrder ) * tCMLeaderDiff->dConstdDOF( tDofType ) );

                        tJacFollower += aWStar * (
                                        tSPGhost->val()( 0 )
                                        * trans( tFollowerFI->dnNdxn( iOrder ) ) * trans( tNormalMatrix )
                                        * tNormalMatrix * tLeaderFI->gradx( iOrder ) * tCMLeaderDiff->dConstdDOF( tDofType ) );
                    }
                }

                // compute the Jacobian for indirect dof dependencies through follower constitutive models
                uint tFollowerNumDofDependencies = mRequestedFollowerGlobalDofTypes.size();
                for( uint iDOF = 0; iDOF < tFollowerNumDofDependencies; iDOF++ )
                {
                    // get dof type
                    const Vector< MSI::Dof_Type > & tDofType = mRequestedFollowerGlobalDofTypes( iDOF );

                    // get index for the dof type
                    sint tIndexDep      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::FOLLOWER );
                    uint tDepStartIndex = mSet->get_jac_dof_assembly_map()( tDofIndexFollower )( tIndexDep, 0 );
                    uint tDepStopIndex  = mSet->get_jac_dof_assembly_map()( tDofIndexFollower )( tIndexDep, 1 );

                    // get sub-matrices
                    auto tJacLeader = mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tDepStartIndex,       tDepStopIndex } );

                    auto tJacFollower = mSet->get_jacobian()(
                            { tFollowerResStartIndex, tFollowerResStopIndex },
                            { tDepStartIndex,      tDepStopIndex } );

                    // if dependency on the dof type
                    if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                    {
                        // add contribution to Jacobian
                        tJacLeader -= aWStar * (
                                        tSPGhost->val()( 0 )
                                        * trans( tLeaderFI->dnNdxn( iOrder ) ) * trans( tNormalMatrix )
                                        * tNormalMatrix * tFollowerConductivity * tFollowerFI->dnNdxn( iOrder ) );

                        tJacFollower += aWStar * (
                                        tSPGhost->val()( 0 )
                                        * trans( tFollowerFI->dnNdxn( iOrder ) ) * trans( tNormalMatrix )
                                        * tNormalMatrix * tFollowerConductivity * tFollowerFI->dnNdxn( iOrder ) );
                    }

                    // if diffusion CM depends on dof type
                    if( tCMFollowerDiff->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to Jacobian
                        tJacLeader -= aWStar * (
                                tSPGhost->val()( 0 )
                                * trans( tLeaderFI->dnNdxn( iOrder ) ) * trans( tNormalMatrix )
                                * tNormalMatrix * tFollowerFI->gradx( iOrder ) * tCMFollowerDiff->dConstdDOF( tDofType ) );

                        tJacFollower += aWStar * (
                                tSPGhost->val()( 0 )
                                * trans( tFollowerFI->dnNdxn( iOrder ) ) * trans( tNormalMatrix )
                                * tNormalMatrix * tFollowerFI->gradx( iOrder ) * tCMFollowerDiff->dConstdDOF( tDofType ) );
                    }
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Diffusion_Virtual_Work_Ghost::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Diffusion_Virtual_Work_Ghost::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Diffusion_Virtual_Work_Ghost::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Diffusion_Virtual_Work_Ghost::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Diffusion_Virtual_Work_Ghost::compute_dRdp - This function does nothing.");
        }

        //------------------------------------------------------------------------------

        void IWG_Diffusion_Virtual_Work_Ghost::get_flat_normal_matrix(
                Matrix< DDRMat > & aFlatNormal,
                uint               aOrder )
        {
            // get spatial dimensions
            uint tSpaceDim = mNormal.numel();

            // switch on the ghost order
            switch( aOrder )
            {
                case 1 :
                {
                    switch ( tSpaceDim )
                    {
                        case 2 :
                        {
                            aFlatNormal = trans( mNormal );
                            break;
                        }
                        case 3 :
                        {
                            aFlatNormal = trans( mNormal );
                            break;
                        }
                        default:
                        {
                            MORIS_ERROR( false, "IWG_Diffusion_Virtual_Work_Ghost::get_flat_normal_matrix - Spatial dimensions can only be 2, 3." );
                        }
                    }
                    break;
                }

                case 2 :
                {
                    switch ( tSpaceDim )
                    {
                        case 2 :
                        {
                            // set the normal matrix size
                            aFlatNormal.set_size( 2, 3, 0.0 );

                            // fill the normal matrix
                            aFlatNormal( 0, 0 ) = mNormal( 0 );
                            aFlatNormal( 1, 1 ) = mNormal( 1 );

                            aFlatNormal( 0, 2 ) = mNormal( 1 );
                            aFlatNormal( 1, 2 ) = mNormal( 0 );

                            break;
                        }
                        case 3 :
                        {
                            // set the normal matrix size
                            aFlatNormal.set_size( 3, 6, 0.0 );

                            // fill the normal matrix
                            aFlatNormal( 0, 0 ) = mNormal( 0 );
                            aFlatNormal( 1, 1 ) = mNormal( 1 );
                            aFlatNormal( 2, 2 ) = mNormal( 2 );

                            aFlatNormal( 1, 3 ) = mNormal( 2 );
                            aFlatNormal( 2, 3 ) = mNormal( 1 );

                            aFlatNormal( 0, 4 ) = mNormal( 2 );
                            aFlatNormal( 2, 4 ) = mNormal( 0 );

                            aFlatNormal( 0, 5 ) = mNormal( 1 );
                            aFlatNormal( 1, 5 ) = mNormal( 0 );

                            break;
                        }
                        default:
                        {
                            MORIS_ERROR( false, "IWG_Diffusion_Virtual_Work_Ghost::get_flat_normal_matrix - Spatial dimensions can only be 2, 3." );
                        }
                    }
                    break;
                }

                case 3 :
                {
                    switch ( tSpaceDim )
                    {
                        case 2 :
                        {
                            // set the normal matrix size
                            aFlatNormal.set_size( 3, 4, 0.0 );

                            aFlatNormal( 0, 0 ) = mNormal( 0 );
                            aFlatNormal( 1, 1 ) = mNormal( 1 );

                            aFlatNormal( 0, 2 ) = mNormal( 1 );
                            aFlatNormal( 1, 3 ) = mNormal( 0 );

                            real tSqrtOf2 = std::sqrt( 2 );

                            aFlatNormal( 2, 2 ) = tSqrtOf2 * mNormal( 0 );
                            aFlatNormal( 2, 3 ) = tSqrtOf2 * mNormal( 1 );
                            break;
                        }
                        case 3 :
                        {
                            // set the normal matrix size
                            aFlatNormal.set_size( 6, 10, 0.0 );

                            aFlatNormal( 0, 0 ) = mNormal( 0 );
                            aFlatNormal( 1, 1 ) = mNormal( 1 );
                            aFlatNormal( 2, 2 ) = mNormal( 2 );

                            aFlatNormal( 0, 3 ) = mNormal( 1 );
                            aFlatNormal( 0, 4 ) = mNormal( 2 );

                            aFlatNormal( 1, 5 ) = mNormal( 0 );
                            aFlatNormal( 1, 6 ) = mNormal( 2 );

                            aFlatNormal( 2, 7 ) = mNormal( 0 );
                            aFlatNormal( 2, 8 ) = mNormal( 1 );

                            real tSqrtOf2 = std::sqrt( 2 );

                            aFlatNormal( 3, 3 ) = tSqrtOf2 * mNormal( 0 );
                            aFlatNormal( 3, 5 ) = tSqrtOf2 * mNormal( 1 );
                            aFlatNormal( 3, 9 ) = tSqrtOf2 * mNormal( 2 );

                            aFlatNormal( 4, 6 ) = tSqrtOf2 * mNormal( 1 );
                            aFlatNormal( 4, 8 ) = tSqrtOf2 * mNormal( 2 );
                            aFlatNormal( 4, 9 ) = tSqrtOf2 * mNormal( 0 );

                            aFlatNormal( 5, 4 ) = tSqrtOf2 * mNormal( 0 );
                            aFlatNormal( 5, 7 ) = tSqrtOf2 * mNormal( 2 );
                            aFlatNormal( 5, 9 ) = tSqrtOf2 * mNormal( 1 );
                            break;
                        }
                        default:
                        {
                            MORIS_ERROR( false, "IWG_Diffusion_Virtual_Work_Ghost::get_flat_normal_matrix - Spatial dimensions can only be 2, 3." );
                        }
                    }
                    break;
                }

                default:
                {
                    MORIS_ERROR( false, "IWG_Diffusion_Virtual_Work_Ghost::get_flat_normal_matrix - order not supported" );
                }
            }
        }
        //------------------------------------------------------------------------------
}    // namespace moris::fem
