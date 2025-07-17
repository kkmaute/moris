/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Incompressible_NS_Convective_Velocity_Ghost.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Incompressible_NS_Convective_Velocity_Ghost.hpp"
//LINALG/src
#include "fn_trans.hpp"
#include "fn_vectorize.hpp"
#include "fn_isfinite.hpp"

namespace moris::fem
{

    //------------------------------------------------------------------------------

    IWG_Incompressible_NS_Convective_Velocity_Ghost::IWG_Incompressible_NS_Convective_Velocity_Ghost()
    {
        // set ghost flag
        mIsGhost = true;

        // set size for the stabilization parameter pointer cell
        mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

        // populate the stabilization map
        mStabilizationMap[ "ConvectiveGhost" ] = static_cast< uint >( IWG_Stabilization_Type::CONVECTIVE_GHOST );
    }

    //------------------------------------------------------------------------------

    void IWG_Incompressible_NS_Convective_Velocity_Ghost::compute_residual( real aWStar )
    {
        // check leader field interpolators
#ifdef MORIS_HAVE_DEBUG
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get follower index for residual dof type, indices for assembly
            uint tFollowerDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );
            uint tFollowerResStartIndex = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 0 );
            uint tFollowerResStopIndex  = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 1 );

            // get the leader field interpolator for the residual dof type
            Field_Interpolator * tFILeader = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the follower field interpolator for the residual dof type
            Field_Interpolator * tFIFollower  = mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the convective stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPConvective =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::CONVECTIVE_GHOST ) );

            // get flattened derivatives dnNdxn for leader and follower
            Matrix< DDRMat > tLeaderdNdx;
            this->compute_dnNdxn( tLeaderdNdx, mtk::Leader_Follower::LEADER );
            Matrix< DDRMat > tFollowerdNdx;
            this->compute_dnNdxn( tFollowerdNdx, mtk::Leader_Follower::FOLLOWER );

            // premultiply common terms
            Matrix< DDRMat > tConvectivePreMultiply = vectorize(
                    tSPConvective->val()( 0 ) * ( tFILeader->gradx( 1 ) - tFIFollower->gradx( 1 ) ) );

            // compute leader residual
            mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex },
                    { 0, 0 } ) += aWStar * (
                            trans( tLeaderdNdx ) * tConvectivePreMultiply );

            // compute follower residual
            mSet->get_residual()( 0 )(
                    { tFollowerResStartIndex, tFollowerResStopIndex },
                    { 0, 0 } ) -= aWStar * (
                            trans( tFollowerdNdx ) * tConvectivePreMultiply );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Incompressible_NS_Convective_Velocity_Ghost::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Convective_Velocity_Ghost::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get follower index for residual dof type, indices for assembly
            uint tFollowerDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );
            uint tFollowerResStartIndex = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 0 );
            uint tFollowerResStopIndex  = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 1 );

            // get the leader field interpolator for residual dof type
            Field_Interpolator * tFILeader = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the follower field interpolator for residual dof type
            Field_Interpolator * tFIFollower  = mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the convective stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPConvective =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::CONVECTIVE_GHOST ) );

            // get flattened derivatives dnNdxn for leader and follower
            Matrix< DDRMat > tLeaderdNdx;
            this->compute_dnNdxn( tLeaderdNdx, mtk::Leader_Follower::LEADER );
            Matrix< DDRMat > tFollowerdNdx;
            this->compute_dnNdxn( tFollowerdNdx, mtk::Leader_Follower::FOLLOWER );

            // get number of leader dependencies
            uint tLeaderNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();
            uint tFollowerNumDofDependencies  = mRequestedFollowerGlobalDofTypes.size();

            // compute the jacobian for indirect dof dependencies through leader
            for( uint iDOF = 0; iDOF < tLeaderNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Vector< MSI::Dof_Type > & tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // compute jacobian direct dependencies
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // dRM/dM
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * (
                                    tSPConvective->val()( 0 ) * trans( tLeaderdNdx ) * tLeaderdNdx );

                    // dRS/dM
                    mSet->get_jacobian()(
                            { tFollowerResStartIndex,  tFollowerResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -= aWStar * (
                                    tSPConvective->val()( 0 ) * trans( tFollowerdNdx ) * tLeaderdNdx );
                }

                // if stabilization parameter dependency on the dof type
                if ( tSPConvective->check_dof_dependency( tDofType ) )
                {
                    // premultiply common terms
                    Matrix< DDRMat > tConvectivePreMultiply = vectorize(
                            tFILeader->gradx( 1 ) - tFIFollower->gradx( 1 ) ) ;

                    tConvectivePreMultiply = tConvectivePreMultiply * tSPConvective->dSPdLeaderDOF( tDofType );

                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * (
                                    trans( tLeaderdNdx ) * tConvectivePreMultiply );

                    mSet->get_jacobian()(
                            { tFollowerResStartIndex,  tFollowerResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -= aWStar * (
                                    trans( tFollowerdNdx ) * tConvectivePreMultiply );
                }
            }

            // compute the jacobian for indirect dof dependencies through follower
            for( uint iDOF = 0; iDOF < tFollowerNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Vector< MSI::Dof_Type > tDofType = mRequestedFollowerGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex        = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::FOLLOWER );
                uint tFollowerDepStartIndex = mSet->get_jac_dof_assembly_map()( tFollowerDofIndex )( tDofDepIndex, 0 );
                uint tFollowerDepStopIndex  = mSet->get_jac_dof_assembly_map()( tFollowerDofIndex )( tDofDepIndex, 1 );

                // compute jacobian direct dependencies
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // dRM/dS
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tFollowerDepStartIndex,  tFollowerDepStopIndex } ) -= aWStar * (
                                    tSPConvective->val()( 0 ) * trans( tLeaderdNdx ) * tFollowerdNdx );

                    // dRS/dS
                    mSet->get_jacobian()(
                            { tFollowerResStartIndex, tFollowerResStopIndex },
                            { tFollowerDepStartIndex, tFollowerDepStopIndex } ) += aWStar * (
                                    tSPConvective->val()( 0 ) * trans( tFollowerdNdx ) * tFollowerdNdx );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Ghost_Normal_Field::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Convective_Velocity_Ghost::compute_jacobian_and_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Convective_Velocity_Ghost::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Convective_Velocity_Ghost::compute_dRdp( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Convective_Velocity_Ghost::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Convective_Velocity_Ghost::compute_dnNdxn(
                Matrix< DDRMat >  & adNdx,
                mtk::Leader_Follower   aIsLeader )
        {
            // get the field interpolator for residual dof type
            Field_Interpolator * tFIVelocity =
                    this->get_field_interpolator_manager( aIsLeader )->
                    get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // init size for dnNdtn
            uint tNumField = tFIVelocity->get_number_of_fields();
            uint tNumRow = tFIVelocity->dnNdxn( 1 ).n_rows();
            uint tNumCol = tFIVelocity->dnNdxn( 1 ).n_cols();
            adNdx.set_size( tNumField * tNumRow, tNumField * tNumCol , 0.0 );

            // loop over the fields
            for( uint iField = 0; iField < tNumField; iField++ )
            {
                // fill the matrix for each dimension
                adNdx(
                        { iField * tNumRow, ( iField + 1 ) * tNumRow - 1 },
                        { iField * tNumCol, ( iField + 1 ) * tNumCol - 1 } ) =
                                tFIVelocity->dnNdxn( 1 ).matrix_data();
            }
        }

        //------------------------------------------------------------------------------
}    // namespace moris::fem
