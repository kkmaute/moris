/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Ghost_Normal_Field.cpp
 *
 */

#include "cl_FEM_IWG_Ghost_Normal_Field.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IWG_Ghost_Normal_Field::IWG_Ghost_Normal_Field()
        {
            // set ghost flag
            mIsGhost = true;

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "GhostSP" ] = static_cast< uint >( IWG_Stabilization_Type::GHOST_SP );
        }

        //------------------------------------------------------------------------------

        void IWG_Ghost_Normal_Field::compute_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader and follower field interpolators
            this->check_field_interpolators( mtk::Leader_Follower::LEADER );
            this->check_field_interpolators( mtk::Leader_Follower::FOLLOWER );
#endif

            // set interpolation order
            IWG::set_interpolation_order();

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get follower index for residual dof type, indices for assembly
            uint tFollowerDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );
            uint tFollowerResStartIndex = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 0 );
            uint tFollowerResStopIndex  = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 1 );

            // get the leader field interpolator for residual dof type
            Field_Interpolator * tFILeader =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the follower field interpolator for residual dof type
            Field_Interpolator * tFIFollower  =
                    mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the stabilization parameter
            std::shared_ptr< Stabilization_Parameter > & tSP =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::GHOST_SP ) );

            // get the part of the matrix to assemble into
            auto tLeaderRes = mSet->get_residual()( 0 )( { tLeaderResStartIndex, tLeaderResStopIndex }, { 0, 0 } );
            auto tFollowerRes  = mSet->get_residual()( 0 )( { tFollowerResStartIndex,  tFollowerResStopIndex  }, { 0, 0 } );

            // loop over the interpolation order
            for ( uint iOrder = 1; iOrder <= mOrder; iOrder++ )
            {
                // set the order for stabilization parameters
                tSP->set_interpolation_order( iOrder );

                // get normal matrix
                Matrix< DDRMat > tFlatNormal;
                this->get_flat_normal_matrix( tFlatNormal, iOrder );

                // premultiply common terms
                Matrix< DDRMat > tPreMultiply = vectorize(
                        tSP->val()( 0 ) * tFlatNormal * ( tFILeader->gradx( iOrder ) - tFIFollower->gradx( iOrder ) ) );

                // get flattened directional derivatives for leader and follower
                Matrix< DDRMat > tLeaderdNdxFlat;
                this->compute_flat_dnNdxn( tLeaderdNdxFlat, iOrder, mtk::Leader_Follower::LEADER );
                Matrix< DDRMat > tFollowerdNdxFlat;
                this->compute_flat_dnNdxn( tFollowerdNdxFlat, iOrder, mtk::Leader_Follower::FOLLOWER );

                // compute leader residual
                tLeaderRes += aWStar * ( tLeaderdNdxFlat * tPreMultiply );

                // compute follower residual
                tFollowerRes -= aWStar * ( tFollowerdNdxFlat * tPreMultiply );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Ghost_Normal_Field::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Ghost_Normal_Field::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader and follower field interpolators
            this->check_field_interpolators( mtk::Leader_Follower::LEADER );
            this->check_field_interpolators( mtk::Leader_Follower::FOLLOWER );
#endif

            // set interpolation order
            IWG::set_interpolation_order();

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get follower index for residual dof type, indices for assembly
            uint tFollowerDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );
            uint tFollowerResStartIndex = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 0 );
            uint tFollowerResStopIndex  = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 1 );

            // get the leader field interpolator for residual dof type
            Field_Interpolator * tFILeader =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the follower field interpolator for residual dof type
            Field_Interpolator * tFIFollower  =
                    mFollowerFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the stabilization parameter
            std::shared_ptr< Stabilization_Parameter > & tSP =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::GHOST_SP ) );

            // get number of leader and follower dependencies
            uint tLeaderNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();
            uint tFollowerNumDofDependencies  = mRequestedFollowerGlobalDofTypes.size();

            // loop over the interpolation orders
            for ( uint iOrder = 1; iOrder <= mOrder; iOrder++ )
            {
                // set the order for the stabilization parameter
                tSP->set_interpolation_order( iOrder );

                // get flattened normal
                Matrix< DDRMat > tFlatNormal;
                this->get_flat_normal_matrix( tFlatNormal, iOrder );

                // get flattened directional derivatives for leader and follower
                Matrix< DDRMat > tLeaderdNdxFlat;
                this->compute_flat_dnNdxn( tLeaderdNdxFlat, iOrder, mtk::Leader_Follower::LEADER );
                Matrix< DDRMat > tFollowerdNdxFlat;
                this->compute_flat_dnNdxn( tFollowerdNdxFlat, iOrder, mtk::Leader_Follower::FOLLOWER );

                // compute the jacobian for indirect dof dependencies through leader
                for( uint iDOF = 0; iDOF < tLeaderNumDofDependencies; iDOF++ )
                {
                    // get the dof type
                    Cell< MSI::Dof_Type > & tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                    // get the index for the dof type
                    sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                    uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                    uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                    // get the part of the Jacobian to assemble into
                    auto tJacLeaderLeader = mSet->get_jacobian()( { tLeaderResStartIndex, tLeaderResStopIndex },  { tLeaderDepStartIndex, tLeaderDepStopIndex } );
                    auto tJacFollowerLeader  = mSet->get_jacobian()( { tFollowerResStartIndex,  tFollowerResStopIndex  },  { tLeaderDepStartIndex, tLeaderDepStopIndex } );

                    // compute jacobian direct dependencies
                    if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                    {
                        // dRM/dM
                        tJacLeaderLeader += aWStar * ( tLeaderdNdxFlat * tSP->val()( 0 ) * trans( tLeaderdNdxFlat ) );

                        // dRS/dM
                        tJacFollowerLeader -= aWStar * ( tFollowerdNdxFlat * tSP->val()( 0 ) * trans( tLeaderdNdxFlat ) );
                    }

                    // if stabilization parameter dependency on the dof type
                    if ( tSP->check_dof_dependency( tDofType ) )
                    {
                        // premultiply common terms
                        Matrix< DDRMat > tPreMultiply = vectorize(
                                tFlatNormal * ( tFILeader->gradx( iOrder ) - tFIFollower->gradx( iOrder ) ) );

                        tPreMultiply = tPreMultiply * tSP->dSPdLeaderDOF( tDofType );

                        // add contribution to jacobian
                        tJacLeaderLeader += aWStar * ( tLeaderdNdxFlat * tPreMultiply );

                        tJacFollowerLeader -= aWStar * ( tFollowerdNdxFlat  * tPreMultiply );
                    }
                }

                // compute the jacobian for indirect dof dependencies through follower
                for( uint iDOF = 0; iDOF < tFollowerNumDofDependencies; iDOF++ )
                {
                    // get the dof type
                    Cell< MSI::Dof_Type > tDofType = mRequestedFollowerGlobalDofTypes( iDOF );

                    // get the index for the dof type
                    sint tDofDepIndex        = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::FOLLOWER );
                    uint tFollowerDepStartIndex = mSet->get_jac_dof_assembly_map()( tFollowerDofIndex )( tDofDepIndex, 0 );
                    uint tFollowerDepStopIndex  = mSet->get_jac_dof_assembly_map()( tFollowerDofIndex )( tDofDepIndex, 1 );

                    // get the part of the Jacobian to assemble into
                    auto tJacLeaderFollower = mSet->get_jacobian()( { tLeaderResStartIndex, tLeaderResStopIndex },  { tFollowerDepStartIndex,  tFollowerDepStopIndex } );
                    auto tJacFollowerFollower  = mSet->get_jacobian()( { tFollowerResStartIndex,  tFollowerResStopIndex  },  { tFollowerDepStartIndex,  tFollowerDepStopIndex } );

                    // compute jacobian direct dependencies
                    if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                    {
                        // dRM/dS
                        tJacLeaderFollower -= aWStar * ( tLeaderdNdxFlat * tSP->val()( 0 ) * trans( tFollowerdNdxFlat ) );

                        // dRS/dS
                        tJacFollowerFollower += aWStar * ( tFollowerdNdxFlat * tSP->val()( 0 ) * trans( tFollowerdNdxFlat ) );
                    }
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Ghost_Normal_Field::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Ghost_Normal_Field::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Ghost_Normal_Field::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Ghost_Normal_Field::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Ghost_Normal_Field::compute_dRdp - This function does nothing.");
        }

        //------------------------------------------------------------------------------

        void IWG_Ghost_Normal_Field::get_flat_normal_matrix(
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
                            MORIS_ERROR( false, "IWG_Ghost_Normal_Field::get_flat_normal_matrix - Spatial dimensions can only be 2, 3." );
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
                            MORIS_ERROR( false, "IWG_Ghost_Normal_Field::get_flat_normal_matrix - Spatial dimensions can only be 2, 3." );
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
                            MORIS_ERROR( false, "IWG_Ghost_Normal_Field::get_flat_normal_matrix - Spatial dimensions can only be 2, 3." );
                        }
                    }
                    break;
                }

                default:
                {
                    MORIS_ERROR( false, "IWG_Ghost_Normal_Field::get_flat_normal_matrix - order not supported" );
                }
            }
        }

        //------------------------------------------------------------------------------

        void IWG_Ghost_Normal_Field::compute_flat_dnNdxn(
                Matrix< DDRMat >  & aFlatdnNdxn,
                uint                aOrder,
                mtk::Leader_Follower   aIsLeader )
        {
            // get flattened normal
            Matrix< DDRMat > tFlatNormal;
            this->get_flat_normal_matrix( tFlatNormal, aOrder );

            // get the residual dof type FI
            Field_Interpolator * tFIResidual =
                    this->get_field_interpolator_manager( aIsLeader )->
                    get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get number of fields
            uint tNumFields = tFIResidual->get_number_of_fields();

            // get flat dnNdxn (dnNdxn . normal)
            Matrix< DDRMat > tdnNdxn = trans( tFIResidual->dnNdxn( aOrder ) ) * trans( tFlatNormal );
            uint tNumRows = tdnNdxn.n_rows();
            uint tNumCols = tdnNdxn.n_cols();

            // set size for block flat dnNdxn (dnNdxn . normal)
            aFlatdnNdxn.set_size( tNumFields * tNumRows, tNumFields * tNumCols, 0.0 );

            // loop over the number of fields
            for( uint iField = 0; iField < tNumFields; iField++ )
            {
                // fill block flat dnNdxn
                aFlatdnNdxn(
                        { iField * tNumRows, ( iField + 1 ) * tNumRows - 1 },
                        { iField * tNumCols, ( iField + 1 ) * tNumCols - 1 } ) =
                                tdnNdxn.matrix_data();
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

