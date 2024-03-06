/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost.cpp
 *
 */

#include "cl_FEM_IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost.hpp"
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

        IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost::IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost()
        {
            // set ghost flag
            mIsGhost = true;
            init_property( "Thickness", IWG_Property_Type::THICKNESS );
            init_constitutive_model( "ElastLinIso", IWG_Constitutive_Type::ELAST_LIN_ISO );
            init_stabilization_parameter( "GhostVW", IWG_Stabilization_Type::GHOST_VW );
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost::compute_residual( real aWStar )
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

            // get indices for SP, CM and properties
            const std::shared_ptr< Constitutive_Model > &tLeaderCM   = get_leader_constitutive_model( IWG_Constitutive_Type::ELAST_LIN_ISO );
            const std::shared_ptr< Constitutive_Model > &tFollowerCM = get_follower_constitutive_model( IWG_Constitutive_Type::ELAST_LIN_ISO );

            // get the stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > &tSP = get_stabilization_parameter( IWG_Stabilization_Type::GHOST_VW );

            // get thickness property
            const std::shared_ptr< Property > &tPropThickness = get_leader_property( IWG_Property_Type::THICKNESS );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // loop over the order
            for ( uint iOrder = 1; iOrder <= mOrder; iOrder++ )
            {
                // FIXME check to be removed
                MORIS_ERROR( iOrder <= 1, "Not implemented for order higher than one yet" );

                // set the order for the stabilization parameter
                tSP->set_interpolation_order( iOrder );

                // get flattened normal matrix
                Matrix< DDRMat > tNormalMatrix;
                this->get_flat_normal_matrix( tNormalMatrix, 2 );

                // compute the jump in traction
                Matrix< DDRMat > tGradJump =
                        tLeaderCM->traction( get_normal() ) - tFollowerCM->traction( get_normal() );

                // compute the residual
                mSet->get_residual()( 0 )(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { 0, 0 } ) += aWStar * ( tSP->val()( 0 ) * trans( tLeaderCM->testStrain() ) * trans( tNormalMatrix ) * tGradJump );

                mSet->get_residual()( 0 )(
                        { tFollowerResStartIndex, tFollowerResStopIndex },
                        { 0, 0 } ) -= aWStar * ( tSP->val()( 0 ) * trans( tFollowerCM->testStrain() ) * trans( tNormalMatrix ) * tGradJump );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost::compute_jacobian( real aWStar )
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

            // get number of leader dof dependencies
            uint tLeaderNumDofDependencies = get_requested_leader_dof_types().size();

            // get the stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > &tSP = get_stabilization_parameter( IWG_Stabilization_Type::GHOST_VW );

            // get thickness property
            const std::shared_ptr< Property > &tPropThickness = get_leader_property( IWG_Property_Type::THICKNESS );

            const std::shared_ptr< Constitutive_Model > &tLeaderCM   = get_leader_constitutive_model( IWG_Constitutive_Type::ELAST_LIN_ISO );
            const std::shared_ptr< Constitutive_Model > &tFollowerCM = get_follower_constitutive_model( IWG_Constitutive_Type::ELAST_LIN_ISO );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // order 1
            for ( uint iOrder = 1; iOrder <= mOrder; iOrder++ )
            {
                // FIXME check to be removed
                MORIS_ERROR( iOrder <= 1, "Not implemented for order higher than one yet" );

                // set the order for the stabilization parameter
                tSP->set_interpolation_order( iOrder );

                // get normal matrix
                Matrix< DDRMat > tNormalMatrix;
                this->get_flat_normal_matrix( tNormalMatrix, 2 );

                // compute the jacobian for indirect dof dependencies through leader constitutive models
                for ( uint iDOF = 0; iDOF < tLeaderNumDofDependencies; iDOF++ )
                {
                    // get the dof type
                    const Vector< MSI::Dof_Type > &tDofType = get_requested_leader_dof_types()( iDOF );

                    // get the index for the dof type
                    sint tIndexDep      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                    uint tDepStartIndex = mSet->get_jac_dof_assembly_map()( tDofIndexLeader )( tIndexDep, 0 );
                    uint tDepStopIndex  = mSet->get_jac_dof_assembly_map()( tDofIndexLeader )( tIndexDep, 1 );

                    // if dependency on the dof type
                    if ( tLeaderCM->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tDepStartIndex, tDepStopIndex } ) += aWStar * ( tSP->val()( 0 ) * trans( tLeaderCM->testStrain() ) * trans( tNormalMatrix ) * tLeaderCM->dTractiondDOF( tDofType, get_normal() ) );

                        mSet->get_jacobian()(
                                { tFollowerResStartIndex, tFollowerResStopIndex },
                                { tDepStartIndex, tDepStopIndex } ) -= aWStar * ( tSP->val()( 0 ) * trans( tFollowerCM->testStrain() ) * trans( tNormalMatrix ) * tLeaderCM->dTractiondDOF( tDofType, get_normal() ) );
                    }
                }

                // compute the jacobian for indirect dof dependencies through follower constitutive models
                uint tFollowerNumDofDependencies = get_requested_follower_dof_types().size();
                for ( uint iDOF = 0; iDOF < tFollowerNumDofDependencies; iDOF++ )
                {
                    // get dof type
                    Vector< MSI::Dof_Type > tDofType = get_requested_follower_dof_types()( iDOF );

                    // get index for the dof type
                    sint tIndexDep      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::FOLLOWER );
                    uint tDepStartIndex = mSet->get_jac_dof_assembly_map()( tDofIndexFollower )( tIndexDep, 0 );
                    uint tDepStopIndex  = mSet->get_jac_dof_assembly_map()( tDofIndexFollower )( tIndexDep, 1 );

                    // if dependency on the dof type
                    if ( tFollowerCM->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tDepStartIndex, tDepStopIndex } ) -= aWStar * ( tSP->val()( 0 ) * trans( tLeaderCM->testStrain() ) * trans( tNormalMatrix ) * tFollowerCM->dTractiondDOF( tDofType, get_normal() ) );

                        mSet->get_jacobian()(
                                { tFollowerResStartIndex, tFollowerResStopIndex },
                                { tDepStartIndex, tDepStopIndex } ) += aWStar * ( tSP->val()( 0 ) * trans( tFollowerCM->testStrain() ) * trans( tNormalMatrix ) * tFollowerCM->dTractiondDOF( tDofType, get_normal() ) );
                    }
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost::compute_dRdp - This function does nothing." );
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost::get_flat_normal_matrix(
                Matrix< DDRMat > &aFlatNormal,
                uint              aOrder )
        {
            // get spatial dimensions
            uint tSpaceDim = get_normal().numel();

            // get elasticity CM
            const std::shared_ptr< Constitutive_Model > &tCMElasticity = get_leader_constitutive_model( IWG_Constitutive_Type::ELAST_LIN_ISO );

            // switch on the ghost order
            switch ( aOrder )
            {
                case 1:
                {
                    switch ( tSpaceDim )
                    {
                        case 2:
                        {
                            aFlatNormal = trans( get_normal() );
                            break;
                        }
                        case 3:
                        {
                            aFlatNormal = trans( get_normal() );
                            break;
                        }
                        default:
                        {
                            MORIS_ERROR( false, "IWG_Diffusion_Virtual_Work_Ghost::get_flat_normal_matrix - Spatial dimensions can only be 2, 3." );
                        }
                    }
                    break;
                }

                case 2:
                {
                    switch ( tSpaceDim )
                    {
                        case 2:
                        {
                            // set the normal matrix size
                            aFlatNormal.set_size( 2, tCMElasticity->strain().numel(), 0.0 );

                            // fill the normal matrix
                            aFlatNormal( 0, 0 ) = get_normal()( 0 );
                            aFlatNormal( 1, 1 ) = get_normal()( 1 );

                            aFlatNormal( 0, 2 ) = get_normal()( 1 );
                            aFlatNormal( 1, 2 ) = get_normal()( 0 );

                            break;
                        }
                        case 3:
                        {
                            // set the normal matrix size
                            aFlatNormal.set_size( 3, 6, 0.0 );

                            // fill the normal matrix
                            aFlatNormal( 0, 0 ) = get_normal()( 0 );
                            aFlatNormal( 1, 1 ) = get_normal()( 1 );
                            aFlatNormal( 2, 2 ) = get_normal()( 2 );

                            aFlatNormal( 1, 3 ) = get_normal()( 2 );
                            aFlatNormal( 2, 3 ) = get_normal()( 1 );

                            aFlatNormal( 0, 4 ) = get_normal()( 2 );
                            aFlatNormal( 2, 4 ) = get_normal()( 0 );

                            aFlatNormal( 0, 5 ) = get_normal()( 1 );
                            aFlatNormal( 1, 5 ) = get_normal()( 0 );

                            break;
                        }
                        default:
                        {
                            MORIS_ERROR( false, "IWG_Diffusion_Virtual_Work_Ghost::get_flat_normal_matrix - Spatial dimensions can only be 2, 3." );
                        }
                    }
                    break;
                }

                case 3:
                {
                    switch ( tSpaceDim )
                    {
                        case 2:
                        {
                            // set the normal matrix size
                            aFlatNormal.set_size( 3, 4, 0.0 );

                            aFlatNormal( 0, 0 ) = get_normal()( 0 );
                            aFlatNormal( 1, 1 ) = get_normal()( 1 );

                            aFlatNormal( 0, 2 ) = get_normal()( 1 );
                            aFlatNormal( 1, 3 ) = get_normal()( 0 );

                            real tSqrtOf2 = std::sqrt( 2 );

                            aFlatNormal( 2, 2 ) = tSqrtOf2 * get_normal()( 0 );
                            aFlatNormal( 2, 3 ) = tSqrtOf2 * get_normal()( 1 );
                            break;
                        }
                        case 3:
                        {
                            // set the normal matrix size
                            aFlatNormal.set_size( 6, 10, 0.0 );

                            aFlatNormal( 0, 0 ) = get_normal()( 0 );
                            aFlatNormal( 1, 1 ) = get_normal()( 1 );
                            aFlatNormal( 2, 2 ) = get_normal()( 2 );

                            aFlatNormal( 0, 3 ) = get_normal()( 1 );
                            aFlatNormal( 0, 4 ) = get_normal()( 2 );

                            aFlatNormal( 1, 5 ) = get_normal()( 0 );
                            aFlatNormal( 1, 6 ) = get_normal()( 2 );

                            aFlatNormal( 2, 7 ) = get_normal()( 0 );
                            aFlatNormal( 2, 8 ) = get_normal()( 1 );

                            real tSqrtOf2 = std::sqrt( 2 );

                            aFlatNormal( 3, 3 ) = tSqrtOf2 * get_normal()( 0 );
                            aFlatNormal( 3, 5 ) = tSqrtOf2 * get_normal()( 1 );
                            aFlatNormal( 3, 9 ) = tSqrtOf2 * get_normal()( 2 );

                            aFlatNormal( 4, 6 ) = tSqrtOf2 * get_normal()( 1 );
                            aFlatNormal( 4, 8 ) = tSqrtOf2 * get_normal()( 2 );
                            aFlatNormal( 4, 9 ) = tSqrtOf2 * get_normal()( 0 );

                            aFlatNormal( 5, 4 ) = tSqrtOf2 * get_normal()( 0 );
                            aFlatNormal( 5, 7 ) = tSqrtOf2 * get_normal()( 2 );
                            aFlatNormal( 5, 9 ) = tSqrtOf2 * get_normal()( 1 );
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
    } /* namespace fem */
} /* namespace moris */
