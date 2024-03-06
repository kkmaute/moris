/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Gap_Nitsche_Unbiased_Unbiased.cpp
 *
 */

#include "cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Gap_Nitsche_Unbiased.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Stabilization_Parameter.hpp"
#include "fn_assert.hpp"
#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"
#include <memory>

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IWG_Isotropic_Struc_Linear_Contact_Gap_Nitsche_Unbiased::IWG_Isotropic_Struc_Linear_Contact_Gap_Nitsche_Unbiased( sint aBeta )
        {
            // sign for symmetric/unsymmetric Nitsche
            mBeta = aBeta;
            init_property("Thickness", IWG_Property_Type::THICKNESS);
            init_property("Gap", IWG_Property_Type::GAP);
            init_constitutive_model("ElastLinIso", IWG_Constitutive_Type::ELAST_LIN_ISO);
            init_stabilization_parameter("NitscheInterface", IWG_Stabilization_Type::NITSCHE_INTERFACE);
        }

        //------------------------------------------------------------------------------

        void
        IWG_Isotropic_Struc_Linear_Contact_Gap_Nitsche_Unbiased::compute_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader and follower field interpolators
            this->check_field_interpolators( mtk::Leader_Follower::LEADER );
            this->check_field_interpolators( mtk::Leader_Follower::FOLLOWER );
#endif

            // get leader index for residual dof type, indices for assembly
            uint const tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint const tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint const tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get leader field interpolator for the residual dof type
            Field_Interpolator*    tFILeader = get_leader_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );
            Geometry_Interpolator* tGILeader = get_leader_fi_manager()->get_IG_geometry_interpolator();
            // get follower field interpolator for the residual dof type
            Field_Interpolator*    tFIFollower = get_follower_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );
            Geometry_Interpolator* tGIFollower = get_follower_fi_manager()->get_IG_geometry_interpolator();

            // get the elasticity constitutive model
            const std::shared_ptr< Constitutive_Model >& tCMLeaderElasticity = get_leader_constitutive_model(IWG_Constitutive_Type::ELAST_LIN_ISO);

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche = get_stabilization_parameter(IWG_Stabilization_Type::NITSCHE_INTERFACE);

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness = get_leader_property(IWG_Property_Type::THICKNESS);

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            const real tNitsche = tSPNitsche->val()( 0 );

            // normals on leader and follower side (needs to be generalized)
            const Matrix< DDRMat > tNormalLeader = get_normal();

            // normal projection operator
            const Matrix< DDRMat > tNormalProjectorLeader = tNormalLeader * trans( tNormalLeader );

            // initial gap
            const real tGap = norm( tGIFollower->valx() - tGILeader->valx() );

            // compute the jump
            const Matrix< DDRMat > tJumpLeader = tFILeader->val() - tFIFollower->val() - tNormalLeader * tGap;

            // compute projection of displacement jump onto normal
            const real tNormalJumpLeader = dot( tJumpLeader, tNormalLeader );

            // evaluate tractions
            const Matrix< DDRMat > tTractionLeader = tCMLeaderElasticity->traction( tNormalLeader );

            // compute contact pressure
            const real tIfcPressureLeader = dot( tTractionLeader, tNormalLeader );

            // check for contact on leader side
            if ( tIfcPressureLeader - tNitsche * tNormalJumpLeader < 0 )
            {
                // compute leader residual
                mSet->get_residual()( 0 )(
                        { tLeaderResStartIndex, tLeaderResStopIndex } ) +=
                        aWStar * 0.5 * (
                                -tFILeader->N_trans() * tNormalProjectorLeader * tTractionLeader
                                + mBeta * tCMLeaderElasticity->testTraction_trans( tNormalLeader, mResidualDofType( 0 ) ) * tNormalProjectorLeader * tJumpLeader
                                + tNitsche * tFILeader->N_trans() * tNormalProjectorLeader * tJumpLeader );
            }
            else
            {
                mSet->get_residual()( 0 )(
                        { tLeaderResStartIndex, tLeaderResStopIndex } ) +=
                        aWStar * 0.5 * ( -mBeta / tNitsche * tCMLeaderElasticity->testTraction_trans( tNormalLeader, mResidualDofType( 0 ) ) * tNormalProjectorLeader * tTractionLeader );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Isotropic_Struc_Linear_Contact_Gap::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Isotropic_Struc_Linear_Contact_Gap_Nitsche_Unbiased::compute_jacobian( real aWStar )
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
            Field_Interpolator*    tFILeader = get_leader_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );
            Geometry_Interpolator* tGILeader = get_leader_fi_manager()->get_IG_geometry_interpolator();
            // get follower field interpolator for the residual dof type
            Field_Interpolator*    tFIFollower = get_follower_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );
            Geometry_Interpolator* tGIFollower = get_follower_fi_manager()->get_IG_geometry_interpolator();

            real tGap = norm( tGIFollower->valx() - tGILeader->valx() );

            // get the elasticity constitutive model
            const std::shared_ptr< Constitutive_Model >& tCMLeaderElasticity = get_leader_constitutive_model(IWG_Constitutive_Type::ELAST_LIN_ISO);

            const std::shared_ptr< Constitutive_Model > &tCMFollowerElasticity = get_follower_constitutive_model(IWG_Constitutive_Type::ELAST_LIN_ISO);

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche = get_stabilization_parameter(IWG_Stabilization_Type::NITSCHE_INTERFACE);

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness = get_leader_property(IWG_Property_Type::THICKNESS);

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            const real tNitsche      = tSPNitsche->val()( 0 );
            const real tLeaderWeight = tSPNitsche->val()( 1 );

            // normals on leader and follower side (needs to be generalized)
            const Matrix< DDRMat > tNormalLeader   = get_normal();
            const Matrix< DDRMat > tNormalFollower = -1.0 * get_normal();

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
            const uint tLeaderNumDofDependencies = get_requested_leader_dof_types().size();

            // compute the Jacobian for indirect dof dependencies through leader constitutive models
            for ( uint iDOF = 0; iDOF < tLeaderNumDofDependencies; iDOF++ )
            {
                // get the dof type
                const Vector< MSI::Dof_Type >& tDofType = get_requested_leader_dof_types()( iDOF );

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
            uint tFollowerNumDofDependencies = get_requested_follower_dof_types().size();
            for ( uint iDOF = 0; iDOF < tFollowerNumDofDependencies; iDOF++ )
            {
                // get dof type
                const Vector< MSI::Dof_Type >& tDofType = get_requested_follower_dof_types()( iDOF );

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
        IWG_Isotropic_Struc_Linear_Contact_Gap_Nitsche_Unbiased::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Contact_Gap_Nitsche_Unbiased::compute_jacobian_and_residual - This function does nothing." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Isotropic_Struc_Linear_Contact_Gap_Nitsche_Unbiased::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Contact_Gap_Nitsche_Unbiased::compute_dRdp - This function does nothing." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
