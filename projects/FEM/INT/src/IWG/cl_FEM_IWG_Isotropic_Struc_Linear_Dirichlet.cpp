/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Linear_Dirichlet.cpp
 *
 */

#include "cl_FEM_IWG_Isotropic_Struc_Linear_Dirichlet.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Isotropic_Struc_Linear_Dirichlet::IWG_Isotropic_Struc_Linear_Dirichlet( sint aBeta )
        {
            // sign for symmetric/unsymmetric Nitsche
            mBeta = aBeta;
            init_property("Dirichlet", IWG_Property_Type::DIRICHLET);
            init_property("Select", IWG_Property_Type::SELECT);
            init_property("Thickness", IWG_Property_Type::THICKNESS);
            init_constitutive_model("ElastLinIso", IWG_Constitutive_Type::ELAST_LIN_ISO);
            init_stabilization_parameter("DirichletNitsche", IWG_Stabilization_Type::DIRICHLET_NITSCHE);
        }

        //------------------------------------------------------------------------------

        void
        IWG_Isotropic_Struc_Linear_Dirichlet::compute_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here displacement), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get field interpolator for residual dof type
            Field_Interpolator* tFIDispl = get_leader_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the selection matrix property
            const std::shared_ptr< Property >& tPropSelect = get_leader_property(IWG_Property_Type::SELECT);

            // set a default selection matrix if needed
            Matrix< DDRMat > tM;
            if ( tPropSelect == nullptr )
            {
                // get number of fields which should equal spatial dimension
                const uint tSpaceDim = tFIDispl->get_number_of_fields();

                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = tPropSelect->val();

                // skip computing residual if projection matrix is zero
                if ( norm( tM ) < MORIS_REAL_EPS ) return;
            }

            // get the imposed displacement property
            const std::shared_ptr< Property >& tPropDirichlet = get_leader_property(IWG_Property_Type::DIRICHLET);

            // get CM for elasticity
            const std::shared_ptr< Constitutive_Model >& tCMElasticity = get_leader_constitutive_model(IWG_Constitutive_Type::ELAST_LIN_ISO);

            // get SP for Nitsche
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche = get_stabilization_parameter(IWG_Stabilization_Type::DIRICHLET_NITSCHE);

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness = get_leader_property(IWG_Property_Type::THICKNESS);

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // compute jump
            Matrix< DDRMat > tJump = tFIDispl->val() - tPropDirichlet->val();

            // compute the residual
            mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex } ) +=
                    aWStar * ( -tFIDispl->N_trans() * tM * tCMElasticity->traction( get_normal() )                                //
                               + mBeta * tCMElasticity->testTraction_trans( get_normal(), mResidualDofType( 0 ) ) * tM * tJump    //
                               + tSPNitsche->val()( 0 ) * tFIDispl->N_trans() * tM * tJump );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Isotropic_Struc_Linear_Dirichlet::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Isotropic_Struc_Linear_Dirichlet::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here displacement), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get field interpolator for residual dof type
            Field_Interpolator* tFIDispl = get_leader_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the selection matrix property
            const std::shared_ptr< Property >& tPropSelect = get_leader_property(IWG_Property_Type::SELECT);

            // set a default selection matrix if needed
            Matrix< DDRMat > tM;
            if ( tPropSelect == nullptr )
            {
                // get number of fields which should equal spatial dimension
                const uint tSpaceDim = tFIDispl->get_number_of_fields();

                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = tPropSelect->val();

                // skip computing Jacobian if projection matrix is null
                if ( norm( tM ) < MORIS_REAL_EPS ) return;
            }

            // get the imposed displacement property
            const std::shared_ptr< Property >& tPropDirichlet = get_leader_property(IWG_Property_Type::DIRICHLET);

            // get CM for elasticity
            const std::shared_ptr< Constitutive_Model >& tCMElasticity = get_leader_constitutive_model(IWG_Constitutive_Type::ELAST_LIN_ISO);

            // get SP for Nitsche
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche = get_stabilization_parameter(IWG_Stabilization_Type::DIRICHLET_NITSCHE);

            const std::shared_ptr< Property >& tPropThickness = get_leader_property(IWG_Property_Type::THICKNESS);

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // compute jump
            Matrix< DDRMat > tJump = tFIDispl->val() - tPropDirichlet->val();

            // get number of dof dependencies
            uint tNumDofDependencies = get_requested_leader_dof_types().size();
            for ( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the dof type
                const Vector< MSI::Dof_Type >& tDofType = get_requested_leader_dof_types()( iDOF );

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
                    // compute the jacobian for direct dof dependencies
                    tJac += aWStar * (                                                                                                  //
                                    mBeta * tCMElasticity->testTraction_trans( get_normal(), mResidualDofType( 0 ) ) * tM * tFIDispl->N()    //
                                    + tSPNitsche->val()( 0 ) * tFIDispl->N_trans() * tM * tFIDispl->N() );
                }

                // if dependency on the dof type
                if ( tPropDirichlet->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    tJac -= aWStar * (                                                                                                                          //
                                    mBeta * tCMElasticity->testTraction_trans( get_normal(), mResidualDofType( 0 ) ) * tM * tPropDirichlet->dPropdDOF( tDofType )    //
                                    + tSPNitsche->val()( 0 ) * tFIDispl->N_trans() * tM * tPropDirichlet->dPropdDOF( tDofType ) );
                }

                // if dependency on the dof type
                if ( tCMElasticity->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    tJac += aWStar * (                                                                               //
                                    -tFIDispl->N_trans() * tM * tCMElasticity->dTractiondDOF( tDofType, get_normal() )    //
                                    + mBeta * tCMElasticity->dTestTractiondDOF( tDofType, get_normal(), tM * tJump, mResidualDofType( 0 ) ) );
                }

                // if dependency on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    tJac += aWStar * ( tFIDispl->N_trans() * tM * tJump * tSPNitsche->dSPdLeaderDOF( tDofType ) );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Isotropic_Struc_Linear_Dirichlet::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Isotropic_Struc_Linear_Dirichlet::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Dirichlet::compute_jacobian_and_residual - This function does nothing." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Isotropic_Struc_Linear_Dirichlet::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Dirichlet::compute_dRdp - This function does nothing." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

