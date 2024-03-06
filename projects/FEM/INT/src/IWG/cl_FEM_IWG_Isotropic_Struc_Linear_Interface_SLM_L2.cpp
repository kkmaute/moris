/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Linear_Interface_SLM_L2.cpp
 *
 */

#include "cl_FEM_IWG_Isotropic_Struc_Linear_Interface_SLM_L2.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

    	IWG_Isotropic_Struc_Linear_Interface_SLM_L2::IWG_Isotropic_Struc_Linear_Interface_SLM_L2()
        {
            init_property("Normal", IWG_Property_Type::NORMAL);
            init_constitutive_model("ElastLinIso", IWG_Constitutive_Type::ELAST_LIN_ISO);

        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Interface_SLM_L2::compute_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here displacement), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get field interpolator for dof type
            Field_Interpolator* tFILambda = get_leader_fi_manager()->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the property immitating normal
            const std::shared_ptr< Property > & tPropNormal = get_leader_property(IWG_Property_Type::NORMAL);

            // get elasticity CM
            const std::shared_ptr< Constitutive_Model > & tCMElasticity = get_leader_constitutive_model(IWG_Constitutive_Type::ELAST_LIN_ISO);

            // get sub-matrix
            auto tRes = mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex } );

            // compute the residual
            tRes -= aWStar * tFILambda->N_trans() * ( tFILambda->val() - tCMElasticity->traction( tPropNormal->val() ) );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Isotropic_Struc_Linear_Interface_SLM_L2::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Interface_SLM_L2::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here displacement), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get field interpolator for dof type
            Field_Interpolator* tFILambda = get_leader_fi_manager()->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the property immitating normal
            const std::shared_ptr< Property > & tPropNormal = get_leader_property(IWG_Property_Type::NORMAL);

            // get elasticity CM
            const std::shared_ptr< Constitutive_Model > & tCMElasticity = get_leader_constitutive_model(IWG_Constitutive_Type::ELAST_LIN_ISO);

            // get the number of leader dof dependencies
            uint tNumDofDependencies = get_requested_leader_dof_types().size();

            // loop over the leader dof dependencies
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                const Vector< MSI::Dof_Type > & tDofType = get_requested_leader_dof_types()( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // get sub-matrix
                auto tJac = mSet->get_jacobian()(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } );

                // if constitutive model depends on the dof type
                if ( tDofType( 0 ) == MSI::Dof_Type::VX )
                {
                    // compute the contribution to Jacobian
                    tJac -= aWStar  * ( tFILambda->N_trans() * tFILambda->N() );
                }

                // if constitutive model depends on the dof type
                if ( tCMElasticity->check_dof_dependency( tDofType ) )
                {
                    // compute the contribution to Jacobian
                    tJac += aWStar *  ( tFILambda->N_trans() * tCMElasticity->dTractiondDOF( tDofType, tPropNormal->val() ) );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Isotropic_Struc_Linear_Interface_SLM_L2::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Interface_SLM_L2::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Interface_SLM_L2::compute_jacobian_and_residual - This function does nothing.");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Interface_SLM_L2::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Interface_SLM_L2::compute_dRdp - This function does nothing.");
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

