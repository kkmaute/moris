/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Linear_Fluid_Interface.cpp
 *
 */

#include "cl_FEM_IWG_Isotropic_Struc_Linear_Fluid_Interface.hpp"

#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Isotropic_Struc_Linear_Fluid_Interface::IWG_Isotropic_Struc_Linear_Fluid_Interface()
        {
            init_property("Thickness", IWG_Property_Type::THICKNESS);
            init_constitutive_model("ElastLinIso", IWG_Constitutive_Type::ELAST_LIN_ISO);
            init_constitutive_model("Fluid", IWG_Constitutive_Type::FLUID);
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Fluid_Interface::compute_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here displacement), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get field interpolator for the residual dof type
            Field_Interpolator * tFILeader = get_leader_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ) );

            // get follower fluid constitutive model
            const std::shared_ptr< Constitutive_Model > &tCMFollowerFluid = get_follower_constitutive_model(IWG_Constitutive_Type::FLUID);

            // get thickness property
            const std::shared_ptr< Property > & tPropThickness = get_leader_property(IWG_Property_Type::THICKNESS);

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= (tPropThickness!=nullptr) ? tPropThickness->val()(0) : 1;

            // compute the residual
            mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex } ) -= aWStar * (
                            tFILeader->N_trans() * tCMFollowerFluid->traction( get_normal() ) );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Isotropic_Struc_Linear_Fluid_Interface::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Fluid_Interface::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here displacement), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get field interpolator for the residual dof type
            Field_Interpolator * tFILeader = get_leader_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ) );

            // get follower fluid constitutive model
            const std::shared_ptr< Constitutive_Model > &tCMFollowerFluid = get_follower_constitutive_model(IWG_Constitutive_Type::FLUID);

            // get thickness property
            const std::shared_ptr< Property > & tPropThickness = get_leader_property(IWG_Property_Type::THICKNESS);

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= (tPropThickness!=nullptr) ? tPropThickness->val()(0) : 1;

            // get number of leader dof dependencies
            const  uint tNumDofDependencies = get_requested_follower_dof_types().size();

            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get dependent dof type
                const Vector< MSI::Dof_Type > & tDofType = get_requested_follower_dof_types()( iDOF );

                // get the index for the dof type
                const sint tDofDepIndex        = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::FOLLOWER );
                const sint tFollowerDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                const sint tFollowerDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // if dependency on the dof type
                if ( tCMFollowerFluid->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tFollowerDepStartIndex,  tFollowerDepStopIndex  } ) += aWStar * (
                                    - tFILeader->N_trans() * tCMFollowerFluid->dTractiondDOF( tDofType, get_normal() ) );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Isotropic_Struc_Linear_Fluid_Interface::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Fluid_Interface::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Fluid_Interface::compute_jacobian_and_residual - Not implemented.");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Fluid_Interface::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Fluid_Interface::compute_dRdp - Not implemented.");
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

