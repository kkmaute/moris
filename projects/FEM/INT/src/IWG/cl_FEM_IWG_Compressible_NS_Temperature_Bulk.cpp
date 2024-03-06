/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Compressible_NS_Temperature_Bulk.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Compressible_NS_Temperature_Bulk.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Compressible_NS_Temperature_Bulk::IWG_Compressible_NS_Temperature_Bulk()
        {
            init_property( "BodyForce", IWG_Property_Type::BODY_FORCE );
            init_property( "BodyHeatLoad", IWG_Property_Type::BODY_HEAT_LOAD );
            init_constitutive_model( "Fluid", IWG_Constitutive_Type::FLUID );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Temperature_Bulk::compute_residual( real aWStar )
        {
            // check leader field interpolators
#ifdef MORIS_HAVE_DEBUG
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here velocity), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get the Temperature FI
            Field_Interpolator *tFITemp = get_leader_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get Velocity FI
            Field_Interpolator *tFIVelocity = get_leader_fi_manager()->get_field_interpolators_for_type( mDofVelocity );

            // get properties
            std::shared_ptr< Property > tPropBodyForce    = get_leader_property(IWG_Property_Type::BODY_FORCE);
            std::shared_ptr< Property > tPropBodyHeatLoad = get_leader_property(IWG_Property_Type::BODY_HEAT_LOAD);

            // get the compressible fluid constitutive model
            std::shared_ptr< Constitutive_Model > tCMFluid = get_leader_constitutive_model(IWG_Constitutive_Type::FLUID);

            // compute the residual
            mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex },
                    { 0, 0 } ) += aWStar * ( trans( tFITemp->N() ) * tCMFluid->EnergyDot() + trans( tFITemp->dnNdxn( 1 ) ) * ( tCMFluid->flux( CM_Function_Type::WORK ) - tCMFluid->flux( CM_Function_Type::ENERGY ) - tCMFluid->flux( CM_Function_Type::THERMAL ) ) );

            // if there is a body force
            if ( tPropBodyForce != nullptr )
            {
                // get Velocity FI
                Field_Interpolator *tFIDensity = get_leader_fi_manager()->get_field_interpolators_for_type( mDofDensity );

                // add contribution
                mSet->get_residual()( 0 )(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { 0, 0 } ) -= aWStar * ( tFIDensity->val()( 0 ) * trans( tFITemp->N() ) * dot( tPropBodyForce->val(), tFIVelocity->val() ) );
            }

            // if there is a body heat load
            if ( tPropBodyHeatLoad != nullptr )
            {
                // add contribution
                mSet->get_residual()( 0 )(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { 0, 0 } ) -= aWStar * ( trans( tFITemp->N() ) * tPropBodyHeatLoad->val() );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Compressible_NS_Temperature_Bulk::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------
        void IWG_Compressible_NS_Temperature_Bulk::compute_jacobian( real aWStar )
        {
            // check leader field interpolators
#ifdef MORIS_HAVE_DEBUG
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here pressure), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get the field interpolators
            Field_Interpolator *tFITemp     = get_leader_fi_manager()->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );
            Field_Interpolator *tFIVelocity = get_leader_fi_manager()->get_field_interpolators_for_type( mDofVelocity );
            Field_Interpolator *tFIDensity  = get_leader_fi_manager()->get_field_interpolators_for_type( mDofDensity );

            // get the constitutive model
            std::shared_ptr< Constitutive_Model > tCMFluid = get_leader_constitutive_model(IWG_Constitutive_Type::FLUID);

            // get the properties
            std::shared_ptr< Property > tPropBodyForce    = get_leader_property(IWG_Property_Type::BODY_FORCE);
            std::shared_ptr< Property > tPropBodyHeatLoad = get_leader_property(IWG_Property_Type::BODY_HEAT_LOAD);

            // compute the jacobian for dof dependencies
            uint tNumDofDependencies = get_requested_leader_dof_types().size();
            for ( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Vector< MSI::Dof_Type > const &tDofType = get_requested_leader_dof_types()( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // if fluid CM depends on dof type
                if ( tCMFluid->check_dof_dependency( tDofType ) )
                {
                    // add contribution
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * ( trans( tFITemp->N() ) * tCMFluid->dEnergyDotdDOF( tDofType ) + trans( tFITemp->dnNdxn( 1 ) ) * ( tCMFluid->dFluxdDOF( tDofType, CM_Function_Type::WORK ) - tCMFluid->dFluxdDOF( tDofType, CM_Function_Type::ENERGY ) - tCMFluid->dFluxdDOF( tDofType, CM_Function_Type::THERMAL ) ) );
                }

                // if a body force is present
                if ( tPropBodyForce != nullptr )
                {
                    // if the body force depends on the dof type -> indirect dependency
                    if ( tPropBodyForce->check_dof_dependency( tDofType ) )
                    {
                        // compute the jacobian contribution
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -= aWStar * ( tFIDensity->val()( 0 ) * trans( tFITemp->N() ) * trans( tFIVelocity->val() ) * tPropBodyForce->dPropdDOF( tDofType ) );
                    }

                    // if dof type is density and a body force is present
                    if ( tDofType( 0 ) == mDofDensity )
                    {
                        // compute the jacobian contribution
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -= aWStar * ( trans( tFITemp->N() ) * dot( tPropBodyForce->val(), tFIVelocity->val() ) * tFIDensity->N() );
                    }

                    // if dof type is velocity and a body force is present
                    if ( tDofType( 0 ) == mDofVelocity )
                    {
                        // compute the jacobian contribution
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -= aWStar * ( tFIDensity->val()( 0 ) * trans( tFITemp->N() ) * trans( tPropBodyForce->val() ) * tFIVelocity->N() );
                    }
                }

                // if there is a body heat load
                if ( tPropBodyHeatLoad != nullptr )
                {
                    // if the body heat load depends on the dof type -> indirect dependency
                    if ( tPropBodyHeatLoad->check_dof_dependency( tDofType ) )
                    {
                        // compute the jacobian contribution
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex, tLeaderDepStopIndex } ) -= aWStar * ( trans( tFITemp->N() ) * tPropBodyHeatLoad->dPropdDOF( tDofType ) );
                    }
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Compressible_NS_Temperature_Bulk::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Temperature_Bulk::compute_jacobian_and_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Compressible_NS_Temperature_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Temperature_Bulk::compute_dRdp( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Compressible_NS_Temperature_Bulk::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Temperature_Bulk::compute_residual_strong_form(
                Matrix< DDRMat > &aRM,
                real             &aRC )
        {
            MORIS_ERROR( false, "IWG_Compressible_NS_Temperature_Bulk::compute_residual_strong_form - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Temperature_Bulk::compute_jacobian_strong_form(
                Vector< MSI::Dof_Type > aDofTypes,
                Matrix< DDRMat >       &aJM,
                Matrix< DDRMat >       &aJC )
        {
            MORIS_ERROR( false, "IWG_Compressible_NS_Temperature_Bulk::compute_jacobian_strong_form - Not implemented." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
