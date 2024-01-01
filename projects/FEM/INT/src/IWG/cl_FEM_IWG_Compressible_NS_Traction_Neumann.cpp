/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Compressible_NS_Traction_Neumann.cpp
 *
 */

#include "cl_FEM_IWG_Compressible_NS_Traction_Neumann.hpp"
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

        IWG_Compressible_NS_Traction_Neumann::IWG_Compressible_NS_Traction_Neumann()
        {
            // set size for the property pointer cell
            mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Traction" ] = static_cast< uint >( IWG_Property_Type::TRACTION );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Traction_Neumann::compute_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties, constitutive models
            this->check_field_interpolators();
#endif

            // get index for residual dof type, indices for assembly
            uint tDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get field interpolator for residual dof type (velocity)
            Field_Interpolator * tFIVelocity =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get indices for SP, CM, properties
            uint tPropTractionIndex = static_cast< uint >( IWG_Property_Type::TRACTION );

            // compute the residual
            mSet->get_residual()( 0 )(
                    { tResStartIndex, tResStopIndex },
                    { 0, 0 } ) += aWStar * ( -1.0 * trans( tFIVelocity->N() ) * mLeaderProp( tPropTractionIndex )->val() );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Compressible_NS_Traction_Neumann::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Traction_Neumann::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties, constitutive models
            this->check_field_interpolators();
#endif

            // get index for residual dof type, indices for assembly
            uint tDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get field interpolator for residual dof type (velocity)
            Field_Interpolator * tFIVelocity =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get indices for SP, CM, properties
            uint tPropTractionIndex = static_cast< uint >( IWG_Property_Type::TRACTION );

            // compute the jacobian for dof dependencies
            for( uint iDOF = 0; iDOF < mRequestedLeaderGlobalDofTypes.size(); iDOF++ )
            {
                // get dof type
                Vector< MSI::Dof_Type > tDepDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                // get the dof type indices for assembly
                uint tDepDofIndex   = mSet->get_dof_index_for_type( tDepDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tDepStartIndex = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDepDofIndex, 0 );
                uint tDepStopIndex  = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDepDofIndex, 1 );

                // if dependency in the dof type
                if ( mLeaderProp( tPropTractionIndex )->check_dof_dependency( tDepDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tResStartIndex, tResStopIndex },
                            { tDepStartIndex, tDepStopIndex } ) += aWStar * ( -1.0 *
                                    trans( tFIVelocity->N() ) * mLeaderProp( tPropTractionIndex )->dPropdDOF( tDepDofType ) );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Compressible_NS_Traction_Neumann::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Traction_Neumann::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, " IWG_Compressible_NS_Traction_Neumann::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Traction_Neumann::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Compressible_NS_Traction_Neumann::compute_dRdp - Not implemented.");
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

