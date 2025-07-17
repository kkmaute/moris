/*
 * Copyright (c) 2022 University of Colorado
 *Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Diffusion_Convection.cpp
 *
 */

#include "cl_FEM_IWG_Diffusion_Convection.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_isfinite.hpp"

namespace moris::fem
{

    //------------------------------------------------------------------------------

    IWG_Diffusion_Convection::IWG_Diffusion_Convection()
    {
        // set size for the property pointer cell
        mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "HeatTransferCoefficient" ] = static_cast< uint >( IWG_Property_Type::HEAT_TRANSFER_COEFFICIENT );
        mPropertyMap[ "AmbientTemperature" ]      = static_cast< uint >( IWG_Property_Type::AMBIENT_TEMP );
        mPropertyMap[ "Thickness" ]               = static_cast< uint >( IWG_Property_Type::THICKNESS );
    }

    //------------------------------------------------------------------------------

    void
    IWG_Diffusion_Convection::compute_residual( real aWStar )
    {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties, constitutive models
            this->check_field_interpolators();
#endif
            // get heat transfer property
            const std::shared_ptr< Property >& tPropHeatTransfer =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::HEAT_TRANSFER_COEFFICIENT ) );

            // get ambient temperature property
            const std::shared_ptr< Property >& tPropAmbientTemp =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::AMBIENT_TEMP ) );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // get index for residual dof type, indices for assembly
            uint tDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get filed interpolator for residual dof type
            Field_Interpolator* tFI = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // compute the residual
            // N * a * (T - T_ref)
            mSet->get_residual()( 0 )( { tResStartIndex, tResStopIndex } ) +=      //
                    aWStar * ( tFI->N_trans() * tPropHeatTransfer->val()( 0 ) *    //
                               ( tFI->val() - tPropAmbientTemp->val() ) );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Diffusion_Convection::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Convection::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties, constitutive models
            this->check_field_interpolators();
#endif
            // get heat transfer property
            const std::shared_ptr< Property >& tPropHeatTransfer =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::HEAT_TRANSFER_COEFFICIENT ) );

            // get ambient temperature property
            const std::shared_ptr< Property >& tPropAmbientTemp =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::AMBIENT_TEMP ) );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // get index for residual dof type, indices for assembly
            uint tDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get field interpolator for residual dof type
            Field_Interpolator* tFI = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // compute the Jacobian for dof dependencies
            for ( uint iDOF = 0; iDOF < mRequestedLeaderGlobalDofTypes.size(); iDOF++ )
            {
                // get dof type
                const Vector< MSI::Dof_Type >& tDepDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                // get the dof type indices for assembly
                uint tDepDofIndex   = mSet->get_dof_index_for_type( tDepDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tDepStartIndex = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDepDofIndex, 0 );
                uint tDepStopIndex  = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDepDofIndex, 1 );

                // get sub-matrix
                auto tJac = mSet->get_jacobian()(
                        { tResStartIndex, tResStopIndex },
                        { tDepStartIndex, tDepStopIndex } );

                // if dof type is residual dof type
                if ( tDepDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    tJac += aWStar * tPropHeatTransfer->val()( 0 ) * tFI->N_trans() * tFI->N();
                }

                // if dependency of heat transfer coefficient on dof type
                if ( tPropHeatTransfer->check_dof_dependency( tDepDofType ) )
                {
                    // add contribution to Jacobian
                    tJac += aWStar * tFI->N_trans() * ( tFI->val() - tPropAmbientTemp->val() ) *    //
                            tPropHeatTransfer->dPropdDOF( tDepDofType );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Diffusion_Convection::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Convection::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, " IWG_Diffusion_Convection::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Diffusion_Convection::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Diffusion_Convection::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------
}    // namespace moris::fem
