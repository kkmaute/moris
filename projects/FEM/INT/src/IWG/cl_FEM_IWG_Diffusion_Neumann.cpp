/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Diffusion_Neumann.cpp
 *
 */

#include "cl_FEM_IWG_Diffusion_Neumann.hpp"
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

        IWG_Diffusion_Neumann::IWG_Diffusion_Neumann()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Neumann" ] = static_cast< uint >( IWG_Property_Type::NEUMANN );
        }

        //------------------------------------------------------------------------------

        void IWG_Diffusion_Neumann::compute_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check master field interpolators, properties, constitutive models
            this->check_field_interpolators();
#endif

            // get index for residual dof type, indices for assembly
            uint tDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get filed interpolator for residual dof type
            Field_Interpolator * tFI =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get indices for SP, CM, properties
            uint tNeumannIndex = static_cast< uint >( IWG_Property_Type::NEUMANN );

            // compute the residual
            mSet->get_residual()( 0 )(
                    { tResStartIndex, tResStopIndex } ) -=
                            aWStar * ( tFI->N_trans() * mMasterProp( tNeumannIndex )->val() );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Diffusion_Neumann::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Diffusion_Neumann::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check master field interpolators, properties, constitutive models
            this->check_field_interpolators();
#endif

            // get index for residual dof type, indices for assembly
            uint tDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get field interpolator for residual dof type
            Field_Interpolator * tFI =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get indices for SP, CM, properties
            uint tNeumannIndex = static_cast< uint >( IWG_Property_Type::NEUMANN );

            // compute the jacobian for dof dependencies
            for( uint iDOF = 0; iDOF < mRequestedMasterGlobalDofTypes.size(); iDOF++ )
            {
                // get dof type
                const Cell< MSI::Dof_Type > tDepDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the dof type indices for assembly
                uint tDepDofIndex   = mSet->get_dof_index_for_type( tDepDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tDepStartIndex = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDepDofIndex, 0 );
                uint tDepStopIndex  = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDepDofIndex, 1 );

                // if dependency in the dof type
                if ( mMasterProp( tNeumannIndex )->check_dof_dependency( tDepDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tResStartIndex, tResStopIndex },
                            { tDepStartIndex, tDepStopIndex } ) -= aWStar * (
                                    tFI->N_trans() * mMasterProp( tNeumannIndex )->dPropdDOF( tDepDofType ) );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Diffusion_Neumann::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Diffusion_Neumann::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, " IWG_Diffusion_Neumann::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Diffusion_Neumann::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Diffusion_Neumann::compute_dRdp - Not implemented.");
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

