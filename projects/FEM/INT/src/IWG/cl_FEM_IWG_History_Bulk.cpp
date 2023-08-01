/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_History_Bulk.cpp
 *
 */

#include "cl_FEM_IWG_History_Bulk.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_CM_Struc_Linear_Isotropic_Damage.hpp"
// LINALG/src
#include "fn_trans.hpp"
#include "fn_diag_vec.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_History_Bulk::IWG_History_Bulk()
        {
            // set size for the property pointer cell
            mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Weight" ] = static_cast< uint >( IWG_Property_Type::WEIGHT );
            mPropertyMap[ "Lump" ]   = static_cast< uint >( IWG_Property_Type::LUMP );

            // set size for the constitutive model pointer cell
            mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "ElasticDamage" ] = static_cast< uint >( IWG_Constitutive_Type::ELASTIC_DAMAGE );
        }

        //------------------------------------------------------------------------------

        void
        IWG_History_Bulk::compute_residual( real aWStar )
        {
            // check leader field interpolators
#ifdef MORIS_HAVE_DEBUG
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get residual dof type field interpolator
            Field_Interpolator* tFI =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get weight coefficient
            const std::shared_ptr< Property >& tPropWeight =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::WEIGHT ) );

            // get lump coefficient
            const std::shared_ptr< Property >& tPropLump =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::LUMP ) );

            // get the elasticity with damage CM
            const std::shared_ptr< Constitutive_Model >& tCMElasticityDamage =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELASTIC_DAMAGE ) );

            // cast constitutive model base class pointer to elasticity damage constitutive model
            CM_Struc_Linear_Isotropic_Damage* tCMElasticityDamagePtr =
                    dynamic_cast< CM_Struc_Linear_Isotropic_Damage* >( tCMElasticityDamage.get() );

            // if weight property specified
            real tWeight = 1.0;
            if( tPropWeight != nullptr )
            {
                // get weight value from property
                tWeight = tPropWeight->val()( 0 );
            }

            // get sub-matrix
            auto tRes = mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex } );

            // compute the residual
            tRes += aWStar * tWeight * ( tFI->N_trans() * ( tFI->val() - tCMElasticityDamagePtr->history() ) );

            // if lumping
            if ( tPropLump != nullptr )
            {
                // remove consistent contribution
                tRes -= aWStar * tWeight * tFI->N_trans() * tFI->val();

                // get lump parameter
                real tLump = tPropLump->val()( 0 );

                // add scaled consistent contribution
                tRes += aWStar * tWeight * tLump * tFI->N_trans() * tFI->val();

                // get number of space time coefficients
                uint tISize = tFI->get_number_of_space_time_coefficients();

                // grab diagonal of consistent contribution
                Matrix< DDRMat > tDiag    = diag_vec( tFI->N_trans() * tFI->N() );
                real             tDiagSum = sum( tDiag );

                // add lumped contribution
                for ( uint i = 0; i < tISize; i++ )
                {
                    tRes( i ) += aWStar * tWeight * ( 1.0 - tLump ) * tDiag( i ) * tFI->get_coeff()( i ) / tDiagSum;
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_History_Bulk::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_History_Bulk::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get field interpolator for a given dof type
            Field_Interpolator* tFI =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get weight coefficient
            const std::shared_ptr< Property >& tPropWeight =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::WEIGHT ) );

            // get lump coefficient
            const std::shared_ptr< Property >& tPropLump =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::LUMP ) );

            // get the elasticity with damage CM
            const std::shared_ptr< Constitutive_Model >& tCMElasticityDamage =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELASTIC_DAMAGE ) );

            // cast constitutive model base class pointer to elasticity damage constitutive model
            CM_Struc_Linear_Isotropic_Damage* tCMElasticityDamagePtr =
                    dynamic_cast< CM_Struc_Linear_Isotropic_Damage* >( tCMElasticityDamage.get() );

            // if weight property specified
            real tWeight = 1.0;
            if( tPropWeight != nullptr )
            {
                // get weight value from property
                tWeight = tPropWeight->val()( 0 );
            }

            // get the number of leader dof type dependencies
            uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            // loop over leader dof type dependencies
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                const Cell< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // get sub-matrix
                auto tJac = mSet->get_jacobian()(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } );

                // if derivative dof type is residual type
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // compute the Jacobian
                    tJac += aWStar * tWeight * ( tFI->N_trans() * tFI->N() );

                    // if lumping
                    if ( tPropLump != nullptr )
                    {
                        // remove consistent contribution
                        tJac -= aWStar * tWeight * tFI->N_trans() * tFI->N();

                        // get lump parameter
                        real tLump = tPropLump->val()( 0 );

                        // add scaled consistent contribution
                        tJac += aWStar * tWeight * tLump * tFI->N_trans() * tFI->N();

                        // get number of space time coefficients
                        uint tISize = tFI->get_number_of_space_time_coefficients();

                        // grab diagonal of consistent contribution
                        Matrix< DDRMat > tDiag    = diag_vec( tFI->N_trans() * tFI->N() );
                        real             tDiagSum = sum( tDiag );

                        // add lumped contribution
                        for ( uint i = 0; i < tISize; i++ )
                        {
                            tJac( i, i ) += aWStar * tWeight * ( 1.0 - tLump ) * tDiag( i ) / tDiagSum;
                        }
                    }
                }

                // if weight property specified
                if( tPropWeight != nullptr )
                {
                    // if property has dependency on the dof type
                    if ( tPropWeight->check_dof_dependency( tDofType ) )
                    {
                        MORIS_ERROR( false, "IWG_History_Bulk::compute_jacobian - dof dependency of property not supported" );
                    }
                }

                // if constitutive model has dependency on the dof type
                if ( tCMElasticityDamagePtr->check_dof_dependency( tDofType ) )
                {
                    // compute the Jacobian
                    tJac -= aWStar * tWeight * ( tFI->N_trans() * tCMElasticityDamagePtr->dHistorydu( tDofType ) );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_History_Bulk::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_History_Bulk::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_History_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_History_Bulk::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_History_Bulk::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

