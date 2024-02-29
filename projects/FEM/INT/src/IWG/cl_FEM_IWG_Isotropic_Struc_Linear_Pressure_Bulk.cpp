/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Linear_Pressure_Bulk.cpp
 *
 */

#include "cl_FEM_IWG_Isotropic_Struc_Linear_Pressure_Bulk.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_CM_Struc_Linear_Isotropic.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IWG_Isotropic_Struc_Linear_Pressure_Bulk::IWG_Isotropic_Struc_Linear_Pressure_Bulk()
        {
            // set size for the property pointer cell
            mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Thickness" ] = static_cast< uint >( IWG_Property_Type::THICKNESS );

            // set size for the constitutive model pointer cell
            mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "ElastLinIso" ] = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Pressure_Bulk::compute_residual( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here pressure), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get field interpolator for dof type
            Field_Interpolator* tPressureFI     =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));
            Field_Interpolator* tDisplacementFI =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

            // get elasticity CM
            moris::fem::CM_Struc_Linear_Isotropic* tCMElasticity =
                    static_cast < moris::fem::CM_Struc_Linear_Isotropic* > ( mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) ).get());

            // get thickness property
            const std::shared_ptr< Property > & tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= (tPropThickness!=nullptr) ? tPropThickness->val()(0) : 1;

            // compute the residual
            mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex },
                    { 0, 0 } ) += aWStar * (
                    trans( tPressureFI->N() ) * ( tDisplacementFI->div() + tCMElasticity->eval_inv_bulk_modulus() * tPressureFI->val()( 0 ) ) );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Isotropic_Struc_Linear_Pressure_Bulk::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Pressure_Bulk::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here pressure), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get field interpolator for residual dof type
            Field_Interpolator* tPressureFI = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get elasticity CM
            moris::fem::CM_Struc_Linear_Isotropic* tCMElasticity =
                    static_cast < moris::fem::CM_Struc_Linear_Isotropic* > ( mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) ).get());

            // get thickness property
            const std::shared_ptr< Property > & tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= (tPropThickness!=nullptr) ? tPropThickness->val()(0) : 1;

            // compute the jacobian for dof dependencies
            uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Vector< MSI::Dof_Type > & tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // get FI for treated dof type
                Field_Interpolator * tFIDer = mLeaderFIManager->get_field_interpolators_for_type( tDofType( 0 ) );

                // if dof type is pressure
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * (
                                    trans( tPressureFI->N() ) * tCMElasticity->eval_inv_bulk_modulus() * tPressureFI->N() );
                }

                // if dof type is displacement
                if ( tDofType( 0 ) == MSI::Dof_Type::UX )
                {
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * (
                                    trans( tPressureFI->N() ) * tFIDer->div_operator() );
                }

                // if elasticity CM depends on dof type
                if ( tCMElasticity->check_dof_dependency( tDofType ) )
                {
                    mSet->get_jacobian()(
                            { tLeaderResStartIndex, tLeaderResStopIndex },
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * (
                                    trans( tPressureFI->N() ) * tCMElasticity->eval_dInvBulkModulusdDOF( tDofType ) * tPressureFI->val()( 0 ) );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Isotropic_Struc_Linear_Pressure_Bulk::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Pressure_Bulk::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Pressure_Bulk::compute_jacobian_and_residual - Not implemented.");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Pressure_Bulk::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Pressure_Bulk::compute_dRdp - Not implemented.");
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

