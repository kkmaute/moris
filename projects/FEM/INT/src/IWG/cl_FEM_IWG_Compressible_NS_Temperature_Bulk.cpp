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
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "BodyForce" ]     = static_cast< uint >( IWG_Property_Type::BODY_FORCE );
            mPropertyMap[ "BodyHeatLoad" ]  = static_cast< uint >( IWG_Property_Type::BODY_HEAT_LOAD );

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Fluid" ] = static_cast< uint >( IWG_Constitutive_Type::FLUID );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Temperature_Bulk::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here velocity), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the Temperature FI
            Field_Interpolator * tFITemp =  mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get Velocity FI
            Field_Interpolator * tFIVelocity =  mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );

            // get properties
            std::shared_ptr< Property > tPropBodyForce = mMasterProp( static_cast< uint >( IWG_Property_Type::BODY_FORCE ) );
            std::shared_ptr< Property > tPropBodyHeatLoad = mMasterProp( static_cast< uint >( IWG_Property_Type::BODY_HEAT_LOAD ) );

            // get the compressible fluid constitutive model
            std::shared_ptr< Constitutive_Model > tCMFluid = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID ) );

            // compute the residual
            mSet->get_residual()( 0 )(
                    { tMasterResStartIndex, tMasterResStopIndex },
                    { 0, 0 } ) += aWStar * (
                            trans( tFITemp->N() ) * tCMFluid->EnergyDot() +
                            trans( tFITemp->dnNdxn( 1 ) ) *
                            ( tCMFluid->flux( CM_Function_Type::WORK ) -
                              tCMFluid->flux( CM_Function_Type::ENERGY ) -
                              tCMFluid->flux( CM_Function_Type::THERMAL ) ) );

            // if there is a body force
            if ( tPropBodyForce != nullptr )
            {
                // get Velocity FI
                Field_Interpolator * tFIDensity =  mMasterFIManager->get_field_interpolators_for_type( mDofDensity );

                // add contribution
                mSet->get_residual()( 0 )(
                        { tMasterResStartIndex, tMasterResStopIndex },
                        { 0, 0 } ) -= aWStar * (
                                tFIDensity->val()( 0 ) * trans( tFITemp->N() ) * dot( tPropBodyForce->val(), tFIVelocity->val() ) );
            }

            // if there is a body heat load
            if ( tPropBodyHeatLoad != nullptr )
            {
                // add contribution
                mSet->get_residual()( 0 )(
                        { tMasterResStartIndex, tMasterResStopIndex },
                        { 0, 0 } ) -= aWStar * (
                                trans( tFITemp->N() ) * tPropBodyHeatLoad->val() );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Compressible_NS_Temperature_Bulk::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------
        void IWG_Compressible_NS_Temperature_Bulk::compute_jacobian( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here pressure), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the field interpolators
            Field_Interpolator * tFITemp = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));
            Field_Interpolator * tFIVelocity = mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );
            Field_Interpolator * tFIDensity = mMasterFIManager->get_field_interpolators_for_type( mDofDensity );

            // get the constitutive model
            std::shared_ptr< Constitutive_Model > tCMFluid =  mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID ) );

            // get the properties
            std::shared_ptr< Property > tPropBodyForce = mMasterProp( static_cast< uint >( IWG_Property_Type::BODY_FORCE ) );
            std::shared_ptr< Property > tPropBodyHeatLoad = mMasterProp( static_cast< uint >( IWG_Property_Type::BODY_HEAT_LOAD ) );

            // compute the jacobian for dof dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // if fluid CM depends on dof type
                if ( tCMFluid->check_dof_dependency( tDofType ) )
                {
                    // add contribution
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    trans( tFITemp->N() ) * tCMFluid->dEnergyDotdDOF( tDofType ) +
                                    trans( tFITemp->dnNdxn( 1 ) ) *
                                    ( tCMFluid->dFluxdDOF( tDofType, CM_Function_Type::WORK ) -
                                      tCMFluid->dFluxdDOF( tDofType, CM_Function_Type::ENERGY ) -
                                      tCMFluid->dFluxdDOF( tDofType, CM_Function_Type::THERMAL ) ) );
                }

                // if a body force is present
                if ( tPropBodyForce != nullptr )
                {
                    // if the body force depends on the dof type -> indirect dependency
                    if ( tPropBodyForce->check_dof_dependency( tDofType ) )
                    {
                        // compute the jacobian contribution
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                        tFIDensity->val()( 0 ) * trans( tFITemp->N() ) * trans( tFIVelocity->val() ) *
                                        tPropBodyForce->dPropdDOF( tDofType ) );
                    }

                    // if dof type is density and a body force is present
                    if( tDofType( 0 ) == mDofDensity )
                    {
                        // compute the jacobian contribution
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                        trans( tFITemp->N() ) * dot( tPropBodyForce->val(), tFIVelocity->val() ) * tFIDensity->N() );
                    }

                    // if dof type is velocity and a body force is present
                    if( tDofType( 0 ) == mDofVelocity )
                    {
                        // compute the jacobian contribution
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                        tFIDensity->val()( 0 ) * trans( tFITemp->N() ) * trans( tPropBodyForce->val() ) * tFIVelocity->N() );
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
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                        trans( tFITemp->N() ) * tPropBodyHeatLoad->dPropdDOF( tDofType ) );
                    }
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Compressible_NS_Temperature_Bulk::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Temperature_Bulk::compute_jacobian_and_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Compressible_NS_Temperature_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Temperature_Bulk::compute_dRdp( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Compressible_NS_Temperature_Bulk::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Temperature_Bulk::compute_residual_strong_form(
                Matrix< DDRMat > & aRM,
                real             & aRC )
        {
            MORIS_ERROR( false, "IWG_Compressible_NS_Temperature_Bulk::compute_residual_strong_form - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Temperature_Bulk::compute_jacobian_strong_form(
                moris::Cell< MSI::Dof_Type >   aDofTypes,
                Matrix< DDRMat >             & aJM,
                Matrix< DDRMat >             & aJC )
        {
            MORIS_ERROR( false, "IWG_Compressible_NS_Temperature_Bulk::compute_jacobian_strong_form - Not implemented." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

