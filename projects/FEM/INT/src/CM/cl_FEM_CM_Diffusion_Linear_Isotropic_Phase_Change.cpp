/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Diffusion_Linear_Isotropic_Phase_Change.cpp
 *
 */

#include "cl_FEM_CM_Diffusion_Linear_Isotropic_Phase_Change.hpp"
#include "cl_FEM_CM_Diffusion_Linear_Isotropic.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include <iostream>

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"
#include "fn_sum.hpp"
#include "fn_FEM_CM_Phase_State_Functions.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        CM_Diffusion_Linear_Isotropic_Phase_Change::CM_Diffusion_Linear_Isotropic_Phase_Change()
        {
            // set the property pointer cell size
            mProperties.resize( static_cast< uint >( CM_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Conductivity" ]        = static_cast< uint >( CM_Property_Type::CONDUCTIVITY );
            mPropertyMap[ "Density" ]             = static_cast< uint >( CM_Property_Type::DENSITY );
            mPropertyMap[ "HeatCapacity" ]        = static_cast< uint >( CM_Property_Type::HEAT_CAPACITY );
            mPropertyMap[ "LatentHeat" ]          = static_cast< uint >( CM_Property_Type::LATENT_HEAT );
            mPropertyMap[ "PCTemp" ]              = static_cast< uint >( CM_Property_Type::PC_TEMP );
            mPropertyMap[ "PhaseStateFunction" ]  = static_cast< uint >( CM_Property_Type::PHASE_STATE_FUNCTION );
            mPropertyMap[ "PhaseChangeConst" ]    = static_cast< uint >( CM_Property_Type::PHASE_CHANGE_CONST );
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Phase_Change::set_dof_type_list(
                moris::Vector< moris::Vector< MSI::Dof_Type > > aDofTypes,
                moris::Vector< std::string >                  aDofStrings )
        {
            // set dof type list
            Constitutive_Model::set_dof_type_list( aDofTypes );

            // loop over the provided dof types
            for( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
            {
                // get dof type string
                std::string tDofString = aDofStrings( iDof );

                // get dof type
                MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

                // switch on dof type string
                if( tDofString == "Temperature" )
                {
                    mTempDof = tDofType;
                }
                else
                {
                    // error unknown dof string
                    MORIS_ERROR( false,
                            "CM_Diffusion_Linear_Isotropic_Phase_Change::set_dof_type_list - Unknown aDofString : %s \n",
                            tDofString.c_str() );
                }
            }
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Phase_Change::set_local_properties()
        {
            // set the conductivity property
            mPropConductivity = get_property( "Conductivity" );

            // set the heat capacity property
            mPropHeatCapacity = get_property( "HeatCapacity" );

            // set the density property
            mPropDensity = get_property( "Density" );

            // set the latent heat property
            mPropLatentHeat = get_property( "LatentHeat" );

            // set the PC temperature property
            mPropPCTemp = get_property( "PCTemp" );

            // set the PC constant property
            mPropPCConst = get_property( "PhaseChangeConst" );

            // set the PS function property
            mPropPSFunc = get_property( "PhaseStateFunction" );
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_Energy()
        {
            // get the temperature FI
            Field_Interpolator * tFITemp =
                    mFIManager->get_field_interpolators_for_type( mTempDof );

            // compute derivative of Phase State Function
            real tPhaseState = eval_phase_state_function(
                    mPropPCTemp->val()( 0 ),
                    mPropPCConst->val()( 0 ),
                    mPropPSFunc->val()( 0 ),
                    tFITemp );

            // compute enthalpy
            mEnergy = mPropDensity->val()( 0 ) * (
                    + mPropHeatCapacity->val()( 0 ) * tFITemp->val( )
                    + mPropLatentHeat->val()( 0 ) * tPhaseState );
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_EnergyDot()
        {
            // get the temperature FI
            Field_Interpolator * tFITemp =
                    mFIManager->get_field_interpolators_for_type( mTempDof );

            // compute derivative of Phase State Function
            real tdfdT = eval_dFdTemp(
                    mPropPCTemp->val()( 0 ),
                    mPropPCConst->val()( 0 ),
                    mPropPSFunc->val()( 0 ),
                    tFITemp );

            // compute derivative of enthalpy
            mEnergyDot = mPropDensity->val()( 0 ) * (
                    + mPropHeatCapacity->val()( 0 )
                    + mPropLatentHeat->val()( 0 ) * tdfdT ) *
                    tFITemp->gradt( 1 );
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_gradEnergyDot()
        {
            // get the temperature FI
            Field_Interpolator * tFITemp =
                    mFIManager->get_field_interpolators_for_type( mTempDof );

            // compute derivative of Phase State Function
            real tdfdT = eval_dFdTemp(
                    mPropPCTemp->val()( 0 ),
                    mPropPCConst->val()( 0 ),
                    mPropPSFunc->val()( 0 ),
                    tFITemp );

            // compute gradient of
            mGradEnergyDot = mPropDensity->val()( 0 ) *
                    ( mPropHeatCapacity->val()( 0 ) + mPropLatentHeat->val()( 0 ) * tdfdT ) *
                    tFITemp->gradxt();
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_gradEnergy()
        {
            // get the temperature FI
            Field_Interpolator * tFITemp =
                    mFIManager->get_field_interpolators_for_type( mTempDof );

            // compute derivative of Phase State Function
            real tdfdT = eval_dFdTemp(
                    mPropPCTemp->val()( 0 ),
                    mPropPCConst->val()( 0 ),
                    mPropPSFunc->val()( 0 ),
                    tFITemp );

            // compute gradient of
            mGradEnergy = mPropDensity->val()( 0 ) *
                    ( mPropHeatCapacity->val()( 0 ) + mPropLatentHeat->val()( 0 ) * tdfdT ) *
                    tFITemp->gradx( 1 );
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_dEnergydDOF(
                const moris::Vector< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the corresponding FI
            Field_Interpolator * tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the matrix
            mEnergyDof( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // check if density and heat capacity are set
            if ( mPropDensity == nullptr || mPropHeatCapacity == nullptr)
            {
                return;
            }

            // get the temperature FI
            Field_Interpolator * tFITemp =
                    mFIManager->get_field_interpolators_for_type( mTempDof );

            // if the dof type is temperature
            if( aDofTypes( 0 ) == mTempDof )
            {
                // compute derivative of Phase State Function
                const real tdfdT = eval_dFdTemp(
                        mPropPCTemp->val()( 0 ),
                        mPropPCConst->val()( 0 ),
                        mPropPSFunc->val()( 0 ),
                        tFITemp );

                // compute derivative with direct dependency
                mEnergyDof( tDofIndex ) +=
                        mPropDensity->val()( 0 ) * (
                                + mPropHeatCapacity->val()( 0 )
                                + mPropLatentHeat->val()( 0 ) * tdfdT ) *
                                tFITemp->N();
            }

            // if density property depends on the dof type
            if ( mPropDensity->check_dof_dependency( aDofTypes ) )
            {

                // compute derivative of Phase State Function
                real tPhaseState = eval_phase_state_function(
                        mPropPCTemp->val()( 0 ),
                        mPropPCConst->val()( 0 ),
                        mPropPSFunc->val()( 0 ),
                        tFITemp );

                // compute derivative with indirect dependency through properties
                mEnergyDof( tDofIndex ) += (
                        + mPropHeatCapacity->val()( 0 ) * tFITemp->val( )
                        + mPropLatentHeat->val()( 0 ) * tPhaseState ) *
                        mPropDensity->dPropdDOF( aDofTypes );
            }

            // if heat capacity depends on the dof type
            if ( mPropHeatCapacity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mEnergyDof( tDofIndex ) +=
                        mPropDensity->val()( 0 ) *
                        tFITemp->val() *
                        mPropHeatCapacity->dPropdDOF( aDofTypes );
            }

            // if latent heat depends on the dof type
            if ( mPropLatentHeat->check_dof_dependency( aDofTypes ) )
            {
                MORIS_ERROR(false, "CM_Diffusion_Linear_Isotropic_Phase_Change::eval_dEnergydDOF - %s\n",
                        "Dof dependence of Latent heat not implemented.\n");
            }
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_dEnergyDotdDOF(
                const moris::Vector< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the corresponding FI
            Field_Interpolator * tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the matrix
            mEnergyDotDof( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // get the temperature FI
            Field_Interpolator * tFITemp =
                    mFIManager->get_field_interpolators_for_type( mTempDof );

            // compute derivative of Phase State Function
            real tdfdT = eval_dFdTemp(
                    mPropPCTemp->val()( 0 ),
                    mPropPCConst->val()( 0 ),
                    mPropPSFunc->val()( 0 ),
                    tFITemp );

            // if direct dependency on the dof type
            if( aDofTypes( 0 ) == mTempDof )
            {
                const moris::Matrix< DDRMat > dfdDof = eval_dFdTempdDOF(
                        mPropPCTemp->val()( 0 ),
                        mPropPCConst->val()( 0 ),
                        mPropPSFunc->val()( 0 ),
                        tFITemp );

                // compute derivative with direct dependency
                mEnergyDotDof( tDofIndex ) +=
                        mPropDensity->val()( 0 ) * ( mPropHeatCapacity->val()( 0 ) + mPropLatentHeat->val()( 0 ) * tdfdT ) *
                        tFITemp->dnNdtn( 1 )
                        + mPropDensity->val()( 0 ) * mPropLatentHeat->val()( 0 ) * tFITemp->gradt( 1 ) * dfdDof;
            }

            // if density depends on the dof type
            if ( mPropDensity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mEnergyDotDof( tDofIndex ) +=
                        ( mPropHeatCapacity->val()( 0 ) + mPropLatentHeat->val()( 0 ) * tdfdT ) *
                        tFITemp->gradt( 1 ) *
                        mPropDensity->dPropdDOF( aDofTypes );
            }

            // if heat capacity depends on the dof type
            if ( mPropHeatCapacity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mEnergyDotDof( tDofIndex ) +=
                        mPropDensity->val()( 0 ) *
                        tFITemp->gradt( 1 ) *
                        mPropHeatCapacity->dPropdDOF( aDofTypes );
            }

            // if latent heat depends on the dof type
            if ( mPropLatentHeat->check_dof_dependency( aDofTypes ) )
            {
                MORIS_ERROR(false, "CM_Diffusion_Linear_Isotropic_Phase_Change::eval_dEnergyDotdDOF - %s\n",
                        "Dof dependence of Latent heat not implemented.\n");
            }
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_dGradEnergydDOF(
                const moris::Vector< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the corresponding FI
            Field_Interpolator * tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the matrix
            mGradEnergyDof( tDofIndex ).set_size(
                    mSpaceDim,
                    tFIDer->get_number_of_space_time_coefficients(),
                    0.0 );

            // get the temperature FI
            Field_Interpolator * tFITemp =
                    mFIManager->get_field_interpolators_for_type( mTempDof );

            // compute derivative of Phase State Function
            real tdfdT = eval_dFdTemp(
                    mPropPCTemp->val()( 0 ),
                    mPropPCConst->val()( 0 ),
                    mPropPSFunc->val()( 0 ),
                    tFITemp );

            // if direct dependency on the dof type
            if( aDofTypes( 0 ) == mTempDof )
            {
                const moris::Matrix<DDRMat> dfdDof = eval_dFdTempdDOF(
                        mPropPCTemp->val()( 0 ),
                        mPropPCConst->val()( 0 ),
                        mPropPSFunc->val()( 0 ),
                        tFITemp );

                // compute derivative with direct dependency
                mGradEnergyDof( tDofIndex ) +=
                        mPropDensity->val()( 0 ) * ( mPropHeatCapacity->val()(0) + mPropLatentHeat->val()( 0 ) * tdfdT ) *
                        tFITemp->dnNdxn( 1 )
                        + mPropDensity->val()( 0 ) * mPropLatentHeat->val()( 0 ) * tFITemp->gradx( 1 ) * dfdDof;
            }

            // if density depends on the dof type
            if ( mPropDensity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mGradEnergyDof( tDofIndex ) +=
                        ( mPropHeatCapacity->val()( 0 ) + mPropLatentHeat->val()( 0 ) * tdfdT ) *
                        tFITemp->gradx( 1 ) *
                        mPropDensity->dPropdDOF( aDofTypes );
            }

            // if heat capacity depends on the dof type
            if ( mPropHeatCapacity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mGradEnergyDof( tDofIndex ) +=
                        mPropDensity->val()( 0 ) *
                        tFITemp->gradx( 1 ) *
                        mPropHeatCapacity->dPropdDOF( aDofTypes );
            }

            // if latent heat depends on the dof type
            if ( mPropLatentHeat->check_dof_dependency( aDofTypes ) )
            {
                MORIS_ERROR(false, "CM_Diffusion_Linear_Isotropic_Phase_Change::eval_dGradEnergydDOF - %s\n",
                        "Dof dependence of Latent heat not implemented.\n");
            }
        }

        //------------------------------------------------------------------------------

        void CM_Diffusion_Linear_Isotropic_Phase_Change::eval_dGradEnergyDotdDOF(
                const moris::Vector< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the corresponding FI
            Field_Interpolator * tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the matrix
            mGradEnergyDotDof( tDofIndex ).set_size(
                    mSpaceDim,
                    tFIDer->get_number_of_space_time_coefficients(),
                    0.0 );

            // get the temperature FI
            Field_Interpolator * tFITemp =
                    mFIManager->get_field_interpolators_for_type( mTempDof );

            // compute derivative of Phase State Function
            // real tdfdT = this->eval_dFdTemp();
            real tdfdT = eval_dFdTemp(
                    mPropPCTemp->val()( 0 ),
                    mPropPCConst->val()( 0 ),
                    mPropPSFunc->val()( 0 ),
                    tFITemp );

            // if direct dependency on the dof type
            if( aDofTypes( 0 ) == mTempDof )
            {
                const moris::Matrix<DDRMat> dfdDof = eval_dFdTempdDOF(
                        mPropPCTemp->val()( 0 ),
                        mPropPCConst->val()( 0 ),
                        mPropPSFunc->val()( 0 ),
                        tFITemp );

                // compute derivative with direct dependency
                mGradEnergyDotDof( tDofIndex ) +=
                        mPropDensity->val()( 0 ) * ( mPropHeatCapacity->val()(0) + mPropLatentHeat->val()( 0 ) * tdfdT ) *
                        tFITemp->d2Ndxt()
                        + mPropDensity->val()( 0 ) * mPropLatentHeat->val()( 0 ) * tFITemp->gradxt() * dfdDof;
            }

            // if density depends on the dof type
            if ( mPropDensity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mGradEnergyDotDof( tDofIndex ) +=
                        ( mPropHeatCapacity->val()( 0 ) + mPropLatentHeat->val()( 0 ) * tdfdT ) *
                        tFITemp->gradxt() *
                        mPropDensity->dPropdDOF( aDofTypes );
            }

            // if heat capacity depends on the dof type
            if ( mPropHeatCapacity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mGradEnergyDotDof( tDofIndex ) +=
                        mPropDensity->val()( 0 ) *
                        tFITemp->gradxt() *
                        mPropHeatCapacity->dPropdDOF( aDofTypes );
            }

            // if latent heat depends on the dof type
            if ( mPropLatentHeat->check_dof_dependency( aDofTypes ) )
            {
                MORIS_ERROR(false, "CM_Diffusion_Linear_Isotropic_Phase_Change::eval_dGradEnergyDotdDOF - %s\n",
                        "Dof dependence of Latent heat not implemented.\n");
            }
       }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

