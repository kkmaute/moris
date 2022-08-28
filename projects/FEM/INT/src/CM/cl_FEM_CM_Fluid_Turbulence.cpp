/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Fluid_Turbulence.cpp
 *
 */

#include "cl_FEM_CM_Fluid_Turbulence.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"
// LINALG/src
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "op_minus.hpp"

namespace moris
{
    namespace fem
    {

        //--------------------------------------------------------------------------------------------------------------

        CM_Fluid_Turbulence::CM_Fluid_Turbulence()
        {
            // set the property pointer cell size
            mProperties.resize( static_cast< uint >( CM_Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "Density" ]      = static_cast< uint >( CM_Property_Type::DENSITY );
            mPropertyMap[ "Viscosity" ]    = static_cast< uint >( CM_Property_Type::VISCOSITY );
            mPropertyMap[ "KinViscosity" ] = static_cast< uint >( CM_Property_Type::KIN_VISCOSITY );

            // FIXME for now only 1st order allowed
            uint tOrder = 1;

            // init storage for evaluation
            mdChidx.resize( tOrder );
            mdFv1dx.resize( tOrder );
            mdTurbDynViscdx.resize( tOrder );
            mdEffDynViscdx.resize( tOrder );

            // init flag for evaluation
            mdChidxEval.set_size( tOrder, 1, true );
            mdFv1dxEval.set_size( tOrder, 1, true );
            mdTurbDynViscdxEval.set_size( tOrder, 1, true );
            mdEffDynViscdxEval.set_size( tOrder, 1, true );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Fluid_Turbulence::reset_eval_flags()
        {
            // call parent implementation
            Constitutive_Model::reset_eval_flags();

            // reset child specific eval flags for chi
            mChiEval = true;
            mdChiduEval.fill( true );
            mdChidxEval.fill( true );
            mdChidxduEval.fill( true );

            // reset child specific eval flags for fv1
            mFv1Eval = true;
            mdFv1duEval.fill( true );
            mdFv1dxEval.fill( true );
            mdFv1dxduEval.fill( true );

            // reset child specific eval flags for turbulence dynamic viscosity
            mTurbDynViscEval = true;
            mdTurbDynViscduEval.fill( true );
            mdTurbDynViscdxEval.fill( true );
            mdTurbDynViscdxduEval.fill( true );

            // reset child specific eval flags for effective dynamic viscosity
            mEffDynViscEval = true;
            mdEffDynViscduEval.fill( true );
            mdEffDynViscdxEval.fill( true );
            mdEffDynViscdxduEval.fill( true );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Fluid_Turbulence::build_global_dof_type_list()
        {
            // call parent implementation
            Constitutive_Model::build_global_dof_type_list();

            // get number of dof types
            uint tNumGlobalDofTypes = mGlobalDofTypes.size();

            // init child specific eval flags
            mdChiduEval.set_size( tNumGlobalDofTypes, 1, true );
            mdFv1duEval.set_size( tNumGlobalDofTypes, 1, true );
            mdTurbDynViscduEval.set_size( tNumGlobalDofTypes, 1, true );
            mdEffDynViscduEval.set_size( tNumGlobalDofTypes, 1, true );

            // FIXME for now only 1st order allowed
            uint tOrder = 1;
            mdChidxduEval.set_size( tOrder, tNumGlobalDofTypes, true );
            mdFv1dxduEval.set_size( tOrder, tNumGlobalDofTypes, true );
            mdEffDynViscdxduEval.set_size( tOrder, tNumGlobalDofTypes, true );
            mdTurbDynViscdxduEval.set_size( tOrder, tNumGlobalDofTypes, true );

            // init child specific storage
            mdChidu.resize( tNumGlobalDofTypes );
            mdFv1du.resize( tNumGlobalDofTypes );
            mdTurbDynViscdu.resize( tNumGlobalDofTypes );
            mdEffDynViscdu.resize( tNumGlobalDofTypes );

            mdChidxdu.resize( tOrder );
            mdFv1dxdu.resize( tOrder );
            mdEffDynViscdxdu.resize( tOrder );
            mdTurbDynViscdxdu.resize( tOrder );
            for ( uint iOrder = 0; iOrder < tOrder; iOrder++ )
            {
                mdChidxdu( iOrder ).resize( tNumGlobalDofTypes );
                mdFv1dxdu( iOrder ).resize( tNumGlobalDofTypes );
                mdEffDynViscdxdu( iOrder ).resize( tNumGlobalDofTypes );
                mdTurbDynViscdxdu( iOrder ).resize( tNumGlobalDofTypes );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Fluid_Turbulence::set_dof_type_list(
                moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes,
                moris::Cell< std::string >                  aDofStrings )
        {
            // set dof type list
            Constitutive_Model::set_dof_type_list( aDofTypes );

            // loop over the provided dof type
            for ( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
            {
                // get dof type string
                std::string tDofString = aDofStrings( iDof );

                // get dof type
                MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

                // switch on dof type string
                if ( tDofString == "Velocity" )
                {
                    mDofVelocity = tDofType;
                }
                else if ( tDofString == "Pressure" )
                {
                    mDofPressure = tDofType;
                }
                else if ( tDofString == "Viscosity" )
                {
                    mDofViscosity = tDofType;
                }
                else
                {
                    // error unknown dof string
                    MORIS_ERROR( false,
                            "CM_Fluid_Turbulence::set_dof_type_list - Unknown aDofString : %s \n",
                            tDofString.c_str() );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Fluid_Turbulence::set_local_properties()
        {
            // set the density property
            mPropDensity = get_property( "Density" );

            // set the dynamic viscosity property
            mPropViscosity = get_property( "Viscosity" );

            // set the kinematic viscosity property
            mPropKinViscosity = get_property( "KinViscosity" );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Fluid_Turbulence::eval_flux()
        {
            // get the pressure FI
            Field_Interpolator* tPressureFI =
                    mFIManager->get_field_interpolators_for_type( mDofPressure );

            // create identity matrix
            Matrix< DDRMat > tI( mSpaceDim, 1, 1.0 );

            // evaluate pressure contribution to flux
            Matrix< DDRMat > tP( ( mSpaceDim - 1 ) * 3, 1, 0.0 );

            tP( { 0, mSpaceDim - 1 }, { 0, 0 } ) = tI * tPressureFI->val();

            // compute flux
            mFlux = -1.0 * tP + 2.0 * this->effective_dynamic_viscosity()( 0 ) * this->strain();
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Fluid_Turbulence::eval_dFluxdDOF(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof FI
            Field_Interpolator* tFI =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init mdFluxdDof
            mdFluxdDof( tDofIndex ).set_size( ( mSpaceDim - 1 ) * 3, tFI->get_number_of_space_time_coefficients() );

            // if velocity dof
            if ( aDofTypes( 0 ) == mDofVelocity )
            {
                // build dfluxdv
                mdFluxdDof( tDofIndex ) =
                        2.0 * this->effective_dynamic_viscosity()( 0 ) * this->dStraindDOF( aDofTypes );
            }
            // if pressure dof
            else if ( aDofTypes( 0 ) == mDofPressure )
            {
                // create identity matrix
                Matrix< DDRMat > tI( mSpaceDim, 1, 1.0 );
                Matrix< DDRMat > tII( ( mSpaceDim - 1 ) * 3, 1, 0.0 );

                tII( { 0, mSpaceDim - 1 }, { 0, 0 } ) = tI.matrix_data();

                // build the dfluxdp
                mdFluxdDof( tDofIndex ) = -1.0 * tII * tFI->N();
            }
            else
            {
                mdFluxdDof( tDofIndex ).fill( 0.0 );
            }

            // add contribution from the derivative of the effective dynamic viscosity wrt dof
            mdFluxdDof( tDofIndex ) +=
                    2.0 * this->strain() * this->deffdynviscdu( aDofTypes );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Fluid_Turbulence::eval_dfluxdx( uint aOrder )
        {
            // only 1st order supported
            MORIS_ERROR( aOrder == 1, "CM_Fluid_Turbulence::eval_dfluxdx - only 1st order supported." );

            // get the pressure FI
            Field_Interpolator* tPressureFI =
                    mFIManager->get_field_interpolators_for_type( mDofPressure );

            // create identity matrix
            Matrix< DDRMat > tI( mSpaceDim, 1, 1.0 );
            Matrix< DDRMat > tP( ( mSpaceDim - 1 ) * 3, 1, 0.0 );

            tP( { 0, mSpaceDim - 1 }, { 0, 0 } ) = tI.matrix_data();

            // flatten deffdynviscdx
            Matrix< DDRMat > tdEffDynViscdxFlat;
            this->flatten_normal( this->deffdynviscdx( aOrder ), tdEffDynViscdxFlat );

            // evaluate dfluxdx
            mdFluxdx( aOrder - 1 ) =
                    trans( tPressureFI->gradx( aOrder ) * tP ) - 2.0 * this->effective_dynamic_viscosity() * this->dstraindx( aOrder ) + 2.0 * tdEffDynViscdxFlat * this->strain();
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Fluid_Turbulence::eval_divflux()
        {
            // get the pressure FI
            Field_Interpolator* tPressureFI =
                    mFIManager->get_field_interpolators_for_type( mDofPressure );

            // flatten deffdynviscdx
            Matrix< DDRMat > tdEffDynViscdxFlat;
            this->flatten_normal( this->deffdynviscdx( 1 ), tdEffDynViscdxFlat );

            // compute flux
            mDivFlux = -1.0 * tPressureFI->gradx( 1 ) + 2.0 * this->effective_dynamic_viscosity()( 0 ) * this->divstrain() + 2.0 * tdEffDynViscdxFlat * this->strain();
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Fluid_Turbulence::eval_ddivfluxdu( const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the corresponding FI
            Field_Interpolator* tFI =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for ddivflux/du
            mddivfluxdu( tDofIndex ).set_size( mSpaceDim, tFI->get_number_of_space_time_coefficients() );

            // if velocity dof
            if ( aDofTypes( 0 ) == mDofVelocity )
            {
                // flatten deffdynviscdx
                Matrix< DDRMat > tdEffDynViscdxFlat;
                this->flatten_normal( this->deffdynviscdx( 1 ), tdEffDynViscdxFlat );

                // fill ddivstrain/dv
                mddivfluxdu( tDofIndex ) =
                        2.0 * this->effective_dynamic_viscosity()( 0 ) * this->ddivstraindu( aDofTypes ) + 2.0 * tdEffDynViscdxFlat * this->dStraindDOF( aDofTypes );
            }
            // if pressure dof
            else if ( aDofTypes( 0 ) == mDofPressure )
            {
                // fill ddivstrain/dp
                mddivfluxdu( tDofIndex ) = -1.0 * tFI->dnNdxn( 1 );
            }
            else
            {
                mddivfluxdu( tDofIndex ).fill( 0.0 );
            }

            // add contribution from derivative of effective dynamic viscosity wrt dof
            mddivfluxdu( tDofIndex ) +=
                    2.0 * this->divstrain() * this->deffdynviscdu( aDofTypes );

            // get the full strain tensor
            const Matrix< DDRMat >& tStrain = this->strain();
            Matrix< DDRMat >        tStrainFull;
            switch ( mSpaceDim )
            {
                case 2:
                {
                    tStrainFull = {
                        { tStrain( 0 ), tStrain( 2 ) },
                        { tStrain( 2 ), tStrain( 1 ) }
                    };
                    break;
                }

                case 3:
                {
                    tStrainFull = {
                        { tStrain( 0 ), tStrain( 5 ), tStrain( 4 ) },
                        { tStrain( 5 ), tStrain( 1 ), tStrain( 3 ) },
                        { tStrain( 4 ), tStrain( 3 ), tStrain( 2 ) }
                    };
                    break;
                }

                default:
                    MORIS_ERROR( false, "CM_Fluid_Turbulence::eval_ddivfluxdu - only 2 or 3D" );
            }

            // add the derivative of gradx of the effective dynamic viscosity wrt dof
            mddivfluxdu( tDofIndex ) += 2.0 * tStrainFull * this->deffdynviscdxdu( aDofTypes, 1 );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Fluid_Turbulence::eval_traction( const Matrix< DDRMat >& aNormal )
        {
            // flatten the normal
            Matrix< DDRMat > tFlatNormal;
            this->flatten_normal( aNormal, tFlatNormal );

            // compute the traction
            mTraction = tFlatNormal * this->flux();
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Fluid_Turbulence::eval_dTractiondDOF(
                const moris::Cell< MSI::Dof_Type >& aDofTypes,
                const Matrix< DDRMat >&             aNormal )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // flatten normal
            Matrix< DDRMat > tFlatNormal;
            this->flatten_normal( aNormal, tFlatNormal );

            // compute dtractiondu
            mdTractiondDof( tDofIndex ) = tFlatNormal * this->dFluxdDOF( aDofTypes );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Fluid_Turbulence::eval_testTraction(
                const Matrix< DDRMat >&             aNormal,
                const moris::Cell< MSI::Dof_Type >& aTestDofTypes )
        {
            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // flatten the normal
            Matrix< DDRMat > tFlatNormal;
            this->flatten_normal( aNormal, tFlatNormal );

            // get the dof FI
            Field_Interpolator* tFITest =
                    mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) );

            // init mdFluxdDof
            mTestTraction( tTestDofIndex ).set_size( mSpaceDim, tFITest->get_number_of_space_time_coefficients() );

            // add contribution from derivative of effective dynamic viscosity wrt test dof
            mTestTraction( tTestDofIndex ) =
                    2.0 * tFlatNormal * this->strain() * this->deffdynviscdu( aTestDofTypes );

            // if test traction wrt velocity
            if ( aTestDofTypes( 0 ) == mDofVelocity )
            {
                // compute test traction wrt velocity
                mTestTraction( tTestDofIndex ) +=
                        2.0 * this->effective_dynamic_viscosity()( 0 ) * tFlatNormal * this->testStrain();
            }
            // if test traction wrt pressure
            else if ( aTestDofTypes( 0 ) == mDofPressure )
            {
                // create identity matrix
                Matrix< DDRMat > tI( mSpaceDim, 1, 1.0 );
                Matrix< DDRMat > tII( ( mSpaceDim - 1 ) * 3, 1, 0.0 );
                tII( { 0, mSpaceDim - 1 }, { 0, 0 } ) = tI.matrix_data();

                // build the dtesttractiondP
                mTestTraction( tTestDofIndex ) -= tFlatNormal * tII * tFITest->N();
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Fluid_Turbulence::eval_dTestTractiondDOF(
                const moris::Cell< MSI::Dof_Type >& aDofTypes,
                const Matrix< DDRMat >&             aNormal,
                const Matrix< DDRMat >&             aJump,
                const moris::Cell< MSI::Dof_Type >& aTestDofTypes )
        {
            // get test dof type index
            const uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the test dof FI
            Field_Interpolator* tFITest =
                    mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) );

            // get the derivative dof FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init the dTestTractiondDof
            mdTestTractiondDof( tTestDofIndex )( tDofIndex ).set_size( tFITest->get_number_of_space_time_coefficients(), tFIDer->get_number_of_space_time_coefficients() );

            // if viscosity is the test dof
            if ( aTestDofTypes( 0 ) == mDofViscosity && aDofTypes( 0 ) == mDofViscosity )
            {
                // FIXME: Missing second order derivative of effective dynamic viscosity - FD for now
                Constitutive_Model::eval_dtesttractiondu_FD(
                        aDofTypes,
                        aTestDofTypes,
                        mdTestTractiondDof( tTestDofIndex )( tDofIndex ),
                        1e-6,
                        aNormal,
                        aJump,
                        fem::FDScheme_Type::POINT_3_CENTRAL );
            }
            else
            {
                // flatten normal
                Matrix< DDRMat > tFlatNormal;
                this->flatten_normal( aNormal, tFlatNormal );

                // if effective dynamic viscosity depends on test or derivative dof type
                mdTestTractiondDof( tTestDofIndex )( tDofIndex ) =
                        2.0 * trans( tFlatNormal * this->dStraindDOF( aTestDofTypes ) ) * aJump * this->deffdynviscdu( aDofTypes ) + 2.0 * trans( this->deffdynviscdu( aTestDofTypes ) ) * trans( aJump ) * tFlatNormal * this->dStraindDOF( aDofTypes );
            }
        }

        //------------------------------------------------------------------------------

        void
        CM_Fluid_Turbulence::eval_effective_dynamic_viscosity()
        {
            // compute the effective dynamic viscosity
            mEffDynVisc = mPropViscosity->val() + this->turbulent_dynamic_viscosity();
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Fluid_Turbulence::effective_dynamic_viscosity(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Fluid_Turbulence::effective_dynamic_viscosity - Only DEFAULT CM function type known in base class." );

            // if the effective conductivity was not evaluated
            if ( mEffDynViscEval )
            {
                // evaluate the effective conductivity
                this->eval_effective_dynamic_viscosity();

                // set bool for evaluation
                mEffDynViscEval = false;
            }
            // return the effective conductivity value
            return mEffDynVisc;
        }

        //------------------------------------------------------------------------------

        void
        CM_Fluid_Turbulence::eval_deffdynviscdx( uint aOrder )
        {
            // FIXME work only for 1st order
            MORIS_ERROR( aOrder == 1,
                    "CM_Fluid_Turbulence::eval_deffdynviscdx - Works only for 1st order derivative for now." );

            // set matrix size
            mdEffDynViscdx( aOrder - 1 ).set_size( mSpaceDim, 1 );

            // add the derivative of the dynamic viscosity wrt x
            mdEffDynViscdx( aOrder - 1 ) = this->dturbdynviscdx( 1 );

            // if dynamic viscosity depends on space
            if ( mPropViscosity->check_space_dependency( 1 ) )
            {
                // add contribution from dynamic viscosity
                mdEffDynViscdx( aOrder - 1 ) += mPropViscosity->dnPropdxn( 1 );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Fluid_Turbulence::deffdynviscdx(
                uint                  aOrder,
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Fluid_Turbulence::deffdynviscdx - Only DEFAULT CM function type known in base class." );

            MORIS_ERROR( aOrder == 1,
                    "CM_Fluid_Turbulence::deffdynviscdx - Works only for 1st order derivative for now." );

            // if the derivative has not been evaluated yet
            if ( mdEffDynViscdxEval( aOrder - 1 ) )
            {
                // evaluate the derivative
                this->eval_deffdynviscdx( aOrder );

                // set bool for evaluation
                mdEffDynViscdxEval( aOrder - 1 ) = false;
            }

            // return the derivative
            return mdEffDynViscdx( aOrder - 1 );
        }

        //------------------------------------------------------------------------------

        void
        CM_Fluid_Turbulence::eval_deffdynviscdu(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // initialize the matrix for dEffConddu
            mdEffDynViscdu( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // compute the derivative of the effective dynamic viscosity wrt aDofTypes
            mdEffDynViscdu( tDofIndex ) = this->dturbdynviscdu( aDofTypes );

            // if dyn viscosity depends on the dof type
            if ( mPropViscosity->check_dof_dependency( aDofTypes ) )
            {
                mdEffDynViscdu( tDofIndex ) += mPropViscosity->dPropdDOF( aDofTypes );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Fluid_Turbulence::deffdynviscdu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Fluid_Turbulence::deffdynviscdu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Fluid_Turbulence::deffdynviscdu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdEffDynViscduEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_deffdynviscdu( aDofType );

                // set bool for evaluation
                mdEffDynViscduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdEffDynViscdu( tDofIndex );
        }

        //------------------------------------------------------------------------------

        void
        CM_Fluid_Turbulence::eval_deffdynviscdxdu(
                const moris::Cell< MSI::Dof_Type >& aDofTypes,
                uint                                aOrder )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set matrix size
            mdEffDynViscdxdu( aOrder - 1 )( tDofIndex ).set_size( mSpaceDim, tFIDer->get_number_of_space_time_coefficients() );

            // add the derivative of the effective dynamic viscosity wrt aDofTypes and x
            mdEffDynViscdxdu( aOrder - 1 )( tDofIndex ) = this->dturbdynviscdxdu( aDofTypes, 1 );

            // if dynamic viscosity depends on the dof type
            if ( mPropViscosity->check_dof_dependency( aDofTypes ) )
            {
                // if dynamic viscosity depends on space
                if ( mPropViscosity->check_space_dependency( 1 ) )
                {
                    // FIXME dPropdxdu
                    // mdEffDynViscdxdu( aOrder - 1 )( tDofIndex ) += mPropViscosity->dPropdxdDOF( aDofTypes );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Fluid_Turbulence::deffdynviscdxdu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                uint                                aOrder,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Fluid_Turbulence::deffdynviscdxdu - Only DEFAULT CM function type known in base class." );

            MORIS_ERROR( aOrder == 1,
                    "CM_Fluid_Turbulence::deffdynviscdxdu - Works only for 1st order derivative for now." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdEffDynViscdxduEval( aOrder - 1, tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_deffdynviscdxdu( aDofType, aOrder );

                // set bool for evaluation
                mdEffDynViscdxduEval( aOrder - 1, tDofIndex ) = false;
            }

            // return the derivative
            return mdEffDynViscdxdu( aOrder - 1 )( tDofIndex );
        }

        //------------------------------------------------------------------------------

        void
        CM_Fluid_Turbulence::eval_turbulent_dynamic_viscosity()
        {
            // init mTurbDynVisc
            mTurbDynVisc = 0.0;

            // get the viscosity dof type FI
            Field_Interpolator* tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the modified viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if modified viscosity is positive
            if ( tModViscosity >= 0.0 )
            {
                // compute turbulent viscosity
                mTurbDynVisc = mPropDensity->val()( 0 ) * tModViscosity * this->fv1();
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Fluid_Turbulence::turbulent_dynamic_viscosity(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Fluid_Turbulence::turbulent_dynamic_viscosity - Only DEFAULT CM function type known in base class." );

            // if the turbulent dynamic viscosity was not evaluated
            if ( mTurbDynViscEval )
            {
                // evaluate the turbulent dynamic viscosity
                this->eval_turbulent_dynamic_viscosity();

                // set bool for evaluation
                mTurbDynViscEval = false;
            }
            // return the turbulent dynamic viscosity value
            return mTurbDynVisc;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Fluid_Turbulence::eval_dturbdynviscdu(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set matrix size
            mdTurbDynViscdu( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // get the viscosity dof type FI
            Field_Interpolator* tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the modified viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if modified viscosity is positive
            if ( tModViscosity >= 0.0 )
            {
                // add contribution from dfv1du
                mdTurbDynViscdu( tDofIndex ) =
                        mPropDensity->val()( 0 ) * tFIModViscosity->val() * this->dfv1du( aDofTypes );

                // if dof type is viscosity
                if ( aDofTypes( 0 ) == mDofViscosity )
                {
                    // add contribution to dSPdu
                    mdTurbDynViscdu( tDofIndex ) +=
                            mPropDensity->val()( 0 ) * this->fv1() * tFIModViscosity->N();
                }

                // if density depends on dof
                if ( mPropDensity->check_dof_dependency( aDofTypes ) )
                {
                    // add contribution from drhodu
                    mdTurbDynViscdu( tDofIndex ) +=
                            tFIModViscosity->val() * this->fv1() * mPropDensity->dPropdDOF( aDofTypes );
                }
            }
            else
            {
                mdTurbDynViscdu( tDofIndex ).fill( 0.0 );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Fluid_Turbulence::dturbdynviscdu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Fluid_Turbulence::dturbdynviscdu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Fluid_Turbulence::dturbdynviscdu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdTurbDynViscduEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dturbdynviscdu( aDofType );

                // set bool for evaluation
                mdTurbDynViscduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdTurbDynViscdu( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Fluid_Turbulence::eval_dturbdynviscdx( uint aOrder )
        {
            // set matrix size
            mdTurbDynViscdx( aOrder - 1 ).set_size( mSpaceDim, 1 );

            // get the viscosity dof type FI
            Field_Interpolator* tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the modified viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if modified viscosity is positive
            if ( tModViscosity >= 0.0 )
            {
                // compute dTurbDynViscdx
                mdTurbDynViscdx( aOrder - 1 ) =
                        mPropDensity->val()( 0 ) * tFIModViscosity->gradx( 1 ) * this->fv1() + mPropDensity->val()( 0 ) * this->dfv1dx( 1 ) * tModViscosity;

                // if density depends on space
                if ( mPropDensity->check_space_dependency( 1 ) )
                {
                    // add contribution from density space derivative
                    mdTurbDynViscdx( aOrder - 1 ) += this->fv1() * tModViscosity * mPropDensity->dnPropdxn( 1 );
                }
            }
            else
            {
                mdTurbDynViscdx( aOrder - 1 ).fill( 0.0 );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Fluid_Turbulence::dturbdynviscdx(
                uint                  aOrder,
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Fluid_Turbulence::dturbdynviscdx - Only DEFAULT CM function type known in base class." );

            MORIS_ERROR( aOrder == 1,
                    "CM_Fluid_Turbulence::dturbdynviscdx - Works only for 1st order derivative for now." );

            // if the derivative has not been evaluated yet
            if ( mdTurbDynViscdxEval( aOrder - 1 ) )
            {
                // evaluate the derivative
                this->eval_dturbdynviscdx( aOrder );

                // set bool for evaluation
                mdTurbDynViscdxEval( aOrder - 1 ) = false;
            }

            // return the derivative
            return mdTurbDynViscdx( aOrder - 1 );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Fluid_Turbulence::eval_dturbdynviscdxdu(
                const moris::Cell< MSI::Dof_Type >& aDofTypes,
                uint                                aOrder )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set matrix size
            mdTurbDynViscdxdu( aOrder - 1 )( tDofIndex ).set_size( mSpaceDim, tFIDer->get_number_of_space_time_coefficients() );

            // get the viscosity dof type FI
            Field_Interpolator* tFIModViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the modified viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if modified viscosity is positive
            if ( tModViscosity >= 0.0 )
            {
                // add contribution from dfv1du
                mdTurbDynViscdxdu( aOrder - 1 )( tDofIndex ) =
                        mPropDensity->val()( 0 ) * tFIModViscosity->gradx( 1 ) * this->dfv1du( aDofTypes ) + mPropDensity->val()( 0 ) * tFIModViscosity->val()( 0 ) * this->dfv1dxdu( aDofTypes, 1 );

                // if density depends on space
                if ( mPropDensity->check_space_dependency( 1 ) )
                {
                    // add contribution from density space derivative
                    mdTurbDynViscdxdu( aOrder - 1 )( tDofIndex ) +=
                            mPropDensity->dnPropdxn( 1 ) * tFIModViscosity->val()( 0 ) * this->dfv1du( aDofTypes );
                }

                // if dof type is viscosity
                if ( aDofTypes( 0 ) == mDofViscosity )
                {
                    // add contribution to dviscositytdxdu
                    mdTurbDynViscdxdu( aOrder - 1 )( tDofIndex ) +=
                            mPropDensity->val()( 0 ) * this->fv1() * tFIModViscosity->dnNdxn( 1 ) + mPropDensity->val()( 0 ) * this->dfv1dx( 1 ) * tFIModViscosity->N();

                    // if density depends on space
                    if ( mPropDensity->check_space_dependency( 1 ) )
                    {
                        // add contribution from density space derivative
                        mdTurbDynViscdxdu( aOrder - 1 )( tDofIndex ) +=
                                this->fv1() * mPropDensity->dnPropdxn( 1 ) * tFIModViscosity->N();
                    }
                }

                // if density depends on dof
                if ( mPropDensity->check_dof_dependency( aDofTypes ) )
                {
                    mdTurbDynViscdxdu( aOrder - 1 )( tDofIndex ) +=
                            tFIModViscosity->gradx( 1 ) * this->fv1() * mPropDensity->dPropdDOF( aDofTypes ) + this->dfv1dx( 1 ) * tFIModViscosity->val()( 0 ) * mPropDensity->dPropdDOF( aDofTypes );

                    // if density depends on space
                    if ( mPropDensity->check_space_dependency( 1 ) )
                    {
                        // FIXME dPropdxdu
                        // mdTurbDynViscdxdu( aOrder - 1 )( tDofIndex ) +=
                        // tFIModViscosity->val()( 0 ) * tFv1 * mPropDensity->dPropdxdDOF( aDofTypes )
                    }
                }
            }
            else
            {
                mdTurbDynViscdxdu( aOrder - 1 )( tDofIndex ).fill( 0.0 );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Fluid_Turbulence::dturbdynviscdxdu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                uint                                aOrder,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Fluid_Turbulence::dturbdynviscdxdu - Only DEFAULT CM function type known in base class." );

            MORIS_ERROR( aOrder == 1,
                    "CM_Fluid_Turbulence::dturbdynviscdxdu - Works only for 1st order derivative for now." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdTurbDynViscdxduEval( aOrder - 1, tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dturbdynviscdxdu( aDofType, aOrder );

                // set bool for evaluation
                mdTurbDynViscdxduEval( aOrder - 1, tDofIndex ) = false;
            }

            // return the derivative
            return mdTurbDynViscdxdu( aOrder - 1 )( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        CM_Fluid_Turbulence::chi(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Fluid_Turbulence::chi - Only DEFAULT CM function type known in base class." );

            // if the diffusion coefficient was not evaluated
            if ( mChiEval )
            {
                // evaluate chi
                mChi = compute_chi(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity );

                // set bool for evaluation
                mChiEval = false;
            }
            // return the diffusion coefficient
            return mChi;
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Fluid_Turbulence::dchidu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Fluid_Turbulence::dchidu - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Fluid_Turbulence::dchidu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdChiduEval( tDofIndex ) )
            {
                // evaluate the derivative
                compute_dchidu(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity,
                        aDofType,
                        mdChidu( tDofIndex ) );

                // set bool for evaluation
                mdChiduEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdChidu( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Fluid_Turbulence::dchidx(
                uint                  aOrder,
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Fluid_Turbulence::dchidx - Only DEFAULT CM function type known in base class." );

            MORIS_ERROR( aOrder == 1,
                    "CM_Fluid_Turbulence::dchidx - Works only for 1st order derivative for now." );

            // if the derivative has not been evaluated yet
            if ( mdChidxEval( aOrder - 1 ) )
            {
                // evaluate the derivative
                compute_dchidx(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity,
                        mdChidx( aOrder - 1 ) );

                // set bool for evaluation
                mdChidxEval( aOrder - 1 ) = false;
            }

            // return the derivative
            return mdChidx( aOrder - 1 );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Fluid_Turbulence::dchidxdu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                uint                                aOrder,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Fluid_Turbulence::dchidxdu - Only DEFAULT CM function type known in base class." );

            MORIS_ERROR( aOrder == 1,
                    "CM_Fluid_Turbulence::dchidxdu - Works only for 1st order derivative for now." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdChidxduEval( aOrder - 1, tDofIndex ) )
            {
                // evaluate the derivative
                compute_dchidxdu(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity,
                        aDofType,
                        mdChidxdu( aOrder - 1 )( tDofIndex ) );

                // set bool for evaluation
                mdChidxduEval( aOrder - 1, tDofIndex ) = false;
            }

            // return the derivative
            return mdChidxdu( aOrder - 1 )( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        CM_Fluid_Turbulence::fv1(
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Fluid_Turbulence::fv1 - Only DEFAULT CM function type known in base class." );

            // if not evaluated
            if ( mFv1Eval )
            {
                // evaluate
                mFv1 = compute_fv1(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity );

                // set bool for evaluation
                mFv1Eval = false;
            }
            // return
            return mFv1;
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Fluid_Turbulence::dfv1du(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Spalart_Allmaras_Turbulence::dfv1du - Only DEFAULT CM function type known in base class." );

            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "CM_Spalart_Allmaras_Turbulence::dfv1du - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdFv1duEval( tDofIndex ) )
            {
                // evaluate the derivative
                compute_dfv1du(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity,
                        aDofType,
                        mdFv1du( tDofIndex ) );

                // set bool for evaluation
                mdFv1duEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdFv1du( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Fluid_Turbulence::dfv1dx(
                uint                  aOrder,
                enum CM_Function_Type aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Fluid_Turbulence::dfv1dx - Only DEFAULT CM function type known in base class." );

            MORIS_ERROR( aOrder == 1,
                    "CM_Fluid_Turbulence::dfv1dx - Works only for 1st order derivative for now." );

            // if the derivative has not been evaluated yet
            if ( mdFv1dxEval( aOrder - 1 ) )
            {
                // evaluate the derivative
                compute_dfv1dx(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity,
                        mdFv1dx( aOrder - 1 ) );

                // set bool for evaluation
                mdFv1dxEval( aOrder - 1 ) = false;
            }

            // return the derivative
            return mdFv1dx( aOrder - 1 );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        CM_Fluid_Turbulence::dfv1dxdu(
                const moris::Cell< MSI::Dof_Type >& aDofType,
                uint                                aOrder,
                enum CM_Function_Type               aCMFunctionType )
        {
            // check CM function type, base class only supports "DEFAULT"
            MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,
                    "CM_Fluid_Turbulence::dfv1dxdu - Only DEFAULT CM function type known in base class." );

            MORIS_ERROR( aOrder == 1,
                    "CM_Fluid_Turbulence::dfv1dxdu - Works only for 1st order derivative for now." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdFv1dxduEval( aOrder - 1, tDofIndex ) )
            {
                // evaluate the derivative
                compute_dfv1dxdu(
                        { mDofViscosity },
                        mFIManager,
                        mPropKinViscosity,
                        aDofType,
                        mdFv1dxdu( aOrder - 1 )( tDofIndex ) );

                // set bool for evaluation
                mdFv1dxduEval( aOrder - 1, tDofIndex ) = false;
            }

            // return the derivative
            return mdFv1dxdu( aOrder - 1 )( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */

