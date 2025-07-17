/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Fluid_Incompressible_Turbulence.cpp
 *
 */

#include "cl_FEM_CM_Fluid_Incompressible_Turbulence.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"
// LINALG/src
#include "fn_trans.hpp"
#include "fn_dot.hpp"

namespace moris::fem
{

    //--------------------------------------------------------------------------------------------------------------

    CM_Fluid_Incompressible_Turbulence::CM_Fluid_Incompressible_Turbulence()
    {
        // set the property pointer cell size
        mProperties.resize( static_cast< uint >( CM_Property_Type::MAX_ENUM ), nullptr );

        // populate the map
        mPropertyMap[ "Density" ]      = static_cast< uint >( CM_Property_Type::DENSITY );
        mPropertyMap[ "Viscosity" ]    = static_cast< uint >( CM_Property_Type::VISCOSITY );

        // init storage for evaluation
        mdTurbDynViscdx.resize( mMaxSpaceDerOrder );
        mdEffDynViscdx.resize( mMaxSpaceDerOrder );

        // init flag for evaluation
        mdTurbDynViscdxEval.set_size( mMaxSpaceDerOrder, 1, true );
        mdEffDynViscdxEval.set_size( mMaxSpaceDerOrder, 1, true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence::reset_specific_eval_flags()
    {
        // call parent implementation
        CM_Fluid_Incompressible::reset_specific_eval_flags();

        // reset child specific eval flags for turbulence dynamic viscosity
        mTurbDynViscEval = true;
        mdTurbDynViscduEval.fill( true );
        mdTurbDynViscdxEval.fill( true );
        mdTurbDynViscdxduEval.fill( true );
        mTestTurbDynViscEval.fill( true );
        mdTestTurbDynViscduEval.fill( true );

        // reset child specific eval flags for effective dynamic viscosity
        mEffDynViscEval = true;
        mdEffDynViscduEval.fill( true );
        mdEffDynViscdxEval.fill( true );
        mdEffDynViscdxduEval.fill( true );
        mTestEffDynViscEval.fill( true );
        mdTestEffDynViscduEval.fill( true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence::initialize_spec_storage_vars_and_eval_flags()
    {
        // call parent implementation
        CM_Fluid_Incompressible::initialize_spec_storage_vars_and_eval_flags();

        // get number of dof types
        uint tNumGlobalDofTypes = mGlobalDofTypes.size();
        uint tNumDirectDofTypes = mDofTypes.size();

        // init child specific eval flags
        mdTurbDynViscduEval.set_size( tNumGlobalDofTypes, 1, true );
        mdTurbDynViscdxduEval.set_size( mMaxSpaceDerOrder, tNumGlobalDofTypes, true );
        mdTestTurbDynViscduEval.set_size( tNumDirectDofTypes, tNumGlobalDofTypes, true );
        mTestTurbDynViscEval.set_size( tNumDirectDofTypes, 1, true );

        mdEffDynViscduEval.set_size( tNumGlobalDofTypes, 1, true );
        mdEffDynViscdxduEval.set_size( mMaxSpaceDerOrder, tNumGlobalDofTypes, true );
        mdTestEffDynViscduEval.set_size( tNumDirectDofTypes, tNumGlobalDofTypes, true );
        mTestEffDynViscEval.set_size( tNumDirectDofTypes, 1, true );

        // init child specific storage
        mdTurbDynViscdu.resize( tNumGlobalDofTypes );
        mTestTurbDynVisc.resize( tNumDirectDofTypes );
        mdTestTurbDynViscdu.resize( tNumDirectDofTypes );
        mdEffDynViscdu.resize( tNumGlobalDofTypes );
        mTestEffDynVisc.resize( tNumDirectDofTypes );
        mdTestEffDynViscdu.resize( tNumDirectDofTypes );
        for ( uint iDirectDof = 0; iDirectDof < tNumDirectDofTypes; iDirectDof++ )
        {
            mdTestTurbDynViscdu( iDirectDof ).resize( tNumGlobalDofTypes );
            mdTestEffDynViscdu( iDirectDof ).resize( tNumGlobalDofTypes );
        }

        mdEffDynViscdxdu.resize( mMaxSpaceDerOrder );
        mdTurbDynViscdxdu.resize( mMaxSpaceDerOrder );
        for ( uint iOrder = 0; iOrder < mMaxSpaceDerOrder; iOrder++ )
        {
            mdEffDynViscdxdu( iOrder ).resize( tNumGlobalDofTypes );
            mdTurbDynViscdxdu( iOrder ).resize( tNumGlobalDofTypes );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence::set_local_properties()
    {
        // set the density property
        mPropDensity = get_property( "Density" );

        // set the dynamic viscosity property
        mPropViscosity = get_property( "Viscosity" );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence::eval_flux()
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
    CM_Fluid_Incompressible_Turbulence::eval_dFluxdDOF(
            const Vector< MSI::Dof_Type >& aDofTypes )
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
            // derivative wrt velocity
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

            // derivative wrt pressure
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
    CM_Fluid_Incompressible_Turbulence::eval_divflux()
    {
        // get the pressure FI
        Field_Interpolator* tPressureFI =
                mFIManager->get_field_interpolators_for_type( mDofPressure );

        // flatten deffdynviscdx
        Matrix< DDRMat > tdEffDynViscdxFlat;
        this->flatten_normal( this->deffdynviscdx( 1 ), tdEffDynViscdxFlat );

        // compute flux
        mDivFlux = -1.0 * tPressureFI->gradx( 1 )                                        //
                 + 2.0 * this->effective_dynamic_viscosity()( 0 ) * this->divstrain()    //
                 + 2.0 * tdEffDynViscdxFlat * this->strain();
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence::eval_ddivfluxdu( const Vector< MSI::Dof_Type >& aDofTypes )
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
            mddivfluxdu( tDofIndex ) =                                                                  //
                    2.0 * this->effective_dynamic_viscosity()( 0 ) * this->ddivstraindu( aDofTypes )    //
                    + 2.0 * tdEffDynViscdxFlat * this->dStraindDOF( aDofTypes );
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
        mddivfluxdu( tDofIndex ) += 2.0 * this->divstrain() * this->deffdynviscdu( aDofTypes );

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
                MORIS_ERROR( false, "CM_Fluid_Incompressible_Turbulence::eval_ddivfluxdu - only 2 or 3D" );
        }

        // add the derivative of gradx of the effective dynamic viscosity wrt dof
        mddivfluxdu( tDofIndex ) += 2.0 * tStrainFull * this->deffdynviscdxdu( aDofTypes, 1 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence::eval_traction( const Matrix< DDRMat >& aNormal )
    {
        // flatten the normal
        Matrix< DDRMat > tFlatNormal;
        this->flatten_normal( aNormal, tFlatNormal );

        // compute the traction
        mTraction = tFlatNormal * this->flux();
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence::eval_dTractiondDOF(
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Matrix< DDRMat >&        aNormal )
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
    CM_Fluid_Incompressible_Turbulence::eval_testTraction(
            const Matrix< DDRMat >&        aNormal,
            const Vector< MSI::Dof_Type >& aTestDofTypes )
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
    CM_Fluid_Incompressible_Turbulence::eval_dTestTractiondDOF(
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Matrix< DDRMat >&        aNormal,
            const Matrix< DDRMat >&        aJump,
            const Vector< MSI::Dof_Type >& aTestDofTypes )
    {
        // get test dof type index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the test dof FI
        Field_Interpolator* tFITest =
                mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) );

        // get the derivative dof FI
        Field_Interpolator* tFIDer =
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // init the dTestTractiondDof
        mdTestTractiondDof( tTestDofIndex )( tDofIndex ).set_size(    //
                tFITest->get_number_of_space_time_coefficients(),     //
                tFIDer->get_number_of_space_time_coefficients() );

        // flatten normal
        Matrix< DDRMat > tFlatNormal;
        this->flatten_normal( aNormal, tFlatNormal );

        // if effective dynamic viscosity depends on test or derivative dof type
        mdTestTractiondDof( tTestDofIndex )( tDofIndex ) =                                                                      //
                2.0 * (                                                                                                         //
                        trans( tFlatNormal * this->dStraindDOF( aTestDofTypes ) ) * aJump * this->deffdynviscdu( aDofTypes )    //
                        + dot( aJump, tFlatNormal * this->strain() ) * this->dtesteffdynviscdu( aDofTypes, aTestDofTypes )      //
                        + trans( this->testeffdynvisc( aTestDofTypes ) ) * trans( aJump ) * tFlatNormal * this->dStraindDOF( aDofTypes ) );
    }

    //------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence::eval_effective_dynamic_viscosity()
    {
        // compute the effective dynamic viscosity
        mEffDynVisc = mPropViscosity->val() + this->turbulent_dynamic_viscosity();
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Incompressible_Turbulence::effective_dynamic_viscosity(
            enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,                     //
                "CM_Fluid_Incompressible_Turbulence::effective_dynamic_viscosity - "    //
                "Only DEFAULT CM function type known in base class." );

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
    CM_Fluid_Incompressible_Turbulence::eval_deffdynviscdx( uint aOrder )
    {
        // set matrix size
        mdEffDynViscdx( aOrder - 1 ).set_size( mSpaceDim, 1 );

        // add the derivative of the dynamic viscosity wrt x
        mdEffDynViscdx( aOrder - 1 ) = this->dturbdynviscdx( 1 );

        // if dynamic viscosity depends on space
        if ( mPropViscosity->check_space_dependency() )
        {
            // assume that density prop does not depend on x
            MORIS_ERROR( false,
                    "CM_Fluid_Incompressible_Turbulence::eval_deffdynviscdxdu -"    //
                    "Dependence of dyn. viscosity on space not accounted for." );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Incompressible_Turbulence::deffdynviscdx(
            uint                  aOrder,
            enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,       //
                "CM_Fluid_Incompressible_Turbulence::deffdynviscdx - "    //
                "Only DEFAULT CM function type known in base class." );

        MORIS_ERROR( aOrder == 1,                                         //
                "CM_Fluid_Incompressible_Turbulence::deffdynviscdx - "    //
                "Works only for 1st order derivative for now." );

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
    CM_Fluid_Incompressible_Turbulence::eval_deffdynviscdu(
            const Vector< MSI::Dof_Type >& aDofTypes )
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
    CM_Fluid_Incompressible_Turbulence::deffdynviscdu(
            const Vector< MSI::Dof_Type >& aDofType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,       //
                "CM_Fluid_Incompressible_Turbulence::deffdynviscdu - "    //
                "Only DEFAULT CM function type known in base class." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR( this->check_dof_dependency( aDofType ),              //
                "CM_Fluid_Incompressible_Turbulence::deffdynviscdu - "    //
                "No dependency on this dof type." );

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
    CM_Fluid_Incompressible_Turbulence::eval_deffdynviscdxdu(
            const Vector< MSI::Dof_Type >& aDofTypes,
            uint                           aOrder )
    {
        // get the dof type as a uint
        uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( tDofType );

        // add the derivative of the effective dynamic viscosity wrt aDofTypes and x
        mdEffDynViscdxdu( aOrder - 1 )( tDofIndex ) = this->dturbdynviscdxdu( aDofTypes, 1 );

        // if dynamic viscosity depends on space
        if ( mPropViscosity->check_space_dependency() )
        {
            // assume that density prop does not depend on x
            MORIS_ERROR( false,
                    "CM_Fluid_Incompressible_Turbulence::eval_deffdynviscdxdu -"    //
                    "Dependence of dyn. viscosity on space not accounted for." );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Incompressible_Turbulence::deffdynviscdxdu(
            const Vector< MSI::Dof_Type >& aDofType,
            uint                           aOrder,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,         //
                "CM_Fluid_Incompressible_Turbulence::deffdynviscdxdu - "    //
                "Only DEFAULT CM function type known in base class." );

        // check space derivative order
        MORIS_ERROR( aOrder == 1,                                           //
                "CM_Fluid_Incompressible_Turbulence::deffdynviscdxdu - "    //
                "Works only for 1st order derivative for now." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR( this->check_dof_dependency( aDofType ),              //
                "CM_Fluid_Incompressible_Turbulence::deffdynviscdu - "    //
                "No dependency on this dof type." );

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

    const Matrix< DDRMat >&
    CM_Fluid_Incompressible_Turbulence::testeffdynvisc(
            const Vector< MSI::Dof_Type >& aTestDofTypes,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,        //
                "CM_Fluid_Incompressible_Turbulence::testeffdynvisc - "    //
                "Only DEFAULT CM function type known in base class." );

        // check for valid test dof type
        MORIS_ERROR( aTestDofTypes( 0 ) == mDofVelocity || aTestDofTypes( 0 ) == mDofPressure,    //
                "CM_Fluid_Incompressible_Turbulence::testeffdynvisc - "                           //
                "Only velocity and pressure are valid test dof types." );

        // get test dof type index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // if not evaluated
        if ( mTestEffDynViscEval( tTestDofIndex ) )
        {
            // evaluate
            this->eval_testeffdynvisc( aTestDofTypes );

            // set bool for evaluation
            mTestEffDynViscEval( tTestDofIndex ) = false;
        }

        // return value
        return mTestEffDynVisc( tTestDofIndex );
    }

    //------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence::eval_testeffdynvisc(
            const Vector< MSI::Dof_Type >& aTestDofTypes )
    {
        // get the test dof index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // test of effective conductivity evaluated from derivative wrt to dof
        mTestEffDynVisc( tTestDofIndex ) = this->deffdynviscdu( aTestDofTypes );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Incompressible_Turbulence::dtesteffdynviscdu(
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Vector< MSI::Dof_Type >& aTestDofTypes,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,           //
                "CM_Fluid_Incompressible_Turbulence::dtesteffdynviscdu - "    //
                "Only DEFAULT CM function type known in base class." );

        // check for valid test dof type
        MORIS_ERROR( aTestDofTypes( 0 ) == mDofVelocity || aTestDofTypes( 0 ) == mDofPressure,    //
                "CM_Fluid_Incompressible_Turbulence::dtesteffdynviscdu - "                        //
                "Only velocity and pressure are valid test dof types." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR( this->check_dof_dependency( aDofTypes ),                 //
                "CM_Fluid_Incompressible_Turbulence::dtesteffdynviscdu - "    //
                "No dependency on this dof type." );

        // get the test dof index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdTestEffDynViscduEval( tTestDofIndex, tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_dtesteffdynviscdu( aDofTypes, aTestDofTypes );

            // set bool for evaluation
            mdTestEffDynViscduEval( tTestDofIndex, tDofIndex ) = false;
        }

        // return the derivative
        return mdTestEffDynViscdu( tTestDofIndex )( tDofIndex );
    }

    //------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence::eval_dtesteffdynviscdu(
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Vector< MSI::Dof_Type >& aTestDofTypes )
    {
        // get the test dof index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // add contribution from the derivative of the test turbulent dyn. viscosity wrt dof type
        mdTestEffDynViscdu( tTestDofIndex )( tDofIndex ) = this->dtestturbdynviscdu( aDofTypes, aTestDofTypes );

        // if dyn. viscosity depends on the dof type
        if ( mPropViscosity->check_dof_dependency( aTestDofTypes ) )
        {
            // assume that dyn. viscosity prop does not depend on dof
            MORIS_ERROR( false,
                    "CM_Fluid_Incompressible_Turbulence::eval_dtesteffconddu -"    //
                    "Dependence of dyn. viscosity on dof not accounted for." );
        }
    }

    //------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence::eval_turbulent_dynamic_viscosity()
    {
        MORIS_ERROR( false,                                                                  //
                "CM_Fluid_Incompressible_Turbulence::eval_turbulent_dynamic_viscosity - "    //
                "Not implemented in parent class." );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Incompressible_Turbulence::turbulent_dynamic_viscosity(
            enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,                     //
                "CM_Fluid_Incompressible_Turbulence::turbulent_dynamic_viscosity - "    //
                "Only DEFAULT CM function type known in base class." );

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
    CM_Fluid_Incompressible_Turbulence::eval_dturbdynviscdu(
            const Vector< MSI::Dof_Type >& aDofTypes )
    {
        MORIS_ERROR( false,                                                     //
                "CM_Fluid_Incompressible_Turbulence::eval_dturbdynviscdu - "    //
                "Not implemented in parent class." );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Incompressible_Turbulence::dturbdynviscdu(
            const Vector< MSI::Dof_Type >& aDofType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,        //
                "CM_Fluid_Incompressible_Turbulence::dturbdynviscdu - "    //
                "Only DEFAULT CM function type known in base class." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR( this->check_dof_dependency( aDofType ),               //
                "CM_Fluid_Incompressible_Turbulence::dturbdynviscdu - "    //
                "No dependency on this dof type." );

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
    CM_Fluid_Incompressible_Turbulence::eval_dturbdynviscdx(
            uint aOrder )
    {
        MORIS_ERROR( false,                                                     //
                "CM_Fluid_Incompressible_Turbulence::eval_dturbdynviscdx - "    //
                "Not implemented in parent class." );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Incompressible_Turbulence::dturbdynviscdx(
            uint                  aOrder,
            enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,        //
                "CM_Fluid_Incompressible_Turbulence::dturbdynviscdx - "    //
                "Only DEFAULT CM function type known in base class." );

        MORIS_ERROR( aOrder == 1,                                          //
                "CM_Fluid_Incompressible_Turbulence::dturbdynviscdx - "    //
                "Works only for 1st order space derivative." );

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
    CM_Fluid_Incompressible_Turbulence::eval_dturbdynviscdxdu(
            const Vector< MSI::Dof_Type >& aDofTypes,
            uint                           aOrder )
    {
        MORIS_ERROR( false,                                                     //
                "CM_Fluid_Incompressible_Turbulence::eval_dturbdynviscdx - "    //
                "Not implemented in parent class." );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Incompressible_Turbulence::dturbdynviscdxdu(
            const Vector< MSI::Dof_Type >& aDofType,
            uint                           aOrder,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,          //
                "CM_Fluid_Incompressible_Turbulence::dturbdynviscdxdu - "    //
                "Only DEFAULT CM function type known in base class." );

        MORIS_ERROR( aOrder == 1,                                            //
                "CM_Fluid_Incompressible_Turbulence::dturbdynviscdxdu - "    //
                "Works only for 1st order space derivative." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR( this->check_dof_dependency( aDofType ),                 //
                "CM_Fluid_Incompressible_Turbulence::dturbdynviscdxdu - "    //
                "No dependency on this dof type." );

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

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Incompressible_Turbulence::testturbdynvisc(
            const Vector< MSI::Dof_Type >& aTestDofTypes,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,         //
                "CM_Fluid_Incompressible_Turbulence::testturbdynvisc - "    //
                "Only DEFAULT CM function type known in base class." );

        // check for valid test dof type
        MORIS_ERROR( aTestDofTypes( 0 ) == mDofVelocity || aTestDofTypes( 0 ) == mDofPressure,    //
                "CM_Fluid_Incompressible_Turbulence::testturbdynvisc - "                          //
                "Only velocity and pressure are valid test dof types." );

        // get test dof type index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // if not evaluated
        if ( mTestTurbDynViscEval( tTestDofIndex ) )
        {
            // evaluate
            this->eval_testturbdynvisc( aTestDofTypes );

            // set bool for evaluation
            mTestTurbDynViscEval( tTestDofIndex ) = false;
        }

        // return value
        return mTestTurbDynVisc( tTestDofIndex );
    }

    //------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence::eval_testturbdynvisc(
            const Vector< MSI::Dof_Type >& aTestDofTypes )
    {
        MORIS_ERROR( false,                                                      //
                "CM_Fluid_Incompressible_Turbulence::eval_testturbdynvisc - "    //
                "Not implemented in parent class." );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Incompressible_Turbulence::dtestturbdynviscdu(
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Vector< MSI::Dof_Type >& aTestDofTypes,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,            //
                "CM_Fluid_Incompressible_Turbulence::dtestturbdynviscdu - "    //
                "Only DEFAULT CM function type known in base class." );

        // check for valid test dof type
        MORIS_ERROR( aTestDofTypes( 0 ) == mDofVelocity || aTestDofTypes( 0 ) == mDofPressure,    //
                "CM_Fluid_Incompressible_Turbulence::testturbdynvisc - "                          //
                "Only velocity and pressure are valid test dof types." );

        // if aDofType is not an active dof type for the property
        MORIS_ERROR( this->check_dof_dependency( aDofTypes ),
                "CM_Fluid_Incompressible_Turbulence::dtestturbdynviscdu - "    //
                "No dependency on dof type." );

        // get the test dof index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdTestTurbDynViscduEval( tTestDofIndex, tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_dtestturbdynviscdu( aDofTypes, aTestDofTypes );

            // set bool for evaluation
            mdTestTurbDynViscduEval( tTestDofIndex, tDofIndex ) = false;
        }

        // return the derivative
        return mdTestTurbDynViscdu( tTestDofIndex )( tDofIndex );
    }
    //------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence::eval_dtestturbdynviscdu(
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Vector< MSI::Dof_Type >& aTestDofTypes )
    {
        MORIS_ERROR( false,                                                         //
                "CM_Fluid_Incompressible_Turbulence::eval_dtestturbdynviscdu - "    //
                "Not implemented in parent class." );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Incompressible_Turbulence::select_derivative_FD(
            enum CM_Request_Type           aCMRequestType,
            const Vector< MSI::Dof_Type >& aTestDofTypes,
            const Matrix< DDRMat >&        aNormal,
            const Matrix< DDRMat >&        aJump,
            enum CM_Function_Type          aCMFunctionType )
    {
        switch ( aCMRequestType )
        {
            case CM_Request_Type::STRAIN:
            {
                return this->strain( aCMFunctionType );
                break;
            }
            case CM_Request_Type::STRAIN_SPACE_DER:
            {
                mStrain = trans( this->dstraindx( 1, aCMFunctionType ) ) * aJump;
                return mStrain;
                break;
            }
            case CM_Request_Type::FLUX:
            {
                return this->flux( aCMFunctionType );
                break;
            }
            case CM_Request_Type::TRACTION:
            {
                return this->traction( aNormal, aCMFunctionType );
                break;
            }
            case CM_Request_Type::TEST_TRACTION:
            {
                mTraction = this->testTraction_trans( aNormal, aTestDofTypes, aCMFunctionType ) * aJump;
                return mTraction;
                break;
            }
            case CM_Request_Type::DIV_FLUX:
            {
                return this->divflux();
                break;
            }
            case CM_Request_Type::DIV_STRAIN:
            {
                return this->divstrain();
                break;
            }
            case CM_Request_Type::TURB_DYN_VISC:
            {
                return this->turbulent_dynamic_viscosity();
                break;
            }
            case CM_Request_Type::TEST_TURB_DYN_VISC:
            {
                //  get the test dof index
                uint tTestDofIndex                = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );
                mTestTurbDynVisc( tTestDofIndex ) = trans( this->testturbdynvisc( aTestDofTypes ) );
                return mTestTurbDynVisc( tTestDofIndex );
                break;
            }
            case CM_Request_Type::TURB_DYN_VISC_SPACE_DER:
            {
                return this->dturbdynviscdx( 1 );
                break;
            }
            case CM_Request_Type::EFF_DYN_VISC:
            {
                mEffDynVisc = this->effective_dynamic_viscosity();
                return mEffDynVisc;
                break;
            }
            default:
                MORIS_ERROR( false,                                                      //
                        "CM_Fluid_Incompressible_Turbulence::select_derivative_FD - "    //
                        "aCMRequestType undefined" );
                return this->strain( aCMFunctionType );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Incompressible_Turbulence::set_derivative_FD(
            enum CM_Request_Type           aCMRequestType,
            Matrix< DDRMat >&              aDerivativeFD,
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Vector< MSI::Dof_Type >& aTestDofTypes,
            enum CM_Function_Type          aCMFunctionType )
    {
        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        switch ( aCMRequestType )
        {
            case CM_Request_Type::STRAIN:
            {
                // set value to storage
                mdStraindDof( tDofIndex ) = aDerivativeFD;
                break;
            }
            case CM_Request_Type::STRAIN_SPACE_DER:
            {
                // set value to storage
                mdStraindxdu( 0 )( tDofIndex ) = aDerivativeFD;
                break;
            }
            case CM_Request_Type::FLUX:
            {
                // set value to storage
                mdFluxdDof( tDofIndex ) = aDerivativeFD;
                break;
            }
            case CM_Request_Type::TRACTION:
            {
                // set value to storage
                mdTractiondDof( tDofIndex ) = aDerivativeFD;
                break;
            }
            case CM_Request_Type::TEST_TRACTION:
            {
                // get the test dof index
                uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

                // set value to storage
                mdTestTractiondDof( tTestDofIndex )( tDofIndex ) = aDerivativeFD;
                break;
            }
            case CM_Request_Type::DIV_FLUX:
            {
                // set value to storage
                mddivfluxdu( tDofIndex ) = aDerivativeFD;
                break;
            }
            case CM_Request_Type::DIV_STRAIN:
            {
                // set value to storage
                mddivstraindu( tDofIndex ) = aDerivativeFD;
                break;
            }
            case CM_Request_Type::TURB_DYN_VISC:
            {
                // set value to storage
                mdTurbDynViscdu( tDofIndex ) = aDerivativeFD;
                break;
            }
            case CM_Request_Type::TEST_TURB_DYN_VISC:
            {
                // get the test dof index
                uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

                // set value to storage
                mdTestTurbDynViscdu( tTestDofIndex )( tDofIndex ) = aDerivativeFD;
                break;
            }
            case CM_Request_Type::TURB_DYN_VISC_SPACE_DER:
            {
                // set value to storage
                mdTurbDynViscdxdu( 0 )( tDofIndex ) = aDerivativeFD;
                break;
            }
            case CM_Request_Type::EFF_DYN_VISC:
            {
                // set value to storage
                mdEffDynViscdu( tDofIndex ) = aDerivativeFD;
                break;
            }
            default:
                MORIS_ERROR( false,                                                   //
                        "CM_Fluid_Incompressible_Turbulence::set_derivative_FD - "    //
                        "aCMRequestType undefined" );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::fem
