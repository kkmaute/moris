/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Diffusion_Linear_Isotropic_Turbulence.cpp
 *
 */

#include "cl_FEM_CM_Diffusion_Linear_Isotropic_Turbulence.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"
#include "fn_sum.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    CM_Diffusion_Linear_Isotropic_Turbulence::CM_Diffusion_Linear_Isotropic_Turbulence()
    {
        // set the property pointer cell size
        mProperties.resize( static_cast< uint >( CM_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "Conductivity" ] = static_cast< uint >( CM_Property_Type::CONDUCTIVITY );
        mPropertyMap[ "Density" ]      = static_cast< uint >( CM_Property_Type::DENSITY );
        mPropertyMap[ "HeatCapacity" ] = static_cast< uint >( CM_Property_Type::HEAT_CAPACITY );
        mPropertyMap[ "EigenStrain" ]  = static_cast< uint >( CM_Property_Type::EIGEN_STRAIN );

        mPropertyMap[ "TurbulentPrandtl" ]   = static_cast< uint >( CM_Property_Type::TURBULENT_PRANDTL );

        // init storage for evaluation
        mdTurbDynViscdx.resize( mMaxSpaceDerOrder );
        mdEffConddx.resize( mMaxSpaceDerOrder );

        // init flag for evaluation
        mdTurbDynViscdxEval.set_size( mMaxSpaceDerOrder, 1, true );
        mdEffConddxEval.set_size( mMaxSpaceDerOrder, 1, true );
    }

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence::reset_specific_eval_flags()
    {
        // call parent implementation
        CM_Diffusion_Linear_Isotropic::reset_specific_eval_flags();

        // reset child specific eval flags for turbulence dynamic viscosity
        mTurbDynViscEval = true;
        mdTurbDynViscduEval.fill( true );
        mdTurbDynViscdxEval.fill( true );
        mdTurbDynViscdxduEval.fill( true );
        mTestTurbDynViscEval.fill( true );
        mdTestTurbDynViscduEval.fill( true );

        // reset child specific eval flags for effective conductivity
        mEffCondEval = true;
        mdEffCondduEval.fill( true );
        mdEffConddxEval.fill( true );
        mdEffConddxduEval.fill( true );
        mTestEffCondEval.fill( true );
        mdTestEffCondduEval.fill( true );
    }

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence::initialize_spec_storage_vars_and_eval_flags()
    {
        // call parent implementation
        CM_Diffusion_Linear_Isotropic::initialize_spec_storage_vars_and_eval_flags();

        // get number of dof types
        uint tNumGlobalDofTypes = mGlobalDofTypes.size();
        uint tNumDirectDofTypes = mDofTypes.size();

        // init child specific eval flags
        mdTurbDynViscduEval.set_size( tNumGlobalDofTypes, 1, true );
        mdTurbDynViscdxduEval.set_size( mMaxSpaceDerOrder, tNumGlobalDofTypes, true );
        mdTestTurbDynViscduEval.set_size( tNumDirectDofTypes, tNumGlobalDofTypes, true );
        mTestTurbDynViscEval.set_size( tNumDirectDofTypes, 1, true );
        mdEffCondduEval.set_size( tNumGlobalDofTypes, 1, true );
        mdEffConddxduEval.set_size( mMaxSpaceDerOrder, tNumGlobalDofTypes, true );
        mdTestEffCondduEval.set_size( tNumDirectDofTypes, tNumGlobalDofTypes, true );
        mTestEffCondEval.set_size( tNumDirectDofTypes, 1, true );

        // init child specific storage
        mdEffConddu.resize( tNumGlobalDofTypes );
        mTestEffCond.resize( tNumDirectDofTypes );
        mdTestEffConddu.resize( tNumDirectDofTypes );
        mdTurbDynViscdu.resize( tNumGlobalDofTypes );
        mTestTurbDynVisc.resize( tNumDirectDofTypes );
        mdTestTurbDynViscdu.resize( tNumDirectDofTypes );
        for ( uint iDirectDof = 0; iDirectDof < tNumDirectDofTypes; iDirectDof++ )
        {
            mdTestTurbDynViscdu( iDirectDof ).resize( tNumGlobalDofTypes );
            mdTestEffConddu( iDirectDof ).resize( tNumGlobalDofTypes );
        }

        mdTurbDynViscdxdu.resize( mMaxSpaceDerOrder );
        mdEffConddxdu.resize( mMaxSpaceDerOrder );
        for ( uint iOrder = 0; iOrder < mMaxSpaceDerOrder; iOrder++ )
        {
            mdTurbDynViscdxdu( iOrder ).resize( tNumGlobalDofTypes );
            mdEffConddxdu( iOrder ).resize( tNumGlobalDofTypes );
        }
    }

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence::set_local_properties()
    {
        // set the conductivity property
        mPropConductivity = get_property( "Conductivity" );

        // set the heat capacity property
        mPropHeatCapacity = get_property( "HeatCapacity" );

        // set the density property
        mPropDensity = get_property( "Density" );

        // set the eigenstrain property
        mPropEigenStrain = get_property( "EigenStrain" );

        // set the Prandtl turbulence property
        mPropPrandtlT = get_property( "TurbulentPrandtl" );
    }

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence::eval_flux()
    {
        // compute flux
        // Note: this is a numerical flux and not a physical flux which is the negative
        //       of the numerical flux
        mFlux = this->effective_conductivity()( 0 ) * mFIManager->get_field_interpolators_for_type( mTempDof )->gradx( 1 );

        if ( mPropEigenStrain != nullptr )
        {
            MORIS_LOG_INFO( "All theta dof dependencies unchanged" );

            // Get field interpolator for theta
            Field_Interpolator *tFITheta = mFIManager->get_field_interpolators_for_type( mThetaDof );

            // compute normalized gradient of theta
            const Matrix< DDRMat > &tGradTheta = tFITheta->gradx( 1 );

            // compute norm of spatial gradient of theta
            const real tNorm = norm( tGradTheta );

            // add eigen strain contribution
            if ( tNorm > MORIS_REAL_EPS )
            {
                mFlux += mPropConductivity->val()( 0 ) * tGradTheta / tNorm;
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence::eval_dFluxdDOF(
            const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // get the dof type as a uint
        uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( tDofType );

        // get derivative dof type FI
        Field_Interpolator *tFIDer =
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // initialize the matrix
        mdFluxdDof( tDofIndex ).set_size( mSpaceDim, tFIDer->get_number_of_space_time_coefficients() );

        // get temperature FI
        Field_Interpolator *tFITemp =
                mFIManager->get_field_interpolators_for_type( mTempDof );

        // compute derivative with indirect dependency through effective conductivity
        mdFluxdDof( tDofIndex ) = tFITemp->gradx( 1 ) * this->deffconddu( aDofTypes );

        // if direct dependency on the mTempDof type
        if ( aDofTypes( 0 ) == mTempDof )
        {
            // compute derivative with direct dependency
            mdFluxdDof( tDofIndex ) +=
                    this->effective_conductivity()( 0 ) * tFITemp->dnNdxn( 1 );
        }

        // if direct dependency on the mThetaDof type
        if ( aDofTypes( 0 ) == mThetaDof )
        {
            MORIS_LOG_INFO( "All theta dof dependencies unchanged" );

            // get field interpolator for theta
            Field_Interpolator *tFITheta =
                    mFIManager->get_field_interpolators_for_type( mThetaDof );

            // compute normalized gradient of Theta
            const Matrix< DDRMat > &tBTheta    = tFITheta->dnNdxn( 1 );
            const Matrix< DDRMat > &tGradTheta = tFITheta->gradx( 1 );

            // compute norm of spatial gradient of theta
            real tNorm = norm( tGradTheta );

            // add eigen strain contribution
            if ( tNorm > MORIS_REAL_EPS )
            {
                Matrix< DDRMat > tNormGradTheta = tGradTheta / tNorm;
                Matrix< DDRMat > tNormBTheta    = tBTheta / tNorm;

                // compute derivative with direct dependency
                mdFluxdDof( tDofIndex ) +=
                        mPropConductivity->val()( 0 ) * ( tNormBTheta - tNormGradTheta * trans( tNormGradTheta ) * tNormBTheta );
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence::eval_divflux()
    {
        // compute the divergence of the flux
        mDivFlux = this->effective_conductivity() * this->divstrain()    //
                 + trans( this->deffconddx( 1 ) ) * mFIManager->get_field_interpolators_for_type( mTempDof )->gradx( 1 );
    }

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence::eval_ddivfluxdu(
            const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // get the dof type index
        uint tDofIndex =
                mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the corresponding FI
        Field_Interpolator *tFIDer =
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // set size for ddivflux/du
        mddivfluxdu( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

        // add contributions of the derivatives of the effective conductivity and grad effective conductivity
        mddivfluxdu( tDofIndex ) =
                this->effective_conductivity() * this->ddivstraindu( aDofTypes )                                                        //
                + trans( mFIManager->get_field_interpolators_for_type( mTempDof )->gradx( 1 ) ) * this->deffconddxdu( aDofTypes, 1 )    //
                + this->divstrain() * this->deffconddu( aDofTypes );

        // if temperature dof
        if ( aDofTypes( 0 ) == mTempDof )
        {
            // add contributions of the derivatives of the strain and div strain
            mddivfluxdu( tDofIndex ) +=
                    trans( this->deffconddx( 1 ) ) * mFIManager->get_field_interpolators_for_type( mTempDof )->dnNdxn( 1 );
        }
    }

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence::eval_traction(
            const Matrix< DDRMat > &aNormal )
    {
        // compute traction
        mTraction = trans( this->flux() ) * aNormal;
    }

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence::eval_dTractiondDOF(
            const Vector< MSI::Dof_Type > &aDofTypes,
            const Matrix< DDRMat >        &aNormal )
    {
        // get the dof type as a uint
        uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( tDofType );

        // compute derivative
        mdTractiondDof( tDofIndex ) = trans( aNormal ) * this->dFluxdDOF( aDofTypes );
    }

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence::eval_testTraction(
            const Matrix< DDRMat >        &aNormal,
            const Vector< MSI::Dof_Type > &aTestDofTypes )
    {
        // check for acceptable test dof
        MORIS_ERROR( aTestDofTypes( 0 ) == mTempDof,                                  //
                "CM_Diffusion_Linear_Isotropic_Turbulence::eval_dTestTractiondDOF"    //
                " - Only admissible test dof is temperature" );

        // get test dof type index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // add contribution from derivative of effective conductivity wrt test dof
        mTestTraction( tTestDofIndex ) =                                                 //
                trans( aNormal * this->deffconddu( aTestDofTypes ) ) * this->strain()    //
                + this->effective_conductivity()( 0 ) * trans( this->dStraindDOF( aTestDofTypes ) ) * aNormal;
    }

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence::eval_dTestTractiondDOF(
            const Vector< MSI::Dof_Type > &aDofTypes,
            const Matrix< DDRMat >        &aNormal,
            const Vector< MSI::Dof_Type > &aTestDofTypes )
    {
        // get test dof type index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // get the dof type as a uint
        uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( tDofType );

        // if effective dynamic viscosity depends on test or derivative dof type
        mdTestTractiondDof( tTestDofIndex )( tDofIndex ) =                                                   //
                trans( this->dStraindDOF( aTestDofTypes ) ) * aNormal * this->deffconddu( aDofTypes )        //
                + 2.0 * dot( aNormal, this->strain() ) * this->dtesteffconddu( aDofTypes, aTestDofTypes )    //
                + trans( this->deffconddu( aTestDofTypes ) ) * trans( aNormal ) * this->dStraindDOF( aDofTypes );
    }

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence::eval_dTestTractiondDOF(
            const Vector< MSI::Dof_Type > &aDofTypes,
            const Matrix< DDRMat >        &aNormal,
            const Matrix< DDRMat >        &aJump,
            const Vector< MSI::Dof_Type > &aTestDofTypes )
    {
        // get test dof type index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // get the dof type as a uint
        uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( tDofType );

        // if effective dynamic viscosity depends on test or derivative dof type
        mdTestTractiondDof( tTestDofIndex )( tDofIndex ) =                                                   //
                trans( this->testStrain() ) * aNormal * aJump * this->deffconddu( aDofTypes )                //
                + 2.0 * dot( aNormal, this->strain() ) * this->dtesteffconddu( aDofTypes, aTestDofTypes )    //
                + trans( this->testeffcond( aTestDofTypes ) ) * aJump * trans( aNormal ) * this->dStraindDOF( aDofTypes );
    }

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence::eval_const()
    {
        // build an identity matrix
        Matrix< DDRMat > I;
        eye( mSpaceDim, mSpaceDim, I );

        // compute conductivity matrix
        mConst = this->effective_conductivity() * I;
    }

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence::eval_dConstdDOF(
            const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // get the dof type as a uint
        uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( tDofType );

        // get the derivative dof type FI
        Field_Interpolator *tFIDer =
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // reset the matrix
        mdConstdDof( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

        // compute derivative with indirect dependency through properties
        mdConstdDof( tDofIndex ) = this->deffconddu( aDofTypes );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat > &
    CM_Diffusion_Linear_Isotropic_Turbulence::effective_conductivity(
            enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,                      //
                "CM_Diffusion_Linear_Isotropic_Turbulence::effective_conductivity - "    //
                "Only DEFAULT CM function type known in base class." );

        // if the effective conductivity was not evaluated
        if ( mEffCondEval )
        {
            // evaluate the effective conductivity
            this->eval_effective_conductivity();

            // set bool for evaluation
            mEffCondEval = false;
        }
        // return the effective conductivity value
        return mEffCond;
    }

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence::eval_effective_conductivity()
    {
        // compute the effective conductivity
        mEffCond = mPropConductivity->val()( 0 )    //
                 + mPropHeatCapacity->val()( 0 ) * this->turbulent_dynamic_viscosity()( 0 ) / mPropPrandtlT->val()( 0 );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat > &
    CM_Diffusion_Linear_Isotropic_Turbulence::deffconddu(
            const Vector< MSI::Dof_Type > &aDofType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,          //
                "CM_Diffusion_Linear_Isotropic_Turbulence::deffconddu - "    //
                "Only DEFAULT CM function type known in base class." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR( this->check_dof_dependency( aDofType ),                 //
                "CM_Diffusion_Linear_Isotropic_Turbulence::deffconddu - "    //
                "No dependency on this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdEffCondduEval( tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_deffconddu( aDofType );

            // set bool for evaluation
            mdEffCondduEval( tDofIndex ) = false;
        }

        // return the derivative
        return mdEffConddu( tDofIndex );
    }

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence::eval_deffconddu(
            const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // get the dof type as a uint
        uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( tDofType );

        // get derivative dof type FI
        Field_Interpolator *tFIDer =
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // initialize the matrix for dEffConddu
        mdEffConddu( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

        // add contribution from the derivative of the dynamic viscosity wrt aDofTypes
        mdEffConddu( tDofIndex ) =
                mPropHeatCapacity->val()( 0 ) * this->dturbdynviscdu( aDofTypes ) / mPropPrandtlT->val()( 0 );

        // if conductivity depends on the dof type
        if ( mPropConductivity->check_dof_dependency( aDofTypes ) )
        {
            // assume that conductivity prop does not depend on dof
            MORIS_ERROR( false,                                                      //
                    "CM_Diffusion_Linear_Isotropic_Turbulence::eval_deffconddu -"    //
                    "Dependence of conductivity on dof not accounted for." );
        }

        // if heat capacity depends on the dof type
        if ( mPropHeatCapacity->check_dof_dependency( aDofTypes ) )
        {
            // assume that heat capacity prop does not depend on dof
            MORIS_ERROR( false,                                                      //
                    "CM_Diffusion_Linear_Isotropic_Turbulence::eval_deffconddu -"    //
                    "Dependence of heat capacity on dof not accounted for." );
        }

        // if turbulent Prandtl number depends on the dof type
        if ( mPropPrandtlT->check_dof_dependency( aDofTypes ) )
        {
            // assume that urbulent Prandtl  prop does not depend on dof
            MORIS_ERROR( false,                                                      //
                    "CM_Diffusion_Linear_Isotropic_Turbulence::eval_deffconddu -"    //
                    "Dependence of turbulent Prandtl on dof not accounted for." );
        }
    }

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence::eval_deffconddx( uint aOrder )
    {
        // set matrix size
        mdEffConddx( aOrder - 1 ).set_size( mSpaceDim, 1 );

        // add contribution from the derivative of the turbulent dynamic viscosity wrt x
        mdEffConddx( aOrder - 1 ) =
                mPropHeatCapacity->val()( 0 ) * this->dturbdynviscdx( aOrder ) / mPropPrandtlT->val()( 0 );

        // if conductivity depends on space
        if ( mPropConductivity->check_space_dependency() )
        {
            // assume that conductivity prop does not depend on x
            MORIS_ERROR( false,                                                       //
                    "CM_Diffusion_Linear_Isotropic_Turbulence::eval_deffconddx - "    //
                    "Dependence of conductivity on space not accounted for." );
        }

        // if heat capacity depends on space
        if ( mPropHeatCapacity->check_space_dependency() )
        {
            // assume that heat capacity prop does not depend on x
            MORIS_ERROR( false,                                                       //
                    "CM_Diffusion_Linear_Isotropic_Turbulence::eval_deffconddx - "    //
                    "Dependence of heat capacity on space not accounted for." );
        }

        // if turbulent Prandtl depends on space
        if ( mPropPrandtlT->check_space_dependency() )
        {
            // assume that turbulent Prandtl prop does not depend on x
            MORIS_ERROR( false,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::eval_deffconddx - "    //
                    "Dependence of turbulent Prandtl on space not accounted for." );
        }
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat > &
    CM_Diffusion_Linear_Isotropic_Turbulence::deffconddx(
            uint                  aOrder,
            enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,          //
                "CM_Diffusion_Linear_Isotropic_Turbulence::deffconddx - "    //
                " Only DEFAULT CM function type known in base class." );

        // check only first order is asked
        MORIS_ERROR( aOrder == 1,                                            //
                "CM_Diffusion_Linear_Isotropic_Turbulence::deffconddx - "    //
                " Supported only for 1st order space derivative." );

        // if the derivative has not been evaluated yet
        if ( mdEffConddxEval( aOrder - 1 ) )
        {
            // evaluate the derivative
            this->eval_deffconddx( aOrder );

            // set bool for evaluation
            mdEffConddxEval( aOrder - 1 ) = false;
        }

        // return the derivative
        return mdEffConddx( aOrder - 1 );
    }

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence::eval_deffconddxdu(
            const Vector< MSI::Dof_Type > &aDofTypes,
            uint                           aOrder )
    {
        // get the dof type as a uint
        uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( tDofType );

        // get derivative dof type FI
        Field_Interpolator *tFIDer =
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // set matrix size
        mdEffConddxdu( aOrder - 1 )( tDofIndex ).set_size( mSpaceDim, tFIDer->get_number_of_space_time_coefficients() );

        // compute the derivative of the effective conductivity wrt aDofTypes and x
        mdEffConddxdu( aOrder - 1 )( tDofIndex ) =
                mPropHeatCapacity->val()( 0 ) * this->dturbdynviscdxdu( aDofTypes, 1 ) / mPropPrandtlT->val()( 0 );

        // if conductivity depends on space
        if ( mPropConductivity->check_space_dependency() )
        {
            // assume that conductivity prop does not depend on x
            MORIS_ERROR( false,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::eval_deffconddxdu -"    //
                    "Dependence of conductivity on space not accounted for." );
        }

        // if heat capacity depends on the dof type
        if ( mPropHeatCapacity->check_dof_dependency( aDofTypes ) )
        {
            // assume that heat capacity prop does not depend on dof
            MORIS_ERROR( false,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::eval_deffconddxdu -"    //
                    "Dependence of heat capacity on dof not accounted for." );
        }

        // if turbulent Prandtl number depends on the dof type
        if ( mPropPrandtlT->check_dof_dependency( aDofTypes ) )
        {
            // assume that turbulent Prandtl prop does not depend on dof
            MORIS_ERROR( false,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::eval_deffconddxdu -"    //
                    "Dependence of turbulent Prandtl on dof not accounted for." );
        }

        // if heat capacity depends on space
        if ( mPropHeatCapacity->check_space_dependency() )
        {
            // assume that heat capacity prop does not depend on x
            MORIS_ERROR( false,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::eval_deffconddxdu -"    //
                    "Dependence of heat capacity on space not accounted for." );
        }

        // if turbulent Prandtl number depends on space
        if ( mPropPrandtlT->check_space_dependency() )
        {
            // assume that turbulent Prandtl prop does not depend on x
            MORIS_ERROR( false,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::eval_deffconddxdu -"    //
                    "Dependence of turbulent Prandtl on space not accounted for." );
        }
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat > &
    CM_Diffusion_Linear_Isotropic_Turbulence::deffconddxdu(
            const Vector< MSI::Dof_Type > &aDofType,
            uint                           aOrder,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,            //
                "CM_Diffusion_Linear_Isotropic_Turbulence::deffconddxdu - "    //
                "Only DEFAULT CM function type known in base class." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR( this->check_dof_dependency( aDofType ),                   //
                "CM_Diffusion_Linear_Isotropic_Turbulence::deffconddxdu - "    //
                "No dependency on this dof type." );

        // check only first order is asked
        MORIS_ERROR( aOrder == 1,                                              //
                "CM_Diffusion_Linear_Isotropic_Turbulence::deffconddxdu - "    //
                "Supported only for 1st order space derivative." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdEffConddxduEval( aOrder - 1, tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_deffconddxdu( aDofType, aOrder );

            // set bool for evaluation
            mdEffConddxduEval( aOrder - 1, tDofIndex ) = false;
        }

        // return the derivative
        return mdEffConddxdu( aOrder - 1 )( tDofIndex );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat > &
    CM_Diffusion_Linear_Isotropic_Turbulence::testeffcond(
            const Vector< MSI::Dof_Type > &aTestDofTypes,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,           //
                "CM_Diffusion_Linear_Isotropic_Turbulence::testeffcond - "    //
                "Only DEFAULT CM function type known in base class." );

        // check for valid test dof type
        MORIS_ERROR( aTestDofTypes( 0 ) == mTempDof,                          //
                "CM_Diffusion_Linear_Isotropic_Turbulence::testeffcond - "    //
                "Only temperature is a valid test dof type." );

        // get test dof type index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // if not evaluated
        if ( mTestEffCondEval( tTestDofIndex ) )
        {
            // evaluate
            this->eval_testeffcond( aTestDofTypes );

            // set bool for evaluation
            mTestEffCondEval( tTestDofIndex ) = false;
        }

        // return value
        return mTestEffCond( tTestDofIndex );
    }

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence::eval_testeffcond(
            const Vector< MSI::Dof_Type > &aTestDofTypes )
    {
        // get the test dof index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // test of effective conductivity evaluated from derivative wrt to dof
        mTestEffCond( tTestDofIndex ) = this->deffconddu( aTestDofTypes );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat > &
    CM_Diffusion_Linear_Isotropic_Turbulence::dtesteffconddu(
            const Vector< MSI::Dof_Type > &aDofTypes,
            const Vector< MSI::Dof_Type > &aTestDofTypes,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,        //
                "CM_Fluid_Incompressible_Turbulence::dtesteffconddu - "    //
                "Only DEFAULT CM function type known in base class." );

        // check for valid test dof type
        MORIS_ERROR( aTestDofTypes( 0 ) == mTempDof,                             //
                "CM_Diffusion_Linear_Isotropic_Turbulence::dtesteffconddu - "    //
                "Only temperature is a valid test dof type." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR( this->check_dof_dependency( aDofTypes ),                    //
                "CM_Diffusion_Linear_Isotropic_Turbulence::dtesteffconddu - "    //
                "No dependency on this dof type." );

        // get the test dof index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdTestEffCondduEval( tTestDofIndex, tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_dtesteffconddu( aDofTypes, aTestDofTypes );

            // set bool for evaluation
            mdTestEffCondduEval( tTestDofIndex, tDofIndex ) = false;
        }

        // return the derivative
        return mdTestEffConddu( tTestDofIndex )( tDofIndex );
    }

        //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence::eval_dtesteffconddu(
            const Vector< MSI::Dof_Type > &aDofTypes,
            const Vector< MSI::Dof_Type > &aTestDofTypes )
    {
        // get the test dof index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // add contribution from the derivative of the test turbulent dyn. viscosity wrt dof type
        mdTestEffConddu( tTestDofIndex )( tDofIndex ) =    //
                mPropHeatCapacity->val()( 0 ) * this->dtestturbdynviscdu( aDofTypes, aTestDofTypes ) / mPropPrandtlT->val()( 0 );

        // if conductivity depends on the dof type
        if ( mPropConductivity->check_dof_dependency( aTestDofTypes ) )
        {
            // assume that conductivity prop does not depend on dof
            MORIS_ERROR( false,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::eval_dtesteffconddu -"    //
                    "Dependence of conductivity on dof not accounted for." );
        }

        // if heat capacity depends on the dof type
        if ( mPropHeatCapacity->check_dof_dependency( aTestDofTypes ) )
        {
            // assume that heat capacity prop does not depend on dof
            MORIS_ERROR( false,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::eval_dtesteffconddu -"    //
                    "Dependence of heat capacity on dof not accounted for." );
        }

        // if turbulent Prandtl number depends on the dof type
        if ( mPropPrandtlT->check_dof_dependency( aTestDofTypes ) )
        {
            // assume that urbulent Prandtl  prop does not depend on dof
            MORIS_ERROR( false,
                    "CM_Diffusion_Linear_Isotropic_Turbulence::eval_dtesteffconddu -"    //
                    "Dependence of turbulent Prandtl on dof not accounted for." );
        }
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat > &
    CM_Diffusion_Linear_Isotropic_Turbulence::turbulent_dynamic_viscosity(
            enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,                           //
                "CM_Diffusion_Linear_Isotropic_Turbulence::turbulent_dynamic_viscosity - "    //
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

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence::eval_turbulent_dynamic_viscosity()
    {
        MORIS_ERROR( false,                                                                        //
                "CM_Diffusion_Linear_Isotropic_Turbulence::eval_turbulent_dynamic_viscosity - "    //
                "Not implemented in parent class." );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat > &
    CM_Diffusion_Linear_Isotropic_Turbulence::dturbdynviscdu(
            const Vector< MSI::Dof_Type > &aDofType,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,              //
                "CM_Diffusion_Linear_Isotropic_Turbulence::dturbdynviscdu - "    //
                "Only DEFAULT CM function type known in base class." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR( this->check_dof_dependency( aDofType ),                     //
                "CM_Diffusion_Linear_Isotropic_Turbulence::dturbdynviscdu - "    //
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

    //------------------------------------------------------------------------------

    void
    CM_Diffusion_Linear_Isotropic_Turbulence::eval_dturbdynviscdu(
            const Vector< MSI::Dof_Type > &aDofTypes )
    {
        MORIS_ERROR( false,                                                           //
                "CM_Diffusion_Linear_Isotropic_Turbulence::eval_dturbdynviscdu - "    //
                "Not implemented in parent class." );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat > &
    CM_Diffusion_Linear_Isotropic_Turbulence::dturbdynviscdx(
            uint                  aOrder,
            enum CM_Function_Type aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,              //
                "CM_Diffusion_Linear_Isotropic_Turbulence::dturbdynviscdx - "    //
                "Only DEFAULT CM function type known in base class." );

        // check only first order is asked
        MORIS_ERROR( aOrder == 1,                                                 //
                "CM_Diffusion_Linear_Isotropic_Turbulenceh::dturbdynviscdx - "    //
                "Supported only for 1st order space derivative." );

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
    
    //------------------------------------------------------------------------------

       void
       CM_Diffusion_Linear_Isotropic_Turbulence::eval_dturbdynviscdx(
               uint aOrder )
       {
           MORIS_ERROR( false,                                                           //
                   "CM_Diffusion_Linear_Isotropic_Turbulence::eval_dturbdynviscdx - "    //
                   "Not implemented in parent class." );
       }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat > &
    CM_Diffusion_Linear_Isotropic_Turbulence::dturbdynviscdxdu(
            const Vector< MSI::Dof_Type > &aDofType,
            uint                           aOrder,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,                //
                "CM_Diffusion_Linear_Isotropic_Turbulence::dturbdynviscdxdu - "    //
                "Only DEFAULT CM function type known in base class." );

        // check only first order is asked
        MORIS_ERROR( aOrder == 1,                                                  //
                "CM_Diffusion_Linear_Isotropic_Turbulence::dturbdynviscdxdu - "    //
                "Supported only for 1st order space derivative." );

        // if aDofType is not an active dof type for the CM
        MORIS_ERROR( this->check_dof_dependency( aDofType ),                       //
                "CM_Diffusion_Linear_Isotropic_Turbulence::dturbdynviscdxdu - "    //
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

    void
    CM_Diffusion_Linear_Isotropic_Turbulence::eval_dturbdynviscdxdu(
            const Vector< MSI::Dof_Type > &aDofTypes,
            uint                           aOrder )
    {
        MORIS_ERROR( false,                                                             //
                "CM_Diffusion_Linear_Isotropic_Turbulence::eval_dturbdynviscdxdu - "    //
                "Not implemented in parent class." );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat > &
    CM_Diffusion_Linear_Isotropic_Turbulence::testturbdynvisc(
            const Vector< MSI::Dof_Type > &aTestDofTypes,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,               //
                "CM_Diffusion_Linear_Isotropic_Turbulence::testturbdynvisc - "    //
                "Only DEFAULT CM function type known in base class." );

        // check for valid test dof type
        MORIS_ERROR( aTestDofTypes( 0 ) == mTempDof,                              //
                "CM_Diffusion_Linear_Isotropic_Turbulence::testturbdynvisc - "    //
                "Only temperature is a valid test dof type." );

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
    CM_Diffusion_Linear_Isotropic_Turbulence::eval_testturbdynvisc(
            const Vector< MSI::Dof_Type > &aTestDofTypes )
    {
        MORIS_ERROR( false,                                                            //
                "CM_Diffusion_Linear_Isotropic_Turbulence::eval_testturbdynvisc - "    //
                "Not implemented in parent class." );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat > &
    CM_Diffusion_Linear_Isotropic_Turbulence::dtestturbdynviscdu(
            const Vector< MSI::Dof_Type > &aDofTypes,
            const Vector< MSI::Dof_Type > &aTestDofTypes,
            enum CM_Function_Type          aCMFunctionType )
    {
        // check CM function type, base class only supports "DEFAULT"
        MORIS_ASSERT( aCMFunctionType == CM_Function_Type::DEFAULT,                  //
                "CM_Diffusion_Linear_Isotropic_Turbulence::dtestturbdynviscdu - "    //
                "Only DEFAULT CM function type known in base class." );

        // check for valid test dof type
        MORIS_ERROR( aTestDofTypes( 0 ) == mTempDof,                             //
                "CM_Diffusion_Linear_Isotropic_Turbulence::dtestturbdynviscdu - "    //
                "Only temperature is a valid test dof type." );
        
        // if aDofType is not an active dof type for the CM
        MORIS_ERROR( this->check_dof_dependency( aDofTypes ),                        //
                "CM_Diffusion_Linear_Isotropic_Turbulence::dtestturbdynviscdu - "    //
                "No dependency on this dof type." );

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
    CM_Diffusion_Linear_Isotropic_Turbulence::eval_dtestturbdynviscdu(
            const Vector< MSI::Dof_Type > &aDofTypes,
            const Vector< MSI::Dof_Type > &aTestDofTypes )
    {
        MORIS_ERROR( false,                                                               //
                "CM_Diffusion_Linear_Isotropic_Turbulence::eval_dtestturbdynviscdu - "    //
                "Not implemented in parent class." );
    }

    //------------------------------------------------------------------------------

}    // namespace moris::fem
