/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Fluid_Compressible_Ideal.cpp
 *
 */

#include "cl_FEM_CM_Fluid_Compressible_Ideal.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "op_minus.hpp"

namespace moris::fem
{

    //--------------------------------------------------------------------------------------------------------------

    CM_Fluid_Compressible_Ideal::CM_Fluid_Compressible_Ideal()
    {
        // set the property pointer cell size
        mProperties.resize( static_cast< uint >( CM_Property_Type::MAX_ENUM ), nullptr );

        // populate the map
        mPropertyMap[ "IsochoricHeatCapacity" ] = static_cast< uint >( CM_Property_Type::ISOCHORIC_HEAT_CAPACITY );    // constant property
        mPropertyMap[ "SpecificGasConstant" ]   = static_cast< uint >( CM_Property_Type::SPECIFIC_GAS_CONSTANT );      // constant property
        mPropertyMap[ "DynamicViscosity" ]      = static_cast< uint >( CM_Property_Type::DYNAMIC_VISCOSITY );          // may be a fnct. of T
        mPropertyMap[ "ThermalConductivity" ]   = static_cast< uint >( CM_Property_Type::THERMAL_CONDUCTIVITY );       // may be a fnct. of T
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::set_function_pointers()
    {
        switch ( mSpaceDim )
        {
            case ( 2 ):
            {
                m_eval_strain         = &CM_Fluid_Compressible_Ideal::eval_strain_2d;
                m_eval_teststrain     = &CM_Fluid_Compressible_Ideal::eval_teststrain_2d;
                m_eval_velocitymatrix = &CM_Fluid_Compressible_Ideal::eval_velocitymatrix_2d;
                m_unfold_tensor       = &CM_Fluid_Compressible_Ideal::unfold_2d;
                m_flatten_normal      = &CM_Fluid_Compressible_Ideal::flatten_normal_2d;
                mFlatIdentity         = { { 1.0 }, { 1.0 }, { 0.0 } };
                break;
            }
            case ( 3 ):
            {
                m_eval_strain         = &CM_Fluid_Compressible_Ideal::eval_strain_3d;
                m_eval_teststrain     = &CM_Fluid_Compressible_Ideal::eval_teststrain_3d;
                m_eval_velocitymatrix = &CM_Fluid_Compressible_Ideal::eval_velocitymatrix_3d;
                m_unfold_tensor       = &CM_Fluid_Compressible_Ideal::unfold_3d;
                m_flatten_normal      = &CM_Fluid_Compressible_Ideal::flatten_normal_3d;
                mFlatIdentity         = { { 1.0 }, { 1.0 }, { 1.0 }, { 0.0 }, { 0.0 }, { 0.0 } };
                break;
            }
            default:
            {
                MORIS_ERROR( false, "CM_Fluid_Compressible_Ideal::set_function_pointers - this function is currently unused, might be used in the future." );
                break;
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::reset_specific_eval_flags()
    {
        // reset Pressure
        mPressureEval = true;
        mPressureDofEval.fill( true );

        // reset dof derivative of velocity
        mdNveldtEval = true;

        // reset velocity matrix for flattened tensors
        mVelocityMatrixEval = true;

        // reset Thermal Flux ---------------------------------
        mThermalFluxEval = true;
        mThermalFluxDofEval.fill( true );

        // reset work Flux
        mWorkFluxEval = true;
        mWorkFluxDofEval.fill( true );

        // reset energy Flux
        mEnergyFluxEval = true;
        mEnergyFluxDofEval.fill( true );

        // reset mechanical Flux
        // mStressEval = true;
        // mStressDofEval.assign( tNumDofTypes, true );

        // reset Thermal Traction ------------------------------
        mThermalTractionEval = true;
        mThermalTractionDofEval.fill( true );

        // reset work Traction
        mWorkTractionEval = true;
        mWorkTractionDofEval.fill( true );

        // reset energy Traction
        mEnergyTractionEval = true;
        mEnergyTractionDofEval.fill( true );

        // reset Mechanical Traction
        mMechanicalTractionEval = true;
        mMechanicalTractionDofEval.fill( true );

        // reset test tractions --------------------------------

        mThermalTestTractionEval.fill( true );
        mdThermalTestTractiondDofEval.fill( true );

        mMechanicalTestTractionEval.fill( true );
        mdMechanicalTestTractiondDofEval.fill( true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::initialize_spec_storage_vars_and_eval_flags()
    {
        // get number of DoF types
        uint tNumGlobalDofTypes = mGlobalDofTypes.size();
        uint tNumDirectDofTypes = mDofTypes.size();

        // initialize eval flags
        mPressureDofEval.set_size( tNumGlobalDofTypes, 1, true );
        mThermalFluxDofEval.set_size( tNumGlobalDofTypes, 1, true );
        mWorkFluxDofEval.set_size( tNumGlobalDofTypes, 1, true );
        mEnergyFluxDofEval.set_size( tNumGlobalDofTypes, 1, true );
        // mStressDofEval.resize( tNumGlobalDofTypes, true );
        mThermalTractionDofEval.set_size( tNumGlobalDofTypes, 1, true );
        mWorkTractionDofEval.set_size( tNumGlobalDofTypes, 1, true );
        mEnergyTractionDofEval.set_size( tNumGlobalDofTypes, 1, true );
        mMechanicalTractionDofEval.set_size( tNumGlobalDofTypes, 1, true );

        mThermalTestTractionEval.set_size( tNumGlobalDofTypes, 1, true );
        mMechanicalTestTractionEval.set_size( tNumGlobalDofTypes, 1, true );
        mdThermalTestTractiondDofEval.set_size( tNumGlobalDofTypes, 1, true );
        mdMechanicalTestTractiondDofEval.set_size( tNumGlobalDofTypes, 1, true );
        mdThermalTestTractiondDofEval.set_size( tNumDirectDofTypes, tNumGlobalDofTypes, true );
        mdMechanicalTestTractiondDofEval.set_size( tNumDirectDofTypes, tNumGlobalDofTypes, true );

        // initialize storage variables
        mPressureDof.resize( tNumGlobalDofTypes );
        mThermalFluxDof.resize( tNumGlobalDofTypes );
        mWorkFluxDof.resize( tNumGlobalDofTypes );
        mEnergyFluxDof.resize( tNumGlobalDofTypes );
        // mStressDof.resize( tNumGlobalDofTypes );
        mThermalTractionDof.resize( tNumGlobalDofTypes );
        mWorkTractionDof.resize( tNumGlobalDofTypes );
        mEnergyTractionDof.resize( tNumGlobalDofTypes );
        mMechanicalTractionDof.resize( tNumGlobalDofTypes );

        mThermalTestTraction.resize( tNumGlobalDofTypes );
        mMechanicalTestTraction.resize( tNumGlobalDofTypes );
        mdThermalTestTractiondDof.resize( tNumDirectDofTypes );
        mdMechanicalTestTractiondDof.resize( tNumDirectDofTypes );
        for ( uint iDirectDof = 0; iDirectDof < tNumDirectDofTypes; iDirectDof++ )
        {
            mdThermalTestTractiondDof( iDirectDof ).resize( tNumGlobalDofTypes );
            mdMechanicalTestTractiondDof( iDirectDof ).resize( tNumGlobalDofTypes );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::set_dof_type_list(
            const Vector< Vector< MSI::Dof_Type > >& aDofTypes,
            const Vector< std::string >&             aDofStrings )
    {
        // set dof type list
        Constitutive_Model::set_dof_type_list( aDofTypes );

        // loop over the provided dof type
        for ( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
        {
            // get dof type string
            const std::string& tDofString = aDofStrings( iDof );

            // get dof type
            MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

            // switch on dof type string
            if ( tDofString == "Velocity" )
            {
                mDofVelocity = tDofType;
            }
            else if ( tDofString == "Density" )
            {
                mDofDensity = tDofType;
            }
            else if ( tDofString == "Temperature" )
            {
                mDofTemperature = tDofType;
            }
            else
            {
                MORIS_ERROR( false,
                        "CM_Fluid_Compressible_Ideal::set_dof_type_list - Unknown aDofString : %s",
                        tDofString.c_str() );
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::set_local_properties()
    {
        // get the isochoric heat capacity properties
        mPropIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" );

        // get the specific gas constant properties
        mPropSpecificGasConstant = get_property( "SpecificGasConstant" );

        // get the dynamic viscosity properties
        mPropDynamicViscosity = get_property( "DynamicViscosity" );

        // get the thermal conductivity properties
        mPropThermalConductivity = get_property( "ThermalConductivity" );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Ideal::flux( enum CM_Function_Type aCMFunctionType )
    {
        switch ( aCMFunctionType )
        {
            case CM_Function_Type::THERMAL:
                return this->thermal_flux();

            case CM_Function_Type::ENERGY:
                return this->energy_flux();

            case CM_Function_Type::WORK:
                return this->work_flux();

            case CM_Function_Type::MECHANICAL:
                return this->stress();

            case CM_Function_Type::DEFAULT:
                MORIS_ERROR( false, "CM_Fluid_Compressible_Ideal::flux - DEFAULT function type not supported./" );
                return mFlux;

                // unknown CM function type
            default:
                MORIS_ERROR( false, "CM_Fluid_Compressible_Ideal::flux - unknown CM function type for flux." );
                return mFlux;
        }
    }

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Ideal::dFluxdDOF(
            const Vector< MSI::Dof_Type >& aDofTypes,
            enum CM_Function_Type          aCMFunctionType )
    {
        switch ( aCMFunctionType )
        {
            case CM_Function_Type::THERMAL:
                return this->thermal_dFluxdDOF( aDofTypes );

            case CM_Function_Type::ENERGY:
                return this->energy_dFluxdDOF( aDofTypes );

            case CM_Function_Type::WORK:
                return this->work_dFluxdDOF( aDofTypes );

            case CM_Function_Type::MECHANICAL:
                return this->dStressdDOF( aDofTypes );

                // unknown CM function type
            default:
                MORIS_ERROR( false, "CM_Fluid_Compressible_Van_der_Waals::dFluxdDOF - unknown CM function type for flux." );
                return mFlux;
        }
    }

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::eval_thermal_flux()
    {
        // get the temperature FI
        Field_Interpolator* tTempFI = mFIManager->get_field_interpolators_for_type( mDofTemperature );

        // get the thermal conductivity property
        const std::shared_ptr< Property > tThermalConductivity = get_property( "ThermalConductivity" );

        // compute thermal flux q = - k * grad(T)
        mThermalFlux = -1.0 * tThermalConductivity->val()( 0 ) * tTempFI->gradx( 1 );
    }

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Ideal::thermal_flux()
    {
        // if the test strain was not evaluated
        if ( mThermalFluxEval )
        {
            // evaluate the test strain
            this->eval_thermal_flux();

            // set bool for evaluation
            mThermalFluxEval = false;
        }
        // return the test strain value
        return mThermalFlux;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::eval_thermal_dFluxdDOF( const Vector< MSI::Dof_Type >& aDofTypes )
    {
        // get the conductivity property
        const std::shared_ptr< Property > tPropThermalConductivity = get_property( "ThermalConductivity" );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // initialize the matrix
        mThermalFluxDof( tDofIndex ).set_size(    //
                mSpaceDim,                        //
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(),
                0.0 );

        // get temperature FI
        Field_Interpolator* tFITemp =
                mFIManager->get_field_interpolators_for_type( mDofTemperature );

        // if direct dependency on the dof type
        if ( aDofTypes( 0 ) == mDofTemperature )
        {
            // compute derivative with direct dependency
            mThermalFluxDof( tDofIndex ) +=
                    -1.0 * tPropThermalConductivity->val()( 0 ) * tFITemp->dnNdxn( 1 );
        }

        // if indirect dependency of conductivity on the dof type
        if ( tPropThermalConductivity->check_dof_dependency( aDofTypes ) )
        {
            // compute derivative with indirect dependency through properties
            mThermalFluxDof( tDofIndex ) +=
                    -1.0 * tFITemp->gradx( 1 ) * tPropThermalConductivity->dPropdDOF( aDofTypes );
        }
    }

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Ideal::thermal_dFluxdDOF(
            const Vector< MSI::Dof_Type >& aDofTypes )
    {
        // if aDofType is not an active dof type for the CM
        MORIS_ERROR(
                this->check_dof_dependency( aDofTypes ),
                "CM_Fluid_Compressible_Ideal::thermal_dFluxdDOF - no dependency in this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mThermalFluxDofEval( tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_thermal_dFluxdDOF( aDofTypes );

            // set bool for evaluation
            mThermalFluxDofEval( tDofIndex ) = false;
        }

        // return the dof deriv value
        return mThermalFluxDof( tDofIndex );
    }

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::eval_work_flux()
    {
        // compute contribution
        mWorkFlux = this->velocityMatrix() * this->flux( CM_Function_Type::MECHANICAL );
    }

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Ideal::work_flux()
    {
        // if the flux was not evaluated
        if ( mWorkFluxEval )
        {
            // evaluate the flux
            this->eval_work_flux();

            // set bool for evaluation
            mWorkFluxEval = false;
        }
        // return the flux value
        return mWorkFlux;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::eval_work_dFluxdDOF( const Vector< MSI::Dof_Type >& aDofTypes )
    {
        // get the dof type as a uint
        uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( tDofType );

        // get the velocity
        Field_Interpolator* tFIVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity );

        // unfold the flattened stress tensor
        Matrix< DDRMat > tStressTensor;
        this->unfold( this->flux( CM_Function_Type::MECHANICAL ), tStressTensor );

        // initialize the matrix
        mWorkFluxDof( tDofIndex ).set_size(    //
                mSpaceDim,                     //
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(),
                0.0 );

        // direct dependency on the density dof type
        if ( aDofTypes( 0 ) == mDofDensity )
        {
            // compute contribution
            mWorkFluxDof( tDofIndex ) +=
                    this->velocityMatrix() * this->dFluxdDOF( aDofTypes, CM_Function_Type::MECHANICAL );
        }

        // direct dependency on the velocity dof type
        if ( aDofTypes( 0 ) == mDofVelocity )
        {
            // compute contribution
            mWorkFluxDof( tDofIndex ) +=
                    this->velocityMatrix() * this->dFluxdDOF( aDofTypes, CM_Function_Type::MECHANICAL ) +    //
                    tStressTensor * tFIVelocity->N();
        }

        // direct dependency on the temperature dof type
        if ( aDofTypes( 0 ) == mDofTemperature )
        {
            // compute contribution
            mWorkFluxDof( tDofIndex ) +=
                    this->velocityMatrix() * this->dFluxdDOF( aDofTypes, CM_Function_Type::MECHANICAL );
        }
    }

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Ideal::work_dFluxdDOF(
            const Vector< MSI::Dof_Type >& aDofTypes )
    {
        // if aDofType is not an active dof type for the CM
        MORIS_ERROR(
                this->check_dof_dependency( aDofTypes ),
                "CM_Fluid_Compressible_Ideal::work_dFluxdDOF - no dependency in this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mWorkFluxDofEval( tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_work_dFluxdDOF( aDofTypes );

            // set bool for evaluation
            mWorkFluxDofEval( tDofIndex ) = false;
        }

        // return the dof deriv value
        return mWorkFluxDof( tDofIndex );
    }

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::eval_energy_flux()
    {
        // get the velocity
        Matrix< DDRMat > tVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity )->val();

        // compute contribution
        mEnergyFlux = tVelocity * this->Energy();
    }

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Ideal::energy_flux()
    {
        // if the flux was not evaluated
        if ( mEnergyFluxEval )
        {
            // evaluate the flux
            this->eval_energy_flux();

            // set bool for evaluation
            mEnergyFluxEval = false;
        }
        // return the flux value
        return mEnergyFlux;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::eval_energy_dFluxdDOF( const Vector< MSI::Dof_Type >& aDofTypes )
    {
        // get the dof type as a uint
        uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( tDofType );

        // get the velocity
        Field_Interpolator* tFIVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity );

        // initialize the matrix
        mEnergyFluxDof( tDofIndex ).set_size(    //
                mSpaceDim,                       //
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(),
                0.0 );

        // direct dependency on the density dof type
        if ( aDofTypes( 0 ) == mDofDensity )
        {
            // compute contribution
            mEnergyFluxDof( tDofIndex ) +=
                    tFIVelocity->val() * this->dEnergydDOF( aDofTypes );
        }

        // direct dependency on the velocity dof type
        if ( aDofTypes( 0 ) == mDofVelocity )
        {
            // compute contribution
            mEnergyFluxDof( tDofIndex ) +=
                    tFIVelocity->val() * this->dEnergydDOF( aDofTypes ) +    //
                    this->Energy()( 0 ) * tFIVelocity->N();
        }

        // direct dependency on the temperature dof type
        if ( aDofTypes( 0 ) == mDofTemperature )
        {
            // compute contribution
            mEnergyFluxDof( tDofIndex ) +=
                    tFIVelocity->val() * this->dEnergydDOF( aDofTypes );
        }
    }

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Ideal::energy_dFluxdDOF(
            const Vector< MSI::Dof_Type >& aDofTypes )
    {
        // if aDofType is not an active dof type for the CM
        MORIS_ERROR(
                this->check_dof_dependency( aDofTypes ),
                "CM_Fluid_Compressible_Ideal::energy_dFluxdDOF - no dependency in this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mEnergyFluxDofEval( tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_energy_dFluxdDOF( aDofTypes );

            // set bool for evaluation
            mEnergyFluxDofEval( tDofIndex ) = false;
        }

        // return the dof deriv value
        return mEnergyFluxDof( tDofIndex );
    }

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::eval_Energy()
    {
        // get field interpolator values
        Matrix< DDRMat > tDensity     = mFIManager->get_field_interpolators_for_type( mDofDensity )->val();
        Matrix< DDRMat > tVelocity    = mFIManager->get_field_interpolators_for_type( mDofVelocity )->val();
        Matrix< DDRMat > tTemperature = mFIManager->get_field_interpolators_for_type( mDofTemperature )->val();

        // get the heat capacity
        real tIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" )->val()( 0 );

        // compute thermal flux q = - k * grad(T)
        mEnergy = tIsochoricHeatCapacity * tDensity * tTemperature +    //
                  0.5 * trans( tVelocity ) * tVelocity * tDensity;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::eval_dEnergydDOF( const Vector< MSI::Dof_Type >& aDofTypes )
    {
        // get the dof type as a uint
        uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( tDofType );

        // get the FIs
        Field_Interpolator* tFIDensity  = mFIManager->get_field_interpolators_for_type( mDofDensity );
        Field_Interpolator* tFIVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity );
        Field_Interpolator* tFITemp     = mFIManager->get_field_interpolators_for_type( mDofTemperature );

        // get the heat capacity
        real tIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" )->val()( 0 );

        // initialize the matrix
        mEnergyDof( tDofIndex ).set_size(    //
                1,                           //
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(),
                0.0 );

        // direct dependency on the density dof type
        if ( aDofTypes( 0 ) == mDofDensity )
        {
            // compute contribution
            mEnergyDof( tDofIndex ) +=
                    tIsochoricHeatCapacity * tFITemp->val() * tFIDensity->N() +    //
                    0.5 * trans( tFIVelocity->val() ) * tFIVelocity->val() * tFIDensity->N();
        }

        // direct dependency on the velocity dof type
        if ( aDofTypes( 0 ) == mDofVelocity )
        {
            // compute contribution
            mEnergyDof( tDofIndex ) +=
                    tFIDensity->val() * trans( tFIVelocity->val() ) * tFIVelocity->N();
        }

        // direct dependency on the temperature dof type
        if ( aDofTypes( 0 ) == mDofTemperature )
        {
            // compute contribution
            mEnergyDof( tDofIndex ) +=
                    tIsochoricHeatCapacity * tFIDensity->val() * tFITemp->N();
        }
    }

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::CM_Fluid_Compressible_Ideal::eval_EnergyDot()
    {
        // get the FIs
        Field_Interpolator* tFIDensity  = mFIManager->get_field_interpolators_for_type( mDofDensity );
        Field_Interpolator* tFIVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity );
        Field_Interpolator* tFITemp     = mFIManager->get_field_interpolators_for_type( mDofTemperature );

        // get the heat capacity
        real tIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" )->val()( 0 );

        // compute total energy density
        mEnergyDot =
                tIsochoricHeatCapacity * tFIDensity->val() * tFITemp->gradt( 1 ) +                   //
                tIsochoricHeatCapacity * tFITemp->val() * tFIDensity->gradt( 1 ) +                   //
                0.5 * trans( tFIVelocity->val() ) * tFIVelocity->val() * tFIDensity->gradt( 1 ) +    //
                tFIDensity->val() * trans( tFIVelocity->val() ) * trans( tFIVelocity->gradt( 1 ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::eval_dEnergyDotdDOF( const Vector< MSI::Dof_Type >& aDofTypes )
    {
        // get the dof type as a uint
        uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( tDofType );

        // get the FIs
        Field_Interpolator* tFIDensity  = mFIManager->get_field_interpolators_for_type( mDofDensity );
        Field_Interpolator* tFIVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity );
        Field_Interpolator* tFITemp     = mFIManager->get_field_interpolators_for_type( mDofTemperature );

        // get the heat capacity
        real tIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" )->val()( 0 );

        // initialize the matrix
        mEnergyDotDof( tDofIndex ).set_size(    //
                1,                              //
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(),
                0.0 );

        // direct dependency on the density dof type
        if ( aDofTypes( 0 ) == mDofDensity )
        {
            // compute contribution
            mEnergyDotDof( tDofIndex ) +=
                    tIsochoricHeatCapacity * tFITemp->val() * tFIDensity->dnNdtn( 1 ) +                   //
                    tIsochoricHeatCapacity * tFITemp->gradt( 1 ) * tFIDensity->N() +                      //
                    0.5 * trans( tFIVelocity->val() ) * tFIVelocity->val() * tFIDensity->dnNdtn( 1 ) +    //
                    trans( tFIVelocity->val() ) * trans( tFIVelocity->gradt( 1 ) ) * tFIDensity->N();
        }

        // direct dependency on the velocity dof type
        if ( aDofTypes( 0 ) == mDofVelocity )
        {
            // compute contribution
            mEnergyDotDof( tDofIndex ) +=
                    tFIDensity->gradt( 1 ) * trans( tFIVelocity->val() ) * tFIVelocity->N() +    //
                    tFIDensity->val() * tFIVelocity->gradt( 1 ) * tFIVelocity->N() +             //
                    tFIDensity->val() * trans( tFIVelocity->val() ) * this->dNveldt();
        }

        // direct dependency on the temperature dof type
        if ( aDofTypes( 0 ) == mDofTemperature )
        {
            // compute contribution
            mEnergyDotDof( tDofIndex ) +=
                    tIsochoricHeatCapacity * tFIDensity->val() * tFITemp->dnNdtn( 1 ) +    //
                    tIsochoricHeatCapacity * tFIDensity->gradt( 1 ) * tFITemp->N();
        }
    }

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::eval_stress()
    {
        // get velocity FI
        Field_Interpolator* tFIVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity );

        // get the viscosity
        const std::shared_ptr< Property > tPropDynamicViscosity = get_property( "DynamicViscosity" );

        // compute Stress
        mStress = 2.0 * tPropDynamicViscosity->val()( 0 ) *                                          //
                          ( this->strain() - ( 1.0 / 3.0 ) * tFIVelocity->div() * mFlatIdentity )    //
                - mFlatIdentity * this->pressure();
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::eval_dStressdDOF( const Vector< MSI::Dof_Type >& aDofTypes )
    {
        // get the dof type as a uint
        uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( tDofType );

        // get the FIs
        Field_Interpolator* tFIVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity );

        // get the viscosity
        const std::shared_ptr< Property > tPropDynamicViscosity = get_property( "DynamicViscosity" );

        // initialize the matrix
        mdStressdDof( tDofIndex ).set_size(    //
                ( mSpaceDim - 1 ) * 3,         //
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(),
                0.0 );

        // direct dependency on the density dof type
        if ( aDofTypes( 0 ) == mDofDensity )
        {
            // compute contribution
            mdStressdDof( tDofIndex ) -=
                    mFlatIdentity * this->dPressuredDOF( aDofTypes );
        }

        // direct dependency on the velocity dof type
        if ( aDofTypes( 0 ) == mDofVelocity )
        {
            // compute contribution
            mdStressdDof( tDofIndex ) +=                         //
                    2.0 * tPropDynamicViscosity->val()( 0 ) *    //
                    ( this->dStraindDOF( aDofTypes ) - ( 1.0 / 3.0 ) * mFlatIdentity * tFIVelocity->div_operator() );
        }

        // direct dependency on the temperature dof type
        if ( aDofTypes( 0 ) == mDofTemperature )
        {
            // compute contribution
            mdStressdDof( tDofIndex ) -=
                    mFlatIdentity * this->dPressuredDOF( aDofTypes );
        }

        // if indirect dependency of viscosity
        if ( tPropDynamicViscosity->check_dof_dependency( aDofTypes ) )
        {
            // compute derivative with indirect dependency through properties
            mdStressdDof( tDofIndex ) +=                                                               //
                    2.0 * ( this->strain() - ( 1.0 / 3.0 ) * tFIVelocity->div() * mFlatIdentity ) *    //
                    tPropDynamicViscosity->dPropdDOF( aDofTypes );
        }
    }

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Ideal::traction(
            const Matrix< DDRMat >& aNormal,
            enum CM_Function_Type   aCMFunctionType )
    {
        switch ( aCMFunctionType )
        {
            case CM_Function_Type::THERMAL:
                return this->thermal_traction( aNormal );

            case CM_Function_Type::ENERGY:
                return this->energy_traction( aNormal );

            case CM_Function_Type::WORK:
                return this->work_traction( aNormal );

            case CM_Function_Type::MECHANICAL:
                return this->mechanical_traction( aNormal );

                // unknown CM function type
            default:
                MORIS_ERROR( false, "CM_Fluid_Compressible_Ideal::traction - unknown CM function type for traction." );
                return mTraction;
        }
    }

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Ideal::dTractiondDOF(
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Matrix< DDRMat >&        aNormal,
            enum CM_Function_Type          aCMFunctionType )
    {
        switch ( aCMFunctionType )
        {
            case CM_Function_Type::THERMAL:
                return this->thermal_dTractiondDOF( aDofTypes, aNormal );

            case CM_Function_Type::ENERGY:
                return this->energy_dTractiondDOF( aDofTypes, aNormal );

            case CM_Function_Type::WORK:
                return this->work_dTractiondDOF( aDofTypes, aNormal );

            case CM_Function_Type::MECHANICAL:
                return this->mechanical_dTractiondDOF( aDofTypes, aNormal );

                // unknown CM function type
            default:
                MORIS_ERROR( false, "CM_Fluid_Compressible_Ideal::dTractiondDOF - unknown CM function type for traction." );
                return mTraction;
        }
    }

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::eval_thermal_traction( const Matrix< DDRMat >& aNormal )
    {
        // compute the traction
        mThermalTraction = trans( aNormal ) * this->flux( CM_Function_Type::THERMAL );
    }

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Ideal::thermal_traction( const Matrix< DDRMat >& aNormal )
    {
        // if the quantity was not evaluated
        if ( mThermalTractionEval )
        {
            // evaluate the test strain
            this->eval_thermal_traction( aNormal );

            // set bool for evaluation
            mThermalTractionEval = false;
        }
        // return the test strain value
        return mThermalTraction;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::eval_thermal_dTractiondDOF(
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Matrix< DDRMat >&        aNormal )
    {
        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // flatten normal
        // Matrix< DDRMat > tFlatNormal;
        // this->flatten_normal( aNormal, tFlatNormal );

        // direct contribution
        mThermalTractionDof( tDofIndex ) = trans( aNormal ) * this->dFluxdDOF( aDofTypes, CM_Function_Type::THERMAL );
    }

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Ideal::thermal_dTractiondDOF(
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Matrix< DDRMat >&        aNormal )
    {
        // if aDofType is not an active dof type for the CM
        MORIS_ERROR( this->check_dof_dependency( aDofTypes ),
                "CM_Fluid_Compressible_Ideal::thermal_dTractiondDOF - no dependency in this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mThermalTractionDofEval( tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_thermal_dTractiondDOF( aDofTypes, aNormal );

            // set bool for evaluation
            mThermalTractionDofEval( tDofIndex ) = false;
        }

        // return the dof deriv value
        return mThermalTractionDof( tDofIndex );
    }

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::eval_energy_traction( const Matrix< DDRMat >& aNormal )
    {
        // compute the traction
        mEnergyTraction = trans( aNormal ) * this->flux( CM_Function_Type::ENERGY );
    }

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Ideal::energy_traction( const Matrix< DDRMat >& aNormal )
    {
        // if quantity not evaluated
        if ( mEnergyTractionEval )
        {
            // evaluate the test strain
            this->eval_energy_traction( aNormal );

            // set bool for evaluation
            mEnergyTractionEval = false;
        }
        // return the test strain value
        return mEnergyTraction;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::eval_energy_dTractiondDOF(
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Matrix< DDRMat >&        aNormal )
    {
        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // direct contribution
        mEnergyTractionDof( tDofIndex ) = trans( aNormal ) * this->dFluxdDOF( aDofTypes, CM_Function_Type::ENERGY );
    }

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Ideal::energy_dTractiondDOF(
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Matrix< DDRMat >&        aNormal )
    {
        // if aDofType is not an active dof type for the CM
        MORIS_ERROR( this->check_dof_dependency( aDofTypes ),
                "CM_Fluid_Compressible_Ideal::energy_dTractiondDOF - no dependency in this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mEnergyTractionDofEval( tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_energy_dTractiondDOF( aDofTypes, aNormal );

            // set bool for evaluation
            mEnergyTractionDofEval( tDofIndex ) = false;
        }

        // return the dof deriv value
        return mEnergyTractionDof( tDofIndex );
    }

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::eval_work_traction( const Matrix< DDRMat >& aNormal )
    {
        // compute the traction
        mWorkTraction = trans( aNormal ) * this->flux( CM_Function_Type::WORK );
    }

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Ideal::work_traction( const Matrix< DDRMat >& aNormal )
    {
        // if quantity not evaluated
        if ( mWorkTractionEval )
        {
            // evaluate the test strain
            this->eval_work_traction( aNormal );

            // set bool for evaluation
            mWorkTractionEval = false;
        }
        // return the test strain value
        return mWorkTraction;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::eval_work_dTractiondDOF(
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Matrix< DDRMat >&        aNormal )
    {
        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // direct contribution
        mWorkTractionDof( tDofIndex ) = trans( aNormal ) * this->dFluxdDOF( aDofTypes, CM_Function_Type::WORK );
    }

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Ideal::work_dTractiondDOF(
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Matrix< DDRMat >&        aNormal )
    {
        // if aDofType is not an active dof type for the CM
        MORIS_ERROR( this->check_dof_dependency( aDofTypes ),
                "CM_Fluid_Compressible_Ideal::work_dTractiondDOF - no dependency in this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mWorkTractionDofEval( tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_work_dTractiondDOF( aDofTypes, aNormal );

            // set bool for evaluation
            mWorkTractionDofEval( tDofIndex ) = false;
        }

        // return the dof deriv value
        return mWorkTractionDof( tDofIndex );
    }

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::eval_mechanical_traction( const Matrix< DDRMat >& aNormal )
    {
        // flatten the normal
        Matrix< DDRMat > tFlatNormal;
        this->flatten_normal( aNormal, tFlatNormal );

        // compute the traction
        mMechanicalTraction = tFlatNormal * this->flux( CM_Function_Type::MECHANICAL );
    }

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Ideal::mechanical_traction( const Matrix< DDRMat >& aNormal )
    {
        // if quantity not evaluated
        if ( mMechanicalTractionEval )
        {
            // evaluate the test strain
            this->eval_mechanical_traction( aNormal );

            // set bool for evaluation
            mMechanicalTractionEval = false;
        }
        // return the test strain value
        return mMechanicalTraction;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::eval_mechanical_dTractiondDOF(
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Matrix< DDRMat >&        aNormal )
    {
        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // flatten the normal
        Matrix< DDRMat > tFlatNormal;
        this->flatten_normal( aNormal, tFlatNormal );

        // direct contribution
        mMechanicalTractionDof( tDofIndex ) = tFlatNormal * this->dFluxdDOF( aDofTypes, CM_Function_Type::MECHANICAL );
    }

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Ideal::mechanical_dTractiondDOF(
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Matrix< DDRMat >&        aNormal )
    {
        // if aDofType is not an active dof type for the CM
        MORIS_ERROR( this->check_dof_dependency( aDofTypes ),
                "CM_Fluid_Compressible_Ideal::mechanical_dTractiondDOF - no dependency in this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mMechanicalTractionDofEval( tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_mechanical_dTractiondDOF( aDofTypes, aNormal );

            // set bool for evaluation
            mMechanicalTractionDofEval( tDofIndex ) = false;
        }

        // return the dof deriv value
        return mMechanicalTractionDof( tDofIndex );
    }

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::eval_pressure()
    {
        // get field interpolator values
        Matrix< DDRMat > tDensity     = mFIManager->get_field_interpolators_for_type( mDofDensity )->val();
        Matrix< DDRMat > tTemperature = mFIManager->get_field_interpolators_for_type( mDofTemperature )->val();

        // get the specific gas constant
        real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

        // compute the pressure
        mPressure = tDensity * tSpecificGasConstant * tTemperature;
    }

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Ideal::pressure()
    {
        // if the test strain was not evaluated
        if ( mPressureEval )
        {
            // evaluate the test strain
            this->eval_pressure();

            // set bool for evaluation
            mPressureEval = false;
        }

        // return the test strain value
        return mPressure;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::eval_dPressuredDOF( const Vector< MSI::Dof_Type >& aDofTypes )
    {
        // get field interpolators
        Field_Interpolator* tFIDensity = mFIManager->get_field_interpolators_for_type( mDofDensity );
        Field_Interpolator* tFITemp    = mFIManager->get_field_interpolators_for_type( mDofTemperature );

        // get the specific gas constant
        real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // initialize the matrix
        mPressureDof( tDofIndex ).set_size(    //
                1,                             //
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(),
                0.0 );

        // if Density DoF
        if ( aDofTypes( 0 ) == mDofDensity )
        {
            mPressureDof( tDofIndex ) = tSpecificGasConstant * tFITemp->val()( 0 ) * tFIDensity->N();
        }

        // if Temperature DoF
        if ( aDofTypes( 0 ) == mDofTemperature )
        {
            // compute derivative with direct dependency
            mPressureDof( tDofIndex ) = tSpecificGasConstant * tFIDensity->val()( 0 ) * tFITemp->N();
        }
    }

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Ideal::dPressuredDOF(
            const Vector< MSI::Dof_Type >& aDofType )
    {
        // if aDofType is not an active dof type for the CM
        MORIS_ERROR(
                this->check_dof_dependency( aDofType ),
                "CM_Fluid_Compressible_Ideal::dPressuredDOF - no dependency in this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mPressureDofEval( tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_dPressuredDOF( aDofType );

            // set bool for evaluation
            mPressureDofEval( tDofIndex ) = false;
        }

        // return the derivative
        return mPressureDof( tDofIndex );
    }

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::eval_strain_2d()
    {
        // get the velocity spatial gradient from velocity FI
        Matrix< DDRMat > tVelocityGradx = mFIManager->get_field_interpolators_for_type( mDofVelocity )->gradx( 1 );

        // evaluate the strain
        mStrain.set_size( 3, 1, 0.0 );
        mStrain( 0, 0 ) = tVelocityGradx( 0, 0 );
        mStrain( 1, 0 ) = tVelocityGradx( 1, 1 );
        mStrain( 2, 0 ) = 0.5 * ( tVelocityGradx( 1, 0 ) + tVelocityGradx( 0, 1 ) );
    }

    void
    CM_Fluid_Compressible_Ideal::eval_strain_3d()
    {
        // get the velocity spatial gradient from velocity FI
        Matrix< DDRMat > tVelocityGradx = mFIManager->get_field_interpolators_for_type( mDofVelocity )->gradx( 1 );

        // evaluate the strain
        mStrain.set_size( 6, 1, 0.0 );
        mStrain( 0, 0 ) = tVelocityGradx( 0, 0 );
        mStrain( 1, 0 ) = tVelocityGradx( 1, 1 );
        mStrain( 2, 0 ) = tVelocityGradx( 2, 2 );
        mStrain( 3, 0 ) = 0.5 * ( tVelocityGradx( 1, 2 ) + tVelocityGradx( 2, 1 ) );
        mStrain( 4, 0 ) = 0.5 * ( tVelocityGradx( 0, 2 ) + tVelocityGradx( 2, 0 ) );
        mStrain( 5, 0 ) = 0.5 * ( tVelocityGradx( 0, 1 ) + tVelocityGradx( 1, 0 ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::eval_teststrain_2d()
    {
        // compute displacement gradient
        Matrix< DDRMat > tdnNdxn = mFIManager->get_field_interpolators_for_type( mDofVelocity )->dnNdxn( 1 );

        // get number of bases for displacement
        uint tNumBases = mFIManager->get_field_interpolators_for_type( mDofVelocity )->get_number_of_space_time_bases();

        // build the test strain
        mTestStrain.set_size( 3, tNumBases * 2, 0.0 );
        mTestStrain( { 0, 0 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
        mTestStrain( { 2, 2 }, { 0, tNumBases - 1 } ) = 0.5 * tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );

        mTestStrain( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
        mTestStrain( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) = 0.5 * tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
    }

    void
    CM_Fluid_Compressible_Ideal::eval_teststrain_3d()
    {
        // compute displacement gradient
        Matrix< DDRMat > tdnNdxn =
                mFIManager->get_field_interpolators_for_type( mDofVelocity )->dnNdxn( 1 );

        // get number of bases for displacement
        uint tNumBases = mFIManager->get_field_interpolators_for_type( mDofVelocity )->get_number_of_space_time_bases();

        // build the test strain
        mTestStrain.set_size( 6, tNumBases * 3, 0.0 );
        mTestStrain( { 0, 0 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
        mTestStrain( { 4, 4 }, { 0, tNumBases - 1 } ) = 0.5 * tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
        mTestStrain( { 5, 5 }, { 0, tNumBases - 1 } ) = 0.5 * tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );

        mTestStrain( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
        mTestStrain( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } ) = 0.5 * tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
        mTestStrain( { 5, 5 }, { tNumBases, 2 * tNumBases - 1 } ) = 0.5 * tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );

        mTestStrain( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
        mTestStrain( { 3, 3 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 0.5 * tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
        mTestStrain( { 4, 4 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 0.5 * tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::eval_dStraindDOF( const Vector< MSI::Dof_Type >& aDofTypes )
    {
        // get the dof type as a uint
        uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

        // get the dof FI
        Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // get the dof type index
        uint tDofIndex = mGlobalDofTypeMap( tDofType );

        // init mdStraindDof
        mdStraindDof( tDofIndex ).set_size( ( mSpaceDim - 1 ) * 3, tFI->get_number_of_space_time_coefficients(), 0.0 );

        // if velocity dof
        if ( aDofTypes( 0 ) == mDofVelocity )
        {
            // compute derivative
            mdStraindDof( tDofIndex ) += this->testStrain();
        }
    }

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::eval_dNveldt()
    {
        // get the residual dof type FI (here velocity)
        Field_Interpolator* tFIVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity );

        // init size for dnNdtn
        uint tNumRowt = tFIVelocity->get_number_of_fields();
        uint tNumColt = tFIVelocity->dnNdtn( 1 ).n_cols();
        mdNveldt.set_size( tNumRowt, tNumRowt * tNumColt, 0.0 );

        // loop over the fields
        for ( uint iField = 0; iField < tNumRowt; iField++ )
        {
            // fill the matrix for each dimension
            mdNveldt( { iField, iField }, { iField * tNumColt, ( iField + 1 ) * tNumColt - 1 } ) =
                    tFIVelocity->dnNdtn( 1 ).matrix_data();
        }
    }

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Ideal::dNveldt()
    {
        // if the velocity matrix was not evaluated
        if ( mdNveldtEval )
        {
            // evaluate the test strain
            this->eval_dNveldt();

            // set bool for evaluation
            mdNveldtEval = false;
        }

        // return the test strain value
        return mdNveldt;
    }

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::eval_velocitymatrix_2d()
    {
        Matrix< DDRMat > tU = mFIManager->get_field_interpolators_for_type( mDofVelocity )->val();

        // clang-format off
             mVelocityMatrix = {
                    { tU( 0 ),    0.0 , tU( 1 ) },
                    {    0.0 , tU( 1 ), tU( 0 ) } };
        // clang-format on
    }

    void
    CM_Fluid_Compressible_Ideal::eval_velocitymatrix_3d()
    {
        Matrix< DDRMat > tU = mFIManager->get_field_interpolators_for_type( mDofVelocity )->val();

        // clang-format off
            mVelocityMatrix = {
                    { tU( 0 ),    0.0 ,  0.0 ,    0.0 , tU( 2 ), tU( 1 ) },
                    {    0.0 , tU( 1 ),  0.0 , tU( 2 ),    0.0 , tU( 0 ) },
                    {    0.0 ,  0.0 , tU( 2 ), tU( 1 ), tU( 0 ),    0.0  } };
        // clang-format on
    }

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Ideal::velocityMatrix()
    {
        // if the velocity matrix was not evaluated
        if ( mVelocityMatrixEval )
        {
            // evaluate the test strain
            this->eval_velocityMatrix();

            // set bool for evaluation
            mVelocityMatrixEval = false;
        }

        // return the test strain value
        return mVelocityMatrix;
    }

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::unfold_2d(
            const Matrix< DDRMat >& aFlattenedTensor,
            Matrix< DDRMat >&       aExpandedTensor )
    {
        aExpandedTensor = {
            { aFlattenedTensor( 0 ), aFlattenedTensor( 2 ) },
            { aFlattenedTensor( 2 ), aFlattenedTensor( 1 ) }
        };
    }

    void
    CM_Fluid_Compressible_Ideal::unfold_3d(
            const Matrix< DDRMat >& aFlattenedTensor,
            Matrix< DDRMat >&       aExpandedTensor )
    {
        aExpandedTensor = {
            { aFlattenedTensor( 0 ), aFlattenedTensor( 5 ), aFlattenedTensor( 4 ) },
            { aFlattenedTensor( 5 ), aFlattenedTensor( 1 ), aFlattenedTensor( 3 ) },
            { aFlattenedTensor( 4 ), aFlattenedTensor( 3 ), aFlattenedTensor( 2 ) }
        };
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Ideal::flatten_normal_2d(
            const Matrix< DDRMat >& aNormal,
            Matrix< DDRMat >&       aFlatNormal )
    {
        aFlatNormal.set_size( 2, 3, 0.0 );
        aFlatNormal( 0, 0 ) = aNormal( 0, 0 );
        aFlatNormal( 0, 2 ) = aNormal( 1, 0 );
        aFlatNormal( 1, 1 ) = aNormal( 1, 0 );
        aFlatNormal( 1, 2 ) = aNormal( 0, 0 );
    }

    void
    CM_Fluid_Compressible_Ideal::flatten_normal_3d(
            const Matrix< DDRMat >& aNormal,
            Matrix< DDRMat >&       aFlatNormal )
    {
        aFlatNormal.set_size( 3, 6, 0.0 );
        aFlatNormal( 0, 0 ) = aNormal( 0, 0 );
        aFlatNormal( 1, 1 ) = aNormal( 1, 0 );
        aFlatNormal( 2, 2 ) = aNormal( 2, 0 );
        aFlatNormal( 0, 4 ) = aNormal( 2, 0 );
        aFlatNormal( 0, 5 ) = aNormal( 1, 0 );
        aFlatNormal( 1, 3 ) = aNormal( 2, 0 );
        aFlatNormal( 1, 5 ) = aNormal( 0, 0 );
        aFlatNormal( 2, 3 ) = aNormal( 1, 0 );
        aFlatNormal( 2, 4 ) = aNormal( 0, 0 );
    }

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::fem
