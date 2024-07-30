/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Compressible_Newtonian_Fluid.cpp
 *
 */

#include "cl_FEM_CM_Compressible_Newtonian_Fluid.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "op_minus.hpp"

namespace moris
{
    namespace fem
    {

        //--------------------------------------------------------------------------------------------------------------

        CM_Compressible_Newtonian_Fluid::CM_Compressible_Newtonian_Fluid()
        {
            // set the property pointer cell size
            mProperties.resize( static_cast< uint >( CM_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "DynamicViscosity" ]    = static_cast< uint >( CM_Property_Type::DYNAMIC_VISCOSITY );       // may be a fnct. of T
            mPropertyMap[ "ThermalConductivity" ] = static_cast< uint >( CM_Property_Type::THERMAL_CONDUCTIVITY );    // may be a fnct. of T

            // set and populate the material model map
            mMaterialModels.resize( static_cast< uint >( MM_Type::MAX_ENUM ), nullptr );
            mMaterialModelMap[ "ThermodynamicMaterialModel" ] = static_cast< uint >( MM_Type::THERMODYNAMIC_MATERIAL_MODEL );    // constant property
        }

        //--------------------------------------------------------------------------------------------------------------

        Constitutive_Type CM_Compressible_Newtonian_Fluid::get_constitutive_type() const
        {
            return Constitutive_Type::FLUID_COMPRESSIBLE_NEWTONIAN;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Compressible_Newtonian_Fluid::set_function_pointers()
        {
            switch ( mSpaceDim )
            {
                case ( 2 ):
                {
                    m_eval_strain           = &CM_Compressible_Newtonian_Fluid::eval_strain_2d;
                    m_eval_teststrain       = &CM_Compressible_Newtonian_Fluid::eval_teststrain_2d;
                    m_eval_divstrainrate    = &CM_Compressible_Newtonian_Fluid::eval_divstrainrate_2d;
                    m_eval_ddivstrainratedu = &CM_Compressible_Newtonian_Fluid::eval_ddivstrainratedu_2d;
                    m_eval_divDivVel        = &CM_Compressible_Newtonian_Fluid::eval_divDivVel_2d;
                    m_eval_dDivDivVeldu     = &CM_Compressible_Newtonian_Fluid::eval_dDivDivVeldu_2d;
                    m_eval_velocitymatrix   = &CM_Compressible_Newtonian_Fluid::eval_velocitymatrix_2d;
                    m_unfold_tensor         = &CM_Compressible_Newtonian_Fluid::unfold_2d;
                    m_flatten_normal        = &CM_Compressible_Newtonian_Fluid::flatten_normal_2d;
                    mFlatIdentity           = { { 1.0 }, { 1.0 }, { 0.0 } };
                    break;
                }
                case ( 3 ):
                {
                    m_eval_strain           = &CM_Compressible_Newtonian_Fluid::eval_strain_3d;
                    m_eval_teststrain       = &CM_Compressible_Newtonian_Fluid::eval_teststrain_3d;
                    m_eval_divstrainrate    = &CM_Compressible_Newtonian_Fluid::eval_divstrainrate_3d;
                    m_eval_ddivstrainratedu = &CM_Compressible_Newtonian_Fluid::eval_ddivstrainratedu_3d;
                    m_eval_divDivVel        = &CM_Compressible_Newtonian_Fluid::eval_divDivVel_3d;
                    m_eval_dDivDivVeldu     = &CM_Compressible_Newtonian_Fluid::eval_dDivDivVeldu_3d;
                    m_eval_velocitymatrix   = &CM_Compressible_Newtonian_Fluid::eval_velocitymatrix_3d;
                    m_unfold_tensor         = &CM_Compressible_Newtonian_Fluid::unfold_3d;
                    m_flatten_normal        = &CM_Compressible_Newtonian_Fluid::flatten_normal_3d;
                    mFlatIdentity           = { { 1.0 }, { 1.0 }, { 1.0 }, { 0.0 }, { 0.0 }, { 0.0 } };
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "CM_Compressible_Newtonian_Fluid::set_function_pointers - this function is currently unused, might be used in the future." );
                    break;
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Compressible_Newtonian_Fluid::reset_specific_eval_flags()
        {
            // reset dof derivative of velocity
            mdNveldtEval = true;

            // reset velocity matrix for flattened tensors
            mVelocityMatrixEval = true;

            // reset thermal div-strain
            mThermalDivStrainEval = true;
            mThermalDivStrainDofEval.fill( true );

            // reset div-strain-rate
            mDivStrainRateEval = true;
            mDivStrainRateDofEval.fill( true );

            // reset div(div(u)*I)
            mDivDivVelEval = true;
            mDivDivVelDofEval.fill( true );

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

            // reset Thermal Div Flux -----------------------------
            mThermalDivFluxEval = true;
            mThermalDivFluxDofEval.fill( true );

            // reset work Div Flux
            mWorkDivFluxEval = true;
            mWorkDivFluxDofEval.fill( true );

            // reset mechanical Flux
            mMechanicalDivFluxEval = true;
            mMechanicalDivFluxDofEval.fill( true );

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
        CM_Compressible_Newtonian_Fluid::initialize_spec_storage_vars_and_eval_flags()
        {
            // get number of DoF types
            uint tNumGlobalDofTypes = mGlobalDofTypes.size();
            uint tNumDirectDofTypes = mDofTypes.size();

            // initialize eval flags
            mThermalDivStrainDofEval.set_size( tNumGlobalDofTypes, 1, true );
            mDivStrainRateDofEval.set_size( tNumGlobalDofTypes, 1, true );
            mDivDivVelDofEval.set_size( tNumGlobalDofTypes, 1, true );
            mThermalFluxDofEval.set_size( tNumGlobalDofTypes, 1, true );
            mWorkFluxDofEval.set_size( tNumGlobalDofTypes, 1, true );
            mEnergyFluxDofEval.set_size( tNumGlobalDofTypes, 1, true );
            // mStressDofEval.resize( tNumGlobalDofTypes, true );
            mThermalDivFluxDofEval.set_size( tNumGlobalDofTypes, 1, true );
            mWorkDivFluxDofEval.set_size( tNumGlobalDofTypes, 1, true );
            mMechanicalDivFluxDofEval.set_size( tNumGlobalDofTypes, 1, true );
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
            mThermalDivStrainDof.resize( tNumGlobalDofTypes );
            mDivStrainRateDof.resize( tNumGlobalDofTypes );
            mDivDivVelDof.resize( tNumGlobalDofTypes );
            mThermalFluxDof.resize( tNumGlobalDofTypes );
            mWorkFluxDof.resize( tNumGlobalDofTypes );
            mEnergyFluxDof.resize( tNumGlobalDofTypes );
            // mStressDof.resize( tNumGlobalDofTypes );
            mThermalDivFluxDof.resize( tNumGlobalDofTypes );
            mWorkDivFluxDof.resize( tNumGlobalDofTypes );
            mMechanicalDivFluxDof.resize( tNumGlobalDofTypes );
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
        CM_Compressible_Newtonian_Fluid::set_dof_type_list(
                Vector< Vector< MSI::Dof_Type > > aDofTypes,
                Vector< std::string >                  aDofStrings )
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
                else if ( tDofString == "Density" )
                {
                    mDofDensity = tDofType;
                }
                else if ( tDofString == "Pressure" )
                {
                    mDofPressure = tDofType;
                }
                else if ( tDofString == "Temperature" )
                {
                    mDofTemperature = tDofType;
                }
                else
                {
                    MORIS_ERROR( false, "CM_Compressible_Newtonian_Fluid::set_dof_type_list - Unknown aDofString : %s", tDofString.c_str() );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        CM_Compressible_Newtonian_Fluid::set_local_properties()
        {
            // get the dynamic viscosity properties
            mPropDynamicViscosity = get_property( "DynamicViscosity" );

            // get the thermal conductivity properties
            mPropThermalConductivity = get_property( "ThermalConductivity" );
        }

        //------------------------------------------------------------------------------

        void
        CM_Compressible_Newtonian_Fluid::set_local_material_model()
        {
            // get the material model
            mMaterialModel = get_material_model( "ThermodynamicMaterialModel" );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::flux( enum CM_Function_Type aCMFunctionType )
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
                    MORIS_ERROR( false, "CM_Compressible_Newtonian_Fluid::flux - DEFAULT function type not supported./" );
                    return mFlux;

                    // unknown CM function type
                default:
                    MORIS_ERROR( false, "CM_Compressible_Newtonian_Fluid::flux - unknown CM function type for flux." );
                    return mFlux;
            }
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::dFluxdDOF(
                const Vector< MSI::Dof_Type > &aDofTypes,
                enum CM_Function_Type               aCMFunctionType )
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
        CM_Compressible_Newtonian_Fluid::eval_thermal_flux()
        {
            // get the material model
            const std::shared_ptr< Material_Model > tMM = get_material_model( "ThermodynamicMaterialModel" );

            // get the thermal conductivity property
            const std::shared_ptr< Property > tThermalConductivity = get_property( "ThermalConductivity" );

            // compute thermal flux q = - k * grad(T)
            mThermalFlux = -1.0 * tThermalConductivity->val()( 0 ) * tMM->dnTemperaturedxn( 1 );
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::thermal_flux()
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
        CM_Compressible_Newtonian_Fluid::eval_thermal_dFluxdDOF( const Vector< MSI::Dof_Type > &aDofTypes )
        {
            // get the material model
            const std::shared_ptr< Material_Model > tMM = get_material_model( "ThermodynamicMaterialModel" );

            // get the conductivity property
            const std::shared_ptr< Property > tPropThermalConductivity = get_property( "ThermalConductivity" );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // initialize the matrix
            mThermalFluxDof( tDofIndex ).set_size( mSpaceDim, mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(), 0.0 );

            // compute derivative with direct dependency
            if ( tMM->check_dof_dependency( aDofTypes ) )
            {
                mThermalFluxDof( tDofIndex ) +=
                        -1.0 * tPropThermalConductivity->val()( 0 ) * tMM->dnTemperaturedxnDOF( aDofTypes, 1 );
            }

            // if indirect dependency of conductivity on the dof type
            if ( tPropThermalConductivity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mThermalFluxDof( tDofIndex ) +=
                        -1.0 * tMM->dnTemperaturedxn( 1 ) * tPropThermalConductivity->dPropdDOF( aDofTypes );
            }
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::thermal_dFluxdDOF(
                const Vector< MSI::Dof_Type > &aDofTypes )
        {
            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofTypes ),
                    "CM_Compressible_Newtonian_Fluid::thermal_dFluxdDOF - no dependency in this dof type." );

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
        CM_Compressible_Newtonian_Fluid::eval_work_flux()
        {
            // compute contribution
            mWorkFlux = this->velocityMatrix() * this->flux( CM_Function_Type::MECHANICAL );
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::work_flux()
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
        CM_Compressible_Newtonian_Fluid::eval_work_dFluxdDOF( const Vector< MSI::Dof_Type > &aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the velocity
            Field_Interpolator *tFIVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity );

            // unfold the flattened stress tensor
            Matrix< DDRMat > tStressTensor;
            this->unfold( this->flux( CM_Function_Type::MECHANICAL ), tStressTensor );

            // initialize the matrix
            mWorkFluxDof( tDofIndex ).set_size( mSpaceDim, mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(), 0.0 );

            // compute contribution, dependency of flux on dof type
            mWorkFluxDof( tDofIndex ) +=
                    this->velocityMatrix() * this->dFluxdDOF( aDofTypes, CM_Function_Type::MECHANICAL );

            // direct dependency on the velocity dof type
            if ( aDofTypes( 0 ) == mDofVelocity )
            {
                // compute contribution
                mWorkFluxDof( tDofIndex ) += tStressTensor * tFIVelocity->N();
            }
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::work_dFluxdDOF(
                const Vector< MSI::Dof_Type > &aDofTypes )
        {
            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofTypes ),
                    "CM_Compressible_Newtonian_Fluid::work_dFluxdDOF - no dependency in this dof type." );

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
        CM_Compressible_Newtonian_Fluid::eval_energy_flux()
        {
            // get the velocity
            Matrix< DDRMat > tVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity )->val();

            // compute contribution
            mEnergyFlux = tVelocity * this->Energy();
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::energy_flux()
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
        CM_Compressible_Newtonian_Fluid::eval_energy_dFluxdDOF( const Vector< MSI::Dof_Type > &aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the velocity
            Field_Interpolator *tFIVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity );

            // initialize the matrix
            mEnergyFluxDof( tDofIndex ).set_size( mSpaceDim, mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(), 0.0 );

            // direct dependency on the density dof type
            mEnergyFluxDof( tDofIndex ) += tFIVelocity->val() * this->dEnergydDOF( aDofTypes );

            // direct dependency on the velocity dof type
            if ( aDofTypes( 0 ) == mDofVelocity )
            {
                // compute contribution
                mEnergyFluxDof( tDofIndex ) += this->Energy()( 0 ) * tFIVelocity->N();
            }
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::energy_dFluxdDOF(
                const Vector< MSI::Dof_Type > &aDofTypes )
        {
            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofTypes ),
                    "CM_Compressible_Newtonian_Fluid::energy_dFluxdDOF - no dependency in this dof type." );

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

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::divflux( enum CM_Function_Type aCMFunctionType )
        {
            switch ( aCMFunctionType )
            {
                case CM_Function_Type::THERMAL:
                    return this->thermal_divflux();

                case CM_Function_Type::WORK:
                    return this->work_divflux();

                case CM_Function_Type::MECHANICAL:
                    return this->mechanical_divflux();

                case CM_Function_Type::DEFAULT:
                    MORIS_ERROR( false, "CM_Compressible_Newtonian_Fluid::divflux - DEFAULT function type not supported./" );
                    return mDivFlux;

                    // unknown CM function type
                default:
                    MORIS_ERROR( false, "CM_Compressible_Newtonian_Fluid::divflux - unknown CM function type for flux." );
                    return mDivFlux;
            }
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::ddivfluxdu(
                const Vector< MSI::Dof_Type > &aDofTypes,
                enum CM_Function_Type               aCMFunctionType )
        {
            switch ( aCMFunctionType )
            {
                case CM_Function_Type::THERMAL:
                    return this->thermal_ddivfluxdu( aDofTypes );

                case CM_Function_Type::WORK:
                    return this->work_ddivfluxdu( aDofTypes );

                case CM_Function_Type::MECHANICAL:
                    return this->mechanical_ddivfluxdu( aDofTypes );

                case CM_Function_Type::DEFAULT:
                    MORIS_ERROR( false, "CM_Compressible_Newtonian_Fluid::ddivfluxdu - DEFAULT function type not supported./" );
                    return mddivfluxdu( 0 );

                    // unknown CM function type
                default:
                    MORIS_ERROR( false, "CM_Compressible_Newtonian_Fluid::ddivfluxdu - unknown CM function type for flux." );
                    return mddivfluxdu( 0 );
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Compressible_Newtonian_Fluid::eval_mechanical_divflux()
        {
            // get the Dynamic Viscosity
            const std::shared_ptr< Property > tPropDynamicViscosity = get_property( "DynamicViscosity" );

            // compute divflux
            mMechanicalDivFlux = 2.0 * tPropDynamicViscosity->val()( 0 ) * ( this->divstrain( CM_Function_Type::MECHANICAL ) - ( 1.0 / 3.0 ) * this->divDivVel() );

            // FIXME assume that viscosity prop does not depend on x
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::mechanical_divflux()
        {
            // if the test strain was not evaluated
            if ( mMechanicalDivFluxEval )
            {
                // evaluate the test strain
                this->eval_mechanical_divflux();

                // set bool for evaluation
                mMechanicalDivFluxEval = false;
            }
            // return the test strain value
            return mMechanicalDivFlux;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Compressible_Newtonian_Fluid::eval_mechanical_ddivfluxdu( const Vector< MSI::Dof_Type > &aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the FI for the dependent DoF type
            Field_Interpolator *tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get the Dynamic Viscosity
            const std::shared_ptr< Property > tPropDynamicViscosity = get_property( "DynamicViscosity" );

            // set size for ddivflux/du
            mMechanicalDivFluxDof( tDofIndex ).set_size( mSpaceDim, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // if velocity dof
            if ( aDofTypes( 0 ) == mDofVelocity )
            {
                // fill
                mMechanicalDivFluxDof( tDofIndex ) += 2.0 * tPropDynamicViscosity->val()( 0 ) * ( this->ddivstraindu( aDofTypes, CM_Function_Type::MECHANICAL ) - ( 1.0 / 3.0 ) * this->dDivDivVeldu( aDofTypes ) );
            }

            // dependency of the viscosity
            if ( tPropDynamicViscosity->check_dof_dependency( aDofTypes ) )
            {
                // fill ddivstrain/du
                mMechanicalDivFluxDof( tDofIndex ) += 2.0 * ( this->divstrain( CM_Function_Type::MECHANICAL ) + this->divDivVel() ) * tPropDynamicViscosity->dPropdDOF( aDofTypes );
            }
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::mechanical_ddivfluxdu(
                const Vector< MSI::Dof_Type > &aDofTypes )
        {
            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofTypes ),
                    "CM_Compressible_Newtonian_Fluid::mechanical_dDivFluxdDOF - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mMechanicalDivFluxDofEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_mechanical_ddivfluxdu( aDofTypes );

                // set bool for evaluation
                mMechanicalDivFluxDofEval( tDofIndex ) = false;
            }

            // return the dof deriv value
            return mMechanicalDivFluxDof( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Compressible_Newtonian_Fluid::eval_thermal_divflux()
        {
            // get the Thermal Conductivity
            const std::shared_ptr< Property > tPropThermalConductivity = get_property( "ThermalConductivity" );

            // compute divflux
            mThermalDivFlux = -1.0 * tPropThermalConductivity->val() * this->divstrain( CM_Function_Type::THERMAL );

            // FIXME assume that thermal conductivity prop does not depend on x
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::thermal_divflux()
        {
            // if the test strain was not evaluated
            if ( mThermalDivFluxEval )
            {
                // evaluate the test strain
                this->eval_thermal_divflux();

                // set bool for evaluation
                mThermalDivFluxEval = false;
            }
            // return the test strain value
            return mThermalDivFlux;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Compressible_Newtonian_Fluid::eval_thermal_ddivfluxdu( const Vector< MSI::Dof_Type > &aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the FI for the dependent DoF type
            Field_Interpolator *tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get the Dynamic Viscosity
            const std::shared_ptr< Property > tPropThermalConductivity = get_property( "ThermalConductivity" );

            // set size for ddivflux/du
            mThermalDivFluxDof( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // if velocity dof
            if ( aDofTypes( 0 ) == mDofTemperature )
            {
                // fill
                mThermalDivFluxDof( tDofIndex ) -=
                        tPropThermalConductivity->val()( 0 ) * this->ddivstraindu( aDofTypes, CM_Function_Type::THERMAL );
            }

            // dependency of the viscosity
            if ( tPropThermalConductivity->check_dof_dependency( aDofTypes ) )
            {
                // fill
                mThermalDivFluxDof( tDofIndex ) -=
                        this->divstrain( CM_Function_Type::THERMAL ) * tPropThermalConductivity->dPropdDOF( aDofTypes );
            }
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::thermal_ddivfluxdu(
                const Vector< MSI::Dof_Type > &aDofTypes )
        {
            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofTypes ),
                    "CM_Compressible_Newtonian_Fluid::thermal_ddivfluxdu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mThermalDivFluxDofEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_thermal_ddivfluxdu( aDofTypes );

                // set bool for evaluation
                mThermalDivFluxDofEval( tDofIndex ) = false;
            }

            // return the dof deriv value
            return mThermalDivFluxDof( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Compressible_Newtonian_Fluid::eval_work_divflux()
        {
            // get the Velocity FI
            Field_Interpolator *tVelFI = mFIManager->get_field_interpolators_for_type( mDofVelocity );

            // compute divflux
            mWorkDivFlux =
                    trans( this->flux( CM_Function_Type::MECHANICAL ) ) * this->MultipMat() * this->strain() + tVelFI->val_trans() * this->divflux( CM_Function_Type::MECHANICAL );
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::work_divflux()
        {
            // if the test strain was not evaluated
            if ( mWorkDivFluxEval )
            {
                // evaluate the test strain
                this->eval_work_divflux();

                // set bool for evaluation
                mWorkDivFluxEval = false;
            }
            // return the test strain value
            return mWorkDivFlux;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Compressible_Newtonian_Fluid::eval_work_ddivfluxdu( const Vector< MSI::Dof_Type > &aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the FI for the dependent DoF type
            Field_Interpolator *tDepFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );
            Field_Interpolator *tVelFI = mFIManager->get_field_interpolators_for_type( mDofVelocity );

            // set size for ddivflux/du
            mWorkDivFluxDof( tDofIndex ).set_size( 1, tDepFI->get_number_of_space_time_coefficients(), 0.0 );

            // fill
            mWorkDivFluxDof( tDofIndex ) +=
                    trans( this->strain() ) * this->MultipMat() * this->dFluxdDOF( aDofTypes, CM_Function_Type::MECHANICAL ) + tVelFI->val_trans() * this->ddivfluxdu( aDofTypes, CM_Function_Type::MECHANICAL );

            // if velocity dof
            if ( aDofTypes( 0 ) == mDofVelocity )
            {
                // fill
                mWorkDivFluxDof( tDofIndex ) +=
                        trans( this->flux( CM_Function_Type::MECHANICAL ) ) * this->MultipMat() * this->dStraindDOF( aDofTypes ) + trans( this->divflux( CM_Function_Type::MECHANICAL ) ) * tVelFI->N();
            }
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::work_ddivfluxdu(
                const Vector< MSI::Dof_Type > &aDofTypes )
        {
            // if aDofType is not an active dof type for the CM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofTypes ),
                    "CM_Compressible_Newtonian_Fluid::work_ddivfluxdu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mWorkDivFluxDofEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_work_ddivfluxdu( aDofTypes );

                // set bool for evaluation
                mWorkDivFluxDofEval( tDofIndex ) = false;
            }

            // return the dof deriv value
            return mWorkDivFluxDof( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Compressible_Newtonian_Fluid::eval_Energy()
        {
            // get the material model
            const std::shared_ptr< Material_Model > tMM = get_material_model( "ThermodynamicMaterialModel" );

            // get field interpolator values
            Matrix< DDRMat > tVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity )->val();

            // compute thermal flux q = - k * grad(T)
            mEnergy = tMM->Eint() * tMM->density() + 0.5 * trans( tVelocity ) * tVelocity * tMM->density();
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Compressible_Newtonian_Fluid::eval_dEnergydDOF( const Vector< MSI::Dof_Type > &aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the material model
            const std::shared_ptr< Material_Model > tMM = get_material_model( "ThermodynamicMaterialModel" );

            // check material model in debug
            MORIS_ASSERT( tMM != nullptr, "CM_Compressible_Newtonian_Fluid::eval_dEnergydDOF - Material Model not set" );

            // get the velocity FI
            Field_Interpolator *tFIVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity );

            // initialize the matrix
            mEnergyDof( tDofIndex ).set_size( 1, mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(), 0.0 );

            // direct dependency of internal energy on the dof type
            if ( tMM->check_dof_dependency( aDofTypes ) )
            {
                mEnergyDof( tDofIndex ) +=
                        tMM->Eint() * tMM->DensityDOF( aDofTypes ) + tMM->density() * tMM->EintDOF( aDofTypes ) + 0.5 * trans( tFIVelocity->val() ) * tFIVelocity->val() * tMM->DensityDOF( aDofTypes );
            }

            // direct dependency on the velocity dof type
            if ( aDofTypes( 0 ) == mDofVelocity )
            {
                // compute contribution
                mEnergyDof( tDofIndex ) +=
                        tMM->density() * trans( tFIVelocity->val() ) * tFIVelocity->N();
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Compressible_Newtonian_Fluid::eval_EnergyDot()
        {
            // get the velocity FI
            Field_Interpolator *tFIVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity );

            // get the material model
            const std::shared_ptr< Material_Model > tMM = get_material_model( "ThermodynamicMaterialModel" );

            // check material model in debug
            MORIS_ASSERT( tMM != nullptr, "CM_Compressible_Newtonian_Fluid::eval_... - Material Model not set" );

            // compute total energy density
            mEnergyDot =
                    tMM->DensityDot()( 0 ) * ( tMM->Eint() + 0.5 * trans( tFIVelocity->val() ) * tFIVelocity->val() ) + tMM->density()( 0 ) * ( tMM->EintDot() + trans( tFIVelocity->val() ) * trans( tFIVelocity->gradt( 1 ) ) );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Compressible_Newtonian_Fluid::eval_dEnergyDotdDOF( const Vector< MSI::Dof_Type > &aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the material model
            const std::shared_ptr< Material_Model > tMM = get_material_model( "ThermodynamicMaterialModel" );

            // check material model in debug
            MORIS_ASSERT( tMM != nullptr, "CM_Compressible_Newtonian_Fluid::eval_... - Material Model not set" );

            // get the velocity FI
            Field_Interpolator *tFIVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity );

            // initialize the matrix
            mEnergyDotDof( tDofIndex ).set_size( 1, mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(), 0.0 );

            // direct dependency on the density dof type
            if ( tMM->check_dof_dependency( aDofTypes ) )
            {
                // compute contribution from Density (dependent) DoF
                mEnergyDotDof( tDofIndex ) +=
                        ( tMM->Eint() + 0.5 * trans( tFIVelocity->val() ) * tFIVelocity->val() ) * tMM->DensityDotDOF( aDofTypes ) + ( tMM->EintDot() + trans( tFIVelocity->val() ) * trans( tFIVelocity->gradt( 1 ) ) ) * tMM->DensityDOF( aDofTypes );

                // contribution from internal energy
                mEnergyDotDof( tDofIndex ) +=
                        tMM->density()( 0 ) * tMM->EintDotDOF( aDofTypes ) + tMM->DensityDot()( 0 ) * tMM->EintDOF( aDofTypes );
            }

            // direct dependency on the velocity dof type
            if ( aDofTypes( 0 ) == mDofVelocity )
            {
                // compute contribution
                mEnergyDotDof( tDofIndex ) +=
                        tMM->DensityDot()( 0 ) * trans( tFIVelocity->val() ) * tFIVelocity->N() + tMM->density()( 0 ) * tFIVelocity->gradt( 1 ) * tFIVelocity->N() + tMM->density()( 0 ) * trans( tFIVelocity->val() ) * this->dNveldt();
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Compressible_Newtonian_Fluid::eval_stress()
        {
            // get velocity FI
            Field_Interpolator *tFIVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity );

            // get the viscosity
            const std::shared_ptr< Property > tPropDynamicViscosity = get_property( "DynamicViscosity" );

            // compute Stress
            mStress = 2.0 * tPropDynamicViscosity->val()( 0 ) * ( this->strain() - ( 1.0 / 3.0 ) * tFIVelocity->div() * mFlatIdentity );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Compressible_Newtonian_Fluid::eval_dStressdDOF( const Vector< MSI::Dof_Type > &aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the FIs
            Field_Interpolator *tFIVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity );

            // get the viscosity
            const std::shared_ptr< Property > tPropDynamicViscosity = get_property( "DynamicViscosity" );

            // initialize the matrix
            mdStressdDof( tDofIndex ).set_size( ( mSpaceDim - 1 ) * 3, mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(), 0.0 );

            // direct dependency on the velocity dof type
            if ( aDofTypes( 0 ) == mDofVelocity )
            {
                // compute contribution
                mdStressdDof( tDofIndex ) +=
                        2.0 * tPropDynamicViscosity->val()( 0 ) * ( this->dStraindDOF( aDofTypes ) - ( 1.0 / 3.0 ) * mFlatIdentity * tFIVelocity->div_operator() );
            }

            // if indirect dependency of viscosity
            if ( tPropDynamicViscosity->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mdStressdDof( tDofIndex ) +=
                        2.0 * ( this->strain() - ( 1.0 / 3.0 ) * tFIVelocity->div() * mFlatIdentity ) * tPropDynamicViscosity->dPropdDOF( aDofTypes );
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::traction(
                const Matrix< DDRMat > &aNormal,
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
                    MORIS_ERROR( false, "CM_Compressible_Newtonian_Fluid::traction - unknown CM function type for traction." );
                    return mTraction;
            }
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::dTractiondDOF(
                const Vector< MSI::Dof_Type > &aDofTypes,
                const Matrix< DDRMat >             &aNormal,
                enum CM_Function_Type               aCMFunctionType )
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
                    MORIS_ERROR( false, "CM_Compressible_Newtonian_Fluid::dTractiondDOF - unknown CM function type for traction." );
                    return mTraction;
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Compressible_Newtonian_Fluid::eval_thermal_traction( const Matrix< DDRMat > &aNormal )
        {
            // compute the traction
            mThermalTraction = trans( aNormal ) * this->flux( CM_Function_Type::THERMAL );
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::thermal_traction( const Matrix< DDRMat > &aNormal )
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
        CM_Compressible_Newtonian_Fluid::eval_thermal_dTractiondDOF(
                const Vector< MSI::Dof_Type > &aDofTypes,
                const Matrix< DDRMat >             &aNormal )
        {
            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // flatten normal
            // Matrix< DDRMat > tFlatNormal;
            // this->flatten_normal( aNormal, tFlatNormal );

            // direct contribution
            mThermalTractionDof( tDofIndex ) = trans( aNormal ) * this->dFluxdDOF( aDofTypes, CM_Function_Type::THERMAL );
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::thermal_dTractiondDOF(
                const Vector< MSI::Dof_Type > &aDofTypes,
                const Matrix< DDRMat >             &aNormal )
        {
            // if aDofType is not an active dof type for the CM
            MORIS_ERROR( this->check_dof_dependency( aDofTypes ),
                    "CM_Compressible_Newtonian_Fluid::thermal_dTractiondDOF - no dependency in this dof type." );

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
        CM_Compressible_Newtonian_Fluid::eval_energy_traction( const Matrix< DDRMat > &aNormal )
        {
            // compute the traction
            mEnergyTraction = trans( aNormal ) * this->flux( CM_Function_Type::ENERGY );
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::energy_traction( const Matrix< DDRMat > &aNormal )
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
        CM_Compressible_Newtonian_Fluid::eval_energy_dTractiondDOF(
                const Vector< MSI::Dof_Type > &aDofTypes,
                const Matrix< DDRMat >             &aNormal )
        {
            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // direct contribution
            mEnergyTractionDof( tDofIndex ) = trans( aNormal ) * this->dFluxdDOF( aDofTypes, CM_Function_Type::ENERGY );
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::energy_dTractiondDOF(
                const Vector< MSI::Dof_Type > &aDofTypes,
                const Matrix< DDRMat >             &aNormal )
        {
            // if aDofType is not an active dof type for the CM
            MORIS_ERROR( this->check_dof_dependency( aDofTypes ),
                    "CM_Compressible_Newtonian_Fluid::energy_dTractiondDOF - no dependency in this dof type." );

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
        CM_Compressible_Newtonian_Fluid::eval_work_traction( const Matrix< DDRMat > &aNormal )
        {
            // compute the traction
            mWorkTraction = trans( aNormal ) * this->flux( CM_Function_Type::WORK );
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::work_traction( const Matrix< DDRMat > &aNormal )
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
        CM_Compressible_Newtonian_Fluid::eval_work_dTractiondDOF(
                const Vector< MSI::Dof_Type > &aDofTypes,
                const Matrix< DDRMat >             &aNormal )
        {
            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // direct contribution
            mWorkTractionDof( tDofIndex ) = trans( aNormal ) * this->dFluxdDOF( aDofTypes, CM_Function_Type::WORK );
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::work_dTractiondDOF(
                const Vector< MSI::Dof_Type > &aDofTypes,
                const Matrix< DDRMat >             &aNormal )
        {
            // if aDofType is not an active dof type for the CM
            MORIS_ERROR( this->check_dof_dependency( aDofTypes ),
                    "CM_Compressible_Newtonian_Fluid::work_dTractiondDOF - no dependency in this dof type." );

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
        CM_Compressible_Newtonian_Fluid::eval_mechanical_traction( const Matrix< DDRMat > &aNormal )
        {
            // flatten the normal
            Matrix< DDRMat > tFlatNormal;
            this->flatten_normal( aNormal, tFlatNormal );

            // compute the traction
            mMechanicalTraction = tFlatNormal * this->flux( CM_Function_Type::MECHANICAL );
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::mechanical_traction( const Matrix< DDRMat > &aNormal )
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
        CM_Compressible_Newtonian_Fluid::eval_mechanical_dTractiondDOF(
                const Vector< MSI::Dof_Type > &aDofTypes,
                const Matrix< DDRMat >             &aNormal )
        {
            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // flatten the normal
            Matrix< DDRMat > tFlatNormal;
            this->flatten_normal( aNormal, tFlatNormal );

            // direct contribution
            mMechanicalTractionDof( tDofIndex ) = tFlatNormal * this->dFluxdDOF( aDofTypes, CM_Function_Type::MECHANICAL );
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::mechanical_dTractiondDOF(
                const Vector< MSI::Dof_Type > &aDofTypes,
                const Matrix< DDRMat >             &aNormal )
        {
            // if aDofType is not an active dof type for the CM
            MORIS_ERROR( this->check_dof_dependency( aDofTypes ),
                    "CM_Compressible_Newtonian_Fluid::mechanical_dTractiondDOF - no dependency in this dof type." );

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
        CM_Compressible_Newtonian_Fluid::eval_strain_2d()
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
        CM_Compressible_Newtonian_Fluid::eval_strain_3d()
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
        CM_Compressible_Newtonian_Fluid::eval_teststrain_2d()
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
        CM_Compressible_Newtonian_Fluid::eval_teststrain_3d()
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
        CM_Compressible_Newtonian_Fluid::eval_dStraindDOF( const Vector< MSI::Dof_Type > &aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof FI
            Field_Interpolator *tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

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

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::divstrain( enum CM_Function_Type aCMFunctionType )
        {
            switch ( aCMFunctionType )
            {
                case CM_Function_Type::THERMAL:
                    return this->thermal_divstrain();

                case CM_Function_Type::MECHANICAL:
                    return this->divstrainrate();

                case CM_Function_Type::DEFAULT:
                    MORIS_ERROR( false, "CM_Compressible_Newtonian_Fluid::divstrain - DEFAULT function type not supported./" );
                    return mDivStrain;

                    // unknown CM function type
                default:
                    MORIS_ERROR( false, "CM_Compressible_Newtonian_Fluid::divstrain - unknown CM function type for flux." );
                    return mDivStrain;
            }
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::ddivstraindu(
                const Vector< MSI::Dof_Type > &aDofTypes,
                enum CM_Function_Type               aCMFunctionType )
        {
            switch ( aCMFunctionType )
            {
                case CM_Function_Type::THERMAL:
                    return this->thermal_ddivstraindu( aDofTypes );

                case CM_Function_Type::MECHANICAL:
                    return this->ddivstrainratedu( aDofTypes );

                case CM_Function_Type::DEFAULT:
                    MORIS_ERROR( false, "CM_Compressible_Newtonian_Fluid::ddivstraindu - DEFAULT function type not supported./" );
                    return mddivstraindu( 0 );

                    // unknown CM function type
                default:
                    MORIS_ERROR( false, "CM_Compressible_Newtonian_Fluid::ddivstraindu - unknown CM function type for flux." );
                    return mddivstraindu( 0 );
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Compressible_Newtonian_Fluid::eval_thermal_divstrain()
        {
            // get the temperature gradient
            Matrix< DDRMat > td2Tempdx2 =
                    mFIManager->get_field_interpolators_for_type( mDofTemperature )->gradx( 2 );

            // evaluate the divergence of the strain
            mThermalDivStrain = sum( td2Tempdx2( { 0, mSpaceDim - 1 }, { 0, 0 } ) );
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::thermal_divstrain()
        {
            // if the velocity matrix was not evaluated
            if ( mThermalDivStrainEval )
            {
                // evaluate the test strain
                this->eval_thermal_divstrain();

                // set bool for evaluation
                mThermalDivStrainEval = false;
            }

            // return the test strain value
            return mThermalDivStrain;
        }

        //------------------------------------------------------------------------------

        void
        CM_Compressible_Newtonian_Fluid::eval_thermal_ddivstraindu(
                const Vector< MSI::Dof_Type > &aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the temperature FI
            Field_Interpolator *tTempFI = mFIManager->get_field_interpolators_for_type( mDofTemperature );

            // set size for ddivstrain/du
            mThermalDivStrainDof( tDofIndex ).set_size( 1, tTempFI->get_number_of_space_time_coefficients(), 0.0 );

            // if temperature dof type
            if ( aDofTypes( 0 ) == mDofTemperature )
            {
                // get the 2nd order derivative of the shape functions d2Ndx2
                Matrix< DDRMat > tTempd2Ndx2 = tTempFI->dnNdxn( 2 );

                // fill ddivstrain/du
                mThermalDivStrainDof( tDofIndex ) = tTempd2Ndx2.get_row( 0 ) + tTempd2Ndx2.get_row( 1 );

                if ( tTempd2Ndx2.n_rows() == 6 )
                {
                    mThermalDivStrainDof( tDofIndex ) += tTempd2Ndx2.get_row( 2 );
                }
            }
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::thermal_ddivstraindu(
                const Vector< MSI::Dof_Type > &aDofTypes )
        {
            // if aDofType is not an active dof type for the CM
            MORIS_ERROR( this->check_dof_dependency( aDofTypes ),
                    "CM_Compressible_Newtonian_Fluid::thermal_ddivstraindu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mThermalDivStrainDofEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_thermal_ddivstraindu( aDofTypes );

                // set bool for evaluation
                mThermalDivStrainDofEval( tDofIndex ) = false;
            }

            // return the dof deriv value
            return mThermalDivStrainDof( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Compressible_Newtonian_Fluid::eval_divstrainrate_2d()
        {
            // set size for div strain
            mDivStrainRate.set_size( 2, 1, 0.0 );

            // get the velocity gradient
            Matrix< DDRMat > tVelocityGrad =
                    mFIManager->get_field_interpolators_for_type( mDofVelocity )->gradx( 2 );

            // fill div strain
            mDivStrainRate( 0 ) = tVelocityGrad( 0, 0 ) + 0.5 * ( tVelocityGrad( 1, 0 ) + tVelocityGrad( 2, 1 ) );
            mDivStrainRate( 1 ) = 0.5 * ( tVelocityGrad( 2, 0 ) + tVelocityGrad( 0, 1 ) ) + tVelocityGrad( 1, 1 );
        }

        void
        CM_Compressible_Newtonian_Fluid::eval_divstrainrate_3d()
        {
            // set size for div strain
            mDivStrainRate.set_size( 3, 1, 0.0 );

            // get the velocity gradient
            Matrix< DDRMat > tVelocityGrad =
                    mFIManager->get_field_interpolators_for_type( mDofVelocity )->gradx( 2 );

            // fill div strain
            mDivStrainRate( 0 ) = tVelocityGrad( 0, 0 )
                                + 0.5 * ( tVelocityGrad( 1, 0 ) + tVelocityGrad( 5, 1 ) )
                                + 0.5 * ( tVelocityGrad( 2, 0 ) + tVelocityGrad( 4, 2 ) );
            mDivStrainRate( 1 ) = 0.5 * ( tVelocityGrad( 5, 0 ) + tVelocityGrad( 0, 1 ) )
                                + tVelocityGrad( 1, 1 )
                                + 0.5 * ( tVelocityGrad( 2, 1 ) + tVelocityGrad( 3, 2 ) );
            mDivStrainRate( 2 ) = 0.5 * ( tVelocityGrad( 4, 0 ) + tVelocityGrad( 0, 2 ) )
                                + 0.5 * ( tVelocityGrad( 3, 1 ) + tVelocityGrad( 1, 2 ) )
                                + tVelocityGrad( 2, 2 );
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::divstrainrate()
        {
            // if the velocity matrix was not evaluated
            if ( mDivStrainRateEval )
            {
                // evaluate the test strain
                this->eval_divstrainrate();

                // set bool for evaluation
                mDivStrainRateEval = false;
            }

            // return the test strain value
            return mDivStrainRate;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Compressible_Newtonian_Fluid::eval_ddivstrainratedu_2d( const Vector< MSI::Dof_Type > &aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator *tFI =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for ddivstrain/du
            mDivStrainRateDof( tDofIndex ).set_size( 2, tFI->get_number_of_space_time_coefficients(), 0.0 );

            if ( aDofTypes( 0 ) == mDofVelocity )
            {
                // get the 2nd order derivative of the shape functions d2Ndx2
                Matrix< DDRMat > tVelocityd2Ndx2 = tFI->dnNdxn( 2 );

                // get number of space time bases
                uint tNumBases = tFI->get_number_of_space_time_bases();

                // fill ddivstrain/du
                mDivStrainRateDof( tDofIndex )( { 0, 0 }, { 0, tNumBases - 1 } )             = tVelocityd2Ndx2.get_row( 0 ) + 0.5 * tVelocityd2Ndx2.get_row( 1 );
                mDivStrainRateDof( tDofIndex )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } ) = 0.5 * tVelocityd2Ndx2.get_row( 2 );
                mDivStrainRateDof( tDofIndex )( { 1, 1 }, { 0, tNumBases - 1 } )             = 0.5 * tVelocityd2Ndx2.get_row( 2 );
                mDivStrainRateDof( tDofIndex )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = 0.5 * tVelocityd2Ndx2.get_row( 0 ) + tVelocityd2Ndx2.get_row( 1 );
            }
        }

        void
        CM_Compressible_Newtonian_Fluid::eval_ddivstrainratedu_3d( const Vector< MSI::Dof_Type > &aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator *tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for ddivstrain/du
            mDivStrainRateDof( tDofIndex ).set_size( 3, tFI->get_number_of_space_time_coefficients(), 0.0 );

            if ( aDofTypes( 0 ) == mDofVelocity )
            {
                // get the 2nd order derivative of the shape functions d2Ndx2
                Matrix< DDRMat > tVelocityd2Ndx2 = tFI->dnNdxn( 2 );

                // get number of space time bases
                uint tNumBases = tFI->get_number_of_space_time_bases();

                // fill ddivstrain/du
                mDivStrainRateDof( tDofIndex )( { 0, 0 }, { 0, tNumBases - 1 } )                 = tVelocityd2Ndx2.get_row( 0 ) + 0.5 * tVelocityd2Ndx2.get_row( 1 ) + 0.5 * tVelocityd2Ndx2.get_row( 2 );
                mDivStrainRateDof( tDofIndex )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } )     = 0.5 * tVelocityd2Ndx2.get_row( 5 );
                mDivStrainRateDof( tDofIndex )( { 0, 0 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 0.5 * tVelocityd2Ndx2.get_row( 4 );

                mDivStrainRateDof( tDofIndex )( { 1, 1 }, { 0, tNumBases - 1 } )                 = 0.5 * tVelocityd2Ndx2.get_row( 5 );
                mDivStrainRateDof( tDofIndex )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } )     = 0.5 * tVelocityd2Ndx2.get_row( 0 ) + tVelocityd2Ndx2.get_row( 1 ) + 0.5 * tVelocityd2Ndx2.get_row( 2 );
                mDivStrainRateDof( tDofIndex )( { 1, 1 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 0.5 * tVelocityd2Ndx2.get_row( 3 );

                mDivStrainRateDof( tDofIndex )( { 2, 2 }, { 0, tNumBases - 1 } )                 = 0.5 * tVelocityd2Ndx2.get_row( 4 );
                mDivStrainRateDof( tDofIndex )( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } )     = 0.5 * tVelocityd2Ndx2.get_row( 3 );
                mDivStrainRateDof( tDofIndex )( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 0.5 * tVelocityd2Ndx2.get_row( 0 ) + 0.5 * tVelocityd2Ndx2.get_row( 1 ) + tVelocityd2Ndx2.get_row( 2 );
            }
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::ddivstrainratedu(
                const Vector< MSI::Dof_Type > &aDofTypes )
        {
            // if aDofType is not an active dof type for the CM
            MORIS_ERROR( this->check_dof_dependency( aDofTypes ),
                    "CM_Compressible_Newtonian_Fluid::ddivstrainratedu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mDivStrainRateDofEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_ddivstrainratedu( aDofTypes );

                // set bool for evaluation
                mDivStrainRateDofEval( tDofIndex ) = false;
            }

            // return the dof deriv value
            return mDivStrainRateDof( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Compressible_Newtonian_Fluid::eval_divDivVel_2d()
        {
            // initialize div(div(u)*I)
            mDivDivVel.set_size( 2, 1, 0.0 );

            // get the second derivative
            Matrix< DDRMat > td2Veldx2 = mFIManager->get_field_interpolators_for_type( mDofVelocity )->gradx( 2 );

            // fill div strain
            mDivDivVel( 0 ) = td2Veldx2( 0, 0 ) + td2Veldx2( 2, 1 );
            mDivDivVel( 1 ) = td2Veldx2( 2, 0 ) + td2Veldx2( 1, 1 );
        }

        void
        CM_Compressible_Newtonian_Fluid::eval_divDivVel_3d()
        {
            // initialize div(div(u)*I)
            mDivDivVel.set_size( 3, 1, 0.0 );

            // get the second derivative
            Matrix< DDRMat > td2Veldx2 = mFIManager->get_field_interpolators_for_type( mDofVelocity )->gradx( 2 );

            // fill div strain
            mDivDivVel( 0 ) = td2Veldx2( 0, 0 ) + td2Veldx2( 5, 1 ) + td2Veldx2( 4, 2 );
            mDivDivVel( 1 ) = td2Veldx2( 5, 0 ) + td2Veldx2( 1, 1 ) + td2Veldx2( 3, 2 );
            mDivDivVel( 2 ) = td2Veldx2( 4, 0 ) + td2Veldx2( 3, 1 ) + td2Veldx2( 2, 2 );
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::divDivVel()
        {
            // if the velocity matrix was not evaluated
            if ( mDivDivVelEval )
            {
                // evaluate the test strain
                this->eval_divDivVel();

                // set bool for evaluation
                mDivDivVelEval = false;
            }

            // return the test strain value
            return mDivDivVel;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Compressible_Newtonian_Fluid::eval_dDivDivVeldu_2d( const Vector< MSI::Dof_Type > &aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the velocity FI
            Field_Interpolator *tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for ddivstrain/du
            mDivDivVelDof( tDofIndex ).set_size( 2, tFI->get_number_of_space_time_coefficients(), 0.0 );

            if ( aDofTypes( 0 ) == mDofVelocity )
            {
                // get the 2nd order derivative of the shape functions d2Ndx2
                Matrix< DDRMat > tVeld2Ndx2 = tFI->dnNdxn( 2 );

                // get number of space time bases
                uint tNumBases = tFI->get_number_of_space_time_bases();

                mDivDivVelDof( tDofIndex )( { 0, 0 }, { 0, tNumBases - 1 } )             = tVeld2Ndx2.get_row( 0 );
                mDivDivVelDof( tDofIndex )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } ) = tVeld2Ndx2.get_row( 2 );

                mDivDivVelDof( tDofIndex )( { 1, 1 }, { 0, tNumBases - 1 } )             = tVeld2Ndx2.get_row( 2 );
                mDivDivVelDof( tDofIndex )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = tVeld2Ndx2.get_row( 1 );
            }
        }

        void
        CM_Compressible_Newtonian_Fluid::eval_dDivDivVeldu_3d( const Vector< MSI::Dof_Type > &aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the velocity FI
            Field_Interpolator *tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for ddivstrain/du
            mDivDivVelDof( tDofIndex ).set_size( 3, tFI->get_number_of_space_time_coefficients(), 0.0 );

            if ( aDofTypes( 0 ) == mDofVelocity )
            {
                // get the 2nd order derivative of the shape functions d2Ndx2
                Matrix< DDRMat > tVeld2Ndx2 = tFI->dnNdxn( 2 );

                // get number of space time bases
                uint tNumBases = tFI->get_number_of_space_time_bases();

                mDivDivVelDof( tDofIndex )( { 0, 0 }, { 0, tNumBases - 1 } )                 = tVeld2Ndx2.get_row( 0 );
                mDivDivVelDof( tDofIndex )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } )     = tVeld2Ndx2.get_row( 5 );
                mDivDivVelDof( tDofIndex )( { 0, 0 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tVeld2Ndx2.get_row( 4 );

                mDivDivVelDof( tDofIndex )( { 1, 1 }, { 0, tNumBases - 1 } )                 = tVeld2Ndx2.get_row( 5 );
                mDivDivVelDof( tDofIndex )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } )     = tVeld2Ndx2.get_row( 1 );
                mDivDivVelDof( tDofIndex )( { 1, 1 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tVeld2Ndx2.get_row( 3 );

                mDivDivVelDof( tDofIndex )( { 2, 2 }, { 0, tNumBases - 1 } )                 = tVeld2Ndx2.get_row( 4 );
                mDivDivVelDof( tDofIndex )( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } )     = tVeld2Ndx2.get_row( 3 );
                mDivDivVelDof( tDofIndex )( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tVeld2Ndx2.get_row( 2 );
            }
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::dDivDivVeldu(
                const Vector< MSI::Dof_Type > &aDofTypes )
        {
            // if aDofType is not an active dof type for the CM
            MORIS_ERROR( this->check_dof_dependency( aDofTypes ),
                    "CM_Compressible_Newtonian_Fluid::dDivDivVeldu - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mDivDivVelDofEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dDivDivVeldu( aDofTypes );

                // set bool for evaluation
                mDivDivVelDofEval( tDofIndex ) = false;
            }

            // return the dof deriv value
            return mDivDivVelDof( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Compressible_Newtonian_Fluid::eval_dNveldt()
        {
            // get the residual dof type FI (here velocity)
            Field_Interpolator *tFIVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity );

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

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::dNveldt()
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
        CM_Compressible_Newtonian_Fluid::eval_velocitymatrix_2d()
        {
            Matrix< DDRMat > tU = mFIManager->get_field_interpolators_for_type( mDofVelocity )->val();

            mVelocityMatrix = {
                { tU( 0 ), 0.0, tU( 1 ) },
                { 0.0, tU( 1 ), tU( 0 ) }
            };
        }

        void
        CM_Compressible_Newtonian_Fluid::eval_velocitymatrix_3d()
        {
            Matrix< DDRMat > tU = mFIManager->get_field_interpolators_for_type( mDofVelocity )->val();

            mVelocityMatrix = {
                { tU( 0 ), 0.0, 0.0, 0.0, tU( 2 ), tU( 1 ) },
                { 0.0, tU( 1 ), 0.0, tU( 2 ), 0.0, tU( 0 ) },
                { 0.0, 0.0, tU( 2 ), tU( 1 ), tU( 0 ), 0.0 }
            };
        }

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::velocityMatrix()
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
        CM_Compressible_Newtonian_Fluid::unfold_2d(
                const Matrix< DDRMat > &aFlattenedTensor,
                Matrix< DDRMat >       &aExpandedTensor )
        {
            aExpandedTensor = {
                { aFlattenedTensor( 0 ), aFlattenedTensor( 2 ) },
                { aFlattenedTensor( 2 ), aFlattenedTensor( 1 ) }
            };
        }

        void
        CM_Compressible_Newtonian_Fluid::unfold_3d(
                const Matrix< DDRMat > &aFlattenedTensor,
                Matrix< DDRMat >       &aExpandedTensor )
        {
            aExpandedTensor = {
                { aFlattenedTensor( 0 ), aFlattenedTensor( 5 ), aFlattenedTensor( 4 ) },
                { aFlattenedTensor( 5 ), aFlattenedTensor( 1 ), aFlattenedTensor( 3 ) },
                { aFlattenedTensor( 4 ), aFlattenedTensor( 3 ), aFlattenedTensor( 2 ) }
            };
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Compressible_Newtonian_Fluid::flatten_normal_2d(
                const Matrix< DDRMat > &aNormal,
                Matrix< DDRMat >       &aFlatNormal )
        {
            aFlatNormal.set_size( 2, 3, 0.0 );
            aFlatNormal( 0, 0 ) = aNormal( 0, 0 );
            aFlatNormal( 0, 2 ) = aNormal( 1, 0 );
            aFlatNormal( 1, 1 ) = aNormal( 1, 0 );
            aFlatNormal( 1, 2 ) = aNormal( 0, 0 );
        }

        void
        CM_Compressible_Newtonian_Fluid::flatten_normal_3d(
                const Matrix< DDRMat > &aNormal,
                Matrix< DDRMat >       &aFlatNormal )
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

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > &
        CM_Compressible_Newtonian_Fluid::MultipMat()
        {
            // build multiplication matrix
            // for 2D
            if ( mFIManager->get_field_interpolators_for_type( mDofVelocity )->get_number_of_fields() == 2 )
            {
                return mMultipMat2D;
            }
            // for 3D
            else
            {
                return mMultipMat3D;
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
