/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Fluid_Compressible_Van_der_Waals_Test_Tractions.cpp
 *
 */

#include "cl_FEM_CM_Fluid_Compressible_Van_der_Waals.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_dot.hpp"
#include "op_minus.hpp"

namespace moris::fem
{
    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Van_der_Waals::testTraction(
            const Matrix< DDRMat >&        aNormal,
            const Vector< MSI::Dof_Type >& aTestDofTypes,
            enum CM_Function_Type          aCMFunctionType )
    {
        switch ( aCMFunctionType )
        {
            case CM_Function_Type::THERMAL:
                return this->thermal_testTraction( aNormal, aTestDofTypes );

            case CM_Function_Type::MECHANICAL:
                return this->mechanical_testTraction( aNormal, aTestDofTypes );

                // unknown CM function type
            default:
                MORIS_ERROR( false,
                        "CM_Fluid_Compressible_Van_der_Waals::testTraction - unknown CM function type for test-traction." );
                return mTestTraction( 0 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Van_der_Waals::eval_thermal_testTraction(
            const Matrix< DDRMat >&        aNormal,
            const Vector< MSI::Dof_Type >& aTestDofTypes )
    {
        // get test dof type index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // test traction is transpose of the traction dof derivative
        mThermalTestTraction( tTestDofIndex ) = trans( aNormal ) * this->dFluxdDOF( aTestDofTypes, CM_Function_Type::THERMAL );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Van_der_Waals::thermal_testTraction(
            const Matrix< DDRMat >&        aNormal,
            const Vector< MSI::Dof_Type >& aTestDofTypes )
    {
        // get test dof type index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // if the test traction was not evaluated
        if ( mThermalTestTractionEval( tTestDofIndex ) )
        {
            // evaluate the test traction
            this->eval_thermal_testTraction( aNormal, aTestDofTypes );

            // set bool for evaluation
            mThermalTestTractionEval( tTestDofIndex ) = false;
        }
        // return the test traction value
        return mThermalTestTraction( tTestDofIndex );
    }

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Van_der_Waals::eval_mechanical_testTraction(
            const Matrix< DDRMat >&        aNormal,
            const Vector< MSI::Dof_Type >& aTestDofTypes )
    {
        // get test dof type index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // flatten the normal
        Matrix< DDRMat > tFlatNormal;
        this->flatten_normal( aNormal, tFlatNormal );

        // test traction is transpose of the traction dof derivative
        mMechanicalTestTraction( tTestDofIndex ) = tFlatNormal * this->dFluxdDOF( aTestDofTypes, CM_Function_Type::MECHANICAL );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Van_der_Waals::mechanical_testTraction(
            const Matrix< DDRMat >&        aNormal,
            const Vector< MSI::Dof_Type >& aTestDofTypes )
    {
        // get test dof type index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // if the test traction was not evaluated
        if ( mMechanicalTestTractionEval( tTestDofIndex ) )
        {
            // evaluate the test traction
            this->eval_mechanical_testTraction( aNormal, aTestDofTypes );

            // set bool for evaluation
            mMechanicalTestTractionEval( tTestDofIndex ) = false;
        }
        // return the test traction value
        return mMechanicalTestTraction( tTestDofIndex );
    }

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Van_der_Waals::dTestTractiondDOF(
            const Vector< MSI::Dof_Type >& aDofType,
            const Matrix< DDRMat >&        aNormal,
            const Matrix< DDRMat >&        aJump,
            const Vector< MSI::Dof_Type >& aTestDofTypes,
            enum CM_Function_Type          aCMFunctionType )
    {
        switch ( aCMFunctionType )
        {
            case CM_Function_Type::THERMAL:
                return this->thermal_dTestTractiondDOF( aDofType, aNormal, aJump, aTestDofTypes );

            case CM_Function_Type::MECHANICAL:
                return this->mechanical_dTestTractiondDOF( aDofType, aNormal, aJump, aTestDofTypes );

                // unknown CM function type
            default:
                MORIS_ERROR( false,
                        "CM_Fluid_Compressible_Van_der_Waals::dTestTractiondDOF - unknown CM function type for test-traction dof derivative." );
                return mdTestTractiondDof( 0 )( 0 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Van_der_Waals::eval_thermal_dTestTractiondDOF(
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Matrix< DDRMat >&        aNormal,
            const Matrix< DDRMat >&        aJump,
            const Vector< MSI::Dof_Type >& aTestDofTypes )
    {
        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get test dof type index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // initialize the dTestTractiondDof
        mdThermalTestTractiondDof( tTestDofIndex )( tDofIndex ).set_size(                                                       //
                mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) )->get_number_of_space_time_coefficients(),    //
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(),
                0.0 );

        // get the conductivity property
        std::shared_ptr< Property > tPropThermalConductivity = get_property( "ThermalConductivity" );

        // get temperature FI
        Field_Interpolator* tFITemp = mFIManager->get_field_interpolators_for_type( mDofTemperature );

        // if direct dependency on the dof type
        if ( aDofTypes( 0 ) == mDofDensity )
        {
            if ( tPropThermalConductivity->check_dof_dependency( aTestDofTypes ) )
            {
                // add contribution
                mdThermalTestTractiondDof( tTestDofIndex )( tDofIndex ) -=
                        tPropThermalConductivity->dPropdDOF( aTestDofTypes ) *    //
                        aJump( 0, 0 ) * trans( aNormal ) * tFITemp->dnNdxn( 1 );
            }
        }

        // if direct dependency on the dof type
        if ( aDofTypes( 0 ) == mDofTemperature )
        {
            if ( tPropThermalConductivity->check_dof_dependency( aTestDofTypes ) )
            {
                // add contribution
                mdThermalTestTractiondDof( tTestDofIndex )( tDofIndex ) -=
                        tFITemp->dnNdxn( 1 ) * aNormal * aJump( 0, 0 ) *    //
                        tPropThermalConductivity->dPropdDOF( aTestDofTypes );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Van_der_Waals::thermal_dTestTractiondDOF(
            const Vector< MSI::Dof_Type >& aDofType,
            const Matrix< DDRMat >&        aNormal,
            const Matrix< DDRMat >&        aJump,
            const Vector< MSI::Dof_Type >& aTestDofTypes )
    {
        // if aDofType is not an active dof type for the property
        MORIS_ERROR( this->check_dof_dependency( aDofType ),
                "CM_Fluid_Compressible_Van_der_Waals::thermal_dTestTractiondDOF - no dependency on this dof type." );

        // get the test dof index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdThermalTestTractiondDofEval( tTestDofIndex, tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_thermal_dTestTractiondDOF( aDofType, aNormal, aJump, aTestDofTypes );

            // set bool for evaluation
            mdThermalTestTractiondDofEval( tTestDofIndex, tDofIndex ) = false;
        }

        // return the derivative
        return mdThermalTestTractiondDof( tTestDofIndex )( tDofIndex );
    }

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Fluid_Compressible_Van_der_Waals::eval_mechanical_dTestTractiondDOF(
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Matrix< DDRMat >&        aNormal,
            const Matrix< DDRMat >&        aJump,
            const Vector< MSI::Dof_Type >& aTestDofTypes )
    {
        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get test dof type index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // initialize the dTestTractiondDof
        mdMechanicalTestTractiondDof( tTestDofIndex )( tDofIndex ).set_size(                                                    //
                mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) )->get_number_of_space_time_coefficients(),    //
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(),
                0.0 );

        // get the viscosity
        const std::shared_ptr< Property > tPropDynamicViscosity       = get_property( "DynamicViscosity" );
        const std::shared_ptr< Property > tPropCapillarityCoefficient = get_property( "CapillarityCoefficient" );

        real tCapillarityCoefficient = tPropCapillarityCoefficient->val()( 0 );
        real tSpecificGasConstant    = get_property( "SpecificGasConstant" )->val()( 0 );
        real tFirstVdWconstant       = get_property( "FirstVdWconstant" )->val()( 0 );
        real tSecondVdWconstant      = get_property( "SecondVdWconstant" )->val()( 0 );

        // get FIs
        Field_Interpolator* tFIDensity = mFIManager->get_field_interpolators_for_type( mDofDensity );
        // Field_Interpolator * tFIVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity );
        Field_Interpolator* tFITemp = mFIManager->get_field_interpolators_for_type( mDofTemperature );

        // flatten the normal
        Matrix< DDRMat > tFlatNormal;
        this->flatten_normal( aNormal, tFlatNormal );

        //-------------------------------------------------------------------------------
        // RIGHT DOF: Density
        if ( aDofTypes( 0 ) == mDofDensity )
        {
            // LEFT DOF: Density
            if ( aTestDofTypes( 0 ) == mDofDensity )
            {
                // Viscous stress
                // NOTHING!

                // Pressure
                real tFraction = 2.0 * std::pow( tSecondVdWconstant, 2.0 ) * tSpecificGasConstant * tFITemp->val()( 0 ) /    //
                                 std::pow( tSecondVdWconstant - tFIDensity->val()( 0 ), 3.0 );

                mdMechanicalTestTractiondDof( tTestDofIndex )( tDofIndex ) -=
                        tFIDensity->N_trans() * dot( aJump, aNormal ) * ( tFraction - 2.0 * tFirstVdWconstant ) * tFIDensity->N();

                // Korteweg stress
                mdMechanicalTestTractiondDof( tTestDofIndex )( tDofIndex ) +=
                        tCapillarityCoefficient * dot( aJump, aNormal ) * ( trans( tFIDensity->dnNdxn( 1 ) ) * tFIDensity->dnNdxn( 1 ) +    //
                                                                            trans( this->laplaceDensityDof() ) * tFIDensity->N() + tFIDensity->N_trans() * this->laplaceDensityDof() )
                        - trans( tFIDensity->dnNdxn( 1 ) ) * tCapillarityCoefficient * ( aNormal * trans( aJump ) + aJump * trans( aNormal ) ) * tFIDensity->dnNdxn( 1 );

                // indirect contribution
                if ( tPropCapillarityCoefficient->check_dof_dependency( aDofTypes ) )
                {
                    // compute density derivative of korteweg stress
                    Matrix< DDRMat > tdKortewegStressdDensity =
                            tCapillarityCoefficient * mFlatIdentity * tFIDensity->val() * this->laplaceDensityDof() +                //
                            tCapillarityCoefficient * mFlatIdentity * this->laplaceDensity() * tFIDensity->N() +                     //
                            tCapillarityCoefficient * mFlatIdentity * trans( tFIDensity->gradx( 1 ) ) * tFIDensity->dnNdxn( 1 ) +    //
                            ( -1.0 ) * tCapillarityCoefficient * this->densityStrainDof() +                                          //
                            mFlatIdentity * ( tFIDensity->val()( 0 ) * this->laplaceDensity() +                                      //
                                              0.5 * ( trans( tFIDensity->gradx( 1 ) ) * tFIDensity->gradx( 1 ) ) )                   //
                                    * tPropCapillarityCoefficient->dPropdDOF( aDofTypes )                                            //
                            - this->densityStrain() * tPropCapillarityCoefficient->dPropdDOF( aDofTypes );

                    // add contribution
                    mdMechanicalTestTractiondDof( tTestDofIndex )( tDofIndex ) +=
                            trans( tPropDynamicViscosity->dPropdDOF( aTestDofTypes ) ) * ( 1.0 / tCapillarityCoefficient ) *    //
                            trans( aJump ) * tFlatNormal * tdKortewegStressdDensity;
                }
            }

            // LEFT DOF: Velocity
            if ( aTestDofTypes( 0 ) == mDofVelocity )
            {
                // Viscous stress
                if ( tPropDynamicViscosity->check_dof_dependency( aDofTypes ) )
                {
                    mdMechanicalTestTractiondDof( tTestDofIndex )( tDofIndex ) +=
                            trans( tFlatNormal * this->dStraindDOF( aTestDofTypes ) ) * 2.0 * aJump * tPropDynamicViscosity->dPropdDOF( aDofTypes );
                }

                // Pressure
                // NOTHING!

                // Korteweg stress
                // NOTHING!
            }

            // LEFT DOF: Temperature
            if ( aTestDofTypes( 0 ) == mDofTemperature )
            {
                // Viscous stress
                // NOTHING!

                // Pressure
                real tFraction = std::pow( tSecondVdWconstant, 2.0 ) * tSpecificGasConstant /    //
                                 std::pow( tSecondVdWconstant - tFIDensity->val()( 0 ), 2.0 );

                mdMechanicalTestTractiondDof( tTestDofIndex )( tDofIndex ) -=
                        tFITemp->N_trans() * tFraction * dot( aJump, aNormal ) * tFIDensity->N();

                // Korteweg stress
                if ( tPropCapillarityCoefficient->check_dof_dependency( aTestDofTypes ) )
                {
                    // compute density derivative of korteweg stress
                    Matrix< DDRMat > tdKortewegStressdDensity =
                            tCapillarityCoefficient * mFlatIdentity * tFIDensity->val() * this->laplaceDensityDof() +                //
                            tCapillarityCoefficient * mFlatIdentity * this->laplaceDensity() * tFIDensity->N() +                     //
                            tCapillarityCoefficient * mFlatIdentity * trans( tFIDensity->gradx( 1 ) ) * tFIDensity->dnNdxn( 1 ) +    //
                            ( -1.0 ) * tCapillarityCoefficient * this->densityStrainDof();

                    if ( tPropCapillarityCoefficient->check_dof_dependency( aDofTypes ) )
                    {
                        tdKortewegStressdDensity +=
                                mFlatIdentity * ( tFIDensity->val()( 0 ) * this->laplaceDensity() +                       //
                                                  0.5 * ( trans( tFIDensity->gradx( 1 ) ) * tFIDensity->gradx( 1 ) ) )    //
                                        * tPropCapillarityCoefficient->dPropdDOF( aDofTypes )                             //
                                - this->densityStrain() * tPropCapillarityCoefficient->dPropdDOF( aDofTypes );
                    }

                    // add contribution
                    mdMechanicalTestTractiondDof( tTestDofIndex )( tDofIndex ) +=
                            trans( tPropDynamicViscosity->dPropdDOF( aTestDofTypes ) ) * ( 1.0 / tCapillarityCoefficient ) *    //
                            trans( aJump ) * tFlatNormal * tdKortewegStressdDensity;
                }
            }
        }

        //-------------------------------------------------------------------------------
        // RIGHT DOF: Velocity
        if ( aDofTypes( 0 ) == mDofVelocity )
        {
            // LEFT DOF: Density
            if ( aTestDofTypes( 0 ) == mDofDensity )
            {
                // Viscous stress
                if ( tPropDynamicViscosity->check_dof_dependency( aTestDofTypes ) )
                {
                    mdMechanicalTestTractiondDof( tTestDofIndex )( tDofIndex ) +=
                            trans( tPropDynamicViscosity->dPropdDOF( aTestDofTypes ) ) *    //
                            2.0 * trans( aJump ) * tFlatNormal * this->dStraindDOF( aDofTypes );
                }

                // Pressure
                // NOTHING!

                // Korteweg stress
                // NOTHING!
            }

            // LEFT DOF: Velocity
            if ( aTestDofTypes( 0 ) == mDofVelocity )
            {
                // Viscous stress
                if ( tPropDynamicViscosity->check_dof_dependency( aDofTypes ) )
                {
                    mdMechanicalTestTractiondDof( tTestDofIndex )( tDofIndex ) +=
                            trans( tFlatNormal * this->dStraindDOF( aTestDofTypes ) ) *    //
                                    2.0 * aJump * tFlatNormal * this->dStraindDOF( aDofTypes )
                            + trans( tPropDynamicViscosity->dPropdDOF( aTestDofTypes ) ) *    //
                                      2.0 * trans( aJump ) * tFlatNormal * this->dStraindDOF( aDofTypes );
                }

                // Pressure
                // NOTHING!

                // Korteweg stress
                // NOTHING!
            }

            // LEFT DOF: Temperature
            if ( aTestDofTypes( 0 ) == mDofTemperature )
            {
                // Viscous stress
                if ( tPropDynamicViscosity->check_dof_dependency( aTestDofTypes ) )
                {
                    mdMechanicalTestTractiondDof( tTestDofIndex )( tDofIndex ) +=
                            trans( tPropDynamicViscosity->dPropdDOF( tTestDofIndex ) ) *    //
                            2.0 * trans( aJump ) * tFlatNormal * this->dStraindDOF( tDofIndex );
                }

                // Pressure
                // NOTHING!

                // Korteweg stress
                // NOTHING!
            }
        }

        //-------------------------------------------------------------------------------
        // RIGHT DOF: Temperature
        if ( aDofTypes( 0 ) == mDofTemperature )
        {
            // LEFT DOF: Density
            if ( aTestDofTypes( 0 ) == mDofDensity )
            {
                // Viscous stress
                // NOTHING!

                // Pressure
                real tFraction = std::pow( tSecondVdWconstant, 2.0 ) * tSpecificGasConstant /    //
                                 std::pow( tSecondVdWconstant - tFIDensity->val()( 0 ), 2.0 );

                mdMechanicalTestTractiondDof( tTestDofIndex )( tDofIndex ) -=
                        tFIDensity->N_trans() * tFraction * dot( aJump, aNormal ) * tFITemp->N();

                // Korteweg stress
                if ( tPropCapillarityCoefficient->check_dof_dependency( aDofTypes ) )
                {
                    // compute density derivative of korteweg stress
                    Matrix< DDRMat > tdKortewegStressdDensity =
                            tCapillarityCoefficient * mFlatIdentity * tFIDensity->val() * this->laplaceDensityDof() +                //
                            tCapillarityCoefficient * mFlatIdentity * this->laplaceDensity() * tFIDensity->N() +                     //
                            tCapillarityCoefficient * mFlatIdentity * trans( tFIDensity->gradx( 1 ) ) * tFIDensity->dnNdxn( 1 ) +    //
                            ( -1.0 ) * tCapillarityCoefficient * this->densityStrainDof();

                    if ( tPropCapillarityCoefficient->check_dof_dependency( aTestDofTypes ) )
                    {
                        tdKortewegStressdDensity +=
                                mFlatIdentity * ( tFIDensity->val()( 0 ) * this->laplaceDensity() +                       //
                                                  0.5 * ( trans( tFIDensity->gradx( 1 ) ) * tFIDensity->gradx( 1 ) ) )    //
                                        * tPropCapillarityCoefficient->dPropdDOF( aTestDofTypes )                         //
                                - this->densityStrain() * tPropCapillarityCoefficient->dPropdDOF( aTestDofTypes );
                    }

                    // add contribution
                    mdMechanicalTestTractiondDof( tTestDofIndex )( tDofIndex ) +=
                            trans( tdKortewegStressdDensity ) * trans( tFlatNormal ) *    //
                            aJump * ( 1.0 / tCapillarityCoefficient ) * tPropDynamicViscosity->dPropdDOF( aDofTypes );
                }
            }

            // LEFT DOF: Velocity
            if ( aTestDofTypes( 0 ) == mDofVelocity )
            {
                // Viscous stress
                if ( tPropDynamicViscosity->check_dof_dependency( aTestDofTypes ) )
                {
                    mdMechanicalTestTractiondDof( tTestDofIndex )( tDofIndex ) =
                            trans( tFlatNormal * this->dStraindDOF( aTestDofTypes ) ) *    //
                            2.0 * aJump * tPropDynamicViscosity->dPropdDOF( tTestDofIndex );
                }

                // Pressure
                // nothing for now

                // Korteweg stress
                // NOTHING!
            }

            // LEFT DOF: Temperature
            if ( aTestDofTypes( 0 ) == mDofTemperature )
            {
                // Viscous stress
                // NOTHING!

                // Pressure
                // nothing for now

                // Korteweg stress
                // NOTHING!
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    CM_Fluid_Compressible_Van_der_Waals::mechanical_dTestTractiondDOF(
            const Vector< MSI::Dof_Type >& aDofType,
            const Matrix< DDRMat >&        aNormal,
            const Matrix< DDRMat >&        aJump,
            const Vector< MSI::Dof_Type >& aTestDofTypes )
    {
        // if aDofType is not an active dof type for the property
        MORIS_ERROR( this->check_dof_dependency( aDofType ),
                "CM_Fluid_Compressible_Van_der_Waals::mechanical_dTestTractiondDOF - no dependency on this dof type." );

        // get the test dof index
        uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdMechanicalTestTractiondDofEval( tTestDofIndex, tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_mechanical_dTestTractiondDOF( aDofType, aNormal, aJump, aTestDofTypes );

            // set bool for evaluation
            mdMechanicalTestTractiondDofEval( tTestDofIndex, tDofIndex ) = false;
        }

        // return the derivative
        return mdMechanicalTestTractiondDof( tTestDofIndex )( tDofIndex );
    }

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::fem
