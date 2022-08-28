/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Compressible_Newtonian_Fluid_Test_Tractions.cpp
 *
 */

#include "cl_FEM_CM_Compressible_Newtonian_Fluid.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_dot.hpp"
#include "op_minus.hpp"

namespace moris
{
    namespace fem
    {
        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat > & CM_Compressible_Newtonian_Fluid::testTraction(
                const Matrix< DDRMat >             & aNormal,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes,
                enum CM_Function_Type aCMFunctionType )
        {
            switch( aCMFunctionType )
            {
                case CM_Function_Type::THERMAL :
                    return this->thermal_testTraction( aNormal, aTestDofTypes );

                case CM_Function_Type::MECHANICAL :
                    return this->mechanical_testTraction( aNormal, aTestDofTypes );

                    // unknown CM function type
                default :
                    MORIS_ERROR( false ,
                            "CM_Compressible_Newtonian_Fluid::testTraction - unknown CM function type for test-traction." );
                    return mTestTraction( 0 );
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------

        void CM_Compressible_Newtonian_Fluid::eval_thermal_testTraction(
                const Matrix< DDRMat >             & aNormal,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // test traction is transpose of the traction dof derivative
            mThermalTestTraction( tTestDofIndex ) = trans( aNormal ) * this->dFluxdDOF( aTestDofTypes, CM_Function_Type::THERMAL );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat > & CM_Compressible_Newtonian_Fluid::thermal_testTraction(
                const Matrix< DDRMat >             & aNormal,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // if the test traction was not evaluated
            if( mThermalTestTractionEval( tTestDofIndex ) )
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

        void CM_Compressible_Newtonian_Fluid::eval_mechanical_testTraction(
                const Matrix< DDRMat >             & aNormal,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
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

        const Matrix< DDRMat > & CM_Compressible_Newtonian_Fluid::mechanical_testTraction(
                const Matrix< DDRMat >             & aNormal,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // if the test traction was not evaluated
            if( mMechanicalTestTractionEval( tTestDofIndex ) )
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

        const Matrix< DDRMat > & CM_Compressible_Newtonian_Fluid::dTestTractiondDOF(
                const moris::Cell< MSI::Dof_Type > & aDofType,
                const Matrix< DDRMat >             & aNormal,
                const Matrix< DDRMat >             & aJump,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes,
                enum CM_Function_Type                aCMFunctionType )
        {
            switch( aCMFunctionType )
            {
                case CM_Function_Type::THERMAL :
                    return this->thermal_dTestTractiondDOF( aDofType, aNormal, aJump, aTestDofTypes );

                case CM_Function_Type::MECHANICAL :
                    return this->mechanical_dTestTractiondDOF( aDofType, aNormal, aJump, aTestDofTypes );

                    // unknown CM function type
                default :
                    MORIS_ERROR( false ,
                            "CM_Compressible_Newtonian_Fluid::dTestTractiondDOF - unknown CM function type for test-traction dof derivative." );
                    return mdTestTractiondDof( 0 )( 0 );
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------

        void CM_Compressible_Newtonian_Fluid::eval_thermal_dTestTractiondDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                const Matrix< DDRMat >             & aNormal,
                const Matrix< DDRMat >             & aJump,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // initialize the dTestTractiondDof
            mdThermalTestTractiondDof( tTestDofIndex )( tDofIndex ).set_size(
                    mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) )->get_number_of_space_time_coefficients(),
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(), 0.0 );

            // get the conductivity property
            const std::shared_ptr< Property > tPropThermalConductivity = get_property( "ThermalConductivity" );

            // get the material model
            const std::shared_ptr< Material_Model > tMM = get_material_model( "ThermodynamicMaterialModel" );

            // if direct dependency on the dof type
            if( aDofTypes( 0 ) == mDofDensity )
            {
                if( tPropThermalConductivity->check_dof_dependency( aTestDofTypes ) )
                {
                    // add contribution
                    mdThermalTestTractiondDof( tTestDofIndex )( tDofIndex ) -=
                            tPropThermalConductivity->dPropdDOF( aTestDofTypes ) *
                            aJump( 0, 0 ) * trans( aNormal ) * tMM->dnTemperaturedxnDOF( aDofTypes, 1 ) ;
                }
            }

            // if direct dependency on the dof type
            if( aDofTypes( 0 ) == mDofTemperature )
            {
                if( tPropThermalConductivity->check_dof_dependency( aTestDofTypes ) )
                {
                    // add contribution
                    mdThermalTestTractiondDof( tTestDofIndex )( tDofIndex ) -=
                            tMM->dnTemperaturedxnDOF( aDofTypes, 1 ) * aNormal * aJump( 0, 0 ) *
                            tPropThermalConductivity->dPropdDOF( aTestDofTypes );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat > & CM_Compressible_Newtonian_Fluid::thermal_dTestTractiondDOF(
                const moris::Cell< MSI::Dof_Type > & aDofType,
                const Matrix< DDRMat >             & aNormal,
                const Matrix< DDRMat >             & aJump,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // if aDofType is not an active dof type for the property
            MORIS_ERROR( this->check_dof_dependency( aDofType ),
                    "CM_Compressible_Newtonian_Fluid::thermal_dTestTractiondDOF - no dependency on this dof type." );

            // get the test dof index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mdThermalTestTractiondDofEval( tTestDofIndex, tDofIndex ) )
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

        void CM_Compressible_Newtonian_Fluid::eval_mechanical_dTestTractiondDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                const Matrix< DDRMat >             & aNormal,
                const Matrix< DDRMat >             & aJump,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // initialize the dTestTractiondDof
            mdMechanicalTestTractiondDof( tTestDofIndex )( tDofIndex ).set_size(
                    mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) )->get_number_of_space_time_coefficients(),
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(), 0.0 );

            // get the viscosity
            const std::shared_ptr< Property > tPropDynamicViscosity = get_property( "DynamicViscosity" );
            real tSpecificGasConstant = get_property( "SpecificGasConstant" )->val()( 0 );

            // get FIs
            Field_Interpolator * tFIDensity = mFIManager->get_field_interpolators_for_type( mDofDensity );
            //Field_Interpolator * tFIVelocity = mFIManager->get_field_interpolators_for_type( mDofVelocity );
            Field_Interpolator * tFITemp = mFIManager->get_field_interpolators_for_type( mDofTemperature );

            // flatten the normal
            Matrix< DDRMat > tFlatNormal;
            this->flatten_normal( aNormal, tFlatNormal );

            //-------------------------------------------------------------------------------
            // RIGHT DOF: Density
            if( aDofTypes( 0 ) == mDofDensity)
            {
                // LEFT DOF: Density
                if ( aTestDofTypes( 0 ) == mDofDensity )
                {
                    // Viscous stress
                    // NOTHING!

                    // Pressure
                    // NOTHING!
                }

                // LEFT DOF: Velocity
                if ( aTestDofTypes( 0 ) == mDofVelocity )
                {
                    // Viscous stress
                    if( tPropDynamicViscosity->check_dof_dependency( aDofTypes ) )
                    {
                        mdMechanicalTestTractiondDof( tTestDofIndex )( tDofIndex ) +=
                                trans( tFlatNormal * this->dStraindDOF( aTestDofTypes ) ) *
                                2.0 * aJump *
                                tPropDynamicViscosity->dPropdDOF( aDofTypes );
                    }

                    // Pressure
                    // NOTHING!
                }

                // LEFT DOF: Temperature
                if ( aTestDofTypes( 0 ) == mDofTemperature )
                {
                    // Viscous stress
                    // NOTHING!

                    // Pressure
                    mdMechanicalTestTractiondDof( tTestDofIndex )( tDofIndex ) -=
                            tFITemp->N_trans() *
                            tSpecificGasConstant * dot( aJump, aNormal ) *
                            tFIDensity->N();
                }
            }

            //-------------------------------------------------------------------------------
            // RIGHT DOF: Velocity
            if( aDofTypes( 0 ) == mDofVelocity )
            {
                // LEFT DOF: Density
                if ( aTestDofTypes( 0 ) == mDofDensity )
                {
                    // Viscous stress
                    if( tPropDynamicViscosity->check_dof_dependency( aTestDofTypes ) )
                    {
                        mdMechanicalTestTractiondDof( tTestDofIndex )( tDofIndex ) +=
                                trans( tPropDynamicViscosity->dPropdDOF( aTestDofTypes ) ) *
                                2.0 * trans( aJump ) *
                                tFlatNormal * this->dStraindDOF( aDofTypes );
                    }

                    // Pressure
                    // NOTHING!
                }

                // LEFT DOF: Velocity
                if ( aTestDofTypes( 0 ) == mDofVelocity )
                {
                    // Viscous stress
                    if( tPropDynamicViscosity->check_dof_dependency( aDofTypes ) )
                    {
                        mdMechanicalTestTractiondDof( tTestDofIndex )( tDofIndex ) +=
                                trans( tFlatNormal * this->dStraindDOF( aTestDofTypes ) ) *
                                2.0 * aJump *
                                tFlatNormal * this->dStraindDOF( aDofTypes ) +
                                trans( tPropDynamicViscosity->dPropdDOF( aTestDofTypes ) ) *
                                2.0 * trans( aJump ) *
                                tFlatNormal * this->dStraindDOF( aDofTypes );
                    }

                    // Pressure
                    // NOTHING!
                }

                // LEFT DOF: Temperature
                if ( aTestDofTypes( 0 ) == mDofTemperature )
                {
                    // Viscous stress
                    if( tPropDynamicViscosity->check_dof_dependency( aTestDofTypes ) )
                    {
                        mdMechanicalTestTractiondDof( tTestDofIndex )( tDofIndex ) +=
                                trans( tPropDynamicViscosity->dPropdDOF( tTestDofIndex ) ) *
                                2.0 * trans( aJump ) *
                                tFlatNormal * this->dStraindDOF( tDofIndex );
                    }

                    // Pressure
                    // NOTHING!
                }
            }

            //-------------------------------------------------------------------------------
            // RIGHT DOF: Temperature
            if( aDofTypes( 0 ) == mDofTemperature )
            {
                // LEFT DOF: Density
                if ( aTestDofTypes( 0 ) == mDofDensity )
                {
                    // Viscous stress
                    // NOTHING!

                    // Pressure
                    mdMechanicalTestTractiondDof( tTestDofIndex )( tDofIndex ) -=
                            tFIDensity->N_trans() *
                            tSpecificGasConstant * dot( aJump, aNormal ) *
                            tFITemp->N();
                }

                // LEFT DOF: Velocity
                if ( aTestDofTypes( 0 ) == mDofVelocity )
                {
                    // Viscous stress
                    if( tPropDynamicViscosity->check_dof_dependency( aTestDofTypes ) )
                    {
                        mdMechanicalTestTractiondDof( tTestDofIndex )( tDofIndex ) =
                                trans( tFlatNormal * this->dStraindDOF( aTestDofTypes ) ) *
                                2.0 * aJump *
                                tPropDynamicViscosity->dPropdDOF( tTestDofIndex );
                    }

                    // Pressure
                    // nothing for now
                }

                // LEFT DOF: Temperature
                if ( aTestDofTypes( 0 ) == mDofTemperature )
                {
                    // Viscous stress
                    // NOTHING!

                    // Pressure
                    // nothing for now
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat > & CM_Compressible_Newtonian_Fluid::mechanical_dTestTractiondDOF(
                const moris::Cell< MSI::Dof_Type > & aDofType,
                const Matrix< DDRMat >             & aNormal,
                const Matrix< DDRMat >             & aJump,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // if aDofType is not an active dof type for the property
            MORIS_ERROR( this->check_dof_dependency( aDofType ),
                    "CM_Compressible_Newtonian_Fluid::mechanical_dTestTractiondDOF - no dependency on this dof type." );

            // get the test dof index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mdMechanicalTestTractiondDofEval( tTestDofIndex, tDofIndex ) )
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

    } /* namespace fem */
} /* namespace moris */

