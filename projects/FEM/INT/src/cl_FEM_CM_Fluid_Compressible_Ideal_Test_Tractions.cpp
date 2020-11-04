/*
 * cl_FEM_CM_Fluid_Compressible_Ideal_Test_Tractions.cpp
 *
 *  Created on: Oct 27, 2020
 *  Author: wunsch
 */

#include "cl_FEM_CM_Fluid_Compressible_Ideal.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "op_minus.hpp"

namespace moris
{
    namespace fem
    {
        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat > & CM_Fluid_Compressible_Ideal::testTraction(
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
                            "CM_Fluid_Compressible_Ideal::testTraction - unknown CM function type for test-traction." );
                    return mTestTraction( 0 );
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Compressible_Ideal::eval_thermal_testTraction(
                const Matrix< DDRMat >             & aNormal,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // test traction is transpose of the traction dof derivative
            mThermalTestTraction( tTestDofIndex ) = trans( aNormal ) * this->dFluxdDOF( aTestDofTypes, CM_Function_Type::THERMAL );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat > & CM_Fluid_Compressible_Ideal::thermal_testTraction(
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

        void CM_Fluid_Compressible_Ideal::eval_mechanical_testTraction(
                const Matrix< DDRMat >             & aNormal,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            MORIS_ERROR( false , "CM_Fluid_Compressible_Ideal::eval_mechanical_testTraction - not implemented, yet." );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat > & CM_Fluid_Compressible_Ideal::mechanical_testTraction(
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

        const Matrix< DDRMat > & CM_Fluid_Compressible_Ideal::dTestTractiondDOF(
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
                            "CM_Fluid_Compressible_Ideal::dTestTractiondDOF - unknown CM function type for test-traction dof derivative." );
                    return mdTestTractiondDof( 0 )( 0 );
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Compressible_Ideal::eval_thermal_dTestTractiondDOF(
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
            std::shared_ptr< Property > tPropThermalConductivity = get_property( "ThermalConductivity" );

            // get temperature FI
            Field_Interpolator * tFITemp = mFIManager->get_field_interpolators_for_type( mDofTemperature );

            // if direct dependency on the dof type
            if( aDofTypes( 0 ) == mDofDensity )
            {
                if( tPropThermalConductivity->check_dof_dependency( aTestDofTypes ) )
                {
                    // add contribution
                    mdThermalTestTractiondDof( tTestDofIndex )( tDofIndex ) -=
                            tPropThermalConductivity->dPropdDOF( aTestDofTypes ) *
                            aJump( 0, 0 ) * trans( aNormal ) * tFITemp->dnNdxn( 1 );
                }
            }

            // if direct dependency on the dof type
            if( aDofTypes( 0 ) == mDofTemperature )
            {
                if( tPropThermalConductivity->check_dof_dependency( aTestDofTypes ) )
                {
                    // add contribution
                    mdThermalTestTractiondDof( tTestDofIndex )( tDofIndex ) -=
                            tFITemp->dnNdxn( 1 ) * aNormal * aJump( 0, 0 ) *
                            tPropThermalConductivity->dPropdDOF( aTestDofTypes );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat > & CM_Fluid_Compressible_Ideal::thermal_dTestTractiondDOF(
                const moris::Cell< MSI::Dof_Type > & aDofType,
                const Matrix< DDRMat >             & aNormal,
                const Matrix< DDRMat >             & aJump,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // if aDofType is not an active dof type for the property
            MORIS_ERROR( this->check_dof_dependency( aDofType ),
                    "CM_Fluid_Compressible_Ideal::thermal_dTestTractiondDOF - no dependency on this dof type." );

            // get the test dof index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mdThermalTestTractiondDofEval( tTestDofIndex )( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_thermal_dTestTractiondDOF( aDofType, aNormal, aJump, aTestDofTypes );

                // set bool for evaluation
                mdThermalTestTractiondDofEval( tTestDofIndex )( tDofIndex ) = false;
            }

            // return the derivative
            return mdThermalTestTractiondDof( tTestDofIndex )( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------

        void CM_Fluid_Compressible_Ideal::eval_mechanical_dTestTractiondDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                const Matrix< DDRMat >             & aNormal,
                const Matrix< DDRMat >             & aJump,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            MORIS_ERROR( false , "CM_Fluid_Compressible_Ideal::eval_mechanical_dTestTractiondDOF - not implemented, yet." );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat > & CM_Fluid_Compressible_Ideal::mechanical_dTestTractiondDOF(
                const moris::Cell< MSI::Dof_Type > & aDofType,
                const Matrix< DDRMat >             & aNormal,
                const Matrix< DDRMat >             & aJump,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // if aDofType is not an active dof type for the property
            MORIS_ERROR( this->check_dof_dependency( aDofType ),
                    "CM_Fluid_Compressible_Ideal::mechanical_dTestTractiondDOF - no dependency on this dof type." );

            // get the test dof index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mdMechanicalTestTractiondDofEval( tTestDofIndex )( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_mechanical_dTestTractiondDOF( aDofType, aNormal, aJump, aTestDofTypes );

                // set bool for evaluation
                mdMechanicalTestTractiondDofEval( tTestDofIndex )( tDofIndex ) = false;
            }

            // return the derivative
            return mdMechanicalTestTractiondDof( tTestDofIndex )( tDofIndex );
        }

        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
