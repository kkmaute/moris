/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Material_Model_Temperature_Functions.cpp
 *
 */

#include "cl_FEM_Material_Model.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------
    // RETURN FUNCTIONS FOR TEMPERATURE (SECOND EQUATION OF STATE)
    //------------------------------------------------------------------------------

    // if thermodynamic state variable is dependent compute and retrieve value from storage
    const Matrix< DDRMat > &Material_Model::temperature_dep()
    {
        // if the temperature was not evaluated
        if ( mTemperatureEval )
        {
            // evaluate the temperature
            this->eval_temperature();

            // set bool for evaluation
            mTemperatureEval = false;
        }
        // return the temperature value
        return mTemperature;
    }

    // trivial operation: get value from FI
    const Matrix< DDRMat > &Material_Model::temperature_triv()
    {
        // return the temperature value
        return mFIManager->get_field_interpolators_for_type( mDofTemperature )->val();
    }

    //------------------------------------------------------------------------------

    // if thermodynamic state variable is dependent compute and retrieve value from storage
    const Matrix< DDRMat > &Material_Model::TemperatureDot_dep()
    {
        // if the flux was not evaluated
        if ( mTemperatureDotEval )
        {
            // evaluate the flux
            this->eval_TemperatureDot();

            // set bool for evaluation
            mTemperatureDotEval = false;
        }
        // return the flux value
        return mTemperatureDot;
    }

    // trivial operation: get value from FI
    const Matrix< DDRMat > &Material_Model::TemperatureDot_triv()
    {
        // return the temperature rate of change
        return mFIManager->get_field_interpolators_for_type( mDofTemperature )->gradt( 1 );
    }

    //------------------------------------------------------------------------------

    // if thermodynamic state variable is dependent compute and retrieve value from storage
    const Matrix< DDRMat > &Material_Model::dnTemperaturedxn_dep( uint aOrder )
    {
        switch ( aOrder )
        {
            case 1:    // first derivative
            {
                if ( mdTemperaturedxEval )
                {
                    // evaluate the flux
                    this->eval_dTemperaturedx();

                    // set bool for evaluation
                    mdTemperaturedxEval = false;
                }

                // return the flux value
                return mdTemperaturedx;
            }

            case 2:    // second derivative
            {

                if ( md2Temperaturedx2Eval )
                {
                    // evaluate the flux
                    this->eval_d2Temperaturedx2();

                    // set bool for evaluation
                    md2Temperaturedx2Eval = false;
                }

                // return the flux value
                return md2Temperaturedx2;
            }

            default:
            {
                MORIS_ERROR( false, "Material_Model::dnTemperaturedxn - aOrder unknown, only 1 and 2 supported." );
                return mdTemperaturedx;
            }
        }
    }

    // trivial operation: get values from FI
    const Matrix< DDRMat > &Material_Model::dnTemperaturedxn_triv( uint aOrder )
    {
        // return the temperature rate of change
        return mFIManager->get_field_interpolators_for_type( mDofTemperature )->gradx( aOrder );
    }

    //-----------------------------------------------------------------------------
    //-----------------------------------------------------------------------------

    // if thermodynamic state variable is dependent compute and retrieve values from storage
    const Matrix< DDRMat > &Material_Model::TemperatureDOF_dep( const Vector< MSI::Dof_Type > &aDofType )
    {
        // if aDofType is not an active dof type for the MM
        MORIS_ASSERT(
                this->check_dof_dependency( aDofType ),
                "Material_Model::TemperatureDOF - no dependency in this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mTemperatureDofEval( tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_TemperatureDOF( aDofType );

            // set bool for evaluation
            mTemperatureDofEval( tDofIndex ) = false;
        }

        // return the derivative
        return mTemperatureDof( tDofIndex );
    }

    // trivial operation: get values from FI
    const Matrix< DDRMat > &Material_Model::TemperatureDOF_triv( const Vector< MSI::Dof_Type > &aDofType )
    {
        // check DOF deriv is wrt to own DOF-type is with
        if ( aDofType( 0 ) != mDofTemperature )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            if ( mTemperatureDofEval( tDofIndex ) )
            {
                // initialize output matrix with zeros
                mTemperatureDof( tDofIndex ).set_size( 1,    // mSpaceDim
                        mFIManager->get_field_interpolators_for_type( aDofType( 0 ) )->get_number_of_space_time_coefficients(),
                        0.0 );

                // set flag
                mTemperatureDofEval( tDofIndex ) = false;
            }

            // return zero matrix
            return mTemperatureDof( tDofIndex );
        }
        else
        {
            // return the temperature dof deriv
            return mFIManager->get_field_interpolators_for_type( mDofTemperature )->N();
        }
    }

    //-----------------------------------------------------------------------------

    // if thermodynamic state variable is dependent compute and retrieve values from storage
    const Matrix< DDRMat > &Material_Model::TemperatureDotDOF_dep( const Vector< MSI::Dof_Type > &aDofType )
    {
        // if aDofType is not an active dof type for the MM
        MORIS_ASSERT(
                this->check_dof_dependency( aDofType ),
                "Material_Model::TemperatureDotDOF_dep - no dependency in this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mTemperatureDotDofEval( tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_TemperatureDotDOF( aDofType );

            // set bool for evaluation
            mTemperatureDotDofEval( tDofIndex ) = false;
        }

        // return the derivative
        return mTemperatureDotDof( tDofIndex );
    }

    // trivial operation: get values from FI
    const Matrix< DDRMat > &Material_Model::TemperatureDotDOF_triv( const Vector< MSI::Dof_Type > &aDofType )
    {
        // check DOF deriv is wrt to own DOF-type is with
        if ( aDofType( 0 ) != mDofTemperature )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            if ( mTemperatureDotDofEval( tDofIndex ) )
            {
                // initialize output matrix with zeros
                mTemperatureDotDof( tDofIndex ).set_size( 1,    // mSpaceDim
                        mFIManager->get_field_interpolators_for_type( aDofType( 0 ) )->get_number_of_space_time_coefficients(),
                        0.0 );

                // set flag
                mTemperatureDotDofEval( tDofIndex ) = false;
            }

            // return zero matrix
            return mTemperatureDotDof( tDofIndex );
        }
        else
        {
            // return the temperature rate of change dof deriv
            return mFIManager->get_field_interpolators_for_type( mDofTemperature )->dnNdtn( 1 );
        }
    }

    //-----------------------------------------------------------------------------

    // if thermodynamic state variable is dependent compute and retrieve values from storage
    const Matrix< DDRMat > &Material_Model::dnTemperaturedxnDOF_dep( const Vector< MSI::Dof_Type > &aDofType, uint aOrder )
    {
        // if aDofType is not an active dof type for the MM
        MORIS_ASSERT(
                this->check_dof_dependency( aDofType ),
                "Material_Model::dnTemperaturedxnDOF - no dependency in this dof type." );

        // get the dof index
        uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        switch ( aOrder )
        {
            case 1:    // first derivative
            {
                // if the derivative has not been evaluated yet
                if ( mdTemperaturedxDofEval( tDofIndex ) )
                {
                    // evaluate the derivative
                    this->eval_dTemperaturedxDOF( aDofType );

                    // set bool for evaluation
                    mdTemperaturedxDofEval( tDofIndex ) = false;
                }

                // return the derivative
                return mdTemperaturedxDof( tDofIndex );
            }

            case 2:    // second derivative
            {
                // if the derivative has not been evaluated yet
                if ( md2Temperaturedx2DofEval( tDofIndex ) )
                {
                    // evaluate the derivative
                    this->eval_d2Temperaturedx2DOF( aDofType );

                    // set bool for evaluation
                    md2Temperaturedx2DofEval( tDofIndex ) = false;
                }

                // return the derivative
                return md2Temperaturedx2Dof( tDofIndex );
            }

            default:
            {
                MORIS_ERROR( false, "Material_Model::dnTemperaturedxnDOF_dep - aOrder unknown, only 1 and 2 supported." );
                return mdTemperaturedxDof( 0 );
            }
        }
    }

    // trivial operation: get values from FI
    const Matrix< DDRMat > &Material_Model::dnTemperaturedxnDOF_triv( const Vector< MSI::Dof_Type > &aDofType, uint aOrder )
    {
        // check DOF deriv is wrt to own DOF-type is with
        if ( aDofType( 0 ) != mDofTemperature )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            if ( aOrder == 1 )
            {
                if ( mdTemperaturedxDofEval( tDofIndex ) )
                {
                    // initialize output matrix with zeros
                    mdTemperaturedxDof( tDofIndex ).set_size( mSpaceDim, mFIManager->get_field_interpolators_for_type( aDofType( 0 ) )->get_number_of_space_time_coefficients(), 0.0 );

                    // set flag
                    mdTemperaturedxDofEval( tDofIndex ) = false;
                }

                // return zero matrix
                return mdTemperaturedxDof( tDofIndex );
            }
            else if ( aOrder == 2 )
            {
                if ( md2Temperaturedx2DofEval( tDofIndex ) )
                {
                    // initialize output matrix with zeros
                    md2Temperaturedx2Dof( tDofIndex ).set_size( 3 * mSpaceDim - 3, mFIManager->get_field_interpolators_for_type( aDofType( 0 ) )->get_number_of_space_time_coefficients(), 0.0 );

                    // set flag
                    md2Temperaturedx2DofEval( tDofIndex ) = false;
                }

                // return zero matrix
                return md2Temperaturedx2Dof( tDofIndex );
            }
            else
            {
                MORIS_ERROR( false, "Material_Model::dnTemperaturedxnDOF_triv - only orders 1 and 2 implemented." );
                return mdTemperaturedxDof( 0 );
            }
        }
        else
        {
            // return the temperature gradient dof deriv
            return mFIManager->get_field_interpolators_for_type( mDofTemperature )->dnNdxn( aOrder );
        }
    }

    //-----------------------------------------------------------------------------

}    // namespace moris::fem
