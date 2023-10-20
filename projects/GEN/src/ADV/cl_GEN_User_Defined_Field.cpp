/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_User_Defined_Field.cpp
 *
 */

#include "cl_GEN_User_Defined_Field.hpp"

namespace moris::ge
{
    //--------------------------------------------------------------------------------------------------------------

    User_Defined_Field::User_Defined_Field(
            Matrix< DDRMat > aConstants,
            Field_Function   aFieldFunction )
            : Field_Analytic( aConstants )
            , mFieldVariables( aConstants.length() )
    {
        this->set_user_defined_functions( aFieldFunction, nullptr );
        this->import_advs( nullptr );
    }

    //--------------------------------------------------------------------------------------------------------------

    void User_Defined_Field::import_advs( sol::Dist_Vector* aOwnedADVs )
    {
        // Set new variables for the user-defined function calls
        for ( uint iVariableIndex = 0; iVariableIndex < mFieldVariables.size(); iVariableIndex++ )
        {
            mFieldVariables( iVariableIndex ) = mADVManager.get_variable( iVariableIndex );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    User_Defined_Field::get_field_value( const Matrix< DDRMat >& aCoordinates )
    {
        return this->get_field_value_user_defined( aCoordinates, mFieldVariables );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    User_Defined_Field::get_dfield_dadvs( const Matrix< DDRMat >& aCoordinates )
    {
        this->get_dfield_dadvs_user_defined( aCoordinates, mFieldVariables, mSensitivities );
        return mSensitivities;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    User_Defined_Field::get_dfield_dcoordinates(
            const Matrix< DDRMat >& aCoordinates,
            Matrix< DDRMat >&       aSensitivities )
    {
        MORIS_ASSERT( aCoordinates.numel() == aSensitivities.numel(),
                "User_Defined_Field::get_dfield_dcoordinates - inconsistent dimensions" );

        const real tPertubation = 1e-8;

        Matrix< DDRMat > tCoordinates( aCoordinates );

        for ( uint idim = 0; idim < aCoordinates.numel(); ++idim )
        {
            tCoordinates( idim ) += tPertubation;
            real tFieldValue = this->get_field_value_user_defined( tCoordinates, mFieldVariables );

            tCoordinates( idim ) -= 2.0 * tPertubation;
            tFieldValue -= this->get_field_value_user_defined( tCoordinates, mFieldVariables );

            aSensitivities( idim ) = tFieldValue / 2.0 / tPertubation;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    User_Defined_Field::set_user_defined_functions(
            Field_Function       aFieldFunction,
            Sensitivity_Function aSensitivityFunction )
    {
        // Set field evaluation function
        get_field_value_user_defined = aFieldFunction;

        // Check field evaluation function
        MORIS_ERROR( get_field_value_user_defined, "No field evaluation function was provided to a user-defined field." );
        MORIS_ASSERT( std::isfinite( this->get_field_value_user_defined( { { 0.0, 0.0, 0.0 } }, mFieldVariables ) ),
                "There is an error in a user-defined geometry field (field evaluates to nan/infinity)." );

        // Set sensitivity evaluation function
        if ( aSensitivityFunction == nullptr )
        {
            // Sensitivity function was not provided
            get_dfield_dadvs_user_defined = &( User_Defined_Field::no_sensitivities );
        }
        else
        {
            // Sensitivity function was provided
            get_dfield_dadvs_user_defined = aSensitivityFunction;

            // Check sensitivity function
            this->get_dfield_dadvs_user_defined( { { 0.0, 0.0, 0.0 } }, mFieldVariables, mSensitivities );

            // Check for row vector
            MORIS_ERROR( mSensitivities.n_rows() == 1,
                    "A user-defined geometry must provide a row vector for sensitivities." );

            // Check for size
            MORIS_ERROR( mSensitivities.n_cols() == mFieldVariables.size(),
                    "A user-defined geometry must have a sensitivity vector with a length equal to the total "
                    "number of geometry variables (ADVs + constants). sensitivities: %zu   geom variables: %zu \n",
                    mSensitivities.n_cols(),
                    mFieldVariables.size() );

            // Check for values not nan/infinity
            for ( uint tSensitivityIndex = 0; tSensitivityIndex < mSensitivities.n_cols(); tSensitivityIndex++ )
            {
                MORIS_ASSERT( std::isfinite( mSensitivities( tSensitivityIndex ) ),
                        "There is an error in a user-defined geometry sensitivity (evaluates to nan/infinity)." );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    User_Defined_Field::no_sensitivities(
            const Matrix< DDRMat >& aCoordinates,
            const Cell< real >&    aParameters,
            Matrix< DDRMat >&       aSensitivities )
    {
        MORIS_ERROR( false,
                "A sensitivity evaluation function was not provided to a user-defined geometry. "
                "Please make sure that you provide this function, or that sensitivities are not required." );
    }

    //--------------------------------------------------------------------------------------------------------------

}
