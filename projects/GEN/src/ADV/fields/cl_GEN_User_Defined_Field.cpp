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

namespace moris::gen
{
    //--------------------------------------------------------------------------------------------------------------

    User_Defined_Field::User_Defined_Field(
            Field_Function       aFieldFunction,
            Sensitivity_Function aSensitivityFunction,
            Sensitivity_Function aSpatialGradientSensitivity,
            const Vector< ADV >& aADVs,
            std::string          aName )
            : Field_Analytic< 0 >( aADVs, std::move( aName ) )
            , mFieldVariables( aADVs.size() )
            , get_field_value_user_defined( aFieldFunction )
            , get_dfield_dadvs_user_defined( aSensitivityFunction )
            , get_dfield_dcoordinates_select( aSpatialGradientSensitivity )
    {
        this->import_advs( nullptr );
        this->validate_user_defined_functions();
    }

    //--------------------------------------------------------------------------------------------------------------

    User_Defined_Field::User_Defined_Field(
            Field_Function       aFieldFunction,
            const Vector< ADV >& aADVs,
            std::string          aName )
            : User_Defined_Field( aFieldFunction, nullptr, nullptr, aADVs, std::move( aName ) )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void User_Defined_Field::import_advs( sol::Dist_Vector* aOwnedADVs )
    {
        // Set new variables for the user-defined function calls
        for ( uint iVariableIndex = 0; iVariableIndex < mFieldVariables.size(); iVariableIndex++ )
        {
            mFieldVariables( iVariableIndex ) = mADVHandler.get_variable( iVariableIndex );
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
        // Compute via user-defined function if provided, otherwise use finite difference
        this->get_dfield_dcoordinates_select( aCoordinates, mFieldVariables, aSensitivities );
    }

    //--------------------------------------------------------------------------------------------------------------

    void User_Defined_Field::get_dfield_dcoordinates_FD(
            const Matrix< DDRMat >& aCoordinates,
            const Vector< real >&   aParameters,
            Matrix< DDRMat >&       aSensitivities )
    {
        MORIS_ASSERT( aCoordinates.numel() == aSensitivities.numel(),
                "User_Defined_Field::get_dfield_dcoordinates - inconsistent dimensions" );

        // Finite difference perturbation size
        const real tPertubation = 1e-8;

        // Compute sensitivities using central difference
        Matrix< DDRMat > tCoordinates( aCoordinates );    // Copy coordinates to perturb
        for ( uint idim = 0; idim < aCoordinates.numel(); ++idim )
        {
            // Positive perturbation
            tCoordinates( idim ) += tPertubation;
            real tFieldValue = this->get_field_value_user_defined( tCoordinates, aParameters );

            // Negative perturbation
            tCoordinates( idim ) -= 2.0 * tPertubation;
            tFieldValue -= this->get_field_value_user_defined( tCoordinates, aParameters );

            // Reset coordinates
            tCoordinates( idim ) += tPertubation;

            // Compute sensitivity
            aSensitivities( idim ) = tFieldValue / 2.0 / tPertubation;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    User_Defined_Field::validate_user_defined_functions()
    {
        // Check field evaluation function
        MORIS_ERROR( get_field_value_user_defined, "No field evaluation function was provided to a user-defined field." );
        MORIS_ASSERT( std::isfinite( this->get_field_value_user_defined( { { 0.0, 0.0, 0.0 } }, mFieldVariables ) ),
                "There is an error in a user-defined field (field evaluates to nan/infinity)." );

        // Check if sensitivity function was provided
        if ( get_dfield_dadvs_user_defined == nullptr )
        {
            // Sensitivity function was not provided
            MORIS_ERROR( not this->has_advs(),
                    "A sensitivity evaluation function was not provided to user-defined field %s. "
                    "Please make sure that you provide this functionor specify that sensitivities are not required.",
                    this->get_name().c_str() );
            get_dfield_dadvs_user_defined = &( User_Defined_Field::no_sensitivities );
            mSensitivities.set_size( 1, this->mFieldVariables.size(), 0.0 );
        }
        else
        {
            // Check sensitivity function
            this->get_dfield_dadvs_user_defined( { { 0.0, 0.0, 0.0 } }, mFieldVariables, mSensitivities );

            // Check for row vector
            MORIS_ERROR( mSensitivities.n_rows() == 1,
                    "Geometry %s: A user-defined field must provide a row vector for sensitivities.",
                    this->get_name().c_str() );

            // Check for size
            MORIS_ERROR( mSensitivities.n_cols() == mFieldVariables.size(),
                    "Geometry %s: A user-defined field must have a sensitivity vector with a length equal to the total "
                    "number of field variables (ADVs + constants). sensitivities: %zu   geom variables: %zu \n",
                    this->get_name().c_str(),
                    mSensitivities.n_cols(),
                    mFieldVariables.size() );

            // Check for values not nan/infinity
            for ( uint tSensitivityIndex = 0; tSensitivityIndex < mSensitivities.n_cols(); tSensitivityIndex++ )
            {
                MORIS_ASSERT( std::isfinite( mSensitivities( tSensitivityIndex ) ),
                        "Geometry %s: There is an error in a user-defined field sensitivity (evaluates to nan/infinity).",
                        this->get_name().c_str() );
            }
        }

        // Check if a spatial gradient function was provided
        if ( get_dfield_dcoordinates_select == nullptr )
        {
            // Use finite difference to compute spatial gradient if not
            get_dfield_dcoordinates_select = [ this ](
                                                     const Matrix< DDRMat >& aCoordinates,
                                                     const Vector< real >&   aParameters,
                                                     Matrix< DDRMat >&       aSensitivities ) {
                this->get_dfield_dcoordinates_FD(
                        aCoordinates,
                        aParameters,
                        aSensitivities );
            };
            MORIS_LOG_INFO(
                    "No spatial gradient evaluation function was provided to user-defined field %s. "
                    "Finite difference will be used to compute the spatial gradient.",
                    this->get_name().c_str() );
        }
        else
        {
            // Check spatial gradient function
            // FIXME: There is no way to know what the spatial dimension is in this function. Thus there is no check for this
            Matrix< DDRMat > tSensitivities( 1, 3 );
            this->get_dfield_dcoordinates_select( { { 0.0, 0.0, 0.0 } }, mFieldVariables, tSensitivities );

            // // Check for row vector
            // MORIS_ERROR( tSensitivities.n_rows() == 1,
            //         "Geometry %s: A user-defined field must provide a row vector for spatial gradients.",
            //         this->get_name().c_str() );

            // Check for values not nan/infinity
            for ( uint tSensitivityIndex = 0; tSensitivityIndex < tSensitivities.n_cols(); tSensitivityIndex++ )
            {
                MORIS_ASSERT( std::isfinite( tSensitivities( tSensitivityIndex ) ),
                        "Geometry %s: There is an error in a user-defined field spatial gradient (evaluates to nan/infinity).",
                        this->get_name().c_str() );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void User_Defined_Field::no_sensitivities(
            const Matrix< DDRMat >& aCoordinates,
            const Vector< real >&   aParameters,
            Matrix< DDRMat >&       aSensitivities )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::gen
