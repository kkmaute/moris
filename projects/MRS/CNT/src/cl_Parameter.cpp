/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Parameter.cpp
 *
 */

#include "cl_Parameter.hpp"

namespace moris
{
    //--------------------------------------------------------------------------------------------------------------

    Parameter::Parameter( const Parameter& aParameter )
            : mValue( aParameter.mValue )
    {
        // Note: This will have to change if we want to use custom validators.
        if ( aParameter.mValidator )
        {
            mValidator = aParameter.mValidator->copy();
        }
        else
        {
            mValidator = nullptr;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Parameter::~Parameter()
    {
        delete mValidator;
    }

    //--------------------------------------------------------------------------------------------------------------

    std::string Parameter::get_string() const
    {
        return convert_variant_to_string( mValue );
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Parameter::index() const
    {
        return mValue.index();
    }

    //--------------------------------------------------------------------------------------------------------------

    template<>
    Variant Parameter::make_valid_variant( uint aParameterValue )
    {
        // Make value into a variant
        Variant tParameterVariant = make_variant( aParameterValue );

        // If variant is not valid, retry with a vector
        if ( not mValidator->parameter_is_valid( tParameterVariant ) )
        {
            tParameterVariant = make_variant( Vector< uint >( { aParameterValue } ) );
        }

        // Return variant
        return tParameterVariant;
    }

    //--------------------------------------------------------------------------------------------------------------

    template<>
    Variant Parameter::make_valid_variant( real aParameterValue )
    {
        // Make value into a variant
        Variant tParameterVariant = make_variant( aParameterValue );

        // If variant is not valid, retry with a vector
        if ( not mValidator->parameter_is_valid( tParameterVariant ) )
        {
            tParameterVariant = make_variant( Vector< real >( { aParameterValue } ) );
        }

        // Return variant
        return tParameterVariant;
    }

    //--------------------------------------------------------------------------------------------------------------
}
