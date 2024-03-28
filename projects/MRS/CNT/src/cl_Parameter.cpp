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

    std::string convert_variant_to_string( Variant aVariant )
    {
        std::stringstream tStringStream;

        if ( boost::get< bool >( &aVariant ) != nullptr )
        {
            tStringStream << boost::get< bool >( aVariant );
        }
        else if ( boost::get< uint >( &aVariant ) != nullptr )
        {
            tStringStream << boost::get< uint >( aVariant );
        }
        else if ( boost::get< sint >( &aVariant ) != nullptr )
        {
            tStringStream << boost::get< sint >( aVariant );
        }
        else if ( boost::get< real >( &aVariant ) != nullptr )
        {
            tStringStream << boost::get< real >( aVariant );
        }
        else if ( boost::get< std::string >( &aVariant ) != nullptr )
        {
            tStringStream << boost::get< std::string >( aVariant );
        }
        else if ( boost::get< std::pair< std::string, std::string > >( &aVariant ) != nullptr )
        {
            std::pair< std::string, std::string > tPair = boost::get< std::pair< std::string, std::string > >( aVariant );
            tStringStream << tPair.first << "," << tPair.second;
        }
        else
        {
            MORIS_ERROR( false, "Variant conversion error." );
        }

        return tStringStream.str();
    }

    //--------------------------------------------------------------------------------------------------------------

    Parameter::Parameter( const Parameter& aParameter )
            : mValue( aParameter.mValue )
    {
        // Note: This will have to change if we want to use custom validators.
        mValidator = new Validator( *aParameter.mValidator );
    }

    //--------------------------------------------------------------------------------------------------------------

    template<>
    Variant Parameter::make_variant( std::string aParameter )
    {
        // remove whitespaces from string
        split_trim_string( aParameter, ",;" );
        return aParameter;
    }

    //--------------------------------------------------------------------------------------------------------------

    template<>
    Variant Parameter::make_variant( std::pair< std::string, std::string > aParameterValue )
    {
        // extract elements of pair
        std::string tFirst  = aParameterValue.first;
        std::string tSecond = aParameterValue.second;

        // trim off leading and trailing whitespaces
        split_trim_string( tFirst, ",;" );
        split_trim_string( tSecond, ",;" );

        // Return new pair
        return std::make_pair( tFirst, tSecond );
    }

    //--------------------------------------------------------------------------------------------------------------

    template<>
    Variant Parameter::make_variant( const char* aParameterValue )
    {
        // Convert to a string
        std::string tString( aParameterValue );
        return make_variant( tString );
    }

    //--------------------------------------------------------------------------------------------------------------

}
