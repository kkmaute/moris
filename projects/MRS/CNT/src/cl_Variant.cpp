/*
* Copyright (c) 2022 University of Colorado
* Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
*
*------------------------------------------------------------------------------------
*
* cl_Validator.cpp
*
*/

#include "cl_Variant.hpp"

namespace moris
{
    //--------------------------------------------------------------------------------------------------------------

    std::string get_variant_name( uint aVariantIndex )
    {
        switch ( aVariantIndex )
        {
            case 0:
                return "bool";
            case 1:
                return "uint";
            case 2:
                return "sint";
            case 3:
                return "real";
            case 4:
                return "string";
            case 5:
                return "std::pair<string, string>";
            default:
                return "";
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    std::string convert_variant_to_string( Variant aVariant )
    {
        std::stringstream tStringStream;

        switch ( aVariant.index() )
        {
            case 0:
            {
                tStringStream << std::get< bool >( aVariant );
                break;
            }
            case 1:
            {
                tStringStream << std::get< uint >( aVariant );
                break;
            }
            case 2:
            {
                tStringStream << std::get< sint >( aVariant );
                break;
            }
            case 3:
            {
                tStringStream << std::get< real >( aVariant );
                break;
            }
            case 4:
            {
                tStringStream << std::get< std::string >( aVariant );
                break;
            }
            case 5:
            {
                std::pair< std::string, std::string > tPair = std::get< std::pair< std::string, std::string > >( aVariant );
                tStringStream << tPair.first << "," << tPair.second;
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Invalid variant index." );
            }
        }

        return tStringStream.str();
    }

    //--------------------------------------------------------------------------------------------------------------
    
    template<>
    Variant make_variant( std::string aParameter )
    {
        // remove whitespaces from string
        split_trim_string( aParameter, ",;" );
        return aParameter;
    }

    //--------------------------------------------------------------------------------------------------------------

    template<>
    Variant make_variant( std::pair< std::string, std::string > aParameterValue )
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
    Variant make_variant( const char* aParameterValue )
    {
        // Convert to a string
        std::string tString( aParameterValue );
        return make_variant( tString );
    }
    
    //--------------------------------------------------------------------------------------------------------------
}
