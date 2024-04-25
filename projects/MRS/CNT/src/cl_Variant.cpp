/*
* Copyright (c) 2022 University of Colorado
* Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
*
*------------------------------------------------------------------------------------
*
* cl_Variant.cpp
*
*/

#include "cl_Variant.hpp"

namespace moris
{
    //--------------------------------------------------------------------------------------------------------------

    template< typename T >
    void vector_variant_to_string( std::stringstream& aStringStream, Variant aVariant )
    {
        Vector< T > tVector = std::get< Vector< T > >( aVariant );
        for ( auto iVectorElement : tVector )
        {
            aStringStream << iVectorElement;
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
                tStringStream << "\"" << std::get< std::string >( aVariant ) << "\"";
                break;
            }
            case 5:
            {
                std::pair< std::string, std::string > tPair = std::get< std::pair< std::string, std::string > >( aVariant );
                tStringStream << tPair.first << "," << tPair.second;
                break;
            }
            case 6:
            {
                vector_variant_to_string< uint >( tStringStream, aVariant );
                break;
            }
            case 7:
            {
                vector_variant_to_string< real >( tStringStream, aVariant );
                break;
            }
            case 8:
            {
                vector_variant_to_string< std::string >( tStringStream, aVariant );
                break;
            }
            case 9:
            {
                Design_Variable tDesignVariable = std::get< Design_Variable >( aVariant );
                if ( tDesignVariable.is_constant() )
                {
                    tStringStream << tDesignVariable.get_value();
                }
                else
                {
                    tStringStream << "{" << tDesignVariable.get_lower_bound()
                                  << ", " << tDesignVariable.get_value()
                                  << ", " << tDesignVariable.get_upper_bound() << "}";
                }
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
    Variant make_variant( std::string aParameterValue )
    {
        // remove whitespaces from string
        split_trim_string( aParameterValue, ",;" );
        return aParameterValue;
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

    template<>
    Variant make_variant( Vector< const char* > aParameterValue )
    {
        // Convert to a vector of strings
        Vector< std::string > tStringVector;
        tStringVector.reserve( aParameterValue.size() );
        for ( auto iChar : aParameterValue )
        {
            tStringVector.push_back( iChar );
        }
        return make_variant( tStringVector );
    }
    
    //--------------------------------------------------------------------------------------------------------------
}
