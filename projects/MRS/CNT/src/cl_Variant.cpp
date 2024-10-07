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

#define BOOL 0
#define UINT 1
#define SINT 2
#define REAL 3
#define STRING 4
#define PAIR 5
#define VECTOR_UINT 6
#define VECTOR_SINT 7
#define VECTOR_REAL 8
#define VECTOR_STRING 9
#define DESIGN_VARIABLE 10

namespace moris
{
    //--------------------------------------------------------------------------------------------------------------

    Vector< Variant > split_variant( const Variant& aVariant )
    {
        Vector< Variant > tVectorOfVariants;
        switch ( aVariant.index() )
        {
            case VECTOR_UINT:
            {
                auto tVector = std::get< Vector< uint > >( aVariant );
                tVectorOfVariants.reserve( tVector.size() );
                for ( auto iElement : tVector )
                {
                    tVectorOfVariants.push_back( Variant( iElement ) );
                }
                break;
            }
            case VECTOR_SINT:
            {
                auto tVector = std::get< Vector< sint > >( aVariant );
                tVectorOfVariants.reserve( tVector.size() );
                for ( auto iElement : tVector )
                {
                    tVectorOfVariants.push_back( Variant( iElement ) );
                }
                break;
            }
            case VECTOR_REAL:
            {
                auto tVector = std::get< Vector< real > >( aVariant );
                tVectorOfVariants.reserve( tVector.size() );
                for ( auto iElement : tVector )
                {
                    tVectorOfVariants.push_back( Variant( iElement ) );
                }
                break;
            }
            case VECTOR_STRING:
            {
                auto tVector = std::get< Vector< std::string > >( aVariant );
                tVectorOfVariants.reserve( tVector.size() );
                for ( auto iElement : tVector )
                {
                    tVectorOfVariants.push_back( Variant( iElement ) );
                }
                break;
            }
            default:
            {
                tVectorOfVariants = { aVariant };
                break;
            }
        }

        // Return vector of variants
        return tVectorOfVariants;
    }

    //--------------------------------------------------------------------------------------------------------------

    uint get_size( const Variant& aVariant )
    {
        switch ( aVariant.index() )
        {
            case VECTOR_UINT:
            {
                return std::get< Vector< uint > >( aVariant ).size();
            }
            case VECTOR_SINT:
            {
                return std::get< Vector< uint > >( aVariant ).size();
            }
            case VECTOR_REAL:
            {
                return std::get< Vector< real > >( aVariant ).size();
            }
            case VECTOR_STRING:
            {
                return std::get< Vector< std::string > >( aVariant ).size();
            }
            default:
            {
                return 1;
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    template< typename T >
    bool compare( const Variant& aLeft, const Variant& aRight )
    {
        return std::get< T >( aLeft ) == std::get< T >( aRight );
    }

    //--------------------------------------------------------------------------------------------------------------

    bool operator==( const Variant& aLeft, const Variant& aRight )
    {
        // If variants have different indices, they are not equal
        if ( aLeft.index() not_eq aRight.index() )
        {
            return false;
        }

        // Determine equality based on type
        switch ( aLeft.index() )
        {
            case BOOL:
            {
                return compare< bool >( aLeft, aRight );
            }
            case UINT:
            {
                return compare< uint >( aLeft, aRight );
            }
            case SINT:
            {
                return compare< sint >( aLeft, aRight );
            }
            case REAL:
            {
                return compare< real >( aLeft, aRight );
            }
            case STRING:
            {
                return compare< std::string >( aLeft, aRight );
            }
            case PAIR:
            {
                return compare< std::pair< std::string, std::string > >( aLeft, aRight );
            }
            case VECTOR_UINT:
            {
                return compare< Vector< uint > >( aLeft, aRight );
            }
            case VECTOR_SINT:
            {
                return compare< Vector< sint > >( aLeft, aRight );
            }
            case VECTOR_REAL:
            {
                return compare< Vector< real > >( aLeft, aRight );
            }
            case VECTOR_STRING:
            {
                return compare< Vector< std::string > >( aLeft, aRight );
            }
            case DESIGN_VARIABLE:
            {
                return false;
            }
            default:
            {
                MORIS_ERROR( false, "Invalid variant index: %lu.", aLeft.index() );
                return false;
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    template< typename T >
    void vector_variant_to_string( std::stringstream& aStringStream, Variant aVariant )
    {
        aStringStream << "{";
        Vector< T > tVector = std::get< Vector< T > >( aVariant );
        std::string tDelimiter;
        for ( const auto& iVectorElement : tVector )
        {
            aStringStream << tDelimiter << iVectorElement;
            tDelimiter = ", ";
        }
        aStringStream << "}";
    }

    //--------------------------------------------------------------------------------------------------------------

    std::string convert_variant_to_string( Variant aVariant )
    {
        std::stringstream tStringStream;

        switch ( aVariant.index() )
        {
            case BOOL:
            {
                tStringStream << std::get< bool >( aVariant );
                break;
            }
            case UINT:
            {
                tStringStream << std::get< uint >( aVariant );
                break;
            }
            case SINT:
            {
                tStringStream << std::get< sint >( aVariant );
                break;
            }
            case REAL:
            {
                tStringStream << std::get< real >( aVariant );
                break;
            }
            case STRING:
            {
                tStringStream << "\"" << std::get< std::string >( aVariant ) << "\"";
                break;
            }
            case PAIR:
            {
                std::pair< std::string, std::string > tPair = std::get< std::pair< std::string, std::string > >( aVariant );
                tStringStream << tPair.first << "," << tPair.second;
                break;
            }
            case VECTOR_UINT:
            {
                vector_variant_to_string< uint >( tStringStream, aVariant );
                break;
            }
            case VECTOR_SINT:
            {
                vector_variant_to_string< sint >( tStringStream, aVariant );
                break;
            }
            case VECTOR_REAL:
            {
                vector_variant_to_string< real >( tStringStream, aVariant );
                break;
            }
            case VECTOR_STRING:
            {
                vector_variant_to_string< std::string >( tStringStream, aVariant );
                break;
            }
            case DESIGN_VARIABLE:
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
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Invalid variant index: %lu.", aVariant.index() );
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
}    // namespace moris
