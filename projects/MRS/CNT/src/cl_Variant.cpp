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

    Design_Variable::Design_Variable( real aConstantValue )
            : mValue( aConstantValue )
            , mLowerBound( aConstantValue )
            , mUpperBound( aConstantValue )
            , mIsConstant( true )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Design_Variable::Design_Variable( real aLowerBound, real aInitialValue, real aUpperBound )
            : mValue ( aInitialValue )
            , mLowerBound( aLowerBound )
            , mUpperBound( aUpperBound )
            , mIsConstant( false )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Design_Variable::is_constant()
    {
        return mIsConstant;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Design_Variable::get_value()
    {
        return mValue;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Design_Variable::get_lower_bound()
    {
        return mLowerBound;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Design_Variable::get_upper_bound()
    {
        return mUpperBound;
    }

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
            case 6:
                return "Vector<uint>";
            case 7:
                return "Vector<real>";
            case 8:
                return "design variable";
            default:
                return "";
        }
    }

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
