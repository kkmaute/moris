/*
* Copyright (c) 2022 University of Colorado
* Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
*
*------------------------------------------------------------------------------------
*
* cl_Variant.hpp
*
*/

#include <variant>
#include "moris_typedefs.hpp"
#include "cl_Design_Variable.hpp"
#include "fn_Parsing_Tools.hpp"

#pragma once

#define GET_TYPE_NAME( ... ) template<> inline std::string get_type_name< __VA_ARGS__ >(){ return #__VA_ARGS__; }

namespace moris
{
    // Variant typedef
    typedef std::variant< bool, uint, sint, real, std::string, std::pair< std::string, std::string >,
            Vector< uint >, Vector< real >, Vector< std::string >, Design_Variable > Variant;

    /**
     * Gets the name of a data type stored in the Variant, for printing purposes.
     *
     * @tparam T Type in the variant
     * @return Type name
     */
    template< typename T >
    std::string get_type_name()
    {
        MORIS_ERROR( false, "This type is not currently stored in the Variant." );
        return "";
    }

    // Specializations
    GET_TYPE_NAME( bool )
    GET_TYPE_NAME( uint )
    GET_TYPE_NAME( sint )
    GET_TYPE_NAME( real )
    GET_TYPE_NAME( std::string )
    GET_TYPE_NAME( std::pair< std::string, std::string > )
    GET_TYPE_NAME( Vector< uint > )
    GET_TYPE_NAME( Vector< real > )
    GET_TYPE_NAME( Vector< std::string > )
    GET_TYPE_NAME( Design_Variable )

    /**
     * Gets the index of a type with a type argument instead of an instantiated variant.
     *
     * @tparam T Type in the variant
     * @return Variant index
     */
    template< typename T >
    uint variant_index()
    {
        T tValue;
        Variant tVariant = tValue;
        return tVariant.index();
    }

    /**
     * Splits a variant into a vector of variants, if the original variant itself is a vector.
     *
     * @param aVariant Variant, which may be a vector
     * @return Vector of variants, all with size 1
     */
    Vector< Variant > split_variant( const Variant& aVariant );

    /**
     * Gets the size of the underlying data in the variant.
     *
     * @param aVariant Variant
     * @return Size of underlying vector/matrix, or 1 for primitive data types
     */
    uint get_size( const Variant& aVariant );

    /**
     * Comparison operator determining if two variants are equal, based on their stored values.
     *
     * @param aLeft Left variant argument
     * @param aRight Right variant argument
     * @return If both variants store equal values
     */
    bool operator==( const Variant& aLeft, const Variant& aRight );

    /**
     * Gets a string of the value stored inside of a variant.
     *
     * @param aVariant Input variant
     * @return Value as a string
     */
    std::string convert_variant_to_string( Variant aVariant );

    /**
     * Takes the input value and converts it to a variant type.
     *
     * @tparam T Input parameter type
     * @param aParameterValue Input parameter
     */
    template< typename T >
    Variant make_variant( T aParameterValue )
    {
        return aParameterValue;
    }

    // Declare template specializations making variants
    template<> Variant make_variant( std::string );
    template<> Variant make_variant( std::pair< std::string, std::string > );
    template<> Variant make_variant( const char* );
    template<> Variant make_variant( Vector< const char* > );
}
