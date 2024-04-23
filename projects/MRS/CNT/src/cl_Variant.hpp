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

#define GET_TYPE_NAME( ... ) template<> inline std::string get_type_name< __VA_ARGS__ >(){ return "__VA_ARGS__"; }

namespace moris
{
    // Variant typedef
    typedef std::variant< bool, uint, sint, real, std::string, std::pair< std::string, std::string >,
            Vector< uint >, Vector< real >, Design_Variable > Variant;

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
    GET_TYPE_NAME( Design_Variable )

    /**
     * Gets the index of a type with a type argument instead of an instantiated variant.
     *
     * @tparam T Type in the variant
     * @return Variant index
     */
    template< typename T >
    uint get_variant_index()
    {
        T aValue;
        Variant tVariant = aValue;
        return tVariant.index();
    }

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
    template<> Variant make_variant( std::string aParameter );
    template<> Variant make_variant( std::pair< std::string, std::string > aParameterValue );
    template<> Variant make_variant( const char* aParameterValue );
}
