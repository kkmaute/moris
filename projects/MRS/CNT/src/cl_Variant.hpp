/*
* Copyright (c) 2022 University of Colorado
* Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
*
*------------------------------------------------------------------------------------
*
* cl_Validator.hpp
*
*/

#include <variant>
#include "moris_typedefs.hpp"
#include "fn_Parsing_Tools.hpp"

#pragma once

namespace moris
{
    // Variant typedef
    typedef std::variant< bool, uint, sint, real, std::string, std::pair< std::string, std::string > > Variant;

    /**
     * Gets the name of the type that the boost variant contains, based on the above typedef and an index
     *
     * @param aVariantIndex Which index from the boost variant
     * @return Type name
     */
    std::string get_variant_name( uint aVariantIndex );

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
