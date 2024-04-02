/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Validator.cpp
 *
 */

#include "cl_Validator.hpp"

namespace moris
{

    //--------------------------------------------------------------------------------------------------------------

    /**
     * Gets the name of the type that the boost variant contains, based on the above typedef and an index
     *
     * @param aVariantIndex Which index from the boost variant
     * @return Type name
     */
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

    Validator::Validator( uint aTypeIndex )
            : mTypeIndex( aTypeIndex )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Validator::same_type_index( const Variant& aParameterVariant )
    {
        return ( aParameterVariant.index() == mTypeIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    Type_Validator::Type_Validator( uint aTypeIndex )
            : Validator( aTypeIndex )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Type_Validator::parameter_is_valid( const Variant& aParameterVariant )
    {
        return this->same_type_index( aParameterVariant );
    }

    //--------------------------------------------------------------------------------------------------------------

    std::string Type_Validator::get_valid_values()
    {
        return get_variant_name( mTypeIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    Validator* Type_Validator::copy()
    {
        return new Type_Validator( mTypeIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

}
