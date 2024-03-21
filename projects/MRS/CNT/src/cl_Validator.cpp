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
     * @param aWhichIndex Which index from the boost variant
     * @return Type name
     */
    std::string get_which_variant_name( sint aWhichIndex )
    {
        switch ( aWhichIndex )
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

    Validator::Validator( const Variant& aDefaultParameter )
            : mTypeIndex( aDefaultParameter.which() )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void Validator::check_parameter_type( const std::string& aParameterName, const Variant& aParameterVariant )
    {
        MORIS_ERROR( aParameterVariant.which() == mTypeIndex,
                "Parameter %s requires a %s, but was given a %s.",
                aParameterName.c_str(),
                get_which_variant_name( mTypeIndex ).c_str(),
                get_which_variant_name( aParameterVariant.which() ).c_str() );
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Validator::is_parameter_valid( const Variant& aParameterVariant )
    {
        return true;
    }

    //--------------------------------------------------------------------------------------------------------------

    std::string Validator::get_valid_values()
    {
        return "";
    }

    //--------------------------------------------------------------------------------------------------------------

}
