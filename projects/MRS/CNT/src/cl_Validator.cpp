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
