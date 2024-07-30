/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_enum_macros.hpp
 *
 */

#pragma once

#include "cl_Vector.hpp"

namespace moris
{
#define CAT( a, ... ) PRIMITIVE_CAT( a, __VA_ARGS__ )
#define PRIMITIVE_CAT( a, ... ) a##__VA_ARGS__

#define ENUM_MACRO( name, ... ) CAT( ENUM_MACRO_, NARGS( __VA_ARGS__ ) )( name, __VA_ARGS__ )

    // #define NARGS_SEQ( _1, _2, _3, _4, _5, _6, _7, _8, _9, 10, N, ... ) N
    // #define NARGS( ... ) NARGS_SEQ( __VA_ARGS__, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1 )

#include "../../../PRM/ENM/src/fn_enum_macros_generated.hpp"
}    // namespace moris
