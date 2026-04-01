/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_STK_Enums.hpp
 *
 */

#pragma once
#include "fn_enum_macros.hpp"

// FIXME: stk namespace is taken by the library, so need to use a different namespace for these enums
// The SQI_Type should be generalized to just be a QI type, thus allowing us to remove this file entirely
namespace moris::sqi
{
    ENUM_MACRO( SQI_Type,
            VOLUME,
            RAYCAST_SHAPE_DIAMETER,
            INSCRIBED_CIRCLE_SHAPE_DIAMETER,
            SHORTEST_DISTANCE_SHAPE_DIAMETER )
}    // namespace moris::sqi
