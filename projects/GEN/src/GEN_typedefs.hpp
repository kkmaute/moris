/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * GEN_Data_Types.hpp
 *
 */

#ifndef MORIS_GEN_TYPEDEFS_HPP
#define MORIS_GEN_TYPEDEFS_HPP

#include "cl_Bitset.hpp"

// Constants
#define MAX_GEOMETRIES 500

namespace moris::ge
{
    typedef Bitset<MAX_GEOMETRIES> Geometry_Bitset;
}

#endif //MORIS_GEN_TYPEDEFS_HPP

