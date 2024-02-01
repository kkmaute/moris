/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Pdv_Enums.hpp
 *
 */

#pragma once

#include "cl_Bitset.hpp"
#include "cl_Map.hpp"

// Constants
#define MAX_GEOMETRIES 500

namespace moris::ge
{
    // Typdefs
    typedef Bitset< MAX_GEOMETRIES > Geometry_Bitset;

    // Enums
    enum class PDV_Type
    {
        X_COORDINATE,
        Y_COORDINATE,
        Z_COORDINATE,
        DENSITY,
        TEMPERATURE,
        ELASTIC_MODULUS,
        LS1,
        LS2,
        UNDEFINED
    };

    /**
     * Gets a map going from a std::string to a PDV_Type enum. Used to convert parameter list arguments.
     *
     * @return PDV_Type map
     */
    moris::map< std::string, PDV_Type > get_pdv_type_map();
}
