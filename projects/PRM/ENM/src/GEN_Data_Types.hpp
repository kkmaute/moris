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
#include "fn_enum_macros.hpp"
#include "assert.hpp"
#include "cl_Vector.hpp"

// Constants
#define MAX_GEOMETRIES 500

namespace moris::gen
{
    // Typdefs
    typedef Bitset< MAX_GEOMETRIES > Geometry_Bitset;

    // Geometry type enum

    ENUM_MACRO( Geometry_Type,
            LEVEL_SET,
            SURFACE_MESH,
            VOXEL )

    // Field type enum
    ENUM_MACRO( Field_Type,
            NONE,
            CONSTANT,
            LINE,
            CIRCLE,
            SUPERELLIPSE,
            PLANE,
            SPHERE,
            SUPERELLIPSOID,
            SCALED_FIELD,
            COMBINED_FIELDS,
            NODAL_FROM_FILE,
            SIGNED_DISTANCE_OBJECT,
            SIGNED_DISTANCE_IMAGE,
            USER_DEFINED )

    // PDV type enum
    ENUM_MACRO( PDV_Type,
            X_COORDINATE,
            Y_COORDINATE,
            Z_COORDINATE,
            DENSITY,
            TEMPERATURE,
            ELASTIC_MODULUS,
            LS1,
            LS2,
            UNDEFINED )

    // Geometric quantities of interest (GQI) enum
    ENUM_MACRO( GQI_Type,
            CURVATURE,
            SHAPE_DIAMETER,
            FEATURE_SIZE,
            VOLUME,
            SURFACE_AREA )

    // Surface mesh regularization type enum
    ENUM_MACRO( Regularization_Type,
            NONE,
            ISOTROPIC_LAPLACIAN,
            ANSIOTROPIC_LAPLACIAN,
            TAUBIN,
            USER_DEFINED )

    /**
     * Gets a map going from a std::string to a PDV_Type enum. Used to convert parameter list arguments.
     *
     * @return PDV_Type map
     */
    moris::map< std::string, PDV_Type > get_pdv_type_map();
}    // namespace moris::gen
