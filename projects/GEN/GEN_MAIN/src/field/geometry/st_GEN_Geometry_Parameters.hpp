/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * st_GEN_Geometry_Parameters.hpp
 *
 */

#ifndef MORIS_ST_GEN_GEOMETRY_PARAMETERS_HPP
#define MORIS_ST_GEN_GEOMETRY_PARAMETERS_HPP

#include "st_GEN_Field_Parameters.hpp"

namespace moris
{
    namespace ge
    {
        enum class Intersection_Interpolation
        {
            LINEAR,
            MULTILINEAR
        };

        /**
         * This struct contains additional parameters that are used by geometries.
         */
        struct Geometry_Parameters
        {
            Intersection_Interpolation mIntersectionInterpolation = Intersection_Interpolation::LINEAR; //! Intersection type
        };

        /**
         * This is a struct used to simplify \ref moris::ge::Geometry constructors. It contains all geometry parameters.
         */
        struct Geometry_Field_Parameters : Field_Parameters, Geometry_Parameters
        {
        };
    }
}

#endif //MORIS_ST_GEN_GEOMETRY_PARAMETERS_HPP

