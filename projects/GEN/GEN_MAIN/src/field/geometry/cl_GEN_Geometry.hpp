/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Geometry.hpp
 *
 */

#ifndef MORIS_CL_GEN_GEOMETRY_HPP
#define MORIS_CL_GEN_GEOMETRY_HPP

#include "cl_GEN_Field.hpp"
#include "st_GEN_Geometry_Parameters.hpp"

namespace moris
{
    namespace ge
    {
        class Geometry : virtual public Field
        {
        private:
            Geometry_Field_Parameters mParameters;

        public:

            /**
             * Constructor
             */
            Geometry(Geometry_Field_Parameters aParameters);

            /**
             * Copy constructor
             *
             * @param aGeometry Geometry to copy
             */
            Geometry(std::shared_ptr<Geometry> aGeometry);

            /**
             * Sets the intersection interpolation type for this geometry.
             *
             * @param aInterpolationName Intersection interpolation name
             */
            void set_intersection_interpolation(std::string aInterpolationName);

            /**
             * Gets the intersection interpolation type for this geometry.
             *
             * @return Intersection interpolation
             */
            Intersection_Interpolation get_intersection_interpolation();

        };
    }
}

#endif /*MORIS_CL_GEN_GEOMETRY_HPP*/

