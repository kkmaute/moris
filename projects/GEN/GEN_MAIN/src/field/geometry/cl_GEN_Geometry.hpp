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
            Geometry( Geometry_Field_Parameters aParameters );

            /**
             * Copy constructor
             *
             * @param aGeometry Geometry to copy
             */
            Geometry( std::shared_ptr<Geometry> aGeometry );

            /**
             * Sets the intersection interpolation type for this geometry.
             *
             * @param aInterpolationName Intersection interpolation name
             */
            void set_intersection_interpolation( std::string aInterpolationName );

            /**
             * Gets the intersection interpolation type for this geometry.
             *
             * @return Intersection interpolation
             */
            Intersection_Interpolation get_intersection_interpolation();

            /**
             * Sets the intersection mode for this geometry
             * 
             * @param aIntersectionMode the method for which intersections are computed. Modes found in cl_GEN_Pdv_Enums.hpp
             */
            void set_intersection_mode( Intersection_Mode aIntersectionMode );

            /**
             * Gets the mode of intersection used for this geometry
             * 
             * @return Intersection_Mode enum
             */
            Intersection_Mode get_intersection_mode();

            /**
             * Sets the isocontour level for this geometry
             * 
             * @param aIsocontourThreshold The isocontour level that determines the geometry interface
             */
            void set_isocontour_threshold( real aIsocontourThreshold );

            /**
             * Accesses the isocontour level that determines the interface for this geometry
             * 
             * @return the isocontour level that determines the geometry interface
             */
            real get_isocontour_threshold();
            
            /**
             * Sets the isocontour tolerance for this geometry
             *
             * @param aIsocontourTolerance The desired isocontour line tolerance. Should be larger than 1e-14.
             */
            void set_isocontour_tolerance( real aIsocontourTolerance );

            /**
             * Acccesses the isocontour tolerance for this geometry
             *
             * @return isocontour tolerance
             */
            real get_isocontour_tolerance();

            /**
             * @brief sets the intersection tolerance for this geometry
             * 
             * @param aIntersectionTolerance The desired intersection tolerance. Should be larger than 1e-14
             */
            void set_intersection_tolerance( real aIntersectionTolerance );

            /**
             * Accesses the intersection tolerance for this geometry
             * 
             * @return The real value of the intersection tolerance
             */
            real get_intersection_tolerance();

        };
    }
}

#endif /*MORIS_CL_GEN_GEOMETRY_HPP*/

