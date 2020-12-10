#ifndef MORIS_CL_GEN_GEOMETRY_HPP
#define MORIS_CL_GEN_GEOMETRY_HPP

#include "cl_GEN_Field.hpp"

namespace moris
{
    namespace ge
    {
        enum class Intersection_Interpolation
        {
            LINEAR,
            MULTILINEAR
        };

        class Geometry : virtual public Field
        {

        private:
            Intersection_Interpolation mIntersectionInterpolation;

        public:

            /**
             * Constructor
             */
            Geometry(Intersection_Interpolation aInterpolation = Intersection_Interpolation::LINEAR);

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
