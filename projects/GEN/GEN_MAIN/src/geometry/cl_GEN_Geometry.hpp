#ifndef MORIS_CL_GEN_GEOMETRY_HPP
#define MORIS_CL_GEN_GEOMETRY_HPP

#include "cl_GEN_Field_Base.hpp"

namespace moris
{
    namespace ge
    {
        class Geometry : virtual public Field
        {
        public:

            /**
             * Trivial constructor, necessary for clean virtual inheritance without default constructor in base class
             */
            Geometry();

            /**
             * Lets the geometry engine know if sensitivities are available, otherwise it will perform finite
             * differencing instead for intersection locations
             *
             * @return If sensitivities are implemented or not (true in base class)
             */
            virtual bool sensitivities_available();
        };
    }
}

#endif /*MORIS_CL_GEN_GEOMETRY_HPP*/
