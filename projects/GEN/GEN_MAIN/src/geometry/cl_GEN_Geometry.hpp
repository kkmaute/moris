#ifndef MORIS_CL_GEN_GEOMETRY_HPP
#define MORIS_CL_GEN_GEOMETRY_HPP

#include "cl_GEN_Field.hpp"

namespace moris
{
    namespace ge
    {
        class Geometry : virtual public Field
        {

        public:

            /**
             * Constructor for a geometry, neeeded for inheritance
             */
            Geometry();

        };
    }
}

#endif /*MORIS_CL_GEN_GEOMETRY_HPP*/
