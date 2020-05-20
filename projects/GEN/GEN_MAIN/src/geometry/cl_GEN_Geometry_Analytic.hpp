#ifndef MORIS_CL_GEN_GEOMETRY_ANALYTIC_HPP
#define MORIS_CL_GEN_GEOMETRY_ANALYTIC_HPP

#include "cl_GEN_Field_Base.hpp"

namespace moris
{
    namespace ge
    {
        class Geometry_Analytic : virtual public Field
        {
        public:

            /**
             * Trivial constructor, necessary for clean virtual inheritance without default constructor in base class
             */
            Geometry_Analytic();
        };
    }
}

#endif /* MORIS_CL_GEN_GEOMETRY_ANALYTIC_HPP */
