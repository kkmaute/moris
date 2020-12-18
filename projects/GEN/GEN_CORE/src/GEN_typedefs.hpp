#ifndef MORIS_GEN_TYPEDEFS_HPP
#define MORIS_GEN_TYPEDEFS_HPP

#include "cl_Bitset.hpp"

// Constants
#define MAX_GEOMETRIES 16

namespace moris
{
    namespace ge
    {
        typedef Bitset<MAX_GEOMETRIES> Geometry_Bitset;
    }
}

#endif //MORIS_GEN_TYPEDEFS_HPP
