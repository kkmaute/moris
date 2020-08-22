#ifndef MORIS_FN_GEN_CREATE_SIMPLE_MESH_HPP
#define MORIS_FN_GEN_CREATE_SIMPLE_MESH_HPP

#include "cl_MTK_Interpolation_Mesh.hpp"

namespace moris
{
    namespace ge
    {
        /**
         * Creates a simple HMR interpolation mesh (2x2 elements) for GEN testing. There should be no cases where a more
         * complex mesh is needed.
         *
         * @return HMR mesh pointer
         */
        mtk::Interpolation_Mesh* create_simple_mesh(uint aNumXElements, uint aNumYElements);
    }
}

#endif //MORIS_FN_GEN_CREATE_SIMPLE_MESH_HPP
