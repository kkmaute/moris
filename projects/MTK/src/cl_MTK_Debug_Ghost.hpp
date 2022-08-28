/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Debug_Ghost.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_DEBUG_GHOST_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_DEBUG_GHOST_HPP_

#include "cl_MTK_Interpolation_Mesh.hpp"

namespace moris
{
    namespace mtk
    {
        class Debug_Ghost
        {
            public:
                Debug_Ghost();

                Debug_Ghost(Interpolation_Mesh*   aIpMesh);

                ~Debug_Ghost();

                bool perform();

            private:
                Interpolation_Mesh*   mIpMesh;
        };
    }
}

#endif /* PROJECTS_MTK_SRC_CL_MTK_DEBUG_GHOST_HPP_ */

