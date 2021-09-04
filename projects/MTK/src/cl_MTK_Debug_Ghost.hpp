/*
 * cl_MTK_Debug_Ghost.hpp
 *
 *  Created on: Jun 8, 2021
 *      Author: momo
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
