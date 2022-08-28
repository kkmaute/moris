/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Interpolation_Mesh.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_INTERPOLATION_MESH_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_INTERPOLATION_MESH_HPP_

#include "cl_MTK_Mesh_Core.hpp"

#include "assert.hpp"
#include "cl_Matrix.hpp"

namespace moris
{
    namespace mtk
    {
        class Interpolation_Mesh: public virtual Mesh
        {
                // Functions only valid for interpolation mIntegrationMeshes
            public:
                Interpolation_Mesh(){};

                ~Interpolation_Mesh()
                {
                };
        };
    }
}

#endif /* PROJECTS_MTK_SRC_CL_MTK_INTERPOLATION_MESH_HPP_ */

