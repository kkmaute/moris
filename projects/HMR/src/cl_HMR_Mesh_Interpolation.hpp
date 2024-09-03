/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Mesh_Interpolation.hpp
 *
 */

#ifndef PROJECTS_HMR_SRC_CL_HMR_MESH_INTERPOLATION_HPP_
#define PROJECTS_HMR_SRC_CL_HMR_MESH_INTERPOLATION_HPP_

#include <utility>

#include "cl_HMR_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"

namespace moris::hmr
{
    class Interpolation_Mesh_HMR : public virtual Mesh
            , public mtk::Interpolation_Mesh
    {
        //-------------------------------------------------------------------------------
      public:
        //-------------------------------------------------------------------------------

        Interpolation_Mesh_HMR(
                std::shared_ptr< Database > aDatabase,
                uint                        aLagrangeMeshIndex )
                : Mesh(
                          std::move( aDatabase ),
                          aLagrangeMeshIndex )
        {
        }

        //-------------------------------------------------------------------------------

        Interpolation_Mesh_HMR(
                std::shared_ptr< Database > aDatabase,
                uint                        aLagrangeOrder,
                uint                        aLagrangePattern )
                : Mesh(
                          std::move( aDatabase ),
                          aLagrangeOrder,
                          aLagrangePattern )
        {
        }

        //-------------------------------------------------------------------------------

        Interpolation_Mesh_HMR(
                std::shared_ptr< Database > aDatabase,
                uint                        aOrder,
                uint                        aLagrangePattern,
                uint                        aBsplinePattern )
                : Mesh(
                          std::move( aDatabase ),
                          aOrder,
                          aLagrangePattern,
                          aBsplinePattern )
        {
        }

        //-------------------------------------------------------------------------------

        Interpolation_Mesh_HMR(
                std::shared_ptr< Database > aDatabase,
                uint                        aLagrangeOrder,
                uint                        aLagrangePattern,
                uint                        aBSplineOrder,
                uint                        aBsplinePattern )
                : Mesh(
                          std::move( aDatabase ),
                          aLagrangeOrder,
                          aLagrangePattern,
                          aBSplineOrder,
                          aBsplinePattern )
        {
        }

        //-------------------------------------------------------------------------------

        ~Interpolation_Mesh_HMR() override
        {
        }

        //-------------------------------------------------------------------------------

    };    // end class: Interpolation_Mesh_HMR

    //-------------------------------------------------------------------------------

}    // namespace moris::hmr
#endif /* PROJECTS_HMR_SRC_CL_HMR_MESH_INTERPOLATION_HPP_ */
