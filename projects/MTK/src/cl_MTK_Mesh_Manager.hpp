/*
 * cl_MTK_Mesh_Manager.hpp
 *
 *  Created on: Apr 15, 2019
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_MESH_MANAGER_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_MESH_MANAGER_HPP_

#include "typedefs.hpp"
#include "cl_Cell.hpp"

namespace moris
{
    namespace mtk
    {
        class Interpolation_Mesh;
        class Integration_Mesh;

        class Mesh_Manager
        {
            private:

                moris::Cell<Interpolation_Mesh*> mInterpolationMesh;
                moris::Cell<Integration_Mesh*>   mIntegrationMesh;

                moris::Cell< bool >              mIsOwned;

            public:

                Mesh_Manager();

                //--------------------------------------------------------------------

                ~Mesh_Manager();

                //--------------------------------------------------------------------

                uint
                register_mesh_pair(
                        Interpolation_Mesh* aInterpMesh,
                        Integration_Mesh*   aIntegrationMesh,
                        bool                aIsOwned = false );

                //--------------------------------------------------------------------

                void
                get_mesh_pair(
                        moris::moris_index     aPairIndex,
                        Interpolation_Mesh * & aInterpMesh,
                        Integration_Mesh   * & aIntegrationMesh);

                //--------------------------------------------------------------------

                Interpolation_Mesh*
                get_interpolation_mesh(moris::moris_index aMeshIndex);

                //--------------------------------------------------------------------

                Integration_Mesh*
                get_integration_mesh(moris::moris_index aMeshIndex);
        };
    }
}
#endif /* PROJECTS_MTK_SRC_CL_MTK_MESH_MANAGER_HPP_ */
