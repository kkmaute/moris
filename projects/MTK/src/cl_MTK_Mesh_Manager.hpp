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
#include "st_MTK_Mesh_Pair.hpp"

namespace moris
{
    namespace mtk
    {
        class Mesh_Manager
        {
        private:
            moris::Cell<Mesh_Pair> mMeshPairs;

        public:

            Mesh_Manager();

            //--------------------------------------------------------------------

            ~Mesh_Manager();

            //--------------------------------------------------------------------

            uint
            register_mesh_pair(
                    Interpolation_Mesh* aInterpolationMesh,
                    Integration_Mesh*   aIntegrationMesh,
                    bool                aIsOwned = false );

            //--------------------------------------------------------------------

            void
            get_mesh_pair(
                    moris_index          aPairIndex,
                    Interpolation_Mesh * & aInterpolationMesh,
                    Integration_Mesh   * & aIntegrationMesh);

            //--------------------------------------------------------------------

            Interpolation_Mesh*
            get_interpolation_mesh(moris_index aMeshIndex);

            //--------------------------------------------------------------------

            Integration_Mesh*
            get_integration_mesh(moris_index aMeshIndex);
            
        };
    }
}
#endif /* PROJECTS_MTK_SRC_CL_MTK_MESH_MANAGER_HPP_ */
