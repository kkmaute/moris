/*
 * cl_MTK_Mesh_Manager.hpp
 *
 *  Created on: Apr 15, 2019
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_MESH_MANAGER_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_MESH_MANAGER_HPP_

#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "typedefs.hpp"
#include "cl_Cell.hpp"

namespace moris
{
  namespace mtk
  {
  class Mesh_Manager
  {
  private:
      moris::Cell<Interpolation_Mesh*> mInterpolationMesh;
      moris::Cell<Integration_Mesh*>   mIntegrationMesh;

  public:

      Mesh_Manager(){};

      uint
      register_mesh_pair(Interpolation_Mesh* aInterpMesh,
                         Integration_Mesh*   aIntegrationMesh)
      {
          mInterpolationMesh.push_back(aInterpMesh);
          mIntegrationMesh.push_back(aIntegrationMesh);

          return mInterpolationMesh.size()-1;
      }

      void
      get_mesh_pair(moris::moris_index aPairIndex,
                    Interpolation_Mesh * & aInterpMesh,
                    Integration_Mesh   * & aIntegrationMesh)
      {
          aInterpMesh = mInterpolationMesh(aPairIndex);
          aIntegrationMesh = mIntegrationMesh(aPairIndex);
      }

      Interpolation_Mesh*
      get_interpolation_mesh(moris::moris_index aMeshIndex)
      {
          return mInterpolationMesh(aMeshIndex);
      }

      Integration_Mesh*
      get_integration_mesh(moris::moris_index aMeshIndex)
      {
          return mIntegrationMesh(aMeshIndex);
      }


  };

}

}



#endif /* PROJECTS_MTK_SRC_CL_MTK_MESH_MANAGER_HPP_ */
