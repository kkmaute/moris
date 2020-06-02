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

      moris::Cell< bool >              mIsOwned;

  public:

      Mesh_Manager(){};

      ~Mesh_Manager()
      {
         uint tCounter = 0;
         for( bool tIsOwned : mIsOwned )
         {
             if( tIsOwned )
             {
                 delete mIntegrationMesh(tCounter);
                 delete mInterpolationMesh(tCounter);
             }
             tCounter++;
         }

         mInterpolationMesh.clear();
         mIntegrationMesh.clear();
      };

      uint
      register_mesh_pair(Interpolation_Mesh* aInterpMesh,
                         Integration_Mesh*   aIntegrationMesh,
                         bool                aIsOwned = false )
      {
          mInterpolationMesh.push_back(aInterpMesh);
          mIntegrationMesh.push_back(aIntegrationMesh);

          // set flag which tells if mesh is owned by manager. If it is owned it ha to be deleted
          mIsOwned.push_back(aIsOwned);

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
