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
#include "cl_Cell.hpp"

namespace moris
{
  namespace mtk
  {


    class MTK_Manager
    {
    public:

    private:

      // interpolation mesh data
      moris::Cell<std::string>        mInterpolationMeshLabels;
      moris::Cell<Interpolation_Mesh> mInterpolationMeshes;

      // integration mesh data
      moris::Cell<std::string>        mIntegrationMeshLabels;
      moris::Cell<Integration_Mesh>   mIntegrationMeshes;

    }
  }

}



#endif /* PROJECTS_MTK_SRC_CL_MTK_MESH_MANAGER_HPP_ */
