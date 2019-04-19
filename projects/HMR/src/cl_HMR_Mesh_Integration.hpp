/*
 * cl_HMR_Mesh_Interpolation.hpp
 *
 *  Created on: Apr 19, 2019
 *      Author: doble
 */

#ifndef PROJECTS_HMR_SRC_CL_HMR_MESH_INTEGRATION_HPP_
#define PROJECTS_HMR_SRC_CL_HMR_MESH_INTEGRATION_HPP_

#include "cl_HMR_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

namespace moris
{
namespace hmr
{

class Integration_Mesh_HMR : public Mesh, public mtk::Integration_Mesh
{
public:
    Integration_Mesh_HMR(std::shared_ptr< Database > aDatabase,
                           const uint & aLagrangeOrder,
                           const uint & aLagrangePattern):
        Mesh( aDatabase, aLagrangeOrder, aLagrangePattern )
    {

    }
};
}
}


#endif /* PROJECTS_HMR_SRC_CL_HMR_MESH_INTEGRATION_HPP_ */
