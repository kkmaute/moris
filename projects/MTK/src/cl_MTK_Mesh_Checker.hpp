/*
 * cl_MTK_Mesh_Checker.hpp
 *
 *  Created on: Jan 15, 2020
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_MESH_CHECKER_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_MESH_CHECKER_HPP_


#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Cluster.hpp"

namespace moris
{

namespace mtk
{

class Mesh_Checker
{
public:
    Mesh_Checker();
    ~Mesh_Checker();

    bool verify_double_side_sets(Integration_Mesh const * aIgMesh);
    bool verify_side_cluster(Cluster const* aCluster, enum Master_Slave aMasterSlave = Master_Slave::MASTER);

};
}

}

#endif /* PROJECTS_MTK_SRC_CL_MTK_MESH_CHECKER_HPP_ */
