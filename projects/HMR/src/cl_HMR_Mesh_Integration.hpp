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
#include "cl_MTK_Cell_Cluster.hpp"
#include "cl_MTK_Side_Cluster.hpp"

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


    //fixme: IMPLEMENT THESE THREE FUNCTIONS
    mtk::Cell_Cluster const &
    get_cell_cluster(mtk::Cell const & aInterpCell) const
    {
        MORIS_ERROR(0,"Cell clusters not implemented in HMR");
        return *mDummyCluster;
    }

    /*
     * Get block set names
     */
    moris::Cell<std::string>
    get_block_set_names() const
    {
        MORIS_ERROR(0,"get_block_set_names not implemented in HMR");
        return moris::Cell<std::string>(0);
    }

    /*
     * Get cell clusters within a block set
     */
    moris::Cell<mtk::Cell_Cluster const *>
    get_cell_clusters_in_set(moris_index aBlockSetOrdinal) const
    {
        MORIS_ERROR(0,"get_block_set_names not implemented in HMR");
        moris::Cell<mtk::Cell_Cluster const *> tCellInCluster(0);
        return tCellInCluster;
    }

    moris::Cell<mtk::Side_Cluster const *>
    get_side_set_cluster(moris_index aSideSetOrdinal) const
    {
        MORIS_ERROR(0,"get_side_set_cluster not implemented in HMR");
        moris::Cell<mtk::Side_Cluster const *> tSideClusters(0);
        return tSideClusters;
    }

private:
    mtk::Cell_Cluster * mDummyCluster     = nullptr;
};
}
}


#endif /* PROJECTS_HMR_SRC_CL_HMR_MESH_INTEGRATION_HPP_ */
