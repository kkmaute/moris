/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_XTK_Cell_Cluster_Group.cpp
 *
 */

#include "cl_XTK_Cell_Cluster_Group.hpp"
#include "cl_MTK_Cluster.hpp"

using namespace moris;

namespace moris::xtk
{
    //------------------------------------------------------------------------------

    Cell_Cluster_Group::Cell_Cluster_Group(
            const moris_index                         aDiscretizationMeshIndex,
            Vector< std::shared_ptr< mtk::Cluster > > aClusters )
            : mtk::Cell_Cluster_Group( aDiscretizationMeshIndex )
            , mClusters( aClusters )
    {
        // only initialize member variables
    }

    //------------------------------------------------------------------------------

    const Vector< mtk::Cluster const * >
    Cell_Cluster_Group::get_clusters_in_group() const
    {
        // get the number of cluster and initialize a list of raw pointers with this
        uint                           tNumClustersInGroup = mClusters.size();
        Vector< mtk::Cluster const * > tClusters( tNumClustersInGroup );

        // fill list with raw cluster pointers
        for ( uint iCluster = 0; iCluster < tNumClustersInGroup; iCluster++ )
        {
            tClusters( iCluster ) = mClusters( iCluster ).get();
        }

        // return list with raw cluster pointers
        return tClusters;
    }

    //------------------------------------------------------------------------------

}    // namespace moris::xtk
