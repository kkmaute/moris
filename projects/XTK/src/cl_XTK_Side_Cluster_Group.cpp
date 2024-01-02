/*
 * Copyright (c) 2022 University of Colorado 
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. 
 * 
 * ------------------------------------------------------------------------------------ 
 * 
 * cl_ XTK_Side_Cluster_Group.cpp  
 * 
 */

#include "cl_XTK_Side_Cluster_Group.hpp"
#include "cl_MTK_Cluster.hpp"

using namespace moris;

namespace xtk
{
    //------------------------------------------------------------------------------

    Side_Cluster_Group::Side_Cluster_Group( 
            const moris_index                              aDiscretizationMeshIndex,
            Vector< std::shared_ptr< mtk::Cluster > > aClusters,
            std::shared_ptr< mtk::Cluster_Group >          aAssociatedCellClusterGroup )
            : mtk::Side_Cluster_Group( aDiscretizationMeshIndex )
            , mSideClusters( aClusters )
    {
        // check tath the pointer passed is valid
        MORIS_ASSERT( aAssociatedCellClusterGroup.get() != nullptr, 
            "xtk::Side_Cluster_Group::Side_Cluster_Group() - nullptr passed for the cell cluster group" );

        // assign pointer to the associated cell cluster group
        mAssociatedCellClusterGroup = aAssociatedCellClusterGroup;
    }

    //------------------------------------------------------------------------------

    const Vector< mtk::Cluster const* >
    Side_Cluster_Group::get_side_clusters_in_group() const
    {
        // get the number of cluster and initialize a list of raw pointers with this
        uint tNumClustersInGroup = mSideClusters.size();
        Vector< mtk::Cluster const* > tClusters( tNumClustersInGroup );
        
        // fill list with raw cluster pointers
        for( uint iCluster = 0; iCluster < tNumClustersInGroup; iCluster++ )
        {
            tClusters( iCluster ) = mSideClusters( iCluster ).get();
        }
        
        // return list with raw cluster pointers
        return tClusters;
    }

    //------------------------------------------------------------------------------

    mtk::Cluster_Group const* 
    Side_Cluster_Group::get_associated_cell_cluster_group() const
    {
        return mAssociatedCellClusterGroup.get();
    }

    //------------------------------------------------------------------------------

} // namespace xtk
