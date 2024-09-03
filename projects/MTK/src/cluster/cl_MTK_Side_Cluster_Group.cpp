/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_Side_Cluster_Group.cpp
 *
 */

#include "cl_MTK_Side_Cluster_Group.hpp"
#include "cl_MTK_Cluster.hpp"

namespace moris::mtk
{
    //------------------------------------------------------------------------------

    Side_Cluster_Group::Side_Cluster_Group( const moris_index aDiscretizationMeshIndex )
            : Cluster_Group( aDiscretizationMeshIndex )
    {
        mClusterType = mtk::ClusterType::SIDE;
    }

    //------------------------------------------------------------------------------

    moris::real
    Side_Cluster_Group::compute_cluster_group_cell_measure(
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        // access associated cell cluster group
        mtk::Cluster_Group const * tAssociatedCellClusterGroup = this->get_associated_cell_cluster_group();

        // compute the volume derivative on the associated cell cluster group
        return tAssociatedCellClusterGroup->compute_cluster_group_cell_measure( aPrimaryOrVoid, aIsLeader );
    }

    //------------------------------------------------------------------------------

    moris::real
    Side_Cluster_Group::compute_cluster_group_cell_measure_derivative(
            const Matrix< DDRMat >&    aPerturbedVertexCoords,
            uint                       aDirection,
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        // access associated cell cluster group
        mtk::Cluster_Group const * tAssociatedCellClusterGroup = this->get_associated_cell_cluster_group();

        // compute the volume derivative on the associated cell cluster group
        return tAssociatedCellClusterGroup->compute_cluster_group_cell_measure_derivative(
                aPerturbedVertexCoords,
                aDirection,
                aPrimaryOrVoid,
                aIsLeader );
    }

    //------------------------------------------------------------------------------

    moris::real
    Side_Cluster_Group::compute_cluster_group_side_measure(
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        // initialize interface area/length
        real tInterfaceArea = 0.0;

        // get the clusters in this group
        Vector< mtk::Cluster const * > tMySideClusters = this->get_side_clusters_in_group();

        // sum up the interface area/length over all clusters
        for ( uint iCluster = 0; iCluster < tMySideClusters.size(); iCluster++ )
        {
            tInterfaceArea += tMySideClusters( iCluster )->compute_cluster_cell_side_measure( aPrimaryOrVoid, aIsLeader );
        }

        // return sum
        return tInterfaceArea;
    }

    //------------------------------------------------------------------------------

    moris::real
    Side_Cluster_Group::compute_cluster_group_side_measure_derivative(
            const Matrix< DDRMat >&    aPerturbedVertexCoords,
            uint                       aDirection,
            const mtk::Primary_Void    aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        // TODO
        MORIS_ERROR( false, "Side_Cluster_Group::compute_cluster_group_side_measure_derivative() - Not implemented yet." );
        return 0.0;
    }

    //------------------------------------------------------------------------------

}    // namespace moris::mtk
