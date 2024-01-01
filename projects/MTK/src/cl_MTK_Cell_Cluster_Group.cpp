/*
 * Copyright (c) 2022 University of Colorado 
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. 
 * 
 * ------------------------------------------------------------------------------------ 
 * 
 * cl_MTK_Cell_Cluster_Group.cpp  
 * 
 */

#include "cl_MTK_Cell_Cluster_Group.hpp"
#include "cl_MTK_Cluster.hpp"

namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------

        Cell_Cluster_Group::Cell_Cluster_Group( const moris_index aDiscretizationMeshIndex )
                : Cluster_Group( aDiscretizationMeshIndex )
        {
            mClusterType = mtk::ClusterType::CELL;
        }

        //------------------------------------------------------------------------------

        moris::real
        Cell_Cluster_Group::compute_cluster_group_cell_measure(
                const mtk::Primary_Void aPrimaryOrVoid,
                const mtk::Leader_Follower aIsLeader ) const
        {
            // initialize volume measure
            real tVolume = 0.0;

            // get the clusters in this group
            moris::Vector< mtk::Cluster const* > const tMyClusters = this->get_clusters_in_group();

            // sum up volume over all clusters 
            for( uint iCluster = 0; iCluster < tMyClusters.size(); iCluster++ )
            {
                tVolume += tMyClusters( iCluster )->compute_cluster_cell_measure( aPrimaryOrVoid, aIsLeader );
            }

            // return sum
            return tVolume;
        }

        //------------------------------------------------------------------------------
        
        moris::real
        Cell_Cluster_Group::compute_cluster_group_cell_measure_derivative(
                const Matrix< DDRMat > & aPerturbedVertexCoords,
                uint aDirection,
                const mtk::Primary_Void aPrimaryOrVoid,
                const mtk::Leader_Follower aIsLeader ) const
        {
            // TODO
            MORIS_ERROR( false, "Cluster_Group_DataBase::compute_cluster_group_cell_measure_derivative() - not implemented yet." );
            return 0.0;
        }

        //------------------------------------------------------------------------------
        
        moris::real
        Cell_Cluster_Group::compute_cluster_group_side_measure(
                const mtk::Primary_Void aPrimaryOrVoid,
                const mtk::Leader_Follower aIsLeader ) const
        {
            // this function is not implemented in this 
            MORIS_ERROR( false, 
                "Cell_Cluster_Group::compute_cluster_group_side_measure() - "
                "This is a Cell (Bulk) Cluster Group. Side measure not implemented." );

            // return sum
            return 0.0;
        }

        //------------------------------------------------------------------------------
        
        moris::real
        Cell_Cluster_Group::compute_cluster_group_side_measure_derivative(
                const Matrix< DDRMat > & aPerturbedVertexCoords,
                uint aDirection,
                const mtk::Primary_Void aPrimaryOrVoid,
                const mtk::Leader_Follower aIsLeader ) const
        {
            // this function is not implemented in this 
            MORIS_ERROR( false, 
                "Cell_Cluster_Group::compute_cluster_group_side_measure_derivative() - "
                "This is a Cell (Bulk) Cluster Group. Side measure derivative not implemented." );

            // return sum
            return 0.0;
        }

        //------------------------------------------------------------------------------

    } // namespace mtk
}// namespace moris
