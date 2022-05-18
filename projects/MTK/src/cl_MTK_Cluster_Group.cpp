/*
 * Copyright (c) 2022 University of Colorado 
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. 
 * 
 * ------------------------------------------------------------------------------------ 
 * 
 * cl_MTK_Cluster_Group.cpp  
 * 
 */

#include "cl_MTK_Cluster_Group.hpp"
#include "cl_MTK_Cluster.hpp"

namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------

        Cluster_Group::Cluster_Group(
                moris::Cell< std::shared_ptr< Cluster > > aClusters,
                const moris_index              aBsplineMeshListIndex )
                : mClusters( aClusters )
                , mBsplineMeshListIndex( aBsplineMeshListIndex )
        {
            // only initialize the cluster group with the list of clusters, do nothing else in constructor
        }

        //------------------------------------------------------------------------------

        moris::Cell< std::shared_ptr< Cluster > > const&
        Cluster_Group::get_clusters_in_group() const
        {
            return mClusters;
        }

        //------------------------------------------------------------------------------

        moris::real
        Cluster_Group::compute_cluster_group_volume(
                const mtk::Primary_Void aPrimaryOrVoid,
                const mtk::Master_Slave aIsMaster ) const
        {
            // initialize volume measure
            real tVolume = 0.0;

            // sum up volume over all clusters 
            for( uint iCluster = 0; iCluster < mClusters.size(); iCluster++ )
            {
                tVolume += mClusters( iCluster )->compute_cluster_cell_measure( aPrimaryOrVoid, aIsMaster );
            }

            // return sum
            return tVolume;
        }

        //------------------------------------------------------------------------------
        
        moris::real
        Cluster_Group::compute_cluster_group_volume_derivative(
                const Matrix< DDRMat > & aPerturbedVertexCoords,
                uint aDirection,
                const mtk::Primary_Void aPrimaryOrVoid,
                const mtk::Master_Slave aIsMaster ) const
        {
            MORIS_ERROR( false, 
                "Cluster_Group::compute_cluster_group_volume_derivative() - not implemented yet." );

            return 0.0;
        }

        //------------------------------------------------------------------------------
        
        moris::real
        Cluster_Group::compute_cluster_group_side_measure(
                const mtk::Primary_Void aPrimaryOrVoid,
                const mtk::Master_Slave aIsMaster ) const
        {
            // initialize interface/boundary surface/length
            real tInterfaceSize = 0.0;

            // sum up interface/boundary surface/length over all clusters 
            for( uint iCluster = 0; iCluster < mClusters.size(); iCluster++ )
            {
                tInterfaceSize += mClusters( iCluster )->compute_cluster_cell_side_measure( aPrimaryOrVoid, aIsMaster );
            }

            // return sum
            return tInterfaceSize;
        }

        //------------------------------------------------------------------------------
        
        moris::real
        Cluster_Group::compute_cluster_group_side_measure_derivative(
                const Matrix< DDRMat > & aPerturbedVertexCoords,
                uint aDirection,
                const mtk::Primary_Void aPrimaryOrVoid,
                const mtk::Master_Slave aIsMaster ) const
        {
            MORIS_ERROR( false, 
                "Cluster_Group::compute_cluster_group_side_measure_derivative() - not implemented yet." );

            return 0.0;
        }


        //------------------------------------------------------------------------------

    } // namespace moris
} // namespace mtk