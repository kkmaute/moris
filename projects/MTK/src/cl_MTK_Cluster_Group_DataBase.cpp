/*
 * Copyright (c) 2022 University of Colorado 
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. 
 * 
 * ------------------------------------------------------------------------------------ 
 * 
 * cl_MTK_Cluster_Group_DataBase.cpp  
 * 
 */

#include "cl_MTK_Cluster_Group_DataBase.hpp"

namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------

        Cluster_Group_DataBase::Cluster_Group_DataBase(
                const moris_index            aBsplineMeshListIndex,
                const enum mtk::Cluster_Type aClusterType )
                : Cluster_Group( aBsplineMeshListIndex, aClusterType )
        {
            // only initialize member variables
        }

        //------------------------------------------------------------------------------

        moris::Cell< std::shared_ptr< mtk::Cluster > > const&
        Cluster_Group_DataBase::get_clusters_in_group() const
        {
            MORIS_ERROR( false, "Cluster_Group_DataBase::get_clusters_in_group() - Not implemented yet." );
            return mDummy;
        }

        //------------------------------------------------------------------------------

        moris::real
        Cluster_Group_DataBase::compute_cluster_group_volume(
                const mtk::Primary_Void aPrimaryOrVoid,
                const mtk::Master_Slave aIsMaster ) const
        {
            // initialize volume measure
            real tVolume = 0.0;

            // sum up volume over all clusters 
            // for( uint iCluster = 0; iCluster < mClusters.size(); iCluster++ )
            // {
            //     tVolume += mClusters( iCluster )->compute_cluster_cell_measure( aPrimaryOrVoid, aIsMaster );
            // }

            MORIS_ERROR( false, "Cluster_Group_DataBase::compute_cluster_group_volume() - Not implemented yet." );

            // return sum
            return tVolume;
        }

        //------------------------------------------------------------------------------
        
        moris::real
        Cluster_Group_DataBase::compute_cluster_group_volume_derivative(
                const Matrix< DDRMat > & aPerturbedVertexCoords,
                uint aDirection,
                const mtk::Primary_Void aPrimaryOrVoid,
                const mtk::Master_Slave aIsMaster ) const
        {
            MORIS_ERROR( false, "Cluster_Group_DataBase::compute_cluster_group_volume_derivative() - not implemented yet." );
            return 0.0;
        }

        //------------------------------------------------------------------------------
        
        moris::real
        Cluster_Group_DataBase::compute_cluster_group_side_measure(
                const mtk::Primary_Void aPrimaryOrVoid,
                const mtk::Master_Slave aIsMaster ) const
        {
            // initialize interface/boundary surface/length
            real tInterfaceSize = 0.0;

            // sum up interface/boundary surface/length over all clusters 
            // for( uint iCluster = 0; iCluster < mClusters.size(); iCluster++ )
            // {
            //     tInterfaceSize += mClusters( iCluster )->compute_cluster_cell_side_measure( aPrimaryOrVoid, aIsMaster );
            // }

            MORIS_ERROR( false, "Cluster_Group_DataBase::compute_cluster_group_side_measure() - Not implemented yet." );

            // return sum
            return tInterfaceSize;
        }

        //------------------------------------------------------------------------------
        
        moris::real
        Cluster_Group_DataBase::compute_cluster_group_side_measure_derivative(
                const Matrix< DDRMat > & aPerturbedVertexCoords,
                uint aDirection,
                const mtk::Primary_Void aPrimaryOrVoid,
                const mtk::Master_Slave aIsMaster ) const
        {
            MORIS_ERROR( false, "Cluster_Group_DataBase::compute_cluster_group_side_measure_derivative() - not implemented yet." );
            return 0.0;
        }

        //------------------------------------------------------------------------------

    } // namespace mtk
}// namespace moris
