/*
 * Copyright (c) 2022 University of Colorado 
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. 
 * 
 * ------------------------------------------------------------------------------------ 
 * 
 * cl_MTK_Cluster_Group.hpp  
 * 
 */
#ifndef SRC_cl_MTK_Cluster_Group
#define SRC_cl_MTK_Cluster_Group

#include "cl_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_Matrix.hpp"

namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------

        // forward declare the mtk::Cluster
        class Cluster;

        //------------------------------------------------------------------------------

        class Cluster_Group
        {
            //------------------------------------------------------------------------------

          private:

            // list of clusters in group
            moris::Cell< std::shared_ptr< Cluster > > mClusters;

            // index on which B-spline mesh the cluster group is valid
            // NOTE: this is the index of the B-spline mesh wrt. to the list of B-spline meshes enriched using the current Lagrange mesh. 
            // NOTE: This corresponds to the mesh indices used by MSI.
            moris_index mBsplineMeshListIndex;

            //------------------------------------------------------------------------------

          public:
          
            Cluster_Group(
                    moris::Cell< std::shared_ptr< Cluster > > aClusters,
                    const moris_index                          aBsplineMeshListIndex );

            ~Cluster_Group(){};

            //------------------------------------------------------------------------------

            /**
             * @brief Get a the list of clusters in the cluster group
             * 
             * @return moris::Cell< Cluster const* > const& list of clusters in the cluster group
             */
            moris::Cell< std::shared_ptr< Cluster > > const& get_clusters_in_group() const;

            //------------------------------------------------------------------------------

            /**
             * @brief Get the Bspline mesh list index wich the cluster group lives on
             * 
             * @return moris_index Bspline mesh list index
             */
            moris_index get_Bspline_index_for_cluster_group() const;

            //------------------------------------------------------------------------------

            /**
             * @brief Compute the sum of all cluster volumes within the cluster group
             * 
             * @param aPrimaryOrVoid 
             * @param aIsMaster 
             * @return moris::real volume of all clusters 
             */
            moris::real
            compute_cluster_group_volume(
                    const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER ) const;

            //------------------------------------------------------------------------------

            /**
             * @brief compute the derivative of the total cluster group volume
             * 
             * @param aPerturbedVertexCoords 
             * @param aDirection 
             * @param aPrimaryOrVoid 
             * @param aIsMaster 
             * @return moris::real 
             */            
            moris::real
            compute_cluster_group_volume_derivative(
                    const Matrix< DDRMat > & aPerturbedVertexCoords,
                    uint aDirection,
                    const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER ) const;

            //------------------------------------------------------------------------------
            
            /**
             * @brief compute the total cluster group interface/boundary surface/length
             * 
             * @param aPrimaryOrVoid 
             * @param aIsMaster 
             * @return moris::real 
             */
            moris::real
            compute_cluster_group_side_measure(
                    const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER ) const;

            //------------------------------------------------------------------------------
            
            /**
             * @brief compute the derivative of the total cluster group interface/boundary surface/length
             * 
             * @param aPerturbedVertexCoords 
             * @param aDirection 
             * @param aPrimaryOrVoid 
             * @param aIsMaster 
             * @return moris::real 
             */
            moris::real
            compute_cluster_group_side_measure_derivative(
                    const Matrix< DDRMat > & aPerturbedVertexCoords,
                    uint aDirection,
                    const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER ) const;

            //------------------------------------------------------------------------------

        };    // class Cluster

    }    // namespace mtk
}    // namespace moris

#endif /* cl_MTK_Cluster_Group.hpp */