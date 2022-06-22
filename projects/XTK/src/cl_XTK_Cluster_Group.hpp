/*
 * Copyright (c) 2022 University of Colorado 
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. 
 * 
 * ------------------------------------------------------------------------------------ 
 * 
 * cl_XTK_Cluster_Group.hpp  
 * 
 */
#ifndef SRC_cl_XTK_Cluster_Group
#define SRC_cl_XTK_Cluster_Group

#include "cl_MTK_Cluster_Group.hpp"

#include "cl_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_Matrix.hpp"

using namespace moris;

namespace xtk
{
    //------------------------------------------------------------------------------

    class Cluster_Group : public mtk::Cluster_Group
    {
        //------------------------------------------------------------------------------

      protected:

        // list of clusters in group
        moris::Cell< std::shared_ptr< mtk::Cluster > > mClusters;

        // cluster type (mClusterType) and B-spline mesh index (mBsplineMeshListIndex) inherited as member variables from the mtk::Cluster_Group

        //------------------------------------------------------------------------------

      public:
      
        Cluster_Group(
                moris::Cell< std::shared_ptr< mtk::Cluster > > aClusters,
                const moris_index                              aBsplineMeshListIndex,
                const mtk::Cluster_Type                        aClusterType );

        /**
         * @brief default constructor initializing nothing
         * 
         */
        Cluster_Group() = default;

        /**
         * @brief Default Destructor
         * 
         */
        virtual
        ~Cluster_Group() = default;

        //------------------------------------------------------------------------------

        /**
         * @brief Get a the list of clusters in the cluster group
         * 
         * @return moris::Cell< Cluster const* > const& list of clusters in the cluster group
         */
        moris::Cell< std::shared_ptr< mtk::Cluster > > const& get_clusters_in_group() const;

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

    }; // class xtk::Cluster

} // namespace xtk

#endif /* cl_XTK_Cluster_Group.hpp */