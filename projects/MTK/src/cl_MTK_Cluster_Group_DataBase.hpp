/*
 * Copyright (c) 2022 University of Colorado 
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. 
 * 
 * ------------------------------------------------------------------------------------ 
 * 
 * cl_MTK_Cluster_Group_DataBase.hpp
 * 
 */
#ifndef SRC_cl_MTK_Cell_Cluster_DataBase
#define SRC_cl_MTK_Cell_Cluster_DataBase

#include "cl_MTK_Cluster_Group.hpp"

namespace moris
{
    namespace mtk
    {
        // forward declare the mesh for access
        class Mesh;

        //------------------------------------------------------------------------------

        class Cluster_Group_DataBase : public Cluster_Group
        {

            //------------------------------------------------------------------------------
            
        private:

            // pointer to parent mesh holding information
            mtk::Mesh* mMesh;

            // cluster group index in list of cluster groups for its B-spline mesh
            moris_index mClusterGroupIndex;

            // cluster type (mClusterType) and B-spline mesh index (mBsplineMeshListIndex) inherited as member variables from the mtk::Cluster_Group

            // FIXME: temporary placeholder, remove later
            moris::Cell< std::shared_ptr< mtk::Cluster > > mDummy;

            //------------------------------------------------------------------------------

        public:

            //------------------------------------------------------------------------------

            /**
             * @brief Constructor
             * 
             * @param aBsplineMeshListIndex 
             * @param aClusterType 
             */
            Cluster_Group_DataBase(
                    const moris_index       aBsplineMeshListIndex,
                    const enum Cluster_Type aClusterType );

            //------------------------------------------------------------------------------

            /**
             * @brief Default Destructor
             * 
             */
            ~Cluster_Group_DataBase() = default;

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

        }; // class mtk::Cluster_Group_DataBase

    } // namespace moris::mtk
} // namespace moris::mtk

#endif /* cl_MTK_Cell_Cluster_DataBase.hpp */