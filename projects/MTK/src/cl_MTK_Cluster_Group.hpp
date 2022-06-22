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

          protected:

            // index on which B-spline mesh the cluster group is valid
            // NOTE: this is the index of the B-spline mesh wrt. to the list of B-spline meshes enriched using the current Lagrange mesh. 
            // NOTE: This corresponds to the mesh indices used by MSI.
            moris_index mBsplineMeshListIndex;

            // type of Clusters in group
            mtk::Cluster_Type mClusterType;

            //------------------------------------------------------------------------------

          public:

            //------------------------------------------------------------------------------

            /**
             * @brief Constructor
             * 
             * @param aBsplineMeshListIndex 
             * @param aClusterType 
             */
            Cluster_Group(
                    const moris_index  aBsplineMeshListIndex,
                    const Cluster_Type aClusterType );

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
             * @brief Get the Bspline mesh list index wich the cluster group lives on
             * 
             * @return moris_index Bspline mesh list index
             */
            moris_index 
            get_Bspline_index_for_cluster_group() const;

            //------------------------------------------------------------------------------

            /**
             * @brief Get the Bspline mesh list index wich the cluster group lives on
             * 
             * @return moris_index Bspline mesh list index
             */
            mtk::Cluster_Type 
            get_cluster_type_in_group() const;

            //------------------------------------------------------------------------------

            /**
             * @brief Get a the list of clusters in the cluster group
             * 
             * @return moris::Cell< Cluster const* > const& list of clusters in the cluster group
             */
            virtual
            moris::Cell< std::shared_ptr< Cluster > > const& 
            get_clusters_in_group() const = 0;

            //------------------------------------------------------------------------------

            /**
             * @brief Compute the sum of all cluster volumes within the cluster group
             * 
             * @param aPrimaryOrVoid 
             * @param aIsMaster 
             * @return moris::real volume of all clusters 
             */
            virtual
            moris::real
            compute_cluster_group_volume(
                    const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER ) const = 0;

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
            virtual
            moris::real
            compute_cluster_group_volume_derivative(
                    const Matrix< DDRMat > & aPerturbedVertexCoords,
                    uint aDirection,
                    const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER ) const = 0;

            //------------------------------------------------------------------------------
            
            /**
             * @brief compute the total cluster group interface/boundary surface/length
             * 
             * @param aPrimaryOrVoid 
             * @param aIsMaster 
             * @return moris::real 
             */
            virtual
            moris::real
            compute_cluster_group_side_measure(
                    const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER ) const = 0;

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
            virtual
            moris::real
            compute_cluster_group_side_measure_derivative(
                    const Matrix< DDRMat > & aPerturbedVertexCoords,
                    uint aDirection,
                    const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER ) const = 0;

            //------------------------------------------------------------------------------

        };    // class mtk::Cluster_Group

    }    // namespace mtk
}    // namespace moris

#endif /* cl_MTK_Cluster_Group.hpp */