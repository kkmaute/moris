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

            // discretization mesh index
            // NOTE: this is the discretization mesh index as it is used by MSI 
            moris_index mDiscretizationMeshIndex;

            // type of Clusters in group
            mtk::Cluster_Type mClusterType = mtk::Cluster_Type::UNDEFINED;

            //------------------------------------------------------------------------------

          public:

            //------------------------------------------------------------------------------

            /**
             * @brief Constructor
             * 
             * @param aDiscretizationMeshIndex discretization mesh index (in MSI) that the cluster group is associated with
             */
            Cluster_Group( const moris_index aDiscretizationMeshIndex );

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
             * @brief Get the Bspline mesh list index which the cluster group lives on
             * 
             * @return moris_index Bspline mesh list index
             */
            moris_index 
            get_discretization_mesh_index_for_cluster_group() const;

            //------------------------------------------------------------------------------

            /**
             * @brief Get the Bspline mesh list index which the cluster group lives on
             * 
             * @return moris_index Bspline mesh list index
             */
            mtk::Cluster_Type 
            get_cluster_type_in_group() const;

            //------------------------------------------------------------------------------
            // Pure Virtual Functions Handled by this Class's Children
            //------------------------------------------------------------------------------

            /**
             * @brief Compute the sum of all cluster volumes within the cluster group
             * 
             * @param aPrimaryOrVoid 
             * @param aIsLeader 
             * @return moris::real volume of all clusters 
             */
            virtual
            real
            compute_cluster_group_cell_measure(
                    const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const = 0;

            //------------------------------------------------------------------------------

            /**
             * @brief compute the derivative of the total cluster group volume
             * 
             * @param aPerturbedVertexCoords 
             * @param aDirection 
             * @param aPrimaryOrVoid 
             * @param aIsLeader 
             * @return moris::real 
             */            
            virtual
            real
            compute_cluster_group_cell_measure_derivative(
                    const Matrix< DDRMat > & aPerturbedVertexCoords,
                    uint aDirection,
                    const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const = 0;

            //------------------------------------------------------------------------------
            
            /**
             * @brief compute the total cluster group interface/boundary surface/length
             * 
             * @param aPrimaryOrVoid 
             * @param aIsLeader 
             * @return moris::real 
             */
            virtual
            real
            compute_cluster_group_side_measure(
                    const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const = 0;

            //------------------------------------------------------------------------------
            
            /**
             * @brief compute the derivative of the total cluster group interface/boundary surface/length
             * 
             * @param aPerturbedVertexCoords 
             * @param aDirection 
             * @param aPrimaryOrVoid 
             * @param aIsLeader 
             * @return moris::real 
             */
            virtual
            real
            compute_cluster_group_side_measure_derivative(
                    const Matrix< DDRMat > & aPerturbedVertexCoords,
                    uint aDirection,
                    const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const = 0;

            //------------------------------------------------------------------------------

        };    // class mtk::Cluster_Group

    }    // namespace mtk
}    // namespace moris

#endif /* cl_MTK_Cluster_Group.hpp */