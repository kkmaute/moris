/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_Side_Cluster_Group.hpp
 *
 */
#ifndef SRC_cl_MTK_Side_Cluster_Group
#define SRC_cl_MTK_Side_Cluster_Group

#include "cl_MTK_Cluster_Group.hpp"

namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------

        class Cluster;

        //------------------------------------------------------------------------------

        class Side_Cluster_Group : public Cluster_Group
        {
            //------------------------------------------------------------------------------

          protected:
            //------------------------------------------------------------------------------

            /**
             * @brief Get a the list of clusters in the cluster group (how the clusters are accessed is handled by the children)
             *
             * @return Vector< Cluster const* > const& list of clusters in the cluster group
             */
            virtual const Vector< mtk::Cluster const * >
            get_side_clusters_in_group() const = 0;

            //------------------------------------------------------------------------------

            /**
             * @brief Get a the list of clusters in the cluster group (how the cell cluster group is accessed is handled by the children)
             *
             * @return mtk::Cluster_Group const* pointer to the bulk cluster group associated with the current side cluster group
             */
            virtual mtk::Cluster_Group const *
            get_associated_cell_cluster_group() const = 0;

            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------

            /**
             * @brief Constructor
             *
             * @param aDiscretizationMeshIndex discretization mesh index (in MSI) that the cluster group is associated with
             */
            Side_Cluster_Group( const moris_index aDiscretizationMeshIndex );

            //------------------------------------------------------------------------------

            /**
             * @brief default constructor initializing nothing
             *
             */
            Side_Cluster_Group() = default;

            /**
             * @brief Default Destructor
             *
             */
            virtual ~Side_Cluster_Group() = default;

            //------------------------------------------------------------------------------

            /**
             * @brief Compute the sum of all cluster volumes within the cluster group
             *
             * @param aPrimaryOrVoid
             * @param aIsLeader
             * @return moris::real volume of all clusters
             */
            moris::real
            compute_cluster_group_cell_measure(
                    const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const;

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
            moris::real
            compute_cluster_group_cell_measure_derivative(
                    const Matrix< DDRMat >&    aPerturbedVertexCoords,
                    uint                       aDirection,
                    const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const;

            //------------------------------------------------------------------------------

            /**
             * @brief compute the total cluster group interface/boundary surface/length
             *
             * @param aPrimaryOrVoid
             * @param aIsLeader
             * @return moris::real
             */
            moris::real
            compute_cluster_group_side_measure(
                    const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const;

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
            moris::real
            compute_cluster_group_side_measure_derivative(
                    const Matrix< DDRMat >&    aPerturbedVertexCoords,
                    uint                       aDirection,
                    const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const;

            //------------------------------------------------------------------------------

        };    // class mtk::Side_Cluster_Group

    }    // namespace mtk
}    // namespace moris

#endif /* cl_MTK_Side_Cluster_Group.hpp */
