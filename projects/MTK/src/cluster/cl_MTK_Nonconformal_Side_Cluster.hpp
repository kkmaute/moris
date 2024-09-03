//
// Created by frank on 11/30/23.
//

#pragma once

#include "cl_MTK_Double_Side_Cluster.hpp"
#include "cl_MTK_PointPairs.hpp"
#include "cl_Vector.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include <map>
#include <utility>

namespace moris::mtk
{
    class Nonconformal_Side_Cluster : public Double_Side_Cluster
    {
      public:
        Nonconformal_Side_Cluster(
                Cluster const                         *aLeaderSideCluster,
                Cluster const                         *aFollowerSideCluster,
                Vector< IntegrationPointPairs > const &aIntegrationPointPairs,
                Vector< NodalPointPairs > const       &aNodalPointPairs )
                : Double_Side_Cluster( aLeaderSideCluster, aFollowerSideCluster, {} )
                , mIntegrationPointPairs( aIntegrationPointPairs )
                , mNodalPointPairs( aNodalPointPairs ){};

        /**
         * \brief Works similar to the function get_primary_cells_in_cluster() from the base class. The difference is, that
         * the function might return some cells multiple times if they are part of multiple integration point pairs.
         * \example For two leader cell with index 1 and 2 and two follower cells with index 3 and 4. If you have the following
         * integration point pairs:
         *  1 -> 3
         *  1 -> 4
         *  2 -> 3
         *  the function would return a vector with the following cells:
         *    Leader:  1, 1, 2
         *    Followe: 3, 4, 3
         * \param aIsLeader
         * \return A vector with (possibly duplicate) cells that are consistent with the integration pairing.
         */
        Vector< mtk::Cell const * >
        get_nonconforming_primary_cells_in_cluster( const mtk::Leader_Follower aIsLeader ) const;

        /**
         * \brief Similar to the function get_nonconforming_primary_cells_in_cluster(), this function returns cell side ordinals while
         * being consistent with the integration point pairs. I.e. a cell side ordinal might be returned multiple times if it is part of
         * multiple integration point pairs.
         * \param aIsLeader
         * \return
         */
        moris::Matrix< moris::IndexMat >
        get_nonconforming_cell_side_ordinals( const mtk::Leader_Follower aIsLeader ) const;

        [[nodiscard]] Vector< IntegrationPointPairs > const &get_integration_point_pairs() const
        {
            return mIntegrationPointPairs;
        }

        [[nodiscard]] Vector< NodalPointPairs > const &get_nodal_point_pairs() const
        {
            return mNodalPointPairs;
        }

        /**
         * @brief This type is used to describe the cell pairings and the corresponding entries of the integration point pairs and nodal point pairs.
         */
        using NonconformalCellPairing = std::map< std::pair< moris_index, moris_index >, std::pair< moris_index, moris_index > >;

        /**
         * @brief Returns the nonconformal cell pairing and the correct indices of the integration point pairs and nodal point pairs that are used for the pairing.
         * @details Each cell might be paired up with none, one or multiple other cells on the follower side. The pairing between the two cells is given in the
         * first pair (key) of the map in the order <local_leader_cell_index, local_follower_cell_index> (it is given as the <em>local index</em>, i.e. at which
         * position in the list of cells in the cluster, the corresponding entry is stored). The second pair (value) gives the corresponding entries of the
         * integration point pairs and nodal point pairs respectively that are used (reason) for the pairing of the two cells. If either the integration point pairs
         * or the nodal point pairs are not used for the pairing, the corresponding entry is -1.
         * @param aIsLeader
         * @return
         */
        NonconformalCellPairing get_nonconformal_cell_pairing() const;

      private:
        /**
         * @brief The pairs of points that will be used for integration. Each point on the leader side will be paired with a point on the follower side.
         */
        Vector< IntegrationPointPairs > mIntegrationPointPairs;

        /**
         * @brief Maps each node of the cells of the leader side to a corresponding cell and coordinate on the follower side.
         * This is required for the nodal evaluation of IQIs on nonconforming interfaces.
         */
        Vector< NodalPointPairs > mNodalPointPairs;
    };
}    // namespace moris::mtk
