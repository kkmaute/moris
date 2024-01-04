//
// Created by frank on 11/30/23.
//

#ifndef MORIS_CL_MTK_NONCONFORMAL_SIDE_CLUSTER_HPP
#define MORIS_CL_MTK_NONCONFORMAL_SIDE_CLUSTER_HPP

#include "cl_MTK_Double_Side_Cluster.hpp"
#include "cl_MTK_IntegrationPointPairs.hpp"
#include "cl_Vector.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"

namespace moris::mtk
{
    class Nonconformal_Side_Cluster : public Double_Side_Cluster
    {
      public:
        Nonconformal_Side_Cluster(
                Cluster const                         *aLeaderSideCluster,
                Cluster const                         *aFollowerSideCluster,
                Vector< IntegrationPointPairs > const &aIntegrationPointPairs )
                : Double_Side_Cluster( aLeaderSideCluster, aFollowerSideCluster, {} )
                , mIntegrationPointPairs( aIntegrationPointPairs ){};


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

        Matrix< DDRMat > compute_cluster_ig_cell_side_measures( mtk::Primary_Void const aPrimaryOrVoid, mtk::Leader_Follower const aIsLeader ) const override;

        moris::real compute_cluster_cell_side_measure( mtk::Primary_Void const aPrimaryOrVoid, mtk::Leader_Follower const aIsLeader ) const override;

        /**
         * \brief Similar to the function get_nonconforming_primary_cells_in_cluster(), this function returns cell side ordinals while
         * being consistent with the integration point pairs. I.e. a cell side ordinal might be returned multiple times if it is part of
         * multiple integration point pairs.
         * \param aIsLeader
         * \return
         */
        moris::Matrix< moris::IndexMat >
        get_nonconforming_cell_side_ordinals( const mtk::Leader_Follower aIsLeader ) const;

        [[nodiscard]] Vector< IntegrationPointPairs > const &get_integration_point_pairs() const;

        /**
         * \brief This utility method returns a mask of indices that are consistent with the integration point pairs.
         * This can be useful to index into the cells of the cluster in the same order as the integration point pairs.
         * \details E.g. if the leader/follower sides of the integration point pairs are the cells 2, 4, 1, 1 (in this order) and the cells in the cluster are
         * 1, 2, 3, 4 (i.e cell 3 does not get paired with any other cell and 1 is paired twice), the mask would be [1, 3, 0, 0] (corresponding to the 1st, 3rd and
         * 0th cell, thus giving the cells with indices 2, 4, 1 and 1). This works for both leader and follower (primary) cells.
         * \param aIsLeader
         * \return
         */
        Vector< moris_index > get_cell_local_indices( const mtk::Leader_Follower aIsLeader ) const;

      private:
        Vector< IntegrationPointPairs > mIntegrationPointPairs;
    };


}    // namespace moris::mtk


#endif    // MORIS_CL_MTK_NONCONFORMAL_SIDE_CLUSTER_HPP
