/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_XTK_Side_Cluster_Group.hpp
 *
 */
#ifndef SRC_cl_XTK_Side_Cluster_Group
#define SRC_cl_XTK_Side_Cluster_Group

#include "cl_MTK_Side_Cluster_Group.hpp"

namespace moris::mtk
{
    class Cluster;
}

using namespace moris;

namespace moris::xtk
{
    //------------------------------------------------------------------------------

    class Side_Cluster_Group : public mtk::Side_Cluster_Group
    {
        //------------------------------------------------------------------------------

      protected:
        // list of clusters in group
        Vector< std::shared_ptr< mtk::Cluster > > mSideClusters;

        // associated bulk cluster group
        std::shared_ptr< mtk::Cluster_Group > mAssociatedCellClusterGroup;

        //------------------------------------------------------------------------------

        /**
         * @brief Get a the list of clusters in the cluster group (how the clusters are accessed is handled by the children)
         *
         * @return Vector< Cluster const* > const& list of clusters in the cluster group
         */
        Vector< mtk::Cluster const * >
        get_side_clusters_in_group() const override;

        //------------------------------------------------------------------------------

        /**
         * @brief Get a the list of clusters in the cluster group (how the cell cluster group is accessed is handled by the children)
         *
         * @return mtk::Cluster_Group const* pointer to the bulk cluster group associated with the current side cluster group
         */
        mtk::Cluster_Group const *
        get_associated_cell_cluster_group() const override;

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------

        /**
         * @brief Constructor
         *
         * @param aDiscretizationMeshIndex discretization mesh index (in MSI) that the cluster group is associated with
         * @param aClusters
         * @param aAssociatedCellClusterGroup
         */
        Side_Cluster_Group(
                const moris_index                                aDiscretizationMeshIndex,
                const Vector< std::shared_ptr< mtk::Cluster > >& aClusters,
                const std::shared_ptr< mtk::Cluster_Group >&     aAssociatedCellClusterGroup );

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
        ~Side_Cluster_Group() override {};

        //------------------------------------------------------------------------------

    };    // class xtk::Side_Cluster_Group

}    // namespace moris::xtk

#endif /* cl_XTK_Side_Cluster_Group.hpp */
