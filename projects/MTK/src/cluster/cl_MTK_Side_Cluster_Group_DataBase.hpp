/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_Side_Cluster_Group_DataBase.hpp
 *
 */
#ifndef SRC_cl_MTK_Side_Cluster_Group_DataBase
#define SRC_cl_MTK_Side_Cluster_Group_DataBase

#include "cl_MTK_Side_Cluster_Group.hpp"

namespace moris::mtk
{
    //------------------------------------------------------------------------------

    class Mesh;
    class Cluster;

    //------------------------------------------------------------------------------

    class Side_Cluster_Group_DataBase : public Side_Cluster_Group
    {
        //------------------------------------------------------------------------------

      protected:
        // pointer to parent mesh holding information
        mtk::Mesh* mMesh;

        // cluster group index in list of cluster groups for its B-spline mesh
        moris_index mClusterGroupIndex;

        // list of clusters in group
        Vector< mtk::Cluster const * > mSideClusters;

        // associated bulk cluster group
        mtk::Cluster_Group const * mAssociatedCellClusterGroup;

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
        Side_Cluster_Group_DataBase(
                const moris_index                     aDiscretizationMeshIndex,
                const Vector< mtk::Cluster const * >& aClusters,
                mtk::Cluster_Group const *            aAssociatedCellClusterGroup );

        //------------------------------------------------------------------------------

        /**
         * @brief default constructor initializing nothing
         *
         */
        Side_Cluster_Group_DataBase() = default;

        /**
         * @brief Default Destructor
         *
         */
        ~Side_Cluster_Group_DataBase() override = default;

        //------------------------------------------------------------------------------

    };    // class mtk::Side_Cluster_Group_DataBase

}    // namespace moris::mtk

#endif /* cl_MTK_Side_Cluster_Group_DataBase.hpp */
