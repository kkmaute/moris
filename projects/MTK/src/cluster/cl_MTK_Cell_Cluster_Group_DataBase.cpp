/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_Cell_Cluster_Group_DataBase.cpp
 *
 */

#include "cl_MTK_Cell_Cluster_Group_DataBase.hpp"
#include "cl_MTK_Cluster.hpp"

namespace moris::mtk
{
    //------------------------------------------------------------------------------

    Cell_Cluster_Group_DataBase::Cell_Cluster_Group_DataBase(
            const moris_index                     aDiscretizationMeshIndex,
            const Vector< mtk::Cluster const * >& aClusters )
            : mtk::Cell_Cluster_Group( aDiscretizationMeshIndex )
            , mClusters( aClusters )
    {
        // only initialize member variables
    }

    //------------------------------------------------------------------------------

    Vector< mtk::Cluster const * >
    Cell_Cluster_Group_DataBase::get_clusters_in_group() const
    {
        return mClusters;
    }

    //------------------------------------------------------------------------------

}    // namespace moris::mtk
