/*
 * Copyright (c) 2022 University of Colorado 
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. 
 * 
 * ------------------------------------------------------------------------------------ 
 * 
 * cl_ MTK_Side_Cluster_Group_DataBase.cpp  
 * 
 */

#include "cl_MTK_Side_Cluster_Group_DataBase.hpp"
#include "cl_MTK_Cluster.hpp"

namespace moris
{
    namespace mtk
    {
    //------------------------------------------------------------------------------

    Side_Cluster_Group_DataBase::Side_Cluster_Group_DataBase( 
            const moris_index                  aDiscretizationMeshIndex,
            moris::Vector< mtk::Cluster const* > aClusters,
            mtk::Cluster_Group const*          aAssociatedCellClusterGroup )
            : mtk::Side_Cluster_Group( aDiscretizationMeshIndex )
            , mSideClusters( aClusters )
    {
        // check tath the pointer passed is valid
        MORIS_ASSERT( aAssociatedCellClusterGroup != nullptr, 
            "xtk::Side_Cluster_Group_DataBase::Side_Cluster_Group_DataBase() - nullptr passed for the cell cluster group" );

        // assign pointer to the associated cell cluster group
        mAssociatedCellClusterGroup = aAssociatedCellClusterGroup;
    }

    //------------------------------------------------------------------------------

    const moris::Vector< mtk::Cluster const* >
    Side_Cluster_Group_DataBase::get_side_clusters_in_group() const
    {
        return mSideClusters;
    }

    //------------------------------------------------------------------------------

    mtk::Cluster_Group const* 
    Side_Cluster_Group_DataBase::get_associated_cell_cluster_group() const
    {
        return mAssociatedCellClusterGroup;
    }

    //------------------------------------------------------------------------------

    } // namespace mtk
} // namespace moris
