/*
 * Copyright (c) 2022 University of Colorado 
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. 
 * 
 * ------------------------------------------------------------------------------------ 
 * 
 * cl_MTK_Cluster_Group.cpp  
 * 
 */

#include "cl_MTK_Cluster_Group.hpp"
#include "cl_MTK_Cluster.hpp"

namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------

        Cluster_Group::Cluster_Group(
                const moris_index  aBsplineMeshListIndex,
                const Cluster_Type aClusterType )
                : mBsplineMeshListIndex( aBsplineMeshListIndex )
                , mClusterType( aClusterType )
        {
            // only initialize member variables
        }


        //------------------------------------------------------------------------------

        moris_index 
        Cluster_Group::get_Bspline_index_for_cluster_group() const
        {
            return mBsplineMeshListIndex;
        }

        //------------------------------------------------------------------------------

        mtk::Cluster_Type
        Cluster_Group::get_cluster_type_in_group() const
        {
            return mClusterType;
        }

        //------------------------------------------------------------------------------

    } // namespace moris
} // namespace mtk