/*
 * Copyright (c) 2022 University of Colorado 
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. 
 * 
 * ------------------------------------------------------------------------------------ 
 * 
 * cl_XTK_Cell_Cluster_Group.hpp  
 * 
 */
#ifndef SRC_cl_XTK_Cell_Cluster_Group
#define SRC_cl_XTK_Cell_Cluster_Group

#include "cl_MTK_Cell_Cluster_Group.hpp"

namespace moris::mtk
{
    class Cluster;
}

using namespace moris;

namespace xtk
{
    //------------------------------------------------------------------------------

    class Cell_Cluster_Group : public mtk::Cell_Cluster_Group
    {
        //------------------------------------------------------------------------------

    protected:
    
        // list of clusters in group
        moris::Cell< std::shared_ptr< mtk::Cluster > > mClusters;

        //------------------------------------------------------------------------------

        /**
         * @brief Get a the list of clusters in the cluster group (how the clusters are accessed is handled by the children)
         * 
         * @return moris::Cell< Cluster const* > const& list of clusters in the cluster group
         */
        const moris::Cell< mtk::Cluster const* > 
        get_clusters_in_group() const;

        //------------------------------------------------------------------------------

    public:

        //------------------------------------------------------------------------------

        /**
         * @brief Construct a new Cell_Cluster_Group object
         * 
         * @param aDiscretizationMeshIndex discretization mesh index (in MSI) that the cluster group is associated with 
         * @param aClusters cell of pointers to the clusters in the group
         */
        Cell_Cluster_Group( 
                const moris_index                              aDiscretizationMeshIndex,
                moris::Cell< std::shared_ptr< mtk::Cluster > > aClusters );

        //------------------------------------------------------------------------------

        /**
         * @brief default constructor initializing nothing
         * 
         */
        Cell_Cluster_Group() = default;

        /**
         * @brief Default Destructor
         * 
         */
        virtual
        ~Cell_Cluster_Group() = default;


    }; // class xtk::Cell_Cluster_Group

} // namespace xtk

#endif /* cl_XTK_Cell_Cluster_Group.hpp */