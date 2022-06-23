/*
 * Copyright (c) 2022 University of Colorado 
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. 
 * 
 * ------------------------------------------------------------------------------------ 
 * 
 * cl_MTK_Cell_Cluster_Group_DataBase.hpp  
 * 
 */
#ifndef SRC_cl_MTK_Cell_Cluster_Group_DataBase
#define SRC_cl_MTK_Cell_Cluster_Group_DataBase

#include "cl_MTK_Cell_Cluster_Group.hpp"

namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------

        class Mesh;
        class Cluster;

        //------------------------------------------------------------------------------

        class Cell_Cluster_Group_DataBase : public Cell_Cluster_Group
        {
            //------------------------------------------------------------------------------

        protected:
        
            // pointer to parent mesh holding information
            mtk::Mesh* mMesh;

            // cluster group index in list of cluster groups for its B-spline mesh
            moris_index mClusterGroupIndex;

            // list of clusters in group
            moris::Cell< mtk::Cluster const* > mClusters;

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
            Cell_Cluster_Group_DataBase( 
                    const moris_index                  aDiscretizationMeshIndex,
                    moris::Cell< mtk::Cluster const* > aClusters );

            //------------------------------------------------------------------------------

            /**
             * @brief default constructor initializing nothing
             * 
             */
            Cell_Cluster_Group_DataBase() = default;

            /**
             * @brief Default Destructor
             * 
             */
            virtual
            ~Cell_Cluster_Group_DataBase() = default;

            //------------------------------------------------------------------------------

        }; // class mtk::Cell_Cluster_Group_DataBase

    } // namespace xtk
} // namespace moris

#endif /* cl_MTK_Cell_Cluster_Group_DataBase.hpp */